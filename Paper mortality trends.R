#---------------------------#
# MORTALITY TRENDS IN CHILE #
#---------------------------#

required_packages <- c("dplyr", "readr", "haven", 
                       "ggplot2", "gridExtra", 
                       "fitdistrplus", "MASS")

# Install and load required packages
sapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
})

data <- readRDS("ENPG_FULL.RDS") %>% 
  mutate(aux = ifelse(oh1 == "No" & !is.na(oh2) ,1,0)) %>% 
  filter(aux == 0) %>% 
  dplyr::select(-aux) 

data_input <- data %>% 
  mutate(edad_tramo = case_when(between(edad, 15, 29)~1,
                                between(edad, 30,44)~2,
                                between(edad,45,59)~3,
                                between(edad,60,65)~4),
         cvolaj = case_when(oh1 == "No" ~ "ltabs",
                            oh2 == ">30" | oh2 == ">1 año" ~ "fd",
                            sexo == "Mujer" & volajohdia > 0 & volajohdia <= 19.99 ~ "cat1",
                            sexo == "Mujer" & volajohdia >= 20 & volajohdia <= 39.99 ~ "cat2",
                            sexo == "Mujer" & volajohdia >= 40 & volajohdia <= 60  ~ "cat3",
                            sexo == "Mujer" & volajohdia > 60 ~ "cat4",
                            sexo == "Hombre" & volajohdia > 0 & volajohdia <= 39.99 ~ "cat1",
                            sexo == "Hombre" & volajohdia >= 40 & volajohdia <= 59.99 ~ "cat2",
                            sexo == "Hombre" & volajohdia >= 60 & volajohdia<= 100 ~ "cat3",
                            sexo == "Hombre" & volajohdia > 100 ~ "cat4",
                            TRUE ~ NA),
         hed = ifelse(dias_binge > 0,1,0)) %>% 
  dplyr::select(year, sexo, exp, edad_tramo, volajohdia, cvolaj,hed)


input <- data_input %>% 
  filter(!is.na(cvolaj)) %>% 
  group_by(sexo,year, edad_tramo, cvolaj) %>% 
  summarise(weighted_count = sum(exp, na.rm = TRUE)) %>% 
  mutate(prop = round(weighted_count / sum(weighted_count), 2)) %>% 
  dplyr::select(-weighted_count)

data_hed <- data_input %>% 
  filter(!is.na(hed)) %>% 
  group_by(year, sexo, edad_tramo) %>% 
  count(hed) %>% 
  mutate(prop_hed = n/sum(n))

input_male <- input %>% 
  filter(sexo == "Hombre")
input_female <- input %>% 
  filter(sexo == "Mujer")

# Function to fit gamma distribution
# Define the x_vals range (adjust as needed)
x_vals <- seq(0.1, 150, length.out = 1500)

# TRAPEZOIDAL INTEGRATION FUNCTION
trap_int <- function(x, y, rr, prop_abs, rr_form, prop_form) {
  
  # Interval width
  dx <- x[2] - x[1]
  ncgamma <- sum(y[-1] + y[-length(x)]) * dx / 2
  # Normalize the gamma function
  normalized_y <- (1 -  (prop_abs+prop_form)) * y/ncgamma
  
  # Calculate the excess relative risk
  excess_rr <- rr - 1
  
  # Calculate the weighted excess relative risk
  weighted_excess_rr <- normalized_y * excess_rr
  
  # Apply the trapezoidal rule to calculate the numerator
  numerator <- (rr_form-1)*prop_form+sum((weighted_excess_rr[-1] + weighted_excess_rr[-length(x)]) / 2) * dx
  
  # Calculate the denominator
  denominator <- numerator + 1
  
  # Calculate PAF
  paf <- round(numerator / denominator, 3)
  
  return(paf)
}

# LINEAR RELATIVE RISK FUNCTION
rr_linear <- function(x, b){
  exp(x*b)
}

confint_paf <- function(gamma, beta, var_beta, p_abs, rr_form, p_form, rr_function){
  set.seed(145)
  n_sim <- 10000
  simulated_pafs <- numeric(n_sim)
  for (i in 1:n_sim) {
    pca_sim <- rgamma(1000, shape = gamma$estimate["shape"], rate = gamma$estimate["rate"])
    # Calculate the shape and rate parameters from the simulated PCA values
    mean_sim <- mean(pca_sim)
    sd_sim <- sd(pca_sim)
    
    # Calculate shape and rate for the new gamma distribution
    shape_sim <- (mean_sim / sd_sim)^2
    rate_sim <- mean_sim / (sd_sim^2)
    
    # Simulate the gamma distribution based on `k_sim` and `theta_sim`
    y_gamma_sim <- dgamma(x_vals, shape = shape_sim, rate = rate_sim)
    
    # Skip iteration if y_gamma_sim contains NaN values
    if (any(is.nan(y_gamma_sim))) next
    
    # Simulate beta coefficient
    beta_sim <- rnorm(1, beta, sqrt(var_beta))
    rr_sim <- rr_function(x_vals, beta_sim)
    # Simulate proportions of lifetime abstainers and former drinkers for the current age category
    prop_abs_sim <- max(rnorm(1, mean = p_abs, sd = sqrt(p_abs * (1 - p_abs) / 1000)),0.001)
    prop_form_sim <- max(rnorm(1, mean = p_form, sd = sqrt(p_form * (1 - p_form) / 1000)), 0.001)
    
    # Calculate the PAF using the trapezoidal method
    simulated_pafs[i] <- trap_int(x = x_vals, y = y_gamma_sim, rr = rr_sim, 
                                  prop_abs = prop_abs_sim, rr_form = rr_form, prop_form = prop_form_sim)
  }
  
  # Remove NaN values from simulated PAFs
  simulated_pafs <- simulated_pafs[!is.nan(simulated_pafs)]
  
  # Calculate the 95% confidence interval
  paf_lower <- quantile(simulated_pafs, 0.025)
  paf_upper <- quantile(simulated_pafs, 0.975)
  paf_point_estimate <- mean(simulated_pafs)
  
  return(list(
    Point_Estimate = round(paf_point_estimate,3),
    Lower_CI = paf_lower,
    Upper_CI = paf_upper)
  )
}

confint_paf_vcov <- function(gamma, betas, cov_matrix, p_abs, p_form, rr_fd, rr_function) {
  set.seed(145)
  n_sim <- 10000
  simulated_pafs <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    # Simulate PCA mean and SD using the gamma distribution
    pca_sim <- rgamma(1000, shape = gamma$estimate["shape"], rate = gamma$estimate["rate"])
    
    # Calculate the shape and rate parameters from the simulated PCA values
    mean_sim <- mean(pca_sim)
    sd_sim <- sd(pca_sim)
    
    # Calculate shape and rate for the new gamma distribution
    shape_sim <- (mean_sim / sd_sim)^2
    rate_sim <- mean_sim / (sd_sim^2)
    
    # Simulate the gamma distribution
    y_gamma_sim <- dgamma(x_vals, shape = shape_sim, rate = rate_sim)
    
    # Skip iteration if y_gamma_sim contains NaN values
    if (any(is.nan(y_gamma_sim))) next
    
    # Simulate the beta coefficients jointly from a multivariate normal distribution
    beta_sim <- MASS::mvrnorm(1, mu = betas, Sigma = cov_matrix)
    
    rr_sim <- rr_function(x_vals, beta_sim)
    
    # Simulate proportions of lifetime abstainers and former drinkers
    prop_abs_sim <- max(rnorm(1, mean = p_abs, sd = sqrt(p_abs * (1 - p_abs) / 1000)), 0.001)
    prop_form_sim <- max(rnorm(1, mean = p_form, sd = sqrt(p_form * (1 - p_form) / 1000)), 0.001)
    
    # Calculate the PAF using the trapezoidal method
    simulated_pafs[i] <- trap_int(x = x_vals, y = y_gamma_sim, rr = rr_sim, 
                                  prop_abs = prop_abs_sim, rr_form = rr_fd, prop_form = prop_form_sim)
  }
  
  # Remove NaN values from simulated PAFs
  simulated_pafs <- simulated_pafs[!is.nan(simulated_pafs)]
  
  # Calculate the 95% confidence interval
  paf_lower <- quantile(simulated_pafs, 0.025)
  paf_upper <- quantile(simulated_pafs, 0.975)
  paf_point_estimate <- mean(simulated_pafs)
  
  return(list(point_estimate = paf_point_estimate, lower_ci = paf_lower, upper_ci = paf_upper))
}
# Validate with integrate()


# FOR ISCHAEMIC STROKE AND INJURIES WE NEED THE PROPORTION OF CURRENT DRINKERS WHO
# ENGAGE IN HEAVY EPISODIC DRINKING NHED and HED
paf_hed_function <- function(x_60,x_150, y_nhed, y_hed_60, 
                             y_hed_150, rr_nhed, rr_hed_60,
                             rr_hed_150, rr_form,p_abs, 
                             p_form, p_hed) {
  
  trap_int_hed <- function(x, y, rr, prop_abs, rr_form, prop_form, p_hed) {
    
    dx <- x[2] - x[1]
    ncgamma <- sum((y[-1] + y[-length(y)]) / 2) * dx
    
    normalized_y <- ((1 - (prop_abs + prop_form)) * p_hed) * y / ncgamma
    
    excess_rr <- rr - 1
    
    weighted_excess_rr <- normalized_y * excess_rr
    
    numerator <- (rr_form - 1) * prop_form + sum((weighted_excess_rr[-1] + weighted_excess_rr[-length(weighted_excess_rr)]) / 2) * dx
    
    denominator <- numerator + 1
    
    paf <- round(numerator / denominator, 3)
    
    return(paf)
  }
  
  
  int_ri_nhed <- trap_int_hed(
    x = x_60, 
    y = y_nhed, 
    rr = rr_nhed, 
    prop_abs = p_abs, 
    rr_form = rr_form, 
    prop_form = p_form, 
    p_hed = 1 - p_hed
  )
  
  int_ri_hed1 <- trap_int_hed(
    x = x_60, 
    y = y_hed_60, 
    rr = rr_hed_60, 
    prop_abs = p_abs, 
    rr_form = rr_form, 
    prop_form = p_form, 
    p_hed = p_hed
  )
  
  int_ri_hed2 <- trap_int_hed(
    x = x_150, 
    y = y_hed_150, 
    rr = rr_hed_150, 
    prop_abs = p_abs, 
    rr_form = rr_form, 
    prop_form = p_form, 
    p_hed = p_hed
  )
  
  num <- (int_ri_nhed + int_ri_hed1 + int_ri_hed2)
  den <- num + 1
  paf_ri_fem1 <- num / den
  
  return(paf_ri_fem1)
}

confint_paf_hed <- function(gammas, beta, cov_matrix, p_abs, p_form, rr_fd, rr_function_nhed,
                            rr_function_hed, p_hed) {
  set.seed(145)
  n_sim <- 10000
  simulated_pafs <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    # Simulate PCA mean and SD using the gamma distribution
    pca_sim_nhed <- rgamma(1000, shape = gammas[[1]]$estimate["shape"], 
                           rate = gammas[[1]]$estimate["rate"])
    pca_sim_hed <- rgamma(1000, shape = gammas[[2]]$estimate["shape"], 
                          rate = gammas[[2]]$estimate["rate"])
    
    # Calculate the shape and rate parameters from the simulated PCA values
    mean_sim_nhed <- mean(pca_sim_nhed)
    sd_sim_nhed <- sd(pca_sim_nhed)
    mean_sim_hed <- mean(pca_sim_hed)
    sd_sim_hed <- sd(pca_sim_hed)
    
    shape_sim_nhed <- (mean_sim_nhed / sd_sim_nhed)^2
    rate_sim_nhed <- mean_sim_nhed / (sd_sim_nhed^2)
    shape_sim_hed <- (mean_sim_hed / sd_sim_hed)^2
    rate_sim_hed <- mean_sim_hed / (sd_sim_hed^2)
    
    y_gamma_sim_nhed <- dgamma(x_vals_nhed, shape = shape_sim_nhed, rate = rate_sim_nhed)
    y_gamma_sim_hed_60 <- dgamma(x_vals_nhed, shape = shape_sim_hed, rate = rate_sim_hed)
    y_gamma_sim_hed_150 <- dgamma(x_vals_hed, shape = shape_sim_hed, rate = rate_sim_hed)
    
    # Simulate the beta coefficient from a normal distribution
    beta_sim <- rnorm(1, mean = beta, sd = sqrt(cov_matrix))
    
    # Compute the relative risk based on the simulated beta
    rr_sim_nhed <- rr_function_nhed(x = x_vals_nhed, b = beta_sim) # Adjusted to match definition
    rr_sim_hed_60 <- rr_function_hed(x = x_vals_nhed, beta = beta_sim)
    rr_sim_hed_150 <- rr_function_hed(x = x_vals_hed, beta = beta_sim)
    
    # Simulate proportions of lifetime abstainers and former drinkers
    prop_abs_sim <- max(rnorm(1, mean = p_abs, sd = sqrt(p_abs * (1 - p_abs) / 1000)), 0.001)
    prop_form_sim <- max(rnorm(1, mean = p_form, sd = sqrt(p_form * (1 - p_form) / 1000)), 0.001)
    prop_hed_sim <- max(rnorm(1, mean = p_hed, sd = sqrt(p_hed * (1 - p_hed) / 1000)), 0.001)
    
    # Calculate PAF using the customized function
    simulated_pafs[i] <- paf_hed_function(
      x_60 = x_vals_nhed, 
      x_150 = x_vals_hed, 
      y_nhed = y_gamma_sim_nhed, 
      y_hed_60 = y_gamma_sim_hed_60, 
      y_hed_150 = y_gamma_sim_hed_150, 
      rr_nhed = rr_sim_nhed, 
      rr_hed_60 = rr_sim_hed_60, 
      rr_hed_150 = rr_sim_hed_150,
      rr_form = rr_fd, 
      p_abs = prop_abs_sim, 
      p_form = prop_form_sim, 
      p_hed = prop_hed_sim
    )
  }
  
  simulated_pafs <- simulated_pafs[!is.nan(simulated_pafs)]
  if (length(simulated_pafs) == 0) {
    stop("All simulations resulted in NaN values. Please check your input parameters.")
  }
  
  paf_lower <- quantile(simulated_pafs, 0.025)
  paf_upper <- quantile(simulated_pafs, 0.975)
  paf_point_estimate <- mean(simulated_pafs)
  
  return(list(point_estimate = round(paf_point_estimate, 3), 
              lower_ci = round(paf_lower, 3), 
              upper_ci = round(paf_upper, 3)))
}

#---------------------------------#
#------------ 2008 DATA ----------#
#---------------------------------#
cd_fem_list <- list()

# Define the years to process
years <- c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022)

for (year in years) {
  cd_fem_list[[as.character(year)]] <- list()
  
  for (i in 1:4) {  # Iterate over edad_tramo
    cd_fem_list[[as.character(year)]][[i]] <- data_input %>%
      filter(volajohdia > 0, 
             sexo == "Mujer", 
             edad_tramo == i, 
             year == !!year) %>%  # Use the correct value of 'year'
      pull(volajohdia)
  }
}


# numbers represent age group
# 08,10,12,14,16,18,20,22 are years

# Create a list to store gamma fits
g_fem_list <- list()

for (year in names(cd_fem_list)) {
  g_fem_list[[year]] <- lapply(cd_fem_list[[year]], function(data) {
    if (length(data) > 1) {  # Ensure sufficient data to fit
      fitdist(data, "gamma")
    } else {
      NULL  # Handle cases with insufficient data
    }
  })
}

cd_fem_hed_list <- list()

hed_values <- c(0, 1)

# Nested loop to iterate through years, age group, and hed values
for (year in years) {
  cd_fem_hed_list[[as.character(year)]] <- list()
  for (i in 1:4) {
    for (hed_value in hed_values) {
      # Create a descriptive key
      key <- paste0("edad", i, ifelse(hed_value == 0, "_nhed", "_hed"))
      
      # Filter and pull data
      cd_fem_hed_list[[as.character(year)]][[key]] <- data_input %>%
        filter(volajohdia > 0, 
               sexo == "Mujer", 
               edad_tramo == i, 
               hed == !!hed_value, 
               year == !!year) %>%  # Ensure year is correctly interpreted
        pull(volajohdia)
    }
  }
}

# Ensure that gamma fits are organized by "nhed" and "hed"
g_fem_hed_list <- list()

for (year in names(cd_fem_hed_list)) {
  g_fem_hed_list[[year]] <- list()
  
  for (edad_tramo in 1:4) {
    # Extract data for non-heavy drinkers
    nhed_key <- paste0("edad", edad_tramo, "_nhed")
    hed_key <- paste0("edad", edad_tramo, "_hed")
    
    data_nhed <- cd_fem_hed_list[[year]][[nhed_key]]
    data_hed <- cd_fem_hed_list[[year]][[hed_key]]
    
    # Fit gamma distribution if sufficient data exists
    g_fem_hed_list[[year]][[edad_tramo]] <- list(
      nhed = if (length(data_nhed) > 1) fitdist(data_nhed, "gamma") else NULL,
      hed = if (length(data_hed) > 1) fitdist(data_hed, "gamma") else NULL
    )
  }
}

# Filter input_female for rows where cvolaj is "ltabs"
ltabs_fem <- input_female %>%
  filter(cvolaj == "ltabs") %>%
  arrange(year, edad_tramo)

# Create the list for ltabs values by year and age group
p_abs_list_fem <- ltabs_fem %>%
  group_by(year, edad_tramo) %>%
  summarise(prop = list(prop), .groups = "drop") %>%
  split(.$year) %>%
  lapply(function(df) {
    setNames(df$prop, paste0("edad_tramo_", df$edad_tramo))
  })

# Filter input_female for rows where cvolaj is "fd"
fd_fem <- input_female %>%
  filter(cvolaj == "fd") %>%
  arrange(year, edad_tramo)

# Create the list for fd values by year and age group
p_form_list_fem <- fd_fem %>%
  group_by(year, edad_tramo) %>%
  summarise(prop = list(prop), .groups = "drop") %>%
  split(.$year) %>%
  lapply(function(df) {
    setNames(df$prop, paste0("edad_tramo_", df$edad_tramo))
  })

# Filter data for males
data_hed_fem <- data_hed %>%
  filter(sexo == "Mujer", hed == 1) %>%
  arrange(edad_tramo, year)

# Create the list of lists for each edad_tramo
p_hed_list_fem <- data_hed_fem %>%
  group_by(edad_tramo) %>%
  summarise(prop_hed_values = list(prop_hed), .groups = "drop") %>%
  pull(prop_hed_values)


#------------------------------#
#------------MALES-------------#
#------------------------------#

# Create an empty list to store results by year for males
cd_male_list <- list()

# Loop over each year and age group for males
for (year in years) {
  cd_male_list[[as.character(year)]] <- list()
  for (i in 1:4) {
    cd_male_list[[as.character(year)]][[i]] <- data_input %>%
      filter(volajohdia > 0, sexo == "Hombre", edad_tramo == i, 
             year == !!year) %>%
      pull(volajohdia)
  }
}

# Create a list to store gamma fits for males
g_male_list <- list()

# Loop through each year and fit gamma distributions for each edad_tramo
for (year in names(cd_male_list)) {
  g_male_list[[year]] <- lapply(cd_male_list[[year]], function(data) {
    if (length(data) > 1) { # Ensure there is enough data to fit a distribution
      fitdist(data, "gamma")
    } else {
      NULL # Handle cases with insufficient data
    }
  })
}

# PCA of hed and nhed males
cd_male_hed_list <- list()

# Nested loop to iterate through years, age group, and hed values for males
for (year in years) {
  cd_male_hed_list[[as.character(year)]] <- list()
  for (i in 1:4) {
    for (hed_value in hed_values) {
      # Create a descriptive key
      key <- paste0("edad", i, ifelse(hed_value == 0, "_nhed", "_hed"))
      
      # Filter and pull data
      cd_male_hed_list[[as.character(year)]][[key]] <- data_input %>%
        filter(volajohdia > 0, sexo == "Hombre", edad_tramo == i, hed == !!hed_value, year == !!year) %>%
        pull(volajohdia)
    }
  }
}

cd_male_hed_list <- list()

# Nested loop to iterate through years, age group, and hed values for males
for (year in years) {
  cd_male_hed_list[[as.character(year)]] <- list()
  for (i in 1:4) {
    for (hed_value in hed_values) {
      # Create a descriptive key
      key <- paste0("edad", i, ifelse(hed_value == 0, "_nhed", "_hed"))
      
      # Filter and pull data
      cd_male_hed_list[[as.character(year)]][[key]] <- data_input %>%
        filter(volajohdia > 0, sexo == "Hombre", edad_tramo == i, hed == !!hed_value, year == !!year) %>%
        pull(volajohdia)
    }
  }
}

# Create a list to store gamma fits for hed and nhed for males
g_male_hed_list <- list()

# Loop through each year in cd_male_hed_list
for (year in names(cd_male_hed_list)) {
  g_male_hed_list[[year]] <- list()
  
  # Loop through each edad_tramo and hed value
  for (key in names(cd_male_hed_list[[year]])) {
    data <- cd_male_hed_list[[year]][[key]]
    
    # Fit gamma distribution only if there is enough data
    g_male_hed_list[[year]][[key]] <- if (length(data) > 1) {
      fitdist(data, "gamma")
    } else {
      NULL # Handle cases with insufficient data
    }
  }
}

g_male_hed_list <- list()

for (year in names(cd_male_hed_list)) {
  g_male_hed_list[[year]] <- list()
  
  for (edad_tramo in 1:4) {
    # Extract data for non-heavy drinkers
    nhed_key <- paste0("edad", edad_tramo, "_nhed")
    hed_key <- paste0("edad", edad_tramo, "_hed")
    
    data_nhed <- cd_male_hed_list[[year]][[nhed_key]]
    data_hed <- cd_male_hed_list[[year]][[hed_key]]
    
    # Fit gamma distribution if sufficient data exists
    g_male_hed_list[[year]][[edad_tramo]] <- list(
      nhed = if (length(data_nhed) > 1) fitdist(data_nhed, "gamma") else NULL,
      hed = if (length(data_hed) > 1) fitdist(data_hed, "gamma") else NULL
    )
  }
}
# Filter input_female for rows where cvolaj is "ltabs"
ltabs_male <- input_male %>%
  filter(cvolaj == "ltabs") %>%
  arrange(year, edad_tramo)

# Create the list for ltabs values by year and age group
p_abs_list_male <- ltabs_male %>%
  group_by(year, edad_tramo) %>%
  summarise(prop = list(prop), .groups = "drop") %>%
  split(.$year) %>%
  lapply(function(df) {
    setNames(df$prop, paste0("edad_tramo_", df$edad_tramo))
  })

# Filter input_female for rows where cvolaj is "fd"
fd_male <- input_male %>%
  filter(cvolaj == "fd") %>%
  arrange(year, edad_tramo)

# Create the list for fd values by year and age group
p_form_list_male <- fd_male %>%
  group_by(year, edad_tramo) %>%
  summarise(prop = list(prop), .groups = "drop") %>%
  split(.$year) %>%
  lapply(function(df) {
    setNames(df$prop, paste0("edad_tramo_", df$edad_tramo))
  })


# Filter data for males
data_hed_male <- data_hed %>%
  filter(sexo == "Hombre", hed == 1) %>%
  arrange(edad_tramo, year)

# Create the list of lists for each edad_tramo
p_hed_list_male <- data_hed_male %>%
  group_by(edad_tramo) %>%
  summarise(prop_hed_values = list(prop_hed), .groups = "drop") %>%
  pull(prop_hed_values)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#


####################################
# ESTIMATING AAF FOR BREAST CANCER #
####################################

b1_bcan <- 0.01018
var_bcan <- 0.0000007208

# Create a data frame to store results
bcan_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)

for (i in 1:length(bcan_female$Year)) {
  for (j in 1:4) {  # Iterate over edad_tramo
    # Access the relevant elements from the lists
    gamma_fit <- g_fem_list[[as.character(bcan_female$Year[i])]][[j]]
    p_abs <- p_abs_list_fem[[as.character(bcan_female$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_fem[[as.character(bcan_female$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      result_bcan <- tryCatch({
        confint_paf(
          gamma = gamma_fit,
          beta = b1_bcan,
          var_beta = var_bcan,
          p_abs = p_abs,
          p_form = p_form,
          rr_form = 1,
          rr_function = rr_linear
        )
      }, error = function(e) NULL)
      
      # Update results in bcan_female if the calculation succeeds
      if (!is.null(result_bcan) && !is.null(result_bcan$Point_Estimate)) {
        bcan_female[i, paste0("Fem", j, "_point")] <- result_bcan$Point_Estimate
        bcan_female[i, paste0("Fem", j, "_lower")] <- result_bcan$Lower_CI
        bcan_female[i, paste0("Fem", j, "_upper")] <- result_bcan$Upper_CI
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", bcan_female$Year[i], "Group Fem", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", bcan_female$Year[i], "Group Fem", j, "- Invalid data\n")
    }
  }
}

bcan_female <- bcan_female %>% mutate(disease = "Breast Cancer")


# LIP AND ORAL CAVITY CANCER (BOTH SEXES)
b1_locan <- 0.02474
b2_locan <- -0.00004
rr_locan_fd <- 1.2
cov_matrix_locan <- matrix(c(0.000002953, -0.0000000127,
                             -0.0000000127, 0.000000000102), 
                           nrow = 2, ncol = 2, byrow = TRUE)

rr_locan_fun <- function(x, betas){
  # Example: assuming a quadratic model
  rr <- exp(betas[1] * x + betas[2] * x^2)
  return(rr)
}

# Define a data frame to store results

# Define a data frame to store results
locan_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)

# Loop through each year and each female group
for (i in 1:length(locan_female$Year)) {
  for (j in 1:4) {
    # Access the relevant elements from the lists
    gamma_fit <- g_fem_list[[as.character(locan_female$Year[i])]][[j]]
    p_abs <- p_abs_list_fem[[as.character(locan_female$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_fem[[as.character(locan_female$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      # Call confint_paf_vcov with the required inputs
      result <- tryCatch({
        confint_paf_vcov(
          gamma = gamma_fit,
          betas = c(b1_locan, b2_locan),
          cov_matrix = cov_matrix_locan,
          p_abs = p_abs,
          p_form = p_form,
          rr_fd = rr_locan_fd,
          rr_function = rr_locan_fun
        )
      }, error = function(e) NULL)
      
      # Update results in locan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$point_estimate)) {
        locan_female[i, paste0("Fem", j, "_point")] <- result$point_estimate
        locan_female[i, paste0("Fem", j, "_lower")] <- result$lower_ci
        locan_female[i, paste0("Fem", j, "_upper")] <- result$upper_ci
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", locan_female$Year[i], "Group Fem", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", locan_female$Year[i], "Group Fem", j, "- Invalid data\n")
    }
  }
}
locan_female <- locan_female %>% 
  mutate(disease = "Lip and Oral Cavity Cancer")
# Display the results matrix
locan_female

# Define un data frame para almacenar los resultados para hombres
locan_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

# Loop through each year and each female group
for (i in 1:length(locan_male$Year)) {
  for (j in 1:4) {
    # Access the relevant elements from the lists
    gamma_fit <- g_male_list[[as.character(locan_male$Year[i])]][[j]]
    p_abs <- p_abs_list_male[[as.character(locan_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_male[[as.character(locan_male$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      # Call confint_paf_vcov with the required inputs
      result <- tryCatch({
        confint_paf_vcov(
          gamma = gamma_fit,
          betas = c(b1_locan, b2_locan),
          cov_matrix = cov_matrix_locan,
          p_abs = p_abs,
          p_form = p_form,
          rr_fd = rr_locan_fd,
          rr_function = rr_locan_fun
        )
      }, error = function(e) NULL)
      
      # Update results in locan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$point_estimate)) {
        locan_male[i, paste0("Male", j, "_point")] <- result$point_estimate
        locan_male[i, paste0("Male", j, "_lower")] <- result$lower_ci
        locan_male[i, paste0("Male", j, "_upper")] <- result$upper_ci
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", locan_male$Year[i], "Group Male", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", locan_male$Year[i], "Group Male", j, "- Invalid data\n")
    }
  }
}
locan_male <- locan_male %>% 
  mutate(disease = "Lip and Oral Cavity Cancer")
# Display the results matrix
locan_male

################################################
# ESTIMATING AAF FOR OTHER PHARINGEAL CANCER   #
################################################

# OTHER PHARINGEAL CANCER (BOTH SEXES)
b1_opcan <- 0.02474
b2_opcan <- -0.00004
rr_opcan_fd <- 1.2
cov_matrix_opcan <- matrix(c(0.000002953, -0.0000000127,
                             -0.0000000127, 0.000000000102), 
                           nrow = 2, ncol = 2, byrow = TRUE)


# Define a data frame to store results
opcan_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)


for (i in 1:length(opcan_female$Year)) {
  for (j in 1:4) {
    # Access the relevant elements from the lists
    gamma_fit <- g_fem_list[[as.character(opcan_female$Year[i])]][[j]]
    p_abs <- p_abs_list_fem[[as.character(opcan_female$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_fem[[as.character(opcan_female$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      # Call confint_paf_vcov with the required inputs
      result <- tryCatch({
        confint_paf_vcov(
          gamma = gamma_fit,
          betas = c(b1_locan, b2_locan),
          cov_matrix = cov_matrix_locan,
          p_abs = p_abs,
          p_form = p_form,
          rr_fd = rr_locan_fd,
          rr_function = rr_locan_fun
        )
      }, error = function(e) NULL)
      
      # Update results in locan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$point_estimate)) {
        opcan_female[i, paste0("Fem", j, "_point")] <- result$point_estimate
        opcan_female[i, paste0("Fem", j, "_lower")] <- result$lower_ci
        opcan_female[i, paste0("Fem", j, "_upper")] <- result$upper_ci
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", opcan_female$Year[i], "Group Fem", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", opcan_female$Year[i], "Group Fem", j, "- Invalid data\n")
    }
  }
}
opcan_female <- opcan_female %>% 
  mutate(disease = "Other Pharingeal Cancer")

# Display the results matrix
opcan_female

# Define un data frame para almacenar los resultados para hombres
opcan_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)


# Loop through each year and each male group
for (i in 1:length(opcan_male$Year)) {
  for (j in 1:4) {
    # Access the relevant elements from the lists
    gamma_fit <- g_male_list[[as.character(opcan_male$Year[i])]][[j]]
    p_abs <- p_abs_list_male[[as.character(opcan_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_male[[as.character(opcan_male$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      # Call confint_paf_vcov with the required inputs
      result <- tryCatch({
        confint_paf_vcov(
          gamma = gamma_fit,
          betas = c(b1_locan, b2_locan),
          cov_matrix = cov_matrix_locan,
          p_abs = p_abs,
          p_form = p_form,
          rr_fd = rr_locan_fd,
          rr_function = rr_locan_fun
        )
      }, error = function(e) NULL)
      
      # Update results in locan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$point_estimate)) {
        opcan_male[i, paste0("Male", j, "_point")] <- result$point_estimate
        opcan_male[i, paste0("Male", j, "_lower")] <- result$lower_ci
        opcan_male[i, paste0("Male", j, "_upper")] <- result$upper_ci
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", opcan_male$Year[i], "Group Male", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", opcan_male$Year[i], "Group Male", j, "- Invalid data\n")
    }
  }
}
opcan_male <- opcan_male %>% 
  mutate(disease = "Other Pharingeal Cancer")

##########################################
# ESTIMATING AAF FOR OESOPHAFUS CANCER   #
##########################################
rr_oescan_cal <- function(betas,x){
  exp(betas[1]*x+betas[2]*x**3)
}

b1_oescan <- 0.0132063596418668
b2_oescan <- -4.14801974664481*10**-8
rr_oescan_fd <- 1.16
cov_matrix_oescan <- matrix(c(1.525706*10**-7, -6.885205*10**-13,
                              -0.0000000127, 0.000002953), 
                            nrow = 2, ncol = 2, byrow = TRUE)


# Define a data frame to store results
oescan_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)


for (i in 1:length(oescan_female$Year)) {
  for (j in 1:4) {
    # Access the relevant elements from the lists
    gamma_fit <- g_fem_list[[as.character(oescan_female$Year[i])]][[j]]
    p_abs <- p_abs_list_fem[[as.character(oescan_female$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_fem[[as.character(oescan_female$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      # Call confint_paf_vcov with the required inputs
      result <- tryCatch({
        confint_paf_vcov(
          gamma = gamma_fit,
          betas = c(b1_oescan, b2_oescan),
          cov_matrix = cov_matrix_oescan,
          p_abs = p_abs,
          p_form = p_form,
          rr_fd = rr_oescan_fd,
          rr_function = rr_oescan_cal
        )
      }, error = function(e) NULL)
      
      # Update results in locan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$point_estimate)) {
        oescan_female[i, paste0("Fem", j, "_point")] <- result$point_estimate
        oescan_female[i, paste0("Fem", j, "_lower")] <- result$lower_ci
        oescan_female[i, paste0("Fem", j, "_upper")] <- result$upper_ci
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", oescan_female$Year[i], "Group Fem", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", oescan_female$Year[i], "Group Fem", j, "- Invalid data\n")
    }
  }
}


oescan_female <- oescan_female %>% 
  mutate(disease = "Oesophagus Cancer")
# Display the results matrix
oescan_female

# Define un data frame para almacenar los resultados para hombres
oescan_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)


for (i in 1:length(oescan_male$Year)) {
  for (j in 1:4) {
    # Access the relevant elements from the lists
    gamma_fit <- g_male_list[[as.character(oescan_male$Year[i])]][[j]]
    p_abs <- p_abs_list_male[[as.character(oescan_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_male[[as.character(oescan_male$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      # Call confint_paf_vcov with the required inputs
      result <- tryCatch({
        confint_paf_vcov(
          gamma = gamma_fit,
          betas = c(b1_oescan, b2_oescan),
          cov_matrix = cov_matrix_oescan,
          p_abs = p_abs,
          p_form = p_form,
          rr_fd = rr_oescan_fd,
          rr_function = rr_oescan_cal
        )
      }, error = function(e) NULL)
      
      # Update results in locan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$point_estimate)) {
        oescan_male[i, paste0("Male", j, "_point")] <- result$point_estimate
        oescan_male[i, paste0("Male", j, "_lower")] <- result$lower_ci
        oescan_male[i, paste0("Male", j, "_upper")] <- result$upper_ci
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", oescan_male$Year[i], "Group Male", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", oescan_male$Year[i], "Group Male", j, "- Invalid data\n")
    }
  }
}

oescan_male <- oescan_male %>% 
  mutate(disease = "Oesophagus Cancer")

################################################
# ESTIMATING AAF FOR COLON AND RECTUM CANCER   #
################################################

b1_crcan_male <- 0.006806
b1_crcan_fem <- 0.003020
rr_crcan_fd_fem <- 1.05
var_crcan <- 0.000000907
rr_crcan_fd_male <- 2.19

# data frame for females
crcan_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)

# estimating PAF and confidence intervals
for (i in 1:length(crcan_female$Year)) {
  for (j in 1:4) {  # Iterate over edad_tramo
    # Access the relevant elements from the lists
    gamma_fit <- g_fem_list[[as.character(crcan_female$Year[i])]][[j]]
    p_abs <- p_abs_list_fem[[as.character(crcan_female$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_fem[[as.character(crcan_female$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      result <- tryCatch({
        confint_paf(
          gamma = gamma_fit,
          beta = b1_crcan_fem,
          var_beta = var_crcan,
          p_abs = p_abs,
          p_form = p_form,
          rr_form = rr_crcan_fd_fem,
          rr_function = rr_linear
        )
      }, error = function(e) NULL)
      
      # Update results in bcan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$Point_Estimate)) {
        crcan_female[i, paste0("Fem", j, "_point")] <- result$Point_Estimate
        crcan_female[i, paste0("Fem", j, "_lower")] <- result$Lower_CI
        crcan_female[i, paste0("Fem", j, "_upper")] <- result$Upper_CI
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", crcan_female$Year[i], "Group Fem", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", crcan_female$Year[i], "Group Fem", j, "- Invalid data\n")
    }
  }
}

crcan_female <- crcan_female %>% 
  mutate(disease = "Colon and rectum Cancer")

# data frame for males
crcan_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

# PAF estimation for males
for (i in 1:length(crcan_male$Year)) {
  for (j in 1:4) {  # Iterate over edad_tramo
    # Access the relevant elements from the lists
    gamma_fit <- g_male_list[[as.character(crcan_male$Year[i])]][[j]]
    p_abs <- p_abs_list_male[[as.character(crcan_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_male[[as.character(crcan_male$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      result <- tryCatch({
        confint_paf(
          gamma = gamma_fit,
          beta = b1_crcan_male,
          var_beta = var_crcan,
          p_abs = p_abs,
          p_form = p_form,
          rr_form = rr_crcan_fd_male,
          rr_function = rr_linear
        )
      }, error = function(e) NULL)
      
      # Update results in bcan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$Point_Estimate)) {
        crcan_male[i, paste0("Male", j, "_point")] <- result$Point_Estimate
        crcan_male[i, paste0("Male", j, "_lower")] <- result$Lower_CI
        crcan_male[i, paste0("Male", j, "_upper")] <- result$Upper_CI
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", crcan_male$Year[i], "Group Male", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", crcan_male$Year[i], "Group Male", j, "- Invalid data\n")
    }
  }
}
crcan_male <- crcan_male %>% 
  mutate(disease = "Colon and rectum Cancer")

#####################################
# ESTIMATING AAF FOR LIVER CANCER   #
#####################################
b1_lican <- 0.00742949
b2_lican <- 0.0000148593
cov_matrix_lican <- matrix(c(0.000003097,0,
                            0, 0.000003097), nrow = 2, byrow = T)
rr_lican_fd_male <- 1.54
rr_lican_fd_fem <- 2.28

rr_lican_fun <- function(betas, x){
  exp(betas[1]*x-betas[2]*x**2)
}

# Define a data frame to store results
lican_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)

for (i in 1:length(lican_female$Year)) {
  for (j in 1:4) {
    # Access the relevant elements from the lists
    gamma_fit <- g_fem_list[[as.character(lican_female$Year[i])]][[j]]
    p_abs <- p_abs_list_fem[[as.character(lican_female$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_fem[[as.character(lican_female$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      # Call confint_paf_vcov with the required inputs
      result <- tryCatch({
        confint_paf_vcov(
          gamma = gamma_fit,
          betas = c(b1_lican, b2_lican),
          cov_matrix = cov_matrix_lican,
          p_abs = p_abs,
          p_form = p_form,
          rr_fd = rr_lican_fd_fem,
          rr_function = rr_lican_fun
        )
      }, error = function(e) NULL)
      
      # Update results in locan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$point_estimate)) {
        lican_female[i, paste0("Fem", j, "_point")] <- result$point_estimate
        lican_female[i, paste0("Fem", j, "_lower")] <- result$lower_ci
        lican_female[i, paste0("Fem", j, "_upper")] <- result$upper_ci
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", lican_female$Year[i], "Group Fem", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", lican_female$Year[i], "Group Fem", j, "- Invalid data\n")
    }
  }
}
lican_female <- lican_female %>% 
  mutate(disease = "Liver Cancer")

# Define a data frame to store results
lican_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

for (i in 1:length(lican_male$Year)) {
  for (j in 1:4) {
    # Access the relevant elements from the lists
    gamma_fit <- g_male_list[[as.character(lican_male$Year[i])]][[j]]
    p_abs <- p_abs_list_male[[as.character(lican_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_male[[as.character(lican_male$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      # Call confint_paf_vcov with the required inputs
      result <- tryCatch({
        confint_paf_vcov(
          gamma = gamma_fit,
          betas = c(b1_lican, b2_lican),
          cov_matrix = cov_matrix_lican,
          p_abs = p_abs,
          p_form = p_form,
          rr_fd = rr_lican_fd_male,
          rr_function = rr_lican_fun
        )
      }, error = function(e) NULL)
      
      # Update results in locan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$point_estimate)) {
        lican_male[i, paste0("Male", j, "_point")] <- result$point_estimate
        lican_male[i, paste0("Male", j, "_lower")] <- result$lower_ci
        lican_male[i, paste0("Male", j, "_upper")] <- result$upper_ci
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", lican_male$Year[i], "Group Male", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", lican_male$Year[i], "Group Male", j, "- Invalid data\n")
    }
  }
}
lican_male <- lican_male %>% 
  mutate(disease = "Liver Cancer")
######################################
# ESTIMATING AAF FOR LARYNX CANCER   #
#####################################
b1_lxcan <- 0.01462
b2_lxcan <- -0.00002
cov_matrix_lxcan <- matrix(c(3.585e-06, -1.62e-08,
                             -1.62e-08, 1.26e-10), 
                           nrow = 2, ncol = 2, byrow = TRUE)

rr_fd_lxcan <- 1.18
rr_lxcan_fun <- function(x, betas){
  exp(x * betas[1] - x^2 * betas[2])
}

# Define a data frame to store results
lxcan_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)

# PAF estimation for females
for (i in 1:length(lxcan_female$Year)) {
  for (j in 1:4) {
    # Access the relevant elements from the lists
    gamma_fit <- g_fem_list[[as.character(lxcan_female$Year[i])]][[j]]
    p_abs <- p_abs_list_fem[[as.character(lxcan_female$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_fem[[as.character(lxcan_female$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      # Call confint_paf_vcov with the required inputs
      result <- tryCatch({
        confint_paf_vcov(
          gamma = gamma_fit,
          betas = c(b1_lxcan, b2_lxcan),
          cov_matrix = cov_matrix_lxcan,
          p_abs = p_abs,
          p_form = p_form,
          rr_fd = rr_fd_lxcan,
          rr_function = rr_lxcan_fun
        )
      }, error = function(e) NULL)
      
      # Update results in locan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$point_estimate)) {
        lxcan_female[i, paste0("Fem", j, "_point")] <- result$point_estimate
        lxcan_female[i, paste0("Fem", j, "_lower")] <- result$lower_ci
        lxcan_female[i, paste0("Fem", j, "_upper")] <- result$upper_ci
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", lxcan_female$Year[i], "Group Fem", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", lxcan_female$Year[i], "Group Fem", j, "- Invalid data\n")
    }
  }
}
lxcan_female <- lxcan_female %>% 
  mutate(disease = "Larynx Cancer")

# Define a data frame to store results
lxcan_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)
# PAF estimation for males
for (i in 1:length(lxcan_male$Year)) {
  for (j in 1:4) {
    # Access the relevant elements from the lists
    gamma_fit <- g_male_list[[as.character(lxcan_male$Year[i])]][[j]]
    p_abs <- p_abs_list_male[[as.character(lxcan_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_male[[as.character(lxcan_male$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      # Call confint_paf_vcov with the required inputs
      result <- tryCatch({
        confint_paf_vcov(
          gamma = gamma_fit,
          betas = c(b1_lxcan, b2_lxcan),
          cov_matrix = cov_matrix_lxcan,
          p_abs = p_abs,
          p_form = p_form,
          rr_fd = rr_fd_lxcan,
          rr_function = rr_lxcan_fun
        )
      }, error = function(e) NULL)
      
      # Update results in locan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$point_estimate)) {
        lxcan_male[i, paste0("Male", j, "_point")] <- result$point_estimate
        lxcan_male[i, paste0("Male", j, "_lower")] <- result$lower_ci
        lxcan_male[i, paste0("Male", j, "_upper")] <- result$upper_ci
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", lxcan_male$Year[i], "Group Male", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", lxcan_male$Year[i], "Group Male", j, "- Invalid data\n")
    }
  }
}

lxcan_male <- lxcan_male %>% 
  mutate(disease = "Larynx Cancer")

# ag means age group
aaf_cancer_fem <- bind_rows(lican_female, locan_female,bcan_female, lxcan_female,
                            oescan_female, opcan_female, crcan_female) %>% 
  rename(AAF_ag1 = Fem1_point, LL_ag1 = Fem1_lower, UL_ag1 = Fem1_upper,
         AAF_ag2 = Fem2_point, LL_ag2 = Fem2_lower, UL_ag2 = Fem2_upper,
         AAF_ag3 = Fem3_point, LL_ag3 = Fem3_lower, UL_ag3 = Fem3_upper,
         AAF_ag4 = Fem4_point, LL_ag4 = Fem4_lower, UL_ag4 = Fem4_upper)

aaf_cancer_male <- bind_rows(lican_male, locan_male, lxcan_male,
                            oescan_male, opcan_male, crcan_male) %>% 
  rename(AAF_ag1 = Male1_point, LL_ag1 = Male1_lower, UL_ag1 = Male1_upper,
         AAF_ag2 = Male2_point, LL_ag2 = Male2_lower, UL_ag2 = Male2_upper,
         AAF_ag3 = Male3_point, LL_ag3 = Male3_lower, UL_ag3 = Male3_upper,
         AAF_ag4 = Male4_point, LL_ag4 = Male4_lower, UL_ag4 = Male4_upper)

######################
# NEURO-PSYCHIATRIC #
#####################
# EPILEPSY (BOTH SEXES)
b1_epi <- 1.22861
var_epi <- 0.1391974
rr_epi_fun <- function(x, beta){
  exp(beta * (x+0.5)/100)
}
# Define a data frame to store results
epi_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018,2022),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)

# estimating PAF and confidence intervals
for (i in 1:length(epi_female$Year)) {
  for (j in 1:4) {  # Iterate over edad_tramo
    # Access the relevant elements from the lists
    gamma_fit <- g_fem_list[[as.character(epi_female$Year[i])]][[j]]
    p_abs <- p_abs_list_fem[[as.character(epi_female$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_fem[[as.character(epi_female$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      result <- tryCatch({
        confint_paf(
          gamma = gamma_fit,
          beta = b1_epi,
          var_beta = var_epi,
          p_abs = p_abs,
          p_form = p_form,
          rr_form = 1,
          rr_function = rr_epi_fun
        )
      }, error = function(e) NULL)
      
      # Update results in bcan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$Point_Estimate)) {
        epi_female[i, paste0("Fem", j, "_point")] <- result$Point_Estimate
        epi_female[i, paste0("Fem", j, "_lower")] <- result$Lower_CI
        epi_female[i, paste0("Fem", j, "_upper")] <- result$Upper_CI
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", epi_female$Year[i], "Group Fem", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", epi_female$Year[i], "Group Fem", j, "- Invalid data\n")
    }
  }
}
epi_female <- epi_female %>% 
  mutate(disease = "Epilepsy") %>% 
  rename(AAF_ag1 = Fem1_point, LL_ag1 = Fem1_lower, UL_ag1 = Fem1_upper,
         AAF_ag2 = Fem2_point, LL_ag2 = Fem2_lower, UL_ag2 = Fem2_upper,
         AAF_ag3 = Fem3_point, LL_ag3 = Fem3_lower, UL_ag3 = Fem3_upper,
         AAF_ag4 = Fem4_point, LL_ag4 = Fem4_lower, UL_ag4 = Fem4_upper)

# Define a data frame to store results
epi_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

# estimating PAF and confidence intervals
for (i in 1:length(epi_male$Year)) {
  for (j in 1:4) {  # Iterate over age group
    # Access the relevant elements from the lists
    gamma_fit <- g_male_list[[as.character(epi_male$Year[i])]][[j]]
    p_abs <- p_abs_list_male[[as.character(epi_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_male[[as.character(epi_male$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      result <- tryCatch({
        confint_paf(
          gamma = gamma_fit,
          beta = b1_epi,
          var_beta = var_epi,
          p_abs = p_abs,
          p_form = p_form,
          rr_form = 1,
          rr_function = rr_epi_fun
        )
      }, error = function(e) NULL)
      
      # Update results in bcan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$Point_Estimate)) {
        epi_male[i, paste0("Male", j, "_point")] <- result$Point_Estimate
        epi_male[i, paste0("Male", j, "_lower")] <- result$Lower_CI
        epi_male[i, paste0("Male", j, "_upper")] <- result$Upper_CI
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", epi_male$Year[i], "Group Male", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", epi_male$Year[i], "Group Male", j, "- Invalid data\n")
    }
  }
}

epi_male <- epi_male %>% 
  mutate(disease = "Epilepsy") %>% 
  rename(AAF_ag1 = Male1_point, LL_ag1 = Male1_lower, UL_ag1 = Male1_upper,
                AAF_ag2 = Male2_point, LL_ag2 = Male2_lower, UL_ag2 = Male2_upper,
                AAF_ag3 = Male3_point, LL_ag3 = Male3_lower, UL_ag3 = Male3_upper,
                AAF_ag4 = Male4_point, LL_ag4 = Male4_lower, UL_ag4 = Male4_upper)



################
# OTHER CAUSES #
################

# Diabetes Mellitus (Male)
b1_dm_male <- 0.1763703
b2_dm_male <- -0.0728256
fd_dm <- 1.18
vcov_diabetes_male <- matrix(c(0.1681525, -0.2240129,
                               -0.2240129, 0.29964479),
                             nrow = 2, ncol = 2, byrow = TRUE)


rr_diabetes_male_fun <- function(x, betas) {
  b1 <- betas[1]
  b2 <- betas[2]
  
  # Compute the relative risk
  rr <- exp((x / 100)^2 * b1 + (x / 100)^3 * b2)
  
  return(rr)
}

dm_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)


# PAF estimation for males
for (i in 1:length(dm_male$Year)) {
  for (j in 1:4) {
    # Access the relevant elements from the lists
    gamma_fit <- g_male_list[[as.character(dm_male$Year[i])]][[j]]
    p_abs <- p_abs_list_male[[as.character(dm_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_male[[as.character(dm_male$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      # Call confint_paf_vcov with the required inputs
      result <- tryCatch({
        confint_paf_vcov(
          gamma = gamma_fit,
          betas = c(b1_dm_male, b2_dm_male),
          cov_matrix = vcov_diabetes_male,
          p_abs = p_abs,
          p_form = p_form,
          rr_fd = fd_dm,
          rr_function = rr_diabetes_male_fun
        )
      }, error = function(e) NULL)
      
      # Update results in locan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$point_estimate)) {
        dm_male[i, paste0("Male", j, "_point")] <- result$point_estimate
        dm_male[i, paste0("Male", j, "_lower")] <- result$lower_ci
        dm_male[i, paste0("Male", j, "_upper")] <- result$upper_ci
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", dm_male$Year[i], "Group Male", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", dm_male$Year[i], "Group Male", j, "- Invalid data\n")
    }
  }
}
dm_male <- dm_male %>% 
  mutate(disease = "DM2")

# Diabetes Mellitus (Female)
b1_dm_fem <- -1.3133910
b2_dm_fem <- 1.0142390
vcov_diabetes_female <- matrix(c(0.1681525, -0.2240129,
                                 -0.2240129, 0.29964479),
                               nrow = 2, ncol = 2, byrow = TRUE)

rr_diabetes_female_fun <- function(x, betas) {
  b1 <- betas[1]
  b2 <- betas[2]
  
  # Compute the relative risk
  rr <- exp(sqrt(x / 100) * b1 + (x / 100) * b2)
  
  return(rr)
}

dm_fem <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)


# PAF estimation for females
for (i in 1:length(dm_fem$Year)) {
  for (j in 1:4) {
    # Access the relevant elements from the lists
    gamma_fit <- g_fem_list[[as.character(dm_fem$Year[i])]][[j]]
    p_abs <- p_abs_list_fem[[as.character(dm_fem$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_fem[[as.character(dm_fem$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      # Call confint_paf_vcov with the required inputs
      result <- tryCatch({
        confint_paf_vcov(
          gamma = gamma_fit,
          betas = c(b1_dm_fem, b2_dm_fem),
          cov_matrix = vcov_diabetes_female,
          p_abs = p_abs,
          p_form = p_form,
          rr_fd = 1.14,
          rr_function = rr_diabetes_female_fun
        )
      }, error = function(e) NULL)
      
      # Update results in locan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$point_estimate)) {
        dm_fem[i, paste0("Fem", j, "_point")] <- result$point_estimate
        dm_fem[i, paste0("Fem", j, "_lower")] <- result$lower_ci
        dm_fem[i, paste0("Fem", j, "_upper")] <- result$upper_ci
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", dm_fem$Year[i], "Group Fem", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", dm_fem$Year[i], "Group Fem", j, "- Invalid data\n")
    }
  }
}

dm_fem <- dm_fem %>% 
  mutate(disease = "DM2")
# TUBERCULOSIS (BOTH SEXES)
b_tb <- 0.0179695
var_tb <- 0.007215**2


tb_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)

# estimating PAF and confidence intervals
for (i in 1:length(tb_female$Year)) {
  for (j in 1:4) {  # Iterate over edad_tramo
    # Access the relevant elements from the lists
    gamma_fit <- g_fem_list[[as.character(tb_female$Year[i])]][[j]]
    p_abs <- p_abs_list_fem[[as.character(tb_female$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_fem[[as.character(tb_female$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      result <- tryCatch({
        confint_paf(
          gamma = gamma_fit,
          beta = b_tb,
          var_beta = var_tb,
          p_abs = p_abs,
          p_form = p_form,
          rr_form = 1,
          rr_function = rr_linear
        )
      }, error = function(e) NULL)
      
      # Update results in bcan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$Point_Estimate)) {
        tb_female[i, paste0("Fem", j, "_point")] <- result$Point_Estimate
        tb_female[i, paste0("Fem", j, "_lower")] <- result$Lower_CI
        tb_female[i, paste0("Fem", j, "_upper")] <- result$Upper_CI
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", tb_female$Year[i], "Group Fem", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", tb_female$Year[i], "Group Fem", j, "- Invalid data\n")
    }
  }
}
tb_female <- tb_female %>% 
  mutate(disease = "Tuberculosis")

tb_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

# estimating PAF and confidence intervals
for (i in 1:length(tb_male$Year)) {
  for (j in 1:4) {  # Iterate over age group
    # Access the relevant elements from the lists
    gamma_fit <- g_male_list[[as.character(tb_male$Year[i])]][[j]]
    p_abs <- p_abs_list_male[[as.character(tb_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_male[[as.character(tb_male$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      result <- tryCatch({
        confint_paf(
          gamma = gamma_fit,
          beta = b_tb,
          var_beta = var_tb,
          p_abs = p_abs,
          p_form = p_form,
          rr_form = 1,
          rr_function = rr_linear
        )
      }, error = function(e) NULL)
      
      # Update results in bcan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$Point_Estimate)) {
        tb_male[i, paste0("Male", j, "_point")] <- result$Point_Estimate
        tb_male[i, paste0("Male", j, "_lower")] <- result$Lower_CI
        tb_male[i, paste0("Male", j, "_upper")] <- result$Upper_CI
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", tb_male$Year[i], "Group Male", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", tb_male$Year[i], "Group Male", j, "- Invalid data\n")
    }
  }
}

tb_male <- tb_male %>% 
  mutate(disease = "Tuberculosis")

# HIV/AIDS (Female)
b_hiv <- log(1.54)
var_hiv <- 0.078210772**2

# ONLY FOR PPL OVER 49

rr_hiv_fem <- function(x, beta) {
  b1 <- beta[1]
  
  # Initialize the RR as a vector of the same length as x
  rr <- numeric(length(x))
  
  # Case 1: x <= 61
  rr[x <= 49] <- 1
  
  # Case 2: x > 61
  rr[x > 49] <- exp(b1)
  
  return(rr)
}

hiv_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)

# estimating PAF and confidence intervals
for (i in 1:length(hiv_female$Year)) {
  for (j in 1:4) {  # Iterate over edad_tramo
    # Access the relevant elements from the lists
    gamma_fit <- g_fem_list[[as.character(hiv_female$Year[i])]][[j]]
    p_abs <- p_abs_list_fem[[as.character(hiv_female$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_fem[[as.character(hiv_female$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      result <- tryCatch({
        confint_paf(
          gamma = gamma_fit,
          beta = b_hiv,
          var_beta = var_hiv,
          p_abs = p_abs,
          p_form = p_form,
          rr_form = 1,
          rr_function = rr_hiv_fem
        )
      }, error = function(e) NULL)
      
      # Update results in bcan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$Point_Estimate)) {
        hiv_female[i, paste0("Fem", j, "_point")] <- result$Point_Estimate
        hiv_female[i, paste0("Fem", j, "_lower")] <- result$Lower_CI
        hiv_female[i, paste0("Fem", j, "_upper")] <- result$Upper_CI
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", hiv_female$Year[i], "Group Fem", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", hiv_female$Year[i], "Group Fem", j, "- Invalid data\n")
    }
  }
}
hiv_female <- hiv_female %>% 
  mutate(disease = "HIV")

# HIV/AIDS (Male)
# ONLY FOR PPL OVER 61 GRS
rr_hiv_male <- function(x, beta) {
  b1 <- beta[1]
  
  # Initialize the RR as a vector of the same length as x
  rr <- numeric(length(x))
  
  # Case 1: x <= 61
  rr[x <= 61] <- 1
  
  # Case 2: x > 61
  rr[x > 61] <- exp(b1)
  
  return(rr)
}

hiv_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

# estimating PAF and confidence intervals
for (i in 1:length(hiv_male$Year)) {
  for (j in 1:4) {  # Iterate over age group
    # Access the relevant elements from the lists
    gamma_fit <- g_male_list[[as.character(hiv_male$Year[i])]][[j]]
    p_abs <- p_abs_list_male[[as.character(hiv_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_male[[as.character(hiv_male$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      result <- tryCatch({
        confint_paf(
          gamma = gamma_fit,
          beta = b_hiv,
          var_beta = var_hiv,
          p_abs = p_abs,
          p_form = p_form,
          rr_form = 1,
          rr_function = rr_hiv_male
        )
      }, error = function(e) NULL)
      
      # Update results in bcan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$Point_Estimate)) {
        hiv_male[i, paste0("Male", j, "_point")] <- result$Point_Estimate
        hiv_male[i, paste0("Male", j, "_lower")] <- result$Lower_CI
        hiv_male[i, paste0("Male", j, "_upper")] <- result$Upper_CI
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", hiv_male$Year[i], "Group Male", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", hiv_male$Year[i], "Group Male", j, "- Invalid data\n")
    }
  }
}

hiv_male <- hiv_male %>% 
  mutate(disease = "HIV")

# LOWER RESPORATORY INFECTIONS (BOTH SEXES)
b_lri <- 0.004764038
var_lri <- 0.0019220552^2

lri_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)

# estimating PAF and confidence intervals
for (i in 1:length(lri_female$Year)) {
  for (j in 1:4) {  # Iterate over edad_tramo
    # Access the relevant elements from the lists
    gamma_fit <- g_fem_list[[as.character(lri_female$Year[i])]][[j]]
    p_abs <- p_abs_list_fem[[as.character(lri_female$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_fem[[as.character(lri_female$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      result <- tryCatch({
        confint_paf(
          gamma = gamma_fit,
          beta = b_lri,
          var_beta = var_lri,
          p_abs = p_abs,
          p_form = p_form,
          rr_form = 1,
          rr_function = rr_linear
        )
      }, error = function(e) NULL)
      
      # Update results in bcan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$Point_Estimate)) {
        lri_female[i, paste0("Fem", j, "_point")] <- result$Point_Estimate
        lri_female[i, paste0("Fem", j, "_lower")] <- result$Lower_CI
        lri_female[i, paste0("Fem", j, "_upper")] <- result$Upper_CI
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", lri_female$Year[i], "Group Fem", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", lri_female$Year[i], "Group Fem", j, "- Invalid data\n")
    }
  }
}

lri_female <- lri_female %>% 
  mutate(disease = "Lower Respiratory Infection")

lri_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

# estimating PAF and confidence intervals
for (i in 1:length(lri_male$Year)) {
  for (j in 1:4) {  # Iterate over age group
    # Access the relevant elements from the lists
    gamma_fit <- g_male_list[[as.character(lri_male$Year[i])]][[j]]
    p_abs <- p_abs_list_male[[as.character(lri_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_male[[as.character(lri_male$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      result <- tryCatch({
        confint_paf(
          gamma = gamma_fit,
          beta = b_lri,
          var_beta = var_lri,
          p_abs = p_abs,
          p_form = p_form,
          rr_form = 1,
          rr_function = rr_linear
        )
      }, error = function(e) NULL)
      
      # Update results in bcan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$Point_Estimate)) {
        lri_male[i, paste0("Male", j, "_point")] <- result$Point_Estimate
        lri_male[i, paste0("Male", j, "_lower")] <- result$Lower_CI
        lri_male[i, paste0("Male", j, "_upper")] <- result$Upper_CI
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", lri_male$Year[i], "Group Male", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", lri_male$Year[i], "Group Male", j, "- Invalid data\n")
    }
  }
}

lri_male <- lri_male %>% 
  mutate(disease = "Lower Respiratory Infection")

# LIVER CIRRHOSIS (MALES)
b_lc_male <- 0.02793524
var_lc_male <- 0.003832**2
rr_fd_lc <- 3.26

lc_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)


# estimating PAF and confidence intervals
for (i in 1:length(lc_male$Year)) {
  for (j in 1:4) {  # Iterate over age group
    # Access the relevant elements from the lists
    gamma_fit <- g_male_list[[as.character(lc_male$Year[i])]][[j]]
    p_abs <- p_abs_list_male[[as.character(lc_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_male[[as.character(lc_male$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      result <- tryCatch({
        confint_paf(
          gamma = gamma_fit,
          beta = b_lc_male,
          var_beta = var_lc_male,
          p_abs = p_abs,
          p_form = p_form,
          rr_form = rr_fd_lc,
          rr_function = rr_linear
        )
      }, error = function(e) NULL)
      
      # Update results in bcan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$Point_Estimate)) {
        lc_male[i, paste0("Male", j, "_point")] <- result$Point_Estimate
        lc_male[i, paste0("Male", j, "_lower")] <- result$Lower_CI
        lc_male[i, paste0("Male", j, "_upper")] <- result$Upper_CI
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", lc_male$Year[i], "Group Male", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", lc_male$Year[i], "Group Male", j, "- Invalid data\n")
    }
  }
}
lc_male <- lc_male %>% 
  mutate(disease = "Liver Cirrhosis")

# LIVER CIRRHOSIS (FEMALES)
b_lc_fem <- 0.32520349
var_lc_fem <- 0.033130**2
rr_lc_fem_fun <- function(beta, x){
  rr = exp(beta[1]*sqrt(x))
  return(rr)
}

lc_fem <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)


# estimating PAF and confidence intervals
for (i in 1:length(lc_fem$Year)) {
  for (j in 1:4) {  # Iterate over edad_tramo
    # Access the relevant elements from the lists
    gamma_fit <- g_fem_list[[as.character(lc_fem$Year[i])]][[j]]
    p_abs <- p_abs_list_fem[[as.character(lc_fem$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_fem[[as.character(lc_fem$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      result <- tryCatch({
        confint_paf(
          gamma = gamma_fit,
          beta = b_lc_fem,
          var_beta = var_lc_fem,
          p_abs = p_abs,
          p_form = p_form,
          rr_form = rr_fd_lc,
          rr_function = rr_lc_fem_fun
        )
      }, error = function(e) NULL)
      
      # Update results in bcan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$Point_Estimate)) {
        lc_fem[i, paste0("Fem", j, "_point")] <- result$Point_Estimate
        lc_fem[i, paste0("Fem", j, "_lower")] <- result$Lower_CI
        lc_fem[i, paste0("Fem", j, "_upper")] <- result$Upper_CI
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", lc_fem$Year[i], "Group Fem", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", lc_fem$Year[i], "Group Fem", j, "- Invalid data\n")
    }
  }
}

lc_fem <- lc_fem %>% 
  mutate(disease = "Liver Cirrhosis")

# PANCREATITIS 
b1_panc <- 0.013
var_panc <- 0.003803**2

# for former drinkers (both sexes)
rr_panc_fd <- 2.2

panc_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

# estimating PAF and confidence intervals
for (i in 1:length(panc_male$Year)) {
  for (j in 1:4) {  # Iterate over age group
    # Access the relevant elements from the lists
    gamma_fit <- g_male_list[[as.character(panc_male$Year[i])]][[j]]
    p_abs <- p_abs_list_male[[as.character(panc_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_male[[as.character(panc_male$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      result <- tryCatch({
        confint_paf(
          gamma = gamma_fit,
          beta = b1_panc,
          var_beta = var_panc,
          p_abs = p_abs,
          p_form = p_form,
          rr_form = rr_panc_fd,
          rr_function = rr_linear
        )
      }, error = function(e) NULL)
      
      # Update results in bcan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$Point_Estimate)) {
        panc_male[i, paste0("Male", j, "_point")] <- result$Point_Estimate
        panc_male[i, paste0("Male", j, "_lower")] <- result$Lower_CI
        panc_male[i, paste0("Male", j, "_upper")] <- result$Upper_CI
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", panc_male$Year[i], "Group Male", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", panc_male$Year[i], "Group Male", j, "- Invalid data\n")
    }
  }
}
panc_male <- panc_male %>% 
  mutate(disease = "Acute Pancreatitis")

b1_panc_fem <- -0.0277886
b2_panc_fem <- 0.0611466
cov_panc_fem <- matrix(c(0.01127452, -0.00015821,
-0.00015821, 0.01762052), nrow = 2, byrow = T)
rr_panc_fem <- function(betas, x) {
  # Initialize an empty vector to store results
  rr <- numeric(length(x))
  
  # Apply the piecewise function
  rr[x < 3] <- exp(betas[1] * x[x < 3])
  rr[x >= 3 & x < 15] <- exp(betas[1] * x[x >= 3 & x < 15] + betas[2] * ((x[x >= 3 & x < 15] - 3)^2) / 373)
  rr[x >= 15 & x < 40] <- exp(betas[1] * x[x >= 15 & x < 40] + betas[2] * ((x[x >= 3 & x < 15] - 3)^2) / 373 - ((x[x >= 15 & x < 40] - 15)^2 * 37) / (372 * 25))
  rr[x >= 40 & x < 108] <- exp(betas[1] * x[x >= 40 & x < 108] + betas[2] * ((x[x >= 3 & x < 15] - 3)^2) / 373 - ((x[x >= 15 & x < 40] - 15)^2 * 37) / (372 * 25) - ((x[x >= 40 & x < 108] - 40)^2 * 15) / (108^2))
  rr[x >= 108] <- exp(2.327965)
  
  return(rr)
}

panc_fem <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)


# PAF estimation for females
for (i in 1:length(panc_fem$Year)) {
  for (j in 1:4) {
    # Access the relevant elements from the lists
    gamma_fit <- g_fem_list[[as.character(panc_fem$Year[i])]][[j]]
    p_abs <- p_abs_list_fem[[as.character(panc_fem$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_fem[[as.character(panc_fem$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      # Call confint_paf_vcov with the required inputs
      result <- tryCatch({
        confint_paf_vcov(
          gamma = gamma_fit,
          betas = c(b1_panc_fem, b2_panc_fem),
          cov_matrix = cov_panc_fem,
          p_abs = p_abs,
          p_form = p_form,
          rr_fd = rr_panc_fd,
          rr_function = rr_panc_fem
        )
      }, error = function(e) NULL)
      
      # Update results in locan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$point_estimate)) {
        panc_fem[i, paste0("Fem", j, "_point")] <- result$point_estimate
        panc_fem[i, paste0("Fem", j, "_lower")] <- result$lower_ci
        panc_fem[i, paste0("Fem", j, "_upper")] <- result$upper_ci
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", panc_fem$Year[i], "Group Fem", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", panc_fem$Year[i], "Group Fem", j, "- Invalid data\n")
    }
  }
}

panc_fem <- panc_fem %>% 
  mutate(disease = "Acute Pancreatitis")

aaf_other_fem <- bind_rows(panc_fem, lc_fem,lri_female, hiv_female,
                            dm_fem, tb_female, panc_fem) %>% 
  rename(AAF_ag1 = Fem1_point, LL_ag1 = Fem1_lower, UL_ag1 = Fem1_upper,
         AAF_ag2 = Fem2_point, LL_ag2 = Fem2_lower, UL_ag2 = Fem2_upper,
         AAF_ag3 = Fem3_point, LL_ag3 = Fem3_lower, UL_ag3 = Fem3_upper,
         AAF_ag4 = Fem4_point, LL_ag4 = Fem4_lower, UL_ag4 = Fem4_upper)

aaf_other_male <- bind_rows(panc_male, lc_male,lri_male, hiv_male,
                            dm_male, tb_male, panc_male) %>% 
  rename(AAF_ag1 = Male1_point, LL_ag1 = Male1_lower, UL_ag1 = Male1_upper,
         AAF_ag2 = Male2_point, LL_ag2 = Male2_lower, UL_ag2 = Male2_upper,
         AAF_ag3 = Male3_point, LL_ag3 = Male3_lower, UL_ag3 = Male3_upper,
         AAF_ag4 = Male4_point, LL_ag4 = Male4_lower, UL_ag4 = Male4_upper)

############
# INJURIES #
############
x_vals_nhed <- seq(0.1, 60, length.out = 1500)
x_vals_hed <- seq(60,150, length.out = 1500)
# ROAD INJURIES
b1_ri <- 0.00455
var_b1_ri <- 0.001687**2

ri_hed_fun <- function(beta,x){
  rr = 1.49*exp(beta[1]*x)
  return(rr)
}
rr_linear <- function(x, b) {
  exp(b * x)
}


ri_fem <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)


# Loop over each year and each female age group
for (i in 1:length(ri_fem$Year)) {
  for (j in 1:4) {
    # Extract relevant elements from the lists
    gamma_fit_nhed <- g_fem_hed_list[[as.character(ri_fem$Year[i])]][[j]]$nhed
    gamma_fit_hed <- g_fem_hed_list[[as.character(ri_fem$Year[i])]][[j]]$hed
    p_abs <- p_abs_list_fem[[as.character(ri_fem$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_fem[[as.character(ri_fem$Year[i])]][[paste0("edad_tramo_", j)]]
    p_hed <- p_hed_list_fem[[j]][i]
    
    # Combine gamma fits for input
    gammas <- list(gamma_fit_nhed, gamma_fit_hed)
    
    # Validate inputs
    if (!is.null(gammas[[1]]) && !is.null(gammas[[2]]) &&
        !is.na(p_abs) && !is.na(p_form) && !is.na(p_hed)) {
      
      # Try to compute the PAF
      result <- tryCatch({
        confint_paf_hed(
          gammas = gammas, 
          beta = b1_ri, 
          cov_matrix = var_b1_ri, 
          p_abs = p_abs, 
          p_form = p_form, 
          rr_fd = 1, 
          rr_function_nhed = rr_linear,
          rr_function_hed = ri_hed_fun,
          p_hed = p_hed
        )
      }, error = function(e) {
        cat("Error in Year:", ri_fem$Year[i], "Group Fem", j, "\n")
        NULL
      })
      
      # Store results in the data frame if successful
      if (!is.null(result) && !is.null(result$point_estimate)) {
        ri_fem[i, paste0("Fem", j, "_point")] <- result$point_estimate
        ri_fem[i, paste0("Fem", j, "_lower")] <- result$lower_ci
        ri_fem[i, paste0("Fem", j, "_upper")] <- result$upper_ci
      } else {
        cat("Calculation failed for Year:", ri_fem$Year[i], "Group Fem", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", ri_fem$Year[i], "Group Fem", j, "- Invalid data\n")
    }
  }
}

ri_fem <- ri_fem %>% 
  mutate(disease = "Road Injuries")

# Initialize a data frame to store results for males
ri_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

# Loop over each year and each female age group
for (i in 1:length(ri_male$Year)) {
  for (j in 1:4) {
    # Extract relevant elements from the lists
    gamma_fit_nhed <- g_male_hed_list[[as.character(ri_male$Year[i])]][[j]]$nhed
    gamma_fit_hed <- g_male_hed_list[[as.character(ri_male$Year[i])]][[j]]$hed
    p_abs <- p_abs_list_male[[as.character(ri_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_male[[as.character(ri_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_hed <- p_hed_list_male[[j]][i]
    
    # Combine gamma fits for input
    gammas <- list(gamma_fit_nhed, gamma_fit_hed)
    
    # Validate inputs
    if (!is.null(gammas[[1]]) && !is.null(gammas[[2]]) &&
        !is.na(p_abs) && !is.na(p_form) && !is.na(p_hed)) {
      
      # Try to compute the PAF
      result <- tryCatch({
        confint_paf_hed(
          gammas = gammas, 
          beta = b1_ri, 
          cov_matrix = var_b1_ri, 
          p_abs = p_abs, 
          p_form = p_form, 
          rr_fd = 1, 
          rr_function_nhed = rr_linear,
          rr_function_hed = ri_hed_fun,
          p_hed = p_hed
        )
      }, error = function(e) {
        cat("Error in Year:", ri_male$Year[i], "Group Male", j, "\n")
        NULL
      })
      
      # Store results in the data frame if successful
      if (!is.null(result) && !is.null(result$point_estimate)) {
        ri_male[i, paste0("Male", j, "_point")] <- result$point_estimate
        ri_male[i, paste0("Male", j, "_lower")] <- result$lower_ci
        ri_male[i, paste0("Male", j, "_upper")] <- result$upper_ci
      } else {
        cat("Calculation failed for Year:", ri_male$Year[i], "Group Male", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", ri_male$Year[i], "Group Male", j, "- Invalid data\n")
    }
  }
}


# Display results
print(ri_male)

ri_male <- ri_male %>% 
  mutate(disease = "Road Injuries")
#POISONINGS, FALLS, FIRE, HEAT AND HOT SUBSTANCES, DROWING,
# EXPOSURE TO MECHANICAL FORCES, OTHER UNINTENTIONAL INJURIES

# Initialize a data frame to store results for males
injuries_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

# Loop over each year and each female age group
for (i in 1:length(injuries_male$Year)) {
  for (j in 1:4) {
    # Extract relevant elements from the lists
    gamma_fit_nhed <- g_male_hed_list[[as.character(injuries_male$Year[i])]][[j]]$nhed
    gamma_fit_hed <- g_male_hed_list[[as.character(injuries_male$Year[i])]][[j]]$hed
    p_abs <- p_abs_list_male[[as.character(injuries_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_male[[as.character(injuries_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_hed <- p_hed_list_male[[j]][i]
    
    # Combine gamma fits for input
    gammas <- list(gamma_fit_nhed, gamma_fit_hed)
    
    # Validate inputs
    if (!is.null(gammas[[1]]) && !is.null(gammas[[2]]) &&
        !is.na(p_abs) && !is.na(p_form) && !is.na(p_hed)) {
      
      # Try to compute the PAF
      result <- tryCatch({
        confint_paf_hed(
          gammas = gammas, 
          beta = b1_ri, 
          cov_matrix = var_b1_ri, 
          p_abs = p_abs, 
          p_form = p_form, 
          rr_fd = 1, 
          rr_function_nhed = rr_linear,
          rr_function_hed = ri_hed_fun,
          p_hed = p_hed
        )
      }, error = function(e) {
        cat("Error in Year:", injuries_male$Year[i], "Group Male", j, "\n")
        NULL
      })
      
      # Store results in the data frame if successful
      if (!is.null(result) && !is.null(result$point_estimate)) {
        injuries_male[i, paste0("Male", j, "_point")] <- result$point_estimate
        injuries_male[i, paste0("Male", j, "_lower")] <- result$lower_ci
        injuries_male[i, paste0("Male", j, "_upper")] <- result$upper_ci
      } else {
        cat("Calculation failed for Year:", injuries_male$Year[i], "Group Male", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", injuries_male$Year[i], "Group Male", j, "- Invalid data\n")
    }
  }
}

injuries_male <- injuries_male %>% 
  mutate(disease = "Unintentional Injuries")

injuries_fem <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)
# Loop over each year and each female age group
for (i in 1:length(injuries_fem$Year)) {
  for (j in 1:4) {
    # Extract relevant elements from the lists
    gamma_fit_nhed <- g_fem_hed_list[[as.character(injuries_fem$Year[i])]][[j]]$nhed
    gamma_fit_hed <- g_fem_hed_list[[as.character(injuries_fem$Year[i])]][[j]]$hed
    p_abs <- p_abs_list_fem[[as.character(injuries_fem$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_fem[[as.character(injuries_fem$Year[i])]][[paste0("edad_tramo_", j)]]
    p_hed <- p_hed_list_fem[[j]][i]
    
    # Combine gamma fits for input
    gammas <- list(gamma_fit_nhed, gamma_fit_hed)
    
    # Validate inputs
    if (!is.null(gammas[[1]]) && !is.null(gammas[[2]]) &&
        !is.na(p_abs) && !is.na(p_form) && !is.na(p_hed)) {
      
      # Try to compute the PAF
      result <- tryCatch({
        confint_paf_hed(
          gammas = gammas, 
          beta = b1_ri, 
          cov_matrix = var_b1_ri, 
          p_abs = p_abs, 
          p_form = p_form, 
          rr_fd = 1, 
          rr_function_nhed = rr_linear,
          rr_function_hed = ri_hed_fun,
          p_hed = p_hed
        )
      }, error = function(e) {
        cat("Error in Year:", injuries_fem$Year[i], "Group Fem", j, "\n")
        NULL
      })
      
      # Store results in the data frame if successful
      if (!is.null(result) && !is.null(result$point_estimate)) {
        injuries_fem[i, paste0("Fem", j, "_point")] <- result$point_estimate
        injuries_fem[i, paste0("Fem", j, "_lower")] <- result$lower_ci
        injuries_fem[i, paste0("Fem", j, "_upper")] <- result$upper_ci
      } else {
        cat("Calculation failed for Year:", injuries_fem$Year[i], "Group Fem", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", injuries_fem$Year[i], "Group Fem", j, "- Invalid data\n")
    }
  }
}
injuries_fem <- injuries_fem %>% 
  mutate(disease = "Unintentional Injuries")

# self inflicted injuries and interpersonal violence
inj_hed_fun <- function(beta,x){
  rr = 1.70*exp(beta[1]*x)
  return(rr)
}

# Initialize a data frame to store results for males
violence_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

# Loop over each year and each female age group
for (i in 1:length(violence_male$Year)) {
  for (j in 1:4) {
    # Extract relevant elements from the lists
    gamma_fit_nhed <- g_male_hed_list[[as.character(violence_male$Year[i])]][[j]]$nhed
    gamma_fit_hed <- g_male_hed_list[[as.character(violence_male$Year[i])]][[j]]$hed
    p_abs <- p_abs_list_male[[as.character(violence_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_male[[as.character(violence_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_hed <- p_hed_list_male[[j]][i]
    
    # Combine gamma fits for input
    gammas <- list(gamma_fit_nhed, gamma_fit_hed)
    
    # Validate inputs
    if (!is.null(gammas[[1]]) && !is.null(gammas[[2]]) &&
        !is.na(p_abs) && !is.na(p_form) && !is.na(p_hed)) {
      
      # Try to compute the PAF
      result <- tryCatch({
        confint_paf_hed(
          gammas = gammas, 
          beta = b1_ri, 
          cov_matrix = var_b1_ri, 
          p_abs = p_abs, 
          p_form = p_form, 
          rr_fd = 1, 
          rr_function_nhed = rr_linear,
          rr_function_hed = inj_hed_fun,
          p_hed = p_hed
        )
      }, error = function(e) {
        cat("Error in Year:", violence_male$Year[i], "Group Male", j, "\n")
        NULL
      })
      
      # Store results in the data frame if successful
      if (!is.null(result) && !is.null(result$point_estimate)) {
        violence_male[i, paste0("Male", j, "_point")] <- result$point_estimate
        violence_male[i, paste0("Male", j, "_lower")] <- result$lower_ci
        violence_male[i, paste0("Male", j, "_upper")] <- result$upper_ci
      } else {
        cat("Calculation failed for Year:", violence_male$Year[i], "Group Male", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", violence_male$Year[i], "Group Male", j, "- Invalid data\n")
    }
  }
}
violence_male <- violence_male %>% 
  mutate(disease = "Intentional Injuries")

violence_fem <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)

# Loop over each year and each female age group
for (i in 1:length(violence_fem$Year)) {
  for (j in 1:4) {
    # Extract relevant elements from the lists
    gamma_fit_nhed <- g_fem_hed_list[[as.character(violence_fem$Year[i])]][[j]]$nhed
    gamma_fit_hed <- g_fem_hed_list[[as.character(violence_fem$Year[i])]][[j]]$hed
    p_abs <- p_abs_list_fem[[as.character(violence_fem$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_fem[[as.character(violence_fem$Year[i])]][[paste0("edad_tramo_", j)]]
    p_hed <- p_hed_list_fem[[j]][i]
    
    # Combine gamma fits for input
    gammas <- list(gamma_fit_nhed, gamma_fit_hed)
    
    # Validate inputs
    if (!is.null(gammas[[1]]) && !is.null(gammas[[2]]) &&
        !is.na(p_abs) && !is.na(p_form) && !is.na(p_hed)) {
      
      # Try to compute the PAF
      result <- tryCatch({
        confint_paf_hed(
          gammas = gammas, 
          beta = b1_ri, 
          cov_matrix = var_b1_ri, 
          p_abs = p_abs, 
          p_form = p_form, 
          rr_fd = 1, 
          rr_function_nhed = rr_linear,
          rr_function_hed = inj_hed_fun,
          p_hed = p_hed
        )
      }, error = function(e) {
        cat("Error in Year:", violence_fem$Year[i], "Group Fem", j, "\n")
        NULL
      })
      
      # Store results in the data frame if successful
      if (!is.null(result) && !is.null(result$point_estimate)) {
        violence_fem[i, paste0("Fem", j, "_point")] <- result$point_estimate
        violence_fem[i, paste0("Fem", j, "_lower")] <- result$lower_ci
        violence_fem[i, paste0("Fem", j, "_upper")] <- result$upper_ci
      } else {
        cat("Calculation failed for Year:", violence_fem$Year[i], "Group Fem", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", violence_fem$Year[i], "Group Fem", j, "- Invalid data\n")
    }
  }
}

violence_fem <- violence_fem %>% 
  mutate(disease = "Intentional Injuries")

aaf_inj_fem <- bind_rows(violence_fem, ri_fem,injuries_fem) %>% 
  rename(AAF_ag1 = Fem1_point, LL_ag1 = Fem1_lower, UL_ag1 = Fem1_upper,
         AAF_ag2 = Fem2_point, LL_ag2 = Fem2_lower, UL_ag2 = Fem2_upper,
         AAF_ag3 = Fem3_point, LL_ag3 = Fem3_lower, UL_ag3 = Fem3_upper,
         AAF_ag4 = Fem4_point, LL_ag4 = Fem4_lower, UL_ag4 = Fem4_upper)

aaf_inj_male <- bind_rows(violence_male, ri_male,injuries_male) %>% 
  rename(AAF_ag1 = Male1_point, LL_ag1 = Male1_lower, UL_ag1 = Male1_upper,
         AAF_ag2 = Male2_point, LL_ag2 = Male2_lower, UL_ag2 = Male2_upper,
         AAF_ag3 = Male3_point, LL_ag3 = Male3_lower, UL_ag3 = Male3_upper,
         AAF_ag4 = Male4_point, LL_ag4 = Male4_lower, UL_ag4 = Male4_upper)

##################
# CARDIOVASCULAR #
##################

# HYPERTENSIVE HEART DISEASE (MALE)
b1_hhd_male <- 0.0150537
var_hhd_male <- 0.0024196**2
rr_fd_hhd_male <- 1.03

rr_hhd_male_fun <- function(x, b) {

  result <- numeric(length(x))
  
  result[x >= 0 & x < 21] <- exp(b[1] * x[x >= 0 & x < 21] - b[1] * x[x >= 0 & x < 21]^2 / 75)
  result[x >= 21 & x < 75] <- exp(b[1] * x[x >= 21 & x < 75] - b[1] * (x[x >= 21 & x < 75]^2 - (75 * x[x >= 21 & x < 75] - 21 * 75)) / 75)
  result[x >= 75] <- exp(b[1] * x[x >= 75] - b[1])
  
  return(result)
}



# Define un data frame para almacenar los resultados para hombres
hhd_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

# estimating PAF and confidence intervals
for (i in 1:length(hhd_male$Year)) {
  for (j in 1:4) {  # Iterate over age group
    # Access the relevant elements from the lists
    gamma_fit <- g_male_list[[as.character(hhd_male$Year[i])]][[j]]
    p_abs <- p_abs_list_male[[as.character(hhd_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_male[[as.character(hhd_male$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      result <- tryCatch({
        confint_paf(
          gamma = gamma_fit,
          beta = b1_hhd_male,
          var_beta = var_hhd_male,
          p_abs = p_abs,
          p_form = p_form,
          rr_form = rr_fd_hhd_male,
          rr_function = rr_hhd_male_fun
        )
      }, error = function(e) NULL)
      
      # Update results in bcan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$Point_Estimate)) {
        hhd_male[i, paste0("Male", j, "_point")] <- result$Point_Estimate
        hhd_male[i, paste0("Male", j, "_lower")] <- result$Lower_CI
        hhd_male[i, paste0("Male", j, "_upper")] <- result$Upper_CI
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", hhd_male$Year[i], "Group Male", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", hhd_male$Year[i], "Group Male", j, "- Invalid data\n")
    }
  }
}

hhd_male <- hhd_male %>% 
  mutate(disease = "Hypertensive Heart Disease")

# HYPERTENSIVE HEART DISEASE (FEMALE)
b1_hhd_fem <- 0.0154196
var_hhd_fem <- 0.006459**2

rr_hhd_fem <- function(x, b = NULL) {

  result <- numeric(length(x))
  
  # Below 18.9517
  result[x >= 0 & x < 18.9517] <- exp(1)  # Fixed baseline risk
  
  # Between 18.9517 and 75
  result[x >= 18.9517 & x < 75] <- exp(
    -b[1] * x[x >= 18.9517 & x < 75] + 0.0217586 +
      (x[x >= 18.9517 & x < 75]^3 - 20 * (x[x >= 18.9517 & x < 75] - 10)^3 - 10 * (x[x >= 18.9517 & x < 75] - 20)^3) /
      (20^3)
  )
  
  # Above 75
  result[x >= 75] <- exp(-0.9649937)
  
  return(result)
}

# Define a data frame to store results
hhd_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)

# estimating PAF and confidence intervals
for (i in 1:length(hhd_female$Year)) {
  for (j in 1:4) {  # Iterate over edad_tramo
    # Access the relevant elements from the lists
    gamma_fit <- g_fem_list[[as.character(hhd_female$Year[i])]][[j]]
    p_abs <- p_abs_list_fem[[as.character(hhd_female$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_fem[[as.character(hhd_female$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      result <- tryCatch({
        confint_paf(
          gamma = gamma_fit,
          beta = b1_hhd_fem,
          var_beta = var_hhd_fem,
          p_abs = p_abs,
          p_form = p_form,
          rr_form = 1.05,
          rr_function = rr_hhd_fem
        )
      }, error = function(e) NULL)
      
      # Update results in bcan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$Point_Estimate)) {
        hhd_female[i, paste0("Fem", j, "_point")] <- result$Point_Estimate
        hhd_female[i, paste0("Fem", j, "_lower")] <- result$Lower_CI
        hhd_female[i, paste0("Fem", j, "_upper")] <- result$Upper_CI
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", hhd_female$Year[i], "Group Fem", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", hhd_female$Year[i], "Group Fem", j, "- Invalid data\n")
    }
  }
}

hhd_female <- hhd_female %>% 
  mutate(disease = "Hypertensive Heart Disease")

#------------------------------#
# ISCHAEMIC HEART DISEASE MALE #
#------------------------------#

b_ihd_male <- 0.002211
var_ihd_male <- 0.0024196**2
rr_fd_ihd_male <- 1.25

# creating the data frame for males
ihd_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

# estimating PAF and confidence intervals
for (i in 1:length(ihd_male$Year)) {
  for (j in 1:4) {  # Iterate over age group
    # Access the relevant elements from the lists
    gamma_fit <- g_male_list[[as.character(ihd_male$Year[i])]][[j]]
    p_abs <- p_abs_list_male[[as.character(ihd_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_male[[as.character(ihd_male$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      result <- tryCatch({
        confint_paf(
          gamma = gamma_fit,
          beta = b_ihd_male,
          var_beta = var_ihd_male,
          p_abs = p_abs,
          p_form = p_form,
          rr_form = rr_fd_ihd_male,
          rr_function = rr_linear
        )
      }, error = function(e) NULL)
      
      # Update results in bcan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$Point_Estimate)) {
        ihd_male[i, paste0("Male", j, "_point")] <- result$Point_Estimate
        ihd_male[i, paste0("Male", j, "_lower")] <- result$Lower_CI
        ihd_male[i, paste0("Male", j, "_upper")] <- result$Upper_CI
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", ihd_male$Year[i], "Group Male", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", ihd_male$Year[i], "Group Male", j, "- Invalid data\n")
    }
  }
}

ihd_male <- ihd_male %>% 
  mutate(disease = "Ischaemic Heart Disease")

# ISCHAEMIC HEART DISEASE FEMALE
b1_ihd_fem <- -0.0525288
b2_ihd_fem <- 0.0153856
cov_ihd_fem <- matrix(c(0.032510,0,
                        0,0.007925), nrow = 2, byrow = T)

rr_ihd_fem <- function(betas, x){
exp(-betas[1]*x+betas[2]*x*log(x))
}

ihd_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)
# PAF estimation for females
for (i in 1:length(ihd_female$Year)) {
  for (j in 1:4) {
    # Access the relevant elements from the lists
    gamma_fit <- g_fem_list[[as.character(ihd_female$Year[i])]][[j]]
    p_abs <- p_abs_list_fem[[as.character(ihd_female$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_fem[[as.character(ihd_female$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      # Call confint_paf_vcov with the required inputs
      result <- tryCatch({
        confint_paf_vcov(
          gamma = gamma_fit,
          betas = c(b1_ihd_fem, b2_ihd_fem),
          cov_matrix = cov_ihd_fem,
          p_abs = p_abs,
          p_form = p_form,
          rr_fd = 1.54,
          rr_function = rr_ihd_fem
        )
      }, error = function(e) NULL)
      
      # Update results in locan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$point_estimate)) {
        ihd_female[i, paste0("Fem", j, "_point")] <- result$point_estimate
        ihd_female[i, paste0("Fem", j, "_lower")] <- result$lower_ci
        ihd_female[i, paste0("Fem", j, "_upper")] <- result$upper_ci
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", ihd_female$Year[i], "Group Fem", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", ihd_female$Year[i], "Group Fem", j, "- Invalid data\n")
    }
  }
}

ihd_female <- ihd_female %>% 
  mutate(disease = "Ischaemic Heart Disease")

# INTRACEREBRAL HAEMORRHAGE (MALES)
b1_ich_male <- 0.006898937
var_b1_ich_male <- 0.001142**2 
rr_fd_ich_male <- 1.36

# Define un data frame para almacenar los resultados para hombres
ich_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)

# estimating PAF and confidence intervals
for (i in 1:length(ich_male$Year)) {
  for (j in 1:4) {  # Iterate over age group
    # Access the relevant elements from the lists
    gamma_fit <- g_male_list[[as.character(ich_male$Year[i])]][[j]]
    p_abs <- p_abs_list_male[[as.character(ich_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_male[[as.character(ich_male$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      result <- tryCatch({
        confint_paf(
          gamma = gamma_fit,
          beta = b1_ich_male,
          var_beta = var_b1_ich_male,
          p_abs = p_abs,
          p_form = p_form,
          rr_form = rr_fd_ich_male,
          rr_function = rr_linear
        )
      }, error = function(e) NULL)
      
      # Update results in bcan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$Point_Estimate)) {
        ich_male[i, paste0("Male", j, "_point")] <- result$Point_Estimate
        ich_male[i, paste0("Male", j, "_lower")] <- result$Lower_CI
        ich_male[i, paste0("Male", j, "_upper")] <- result$Upper_CI
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", ich_male$Year[i], "Group Male", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", ich_male$Year[i], "Group Male", j, "- Invalid data\n")
    }
  }
}

ich_male <- ich_male %>% 
  mutate(disease = "Intracerebral Haemorrhage")

# INTRACEREBRAL HAEMORRHAGE (FEMALES)
b1_ich_fem <- 0.01466406
var_b1_ich_fem <- 0.003544172**2

rr_fd_ich_fem <- 1.36

ich_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)

# estimating PAF and confidence intervals
for (i in 1:length(ich_female$Year)) {
  for (j in 1:4) {  # Iterate over edad_tramo
    # Access the relevant elements from the lists
    gamma_fit <- g_fem_list[[as.character(ich_female$Year[i])]][[j]]
    p_abs <- p_abs_list_fem[[as.character(ich_female$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_fem[[as.character(ich_female$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      result <- tryCatch({
        confint_paf(
          gamma = gamma_fit,
          beta = b1_ich_fem,
          var_beta = var_b1_ich_fem,
          p_abs = p_abs,
          p_form = p_form,
          rr_form = rr_fd_ich_fem,
          rr_function = rr_linear
        )
      }, error = function(e) NULL)
      
      # Update results in bcan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$Point_Estimate)) {
        ich_female[i, paste0("Fem", j, "_point")] <- result$Point_Estimate
        ich_female[i, paste0("Fem", j, "_lower")] <- result$Lower_CI
        ich_female[i, paste0("Fem", j, "_upper")] <- result$Upper_CI
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", ich_female$Year[i], "Group Fem", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", ich_female$Year[i], "Group Fem", j, "- Invalid data\n")
    }
  }
}

ich_female <- ich_female %>% 
  mutate(disease = "Intracerebral Haemorrhage")

# ISCHAEMIC STROKE MALE
b1_is_male <- -0.141950 
b2_is_male <- 0.039613 
cov_is_male <- matrix(c(0.012866**2,0,
                        0, 0.001782**2), nrow = 2, byrow = T)

rr_is_male_fun <- function(x, betas) {
  b1 <- betas[1]
  b2 <- betas[2]
  # Initialize the RR as a vector of the same length as x
  rr <- numeric(length(x))
  
 rr = exp(betas[1]*sqrt(x)+betas[2]*sqrt(x)*log(x))
  return(rr)
}

is_male <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020,2022),
  Male1_point = NA, Male1_lower = NA, Male1_upper = NA,
  Male2_point = NA, Male2_lower = NA, Male2_upper = NA,
  Male3_point = NA, Male3_lower = NA, Male3_upper = NA,
  Male4_point = NA, Male4_lower = NA, Male4_upper = NA
)


# PAF estimation for males
for (i in 1:length(is_male$Year)) {
  for (j in 1:4) {
    # Access the relevant elements from the lists
    gamma_fit <- g_male_list[[as.character(is_male$Year[i])]][[j]]
    p_abs <- p_abs_list_male[[as.character(is_male$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_male[[as.character(is_male$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      # Call confint_paf_vcov with the required inputs
      result <- tryCatch({
        confint_paf_vcov(
          gamma = gamma_fit,
          betas = c(b1_is_male, b2_is_male),
          cov_matrix = cov_is_male,
          p_abs = p_abs,
          p_form = p_form,
          rr_fd = 0.97,
          rr_function = rr_is_male_fun
        )
      }, error = function(e) NULL)
      
      # Update results in locan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$point_estimate)) {
        is_male[i, paste0("Male", j, "_point")] <- result$point_estimate
        is_male[i, paste0("Male", j, "_lower")] <- result$lower_ci
        is_male[i, paste0("Male", j, "_upper")] <- result$upper_ci
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", is_male$Year[i], "Group Male", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", is_male$Year[i], "Group Male", j, "- Invalid data\n")
    }
  }
}
is_male <- is_male %>% 
  mutate(disease = "Ischaemic Stroke")

# ISCHAEMIC STROKE FEMALE
b1_is_female <- -0.248768
b2_is_female1 <- 0.03708724
cov_is_fem <- matrix(c(0.019163**2, 0, # that values are for the standard error
                       0, 0.000523**2),2,2)

rr_is_fem_fun <- function(x, betas) {
  b1 <- betas[1]
  b2 <- betas[2]
  # Initialize the RR as a vector of the same length as x
  rr <- numeric(length(x))
  
  rr = exp(b1*sqrt(x)+b2*x)

  return(rr)
}

# Define a data frame to store results
is_female <- data.frame(
  Year = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
  Fem1_point = NA, Fem1_lower = NA, Fem1_upper = NA,
  Fem2_point = NA, Fem2_lower = NA, Fem2_upper = NA,
  Fem3_point = NA, Fem3_lower = NA, Fem3_upper = NA,
  Fem4_point = NA, Fem4_lower = NA, Fem4_upper = NA
)



# PAF estimation for females
for (i in 1:length(is_female$Year)) {
  for (j in 1:4) {
    # Access the relevant elements from the lists
    gamma_fit <- g_fem_list[[as.character(is_female$Year[i])]][[j]]
    p_abs <- p_abs_list_fem[[as.character(is_female$Year[i])]][[paste0("edad_tramo_", j)]]
    p_form <- p_form_list_fem[[as.character(is_female$Year[i])]][[paste0("edad_tramo_", j)]]
    
    # Validate inputs
    if (!is.null(gamma_fit) && !is.na(p_abs) && !is.na(p_form)) {
      # Call confint_paf_vcov with the required inputs
      result <- tryCatch({
        confint_paf_vcov(
          gamma = gamma_fit,
          betas = c(b1_is_female, b2_is_female1),
          cov_matrix = cov_is_fem,
          p_abs = p_abs,
          p_form = p_form,
          rr_fd = 0.97,
          rr_function = rr_is_fem_fun
        )
      }, error = function(e) NULL)
      
      # Update results in locan_female if the calculation succeeds
      if (!is.null(result) && !is.null(result$point_estimate)) {
        is_female[i, paste0("Fem", j, "_point")] <- result$point_estimate
        is_female[i, paste0("Fem", j, "_lower")] <- result$lower_ci
        is_female[i, paste0("Fem", j, "_upper")] <- result$upper_ci
      } else {
        # Log if calculation fails
        cat("Calculation failed for Year:", is_female$Year[i], "Group Fem", j, "\n")
      }
    } else {
      # Log missing or invalid data
      cat("Skipping Year:", is_female$Year[i], "Group Fem", j, "- Invalid data\n")
    }
  }
}

is_female <- is_female %>% 
  mutate(disease = "Ischaemic Stroke")

aaf_cv_fem <- bind_rows(is_female, ich_female,hhd_female, ihd_female) %>% 
  rename(AAF_ag1 = Fem1_point, LL_ag1 = Fem1_lower, UL_ag1 = Fem1_upper,
         AAF_ag2 = Fem2_point, LL_ag2 = Fem2_lower, UL_ag2 = Fem2_upper,
         AAF_ag3 = Fem3_point, LL_ag3 = Fem3_lower, UL_ag3 = Fem3_upper,
         AAF_ag4 = Fem4_point, LL_ag4 = Fem4_lower, UL_ag4 = Fem4_upper)

aaf_cv_male <- bind_rows(is_male, ich_male,hhd_male, ihd_male) %>% 
  rename(AAF_ag1 = Male1_point, LL_ag1 = Male1_lower, UL_ag1 = Male1_upper,
         AAF_ag2 = Male2_point, LL_ag2 = Male2_lower, UL_ag2 = Male2_upper,
         AAF_ag3 = Male3_point, LL_ag3 = Male3_lower, UL_ag3 = Male3_upper,
         AAF_ag4 = Male4_point, LL_ag4 = Male4_lower, UL_ag4 = Male4_upper)

aaf_fem <- bind_rows(aaf_cancer_fem, aaf_cv_fem,aaf_inj_fem,aaf_other_fem, epi_female)
aaf_male <- bind_rows(aaf_cancer_male, aaf_cv_male,aaf_inj_male,aaf_other_male, epi_male)

rio::export(aaf_fem, "AAF FEMALES.xlsx")
rio::export(aaf_male, "AAF MALES.xlsx")

#########################
# MORTALITY PREPARATION #
#########################
rm(list = ls());gc()

aaf_fem <- rio::import("AAF FEMALES.xlsx")
aaf_male <- rio::import("AAF MALES.xlsx")

def <- read_rds("data_mortality.rds")

# POPULATION INFORMATION



###########
# AFF = 1 #
###########

def <- def %>%
  mutate(des_men = as.integer(DIAG1 %in% c("F100", "F101", "F102", "F103", "F104", "F105", "F106", "F107", "F108", "F109")), #Mental and behavioral disorders due to alcohol
         deg_nerv = ifelse(DIAG1 == "G312",1,0),
         polineu = ifelse(DIAG1 == "G621",1,0),
         cardiomio = ifelse(DIAG1 == "I426",1,0),
         pancreati_oh = ifelse(DIAG1 == "K860",1,0),
         gastrit = ifelse(DIAG1 == "K292",1,0),
         enven_acc = as.integer(DIAG2 %in% paste0("X45", 0:9)),
         enven_int = as.integer(DIAG2 %in% paste0("X65", 0:9)),
         enven_indet = as.integer(DIAG2 %in% paste0("Y15", 0:9)))

def <- def %>% 
  mutate(aaf1 = rowSums(def[,6:14]))
         

################################
# PARTIALLY ATTRIBUTABLE TO OH #
################################

# NEUROPSYCHIATRIC

# Epilepsy
epilepsy_codes <- c(paste0("C40", 0:9), paste0("C41", 0:9))

# CARDIOVASCULAR

# Hypertensive disease
hhd_codes <- c(paste0("I10", 0:9), paste0("I11", 0:9), paste0("I12", 0:9),
               paste0("I13", 0:9), paste0("I14", 0:9), paste0("I15", 0:9))
# Ischemic heart disease
ihd_codes <- c(paste0("I20", 0:9), paste0("I21", 0:9), paste0("I22", 0:9),
               paste0("I23", 0:9),paste0("I24", 0:9), paste0("I25", 0:9))
# Hemorrhagic stroke
ich_codes <- c(paste0("I60", 0:9), paste0("I61", 0:9), paste0("I62", 0:9))

# Ischemic stroke
is_codes <- c(paste0("I63", 0:9), paste0("I64", 0:9), paste0("I65", 0:9),
              paste0("I66", 0:9))

# CANCER
# Mouth, oropharynx and laryngeal cancer
opcan_codes <- paste0("C0", sprintf("%02d", 0:140))

# Esophageal cancer
oescan_codes <- paste0("C15", 0:9)

# Colon and rectum cancer
crcan_codes <- c(paste0("C18, 0:9"), "C19X", "C20X")

# Liver cancer
lican_codes <- paste0("C22", 0:9)

# Laryngeal cancer
lxcan_codes <- paste0("C32", 0:9)

# Breast cancer
brcan_codes <- paste0("C50", 0:9)

# Other causes attributable to alcohol

# Diabetes
dm2_codes <- paste0("E11", 0:9)

# Tuberculosis
tb_codes <- paste0("A", sprintf("%02d", rep(15:19, each = 10)), 0:9)

# HIV/Aids
hiv_codes <- c(
  paste0("B20", 0:9), paste0("B21", 0:9), paste0("B22", 0:9), paste0("B23", 0:9), paste0("B24", 0:9)
)

# Lower respiratory infection
lri_codes <- c(paste0("J12", 0:9),paste0("J13", 0:9), paste0("J14", 0:9),
               paste0("J15", 0:9), paste0("J16", 0:9), paste0("J17", 0:9),
               paste0("J18", 0:9))
# Liver cirrhosis
lc_codes <- c(paste0("K70", 0:9),paste0("K74", 0:9))

#Acute pancreatitis
panc_codes <- c(paste0("K85", 0:9), "K861")


# MOTOR VEHICLE
ri_codes <-   c("V021", "V022", "V023", "V024", "V025", "V026", "V027", "V028", "V029",
"V031", "V032", "V033", "V034", "V035", "V036", "V037", "V038", "V039",
"V041", "V042", "V043", "V044", "V045", "V046", "V047", "V048", "V049",
"V092", "V093", 
"V123", "V124", "V125", "V126", "V127", "V128", "V129", 
"V133", "V134", "V135", "V136", "V137", "V138", "V139", 
"V143", "V144", "V145", "V146", "V147", "V148", "V149",
"V194", "V195", "V196",
"V203", "V204", "V205", "V206", "V207", "V208", "V209", 
"V213", "V214", "V215", "V216", "V217", "V218", "V219", 
"V223", "V224", "V225", "V226", "V227", "V228", "V229",
"V233", "V234", "V235", "V236", "V237", "V238", "V239",
"V243", "V244", "V245", "V246", "V247", "V248", "V249",
"V253", "V254", "V255", "V256", "V257", "V258", "V259",
"V263", "V264", "V265", "V266", "V267", "V268", "V269",
"V273", "V274", "V275", "V276", "V277", "V278", "V279",
"V283", "V284", "V285", "V286", "V287", "V288", "V289",
"V294", "V295", "V296", "V297", "V298", "V299",
"V304", "V305", "V306", "V307", "V308", "V309", 
"V314", "V315", "V316", "V317", "V318", "V319", 
"V324", "V325", "V326", "V327", "V328", "V329", 
"V334", "V335", "V336", "V337", "V338", "V339", 
"V344", "V345", "V346", "V347", "V348", "V349", 
"V354", "V355", "V356", "V357", "V358", "V359", 
"V364", "V365", "V366", "V367", "V368", "V369", 
"V374", "V375", "V376", "V377", "V378", "V379", 
"V384", "V385", "V386", "V387", "V388", "V389",
"V394", "V395", "V396", "V397", "V398", "V399",
"V404", "V405", "V406", "V407", "V408", "V409",
"V414", "V415", "V416", "V417", "V418", "V419",
"V424", "V425", "V426", "V427", "V428", "V429",
"V434", "V435", "V436", "V437", "V438", "V439",
"V444", "V445", "V446", "V447", "V448", "V449",
"V454", "V455", "V456", "V457", "V458", "V459",
"V464", "V465", "V466", "V467", "V468", "V469",
"V474", "V475", "V476", "V477", "V478", "V479",
"V484", "V485", "V486", "V487", "V488", "V489",
"V494", "V495", "V496", "V497", "V498", "V499",
"V504", "V505", "V506", "V507", "V508", "V509",
"V514", "V515", "V516", "V517", "V518", "V519",
"V524", "V525", "V526", "V527", "V528", "V529",
"V534", "V535", "V536", "V537", "V538", "V539",
"V544", "V545", "V546", "V547", "V548", "V549",
"V554", "V555", "V556", "V557", "V558", "V559",
"V564", "V565", "V566", "V567", "V568", "V569",
"V574", "V575", "V576", "V577", "V578", "V579",
"V584", "V585", "V586", "V587", "V588", "V589",
"V594", "V595", "V596", "V597", "V598", "V599",
"V604", "V605", "V606", "V607", "V608", "V609",
"V614", "V615", "V616", "V617", "V618", "V619",
"V624", "V625", "V626", "V627", "V628", "V629",
"V634", "V635", "V636", "V637", "V638", "V639",
"V644", "V645", "V646", "V647", "V648", "V649",
"V654", "V655", "V656", "V657", "V658", "V659",
"V664", "V665", "V666", "V667", "V668", "V669",
"V674", "V675", "V676", "V677", "V678", "V679",
"V684", "V685", "V686", "V687", "V688", "V689",
"V694", "V695", "V696", "V697", "V698", "V699",
"V704", "V705", "V706", "V707", "V708", "V709",
"V714", "V715", "V716", "V717", "V718", "V719",
"V724", "V725", "V726", "V727", "V728", "V729",
"V734", "V735", "V736", "V737", "V738", "V739",
"V744", "V745", "V746", "V747", "V748", "V749",
"V754", "V755", "V756", "V757", "V758", "V759",
"V764", "V765", "V766", "V767", "V768", "V769",
"V774", "V775", "V776", "V777", "V778", "V779",
"V784", "V785", "V786", "V787", "V788", "V789",
"V794", "V795", "V796", "V797", "V798", "V799",
"V803", "V804", "V805", "V811", "V821", 
"V830", "V831", "V832", "V833", 
"V840", "V841", "V842", "V843", 
"V850", "V851", "V852", "V853", 
"V860", "V861", "V862", "V863", 
"V870", "V871", "V872", "V873", "V874", "V875", "V876", "V877", "V878", "V892")

# UNINTENTIONAL INJURIES
unint_inj_codes <- c(
  paste0("W00", 0:9), # W000 to W009
  paste0("W01", 0:9), # W010 to W019
  paste0("W02", 0:9), # W020 to W029
  paste0("W03", 0:9), # W030 to W039
  paste0("W04", 0:9), # W040 to W049
  paste0("W05", 0:9), # W050 to W059
  paste0("W06", 0:9), # W060 to W069
  paste0("W07", 0:9), # W070 to W079
  paste0("W08", 0:9), # W080 to W089
  paste0("W09", 0:9), # W090 to W099
  paste0("W10", 0:9), # W100 to W109
  paste0("W11", 0:9), # W110 to W119
  paste0("W12", 0:9), # W120 to W129
  paste0("W13", 0:9), # W130 to W139
  paste0("W14", 0:9), # W140 to W149
  paste0("W15", 0:9), # W150 to W159
  paste0("W16", 0:9), # W160 to W169
  paste0("W17", 0:9), # W170 to W179
  paste0("W18", 0:9), # W180 to W189
  paste0("W19", 0:9),
  paste0("X00", 0:9), # X000 to X009
  paste0("X01", 0:9), # X010 to X019
  paste0("X02", 0:9), # X020 to X029
  paste0("X03", 0:9), # X030 to X039
  paste0("X04", 0:9), # X040 to X049
  paste0("X05", 0:9), # X050 to X059
  paste0("X06", 0:9), # X060 to X069
  paste0("X07", 0:9), # X070 to X079
  paste0("X08", 0:9), # X080 to X089
  paste0("X09", 0:9),  # X090 to X099
  paste0("X40", 0:9), # X400 to X409
  paste0("X46", 0:9), # X460 to X469
  paste0("X47", 0:9), # X470 to X479
  paste0("X48", 0:9), # X480 to X489
  paste0("X49", 0:9),  # X490 to X499
  paste0("W65", 0:9), # W650 to W659
  paste0("W66", 0:9), # W660 to W669
  paste0("W67", 0:9), # W670 to W679
  paste0("W68", 0:9), # W680 to W689
  paste0("W69", 0:9), # W690 to W699
  paste0("W70", 0:9), # W700 to W709
  paste0("W71", 0:9), # W710 to W719
  paste0("W72", 0:9), # W720 to W729
  paste0("W73", 0:9), # W730 to W739
  paste0("W74", 0:9),
  paste0("W20", 0:9), paste0("W21", 0:9), paste0("W22", 0:9), paste0("W23", 0:9), paste0("W24", 0:9), paste0("W25", 0:9),
  paste0("W26", 0:9), paste0("W27", 0:9), paste0("W28", 0:9), paste0("W29", 0:9), paste0("W30", 0:9),
  paste0("W31", 0:9), paste0("W32", 0:9), paste0("W33", 0:9), paste0("W34", 0:9), paste0("W35", 0:9),
  paste0("W36", 0:9), paste0("W37", 0:9), paste0("W38", 0:9), paste0("W39", 0:9), paste0("W40", 0:9),
  paste0("W41", 0:9), paste0("W42", 0:9), paste0("W43", 0:9), paste0("W44", 0:9), paste0("W45", 0:9),
  paste0("W46", 0:9), paste0("W47", 0:9), paste0("W48", 0:9), paste0("W49", 0:9), paste0("W50", 0:9),
  paste0("W51", 0:9), paste0("W52", 0:9), paste0("W53", 0:9), paste0("W54", 0:9), paste0("W55", 0:9),
  paste0("W56", 0:9), paste0("W57", 0:9), paste0("W58", 0:9), paste0("W59", 0:9), paste0("W60", 0:9),
  paste0("W61", 0:9), paste0("W62", 0:9), paste0("W63", 0:9), paste0("W64", 0:9),
  paste0("W75", 0:9), paste0("W76", 0:9), paste0("W77", 0:9), paste0("W78", 0:9), paste0("W79", 0:9),
  paste0("W80", 0:9), paste0("W81", 0:9), paste0("W82", 0:9), paste0("W83", 0:9), paste0("W84", 0:9),
  paste0("W85", 0:9), paste0("W86", 0:9), paste0("W87", 0:9), paste0("W88", 0:9), paste0("W89", 0:9),
  paste0("W90", 0:9), paste0("W91", 0:9), paste0("W92", 0:9), paste0("W93", 0:9), paste0("W94", 0:9),
  paste0("W95", 0:9), paste0("W96", 0:9), paste0("W97", 0:9), paste0("W98", 0:9), paste0("W99", 0:9),
  paste0("X10", 0:9), paste0("X11", 0:9), paste0("X12", 0:9), paste0("X13", 0:9), paste0("X14", 0:9),
  paste0("X15", 0:9), paste0("X16", 0:9), paste0("X17", 0:9), paste0("X18", 0:9), paste0("X19", 0:9),
  paste0("X20", 0:9), paste0("X21", 0:9), paste0("X22", 0:9), paste0("X23", 0:9), paste0("X24", 0:9),
  paste0("X25", 0:9), paste0("X26", 0:9), paste0("X27", 0:9), paste0("X28", 0:9), paste0("X29", 0:9),
  paste0("X30", 0:9), paste0("X31", 0:9), paste0("X32", 0:9), paste0("X33", 0:9), paste0("X34", 0:9),
  paste0("X35", 0:9), paste0("X36", 0:9), paste0("X37", 0:9), paste0("X38", 0:9), paste0("X39", 0:9),
  paste0("X50", 0:9), paste0("X51", 0:9), paste0("X52", 0:9), paste0("X53", 0:9), paste0("X54", 0:9),
  paste0("X55", 0:9), paste0("X56", 0:9), paste0("X57", 0:9), paste0("X58", 0:9), paste0("X59", 0:9),
  paste0("Y40", 0:9), paste0("Y41", 0:9), paste0("Y42", 0:9), paste0("Y43", 0:9), paste0("Y44", 0:9),
  paste0("Y45", 0:9), paste0("Y46", 0:9), paste0("Y47", 0:9), paste0("Y48", 0:9), paste0("Y49", 0:9),
  paste0("Y50", 0:9), paste0("Y51", 0:9), paste0("Y52", 0:9), paste0("Y53", 0:9), paste0("Y54", 0:9),
  paste0("Y55", 0:9), paste0("Y56", 0:9), paste0("Y57", 0:9), paste0("Y58", 0:9), paste0("Y59", 0:9),
  paste0("Y60", 0:9), paste0("Y61", 0:9), paste0("Y62", 0:9), paste0("Y63", 0:9), paste0("Y64", 0:9),
  paste0("Y65", 0:9), paste0("Y66", 0:9), paste0("Y67", 0:9), paste0("Y68", 0:9), paste0("Y69", 0:9),
  paste0("Y70", 0:9), paste0("Y71", 0:9), paste0("Y72", 0:9), paste0("Y73", 0:9), paste0("Y74", 0:9),
  paste0("Y75", 0:9), paste0("Y76", 0:9), paste0("Y77", 0:9), paste0("Y78", 0:9), paste0("Y79", 0:9),
  paste0("Y80", 0:9), paste0("Y81", 0:9), paste0("Y82", 0:9), paste0("Y83", 0:9), paste0("Y84", 0:9),
  paste0("Y85", 0:9), paste0("Y86", 0:9), paste0("Y88",0:9), paste0("Y89", 0:9)
  )

int_inj_codes <- c(paste0("Y35", 0:9), paste0("X85", 0:9), paste0("X86", 0:9), paste0("X87", 0:9), paste0("X88", 0:9), paste0("X89", 0:9),
                   paste0("X90", 0:9), paste0("X91", 0:9), paste0("X92", 0:9), paste0("X93", 0:9), paste0("X94", 0:9),
                   paste0("X95", 0:9), paste0("X96", 0:9), paste0("X97", 0:9), paste0("X98", 0:9), paste0("X99", 0:9),
                   paste0("Y00", 0:9), paste0("Y01", 0:9), paste0("Y02", 0:9), paste0("Y03", 0:9), paste0("Y04", 0:9),
                   paste0("Y05", 0:9), paste0("Y06", 0:9), paste0("Y07", 0:9), paste0("Y08", 0:9), paste0("Y09", 0:9),
                   "Y87.1",paste0("Y35", 0:9))


def <- def %>% 
  mutate(epi = if_else(DIAG1 %in% epilepsy_codes, 1, 0),
         ich = if_else(DIAG1 %in% ich_codes, 1, 0),
         is = if_else(DIAG1 %in% is_codes, 1, 0),
         hhd = if_else(DIAG1 %in% hhd_codes, 1, 0),
         bcan = if_else(DIAG1 %in% brcan_codes, 1, 0),
         crcan = if_else(DIAG1 %in% crcan_codes, 1, 0),
         lxcan = if_else(DIAG1 %in% lxcan_codes, 1, 0),
         lican = if_else(DIAG1 %in% lican_codes, 1, 0),
         oescan = if_else(DIAG1 %in% oescan_codes, 1, 0),
         opcan = if_else(DIAG1 %in% opcan_codes, 1, 0),
         dm2 = if_else(DIAG1 %in% dm2_codes, 1, 0),
         ihd = if_else(DIAG1 %in% ihd_codes, 1, 0),
         lri = if_else(DIAG1 %in% lri_codes, 1, 0),
         tb = if_else(DIAG1 %in% tb_codes, 1, 0),
         panc = if_else(DIAG1 %in% panc_codes, 1, 0),
         lc = if_else(DIAG1 %in% lc_codes, 1, 0),
         unint_inj = if_else(DIAG2 %in% unint_inj_codes | DIAG1 %in% unint_inj_codes, 1, 0),
         ri_inj = if_else(DIAG2 %in% ri_codes | DIAG2 %in% ri_codes, 1, 0),
         int_inj = if_else(DIAG1 %in% int_inj_codes|DIAG2 %in% int_inj_codes, 1, 0),
         hiv = if_else(DIAG1 %in% hiv_codes, 1, 0))


def <- def %>% 
  mutate(age_group = case_when(between(age, 15, 29)~1,
                                between(age, 30,44)~2,
                                between(age,45,59)~3,
                                age >= 60~4)) %>% 
  filter(age >= 15)

# Add gender column during transformation
aaf_fem_long <- aaf_fem %>%
  pivot_longer(
    cols = starts_with("AAF_ag"),
    names_to = "age_group",
    names_prefix = "AAF_ag",
    values_to = "point"
  ) %>%
  mutate(
    age_group = as.integer(age_group),
    lower = case_when(
      age_group == 1 ~ LL_ag1,
      age_group == 2 ~ LL_ag2,
      age_group == 3 ~ LL_ag3,
      age_group == 4 ~ LL_ag4
    ),
    upper = case_when(
      age_group == 1 ~ UL_ag1,
      age_group == 2 ~ UL_ag2,
      age_group == 3 ~ UL_ag3,
      age_group == 4 ~ UL_ag4
    ),
    gender = "Mujer" # Explicitly add gender column
  ) %>%
  dplyr::select(year = Year, age_group, gender, point, lower, upper, disease)

aaf_male_long <- aaf_male %>%
  pivot_longer(
    cols = starts_with("AAF_ag"),
    names_to = "age_group",
    names_prefix = "AAF_ag",
    values_to = "point"
  ) %>%
  mutate(
    age_group = as.integer(age_group),
    lower = case_when(
      age_group == 1 ~ LL_ag1,
      age_group == 2 ~ LL_ag2,
      age_group == 3 ~ LL_ag3,
      age_group == 4 ~ LL_ag4
    ),
    upper = case_when(
      age_group == 1 ~ UL_ag1,
      age_group == 2 ~ UL_ag2,
      age_group == 3 ~ UL_ag3,
      age_group == 4 ~ UL_ag4
    ),
    gender = "Hombre" # Explicitly add gender column
  ) %>%
  dplyr::select(year = Year, age_group, gender, point, lower, upper, disease)

# Combine male and female long-form data into one
aaf_long <- bind_rows(aaf_fem_long, aaf_male_long)

# Define all diseases with their filters and genders
disease_filters <- list(
  "Breast Cancer" = list(filter_col = "bcan", genders = c("Mujer")),
  "Liver Cancer" = list(filter_col = "lican", genders = c("Mujer", "Hombre")),
  "Larynx Cancer" = list(filter_col = "lxcan", genders = c("Mujer", "Hombre")),
  "Oesophagus Cancer" = list(filter_col = "oescan", genders = c("Mujer", "Hombre")),
  "Lip and Oral Cavity Cancer" = list(filter_col = "opcan", genders = c("Mujer", "Hombre")),
  "Colon and rectum Cancer" = list(filter_col = "crcan", genders = c("Mujer", "Hombre")),
  "Acute Pancreatitis" = list(filter_col = "panc", genders = c("Mujer", "Hombre")),
  "Epilepsy" = list(filter_col = "epi", genders = c("Mujer", "Hombre")),
  "DM2" = list(filter_col = "dm2", genders = c("Mujer", "Hombre")),
  "HIV" = list(filter_col = "hiv", genders = c("Mujer", "Hombre")),
  "Hypertensive Heart Disease" = list(filter_col = "hhd", genders = c("Mujer", "Hombre")),
  "Intracerebral Haemorrhage" = list(filter_col = "ich", genders = c("Mujer", "Hombre")),
  "Ischaemic Heart Disease" = list(filter_col = "ihd", genders = c("Mujer", "Hombre")),
  "Ischaemic Stroke" = list(filter_col = "is", genders = c("Mujer", "Hombre")),
  "Liver Cirrhosis" = list(filter_col = "lc", genders = c("Mujer", "Hombre")),
  "Lower Respiratory Infection" = list(filter_col = "lri", genders = c("Mujer", "Hombre")),
  "Road Injuries" = list(filter_col = "ri_inj", genders = c("Mujer", "Hombre")),
  "Tuberculosis" = list(filter_col = "tb", genders = c("Mujer", "Hombre")),
  "Unintentional Injuries" = list(filter_col = "unint_inj", genders = c("Mujer", "Hombre")),
  "Intentional Injuries" = list(filter_col = "int_inj", genders = c("Mujer", "Hombre"))
)

# Initialize an empty list to store the results for each disease
all_mortality_results <- list()

# Iterate over each disease in `disease_filters`
for (disease_name in names(disease_filters)) {
  # Extract filter column and genders for the current disease
  filter_info <- disease_filters[[disease_name]]
  filter_col <- filter_info$filter_col
  genders <- filter_info$genders
  
  # Loop through each gender specified for the disease
  for (gender in genders) {
    # Apply the filter and calculate mortality
    mortality_result <- def %>%
      group_by(year, gender, age_group) %>%
      count(.data[[filter_col]]) %>% 
      filter(.data[[filter_col]] == 1) %>% 
      left_join(
        aaf_long[aaf_long$disease == disease_name, ], 
        by = c("year", "age_group", "gender")
      ) %>%
      mutate(
        mort = point * n,
        ll_mort = lower * n,
        up_mort = upper * n,
        disease = disease_name
      ) %>%
      dplyr::select(year, age_group, gender, disease, mort, ll_mort, up_mort)
    
    # Append results to the list
    all_mortality_results[[paste(disease_name, gender, sep = "_")]] <- mortality_result
  }
}

# Combine all results into a single data frame
mortality_results <- bind_rows(all_mortality_results) %>% 
  filter(!is.na(mort))



rio::export(mortality_results, "Mortality Estimates.xlsx")

#---------#
# RESULTS #
#---------#


# STANDARD POPULATION

spw_male <- rio::import("ine_proyecciones.xlsx", sheet = 1) %>% #standard population weights
  mutate(edad = as.numeric(ifelse(edad == "100+", 100, edad)),
         age_group = case_when(between(edad, 15, 29)~1,
                               between(edad, 30,44)~2,
                               between(edad,45,59)~3,
                               edad >= 60 ~4,
                               edad < 15 ~ 0)
  ) %>%
  pivot_longer(cols = starts_with("ano_"), names_to = "year", values_to = "value") %>%
  mutate(year = as.numeric(gsub("ano_", "", year))) %>% 
  group_by(year, age_group) %>% 
  summarise(tot = sum(value)) %>% 
  ungroup() %>%
  group_by(year) %>% 
  mutate(pop = sum(tot),
         spw = tot/pop,
         gender = "Hombre")

# Extract 2018 weights
spw18_male <- spw_male %>% 
  filter(year == 2018) %>% 
  ungroup() %>% 
  select(age_group, spw)  # Remove the year column to make weights fixed for all years

# Repeat 2018 weights for all years
spw_male <- spw_male %>%
  dplyr::select(-spw) %>% 
  left_join(spw18_male, by = "age_group")%>% 
  filter(age_group > 0)

# female
spw_fem <- rio::import("ine_proyecciones.xlsx", sheet = 2) %>% #standard population weights
  mutate(edad = as.numeric(ifelse(edad == "100+", 100, edad)),
         age_group = case_when(between(edad, 15, 29)~1,
                               between(edad, 30,44)~2,
                               between(edad,45,59)~3,
                               edad >= 60 ~4,
                               edad < 15 ~ 0)
  ) %>%
  pivot_longer(cols = starts_with("ano_"), names_to = "year", values_to = "value") %>%
  mutate(year = as.numeric(gsub("ano_", "", year))) %>% 
  group_by(year, age_group) %>% 
  summarise(tot = sum(value)) %>% 
  ungroup() %>%
  group_by(year) %>% 
  mutate(pop = sum(tot),
         spw = tot/pop,
         gender = "Mujer") 

# Extract 2018 weights
spw18_fem <- spw_fem %>% 
  filter(year == 2018) %>% 
  ungroup() %>% 
  select(age_group, spw)  # Remove the year column to make weights fixed for all years

# Repeat 2018 weights for all years
spw_fem <- spw_fem %>% 
  dplyr::select(-spw) %>% 
  left_join(spw18_fem, by = "age_group") %>% 
  filter(age_group > 0)

# both sexes
spw_tot <- spw_male %>% 
  bind_rows(spw_fem) %>% 
  dplyr::select(-spw) %>% 
  group_by(year,age_group) %>% 
  summarise(tot = sum(tot)) %>% 
  mutate(pop =sum(tot),
         spw = tot/pop) %>% 
  filter(age_group > 0)

# AGE-STANDARDIZE MORTALITY (FIGURE 1)
results <- mortality_results %>% 
  group_by(year, age_group) %>% 
  summarise(
    mort = sum(mort, na.rm = TRUE),
    ll_mort = sum(ll_mort, na.rm = TRUE),
    up_mort = sum(up_mort, na.rm = TRUE),
  ) %>% 
  left_join(spw_tot, by = c("year", "age_group")) %>% 
  mutate(
    mort_rate = (mort / tot) * 100000,
    ll_mort_rate = (ll_mort / tot) * 100000,
    up_mort_rate = (up_mort / tot) * 100000
  ) %>%
  summarise(
    year = first(year),  # Ensure year is included
    crude_mort_rate = sum(mort_rate),
    std_mort_rate = sum(mort_rate * spw),
    ll_std_mort_rate = sum(ll_mort_rate * spw),
    up_std_mort_rate = sum(up_mort_rate * spw),
    .groups = "drop"
  )

results_fem <- mortality_results %>% 
  group_by(year, age_group,gender) %>% 
  summarise(
    mort = sum(mort, na.rm = TRUE),
    ll_mort = sum(ll_mort, na.rm = TRUE),
    up_mort = sum(up_mort, na.rm = TRUE),
  ) %>%
  filter(gender == "Mujer") %>% 
  left_join(spw_fem, by = c("year", "age_group")) %>% 
  mutate(
    mort_rate = (mort / tot) * 100000,
    ll_mort_rate = (ll_mort / tot) * 100000,
    up_mort_rate = (up_mort / tot) * 100000
  ) %>%
  group_by(year) %>% 
  summarise(
    year = first(year),  # Ensure year is included
    crude_mort_rate = sum(mort_rate),
    std_mort_rate = sum(mort_rate * spw),
    ll_std_mort_rate = sum(ll_mort_rate * spw),
    up_std_mort_rate = sum(up_mort_rate * spw),
    .groups = "drop"
  )

results_male <- mortality_results %>% 
  filter(gender == "Hombre") %>%  # Filtrar hombres primero
  group_by(year, age_group) %>%  # Asegurar agrupación por año y grupo de edad
  summarise(
    mort = sum(mort, na.rm = TRUE),
    ll_mort = sum(ll_mort, na.rm = TRUE),
    up_mort = sum(up_mort, na.rm = TRUE),
    .groups = "drop"  # Eliminar agrupación después del resumen
  ) %>% 
  left_join(spw_male, by = c("year", "age_group")) %>%  # Unir con los pesos estándar
  mutate(
    mort_rate = (mort / tot) * 100000,  # Tasa de mortalidad cruda por grupo de edad
    ll_mort_rate = (ll_mort / tot) * 100000,
    up_mort_rate = (up_mort / tot) * 100000
  ) %>% 
  group_by(year) %>%  # Agrupar por año para calcular las tasas anuales
  summarise(
    crude_mort_rate = sum(mort_rate),  # Tasa cruda anual
    std_mort_rate = sum(mort_rate * spw),  # Tasa estandarizada anual
    ll_std_mort_rate = sum(ll_mort_rate * spw),  # Límite inferior de la tasa
    up_std_mort_rate = sum(up_mort_rate * spw),  # Límite superior de la tasa
    .groups = "drop"  # Eliminar agrupación para un resultado plano
  )

results_gender <- bind_rows(
  results_fem %>% mutate(gender = "Female"),
  results_male %>% mutate(gender = "Male")
)

# Add a "Total" category to results for Figure 1
results <- results %>%
  mutate(gender = "Total")

# Combine gender-specific data with total data
combined_results <- bind_rows(
  results_gender,  # Gender-specific results
  results          # Total results
)

#----------#
# FIGURE 1 #
#----------#

# Plot the unified figure
ggplot(combined_results, aes(x = year, y = std_mort_rate, linetype = gender, fill = gender)) +
  # Confidence interval ribbons
  geom_ribbon(
    aes(ymin = ll_std_mort_rate, ymax = up_std_mort_rate, fill = gender),
    alpha = 0.3, color = NA
  ) +
  # Lines for age-standardized mortality
  geom_line(
    aes(linetype = gender),
    size = 1.2, color = "black"
  ) +
  # Points for age-standardized mortality
  geom_point(
    aes(shape = gender),
    size = 2, color = "black", fill = "white"
  ) +
  # Customize scales and labels
  scale_x_continuous(
    breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
    labels = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022)
  ) +
  scale_y_continuous(breaks = seq(40, 175, 15)) +
  coord_cartesian(ylim = c(30, 190)) +
  labs(
    x = "Year",
    y = "Age-Standardized Mortality Rate x 100,000 Persons",
    linetype = "",
    fill = "",
    shape = ""
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.position = "top",
    legend.key = element_blank()  # Remove background in legend
  ) +
  # Customize linetypes, fill colors, and shapes for Total and genders
  scale_linetype_manual(
    values = c("Male" = "dotted", "Female" = "dashed", "Total" = "solid"),
    labels = c("Male" = "Men", "Female" = "Women", "Total" = "Total")
  ) +
  scale_fill_manual(
    values = c("Male" = "gray70", "Female" = "gray70", "Total" = "gray40"),
    labels = c("Male" = "Men", "Female" = "Women", "Total" = "Total")
  ) +
  scale_shape_manual(
    values = c("Male" = 21, "Female" = 22, "Total" = 23),
    labels = c("Male" = "Men", "Female" = "Women", "Total" = "Total")
  )
# FIGURE 2
# Prepare gender-specific and total data
gender_data <- mortality_results %>% 
  group_by(year, gender) %>% 
  summarise(
    mort = sum(mort),
    ll_mort = sum(ll_mort),
    up_mort = sum(up_mort),
    .groups = "drop"
  ) %>% 
  left_join(death_sex, by = c("year", "gender")) %>% 
  mutate(
    prop = mort / n,
    ll_prop = ll_mort / n,
    up_prop = up_mort / n
  )

total_data <- mortality_results %>% 
  group_by(year) %>% 
  summarise(
    mort = sum(mort),
    ll_mort = sum(ll_mort),
    up_mort = sum(up_mort),
    .groups = "drop"
  ) %>% 
  left_join(tot_death, by = "year") %>% 
  mutate(
    prop = mort / n,
    ll_prop = ll_mort / n,
    up_prop = up_mort / n
  ) %>% 
  mutate(gender = "Total")
# Combine gender-specific data with total data
combined_data <- bind_rows(gender_data, total_data)

# Plot the combined data
ggplot(combined_data) +
  # Add CI ribbon with distinct gray scales
  geom_ribbon(
    aes(x = year, ymin = ll_prop, ymax = up_prop, fill = gender),
    alpha = 0.3, color = NA
  ) +
  # Add line for proportions
  geom_line(
    aes(x = year, y = prop, linetype = gender),
    color = "black", size = 1.2
  ) +
  # Add points for proportions
  geom_point(
    aes(x = year, y = prop, shape = gender),
    color = "black", size = 2
  ) +
  # Scales and labels
  scale_x_continuous(
    breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
    labels = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022)
  ) +
  scale_y_continuous(breaks = seq(0, 0.3, 0.1)) +
  coord_cartesian(ylim = c(0, 0.3)) +
  labs(
    x = "Year",
    y = "Alcohol-Attributable Deaths / Total Deaths",
    linetype = "",
    shape = "",
    fill = ""
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.position = "top"
  ) +
  scale_linetype_manual(
    values = c("Hombre" = "dotted", "Mujer" = "dashed", "Total" = "solid"),
    labels = c("Hombre" = "Men", "Mujer" = "Women", "Total" = "Total")
  ) +
  scale_fill_manual(
    values = c("Hombre" = "gray70", "Mujer" = "gray70", "Total" = "gray40"),  # Darker for Men/Women, lighter for Total
    labels = c("Hombre" = "Men", "Mujer" = "Women", "Total" = "Total")
  ) +
  scale_shape_manual(
    values = c("Hombre" = 21, "Mujer" = 22, "Total" = 23),
    labels = c("Hombre" = "Men", "Mujer" = "Women", "Total" = "Total")
  )

# AGE GROUP COMPARISON (FIGURE 3)

mortality_results %>% 
  group_by(year, gender, age_group) %>% 
  summarise(mort = sum(mort, na.rm = T),
            ll_mort = sum(ll_mort, na.rm = T),
            up_mort = sum(up_mort, na.rm = T)) %>% 
  left_join(spw_tot, by = c("year","age_group")) %>% 
  mutate(
    mort_rate = (mort / tot) * 100000,  # Tasa de mortalidad cruda por grupo de edad
    ll_mort_rate = (ll_mort / tot) * 100000,
    up_mort_rate = (up_mort / tot) * 100000
  ) %>% 
  ggplot(., aes(x = year, y = mort_rate, linetype = factor(age_group))) +
  geom_line(size = 1, color = "black") +
  geom_point(size = 2, color = "black") +
  scale_x_continuous(
    breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
    labels = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022)
  ) +
  scale_linetype_manual(
    values = c("1" = "solid", "2" = "dashed", "3" = "dotted", "4" = "dotdash"),
    labels = c("1" = "15-29", "2" = "30-44", "3" = "45-59", "4" = "60+")
  ) +
  facet_wrap(~ gender, labeller = as_labeller(c("Mujer" = "Women", "Hombre" = "Men"))) +
  labs(
    x = "Year",
    y = "Mortality rate x 100,000 persons",
    linetype = "Age Group"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.position = "top"
  )

# CAUSE COMPARAISON (FIGURE 4 AND 5)

mortality_results <- mortality_results %>%
  mutate(category = case_when(
    disease %in% c(
      "Breast Cancer", "Colon and rectum Cancer", "Larynx Cancer",
      "Lip and Oral Cavity Cancer", "Liver Cancer", "Oesophagus Cancer"
    ) ~ "Cancer",
    disease %in% c(
      "Intentional Injuries", "Road Injuries", "Unintentional Injuries"
    ) ~ "Injuries",
    disease %in% c(
      "Hypertensive Heart Disease", "Intracerebral Haemorrhage", 
      "Ischaemic Heart Disease", "Ischaemic Stroke"
    ) ~ "Cardiovascular",
    disease %in% c(
      "DM2", "Liver Cirrhosis", "Lower Respiratory Infection", 
      "Tuberculosis", "Acute Pancreatitis", "HIV"
    ) ~ "Other Causes",
    disease %in% c(
      "Epilepsy"
    ) ~ "Neuropsychiatric"))



# FIGURE 4
mortality_results %>%
  filter(gender == "Hombre") %>%
  group_by(year, age_group, category) %>%
  summarise(
    mort = sum(mort, na.rm = TRUE),
    ll_mort = sum(ll_mort, na.rm = TRUE),
    up_mort = sum(up_mort, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(age_group, year) %>%
  mutate(prop_mort = mort / sum(mort, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = year, y = prop_mort, linetype = category, shape = category, col = category)) +  # Map shape to category
  geom_line(size = 1, color = "black") +
  geom_point(size = 2, color = "black") +
  scale_x_continuous(
    breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
    labels = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022)
  ) +
  scale_color_manual(
    values = c(
      "Cancer" = "black", 
      "Cardiovascular" = "darkgray", 
      "Injuries" = "gray", 
      "Neuropsychiatric" = "lightgray", 
      "Other Causes" = "dimgray"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Cancer" = 21,  # Circle
      "Cardiovascular" = 22,  # Square
      "Injuries" = 23,  # Diamond
      "Neuropsychiatric" = 24,  # Triangle-up
      "Other Causes" = 25  # Triangle-down
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Cancer" = "solid", 
      "Cardiovascular" = "dashed", 
      "Injuries" = "dotted", 
      "Neuropsychiatric" = "dotdash", 
      "Other Causes" = "longdash"
    )
  ) +
  facet_wrap(~ age_group, labeller = as_labeller(c("1" = "15-29", "2" = "30-44", "3" = "45-59", "4" = "60+"))) +
  labs(
    x = "Year",
    y = "Proportion of total alcohol-attributable deaths by cause",
    linetype = "Category",
    shape = "Category"  # Add shape legend title
  ) +
  theme_minimal()

# FIGURE 5
mortality_results %>%
  filter(gender == "Mujer") %>%
  group_by(year, age_group, category) %>%
  summarise(
    mort = sum(mort, na.rm = TRUE),
    ll_mort = sum(ll_mort, na.rm = TRUE),
    up_mort = sum(up_mort, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(age_group, year) %>%
  mutate(prop_mort = mort / sum(mort, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = year, y = prop_mort, linetype = category, shape = category, col = category)) +  # Map shape to category
  geom_line(size = 1, color = "black") +
  geom_point(size = 2, color = "black") +
  scale_x_continuous(
    breaks = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022),
    labels = c(2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022)
  ) +
  scale_color_manual(
    values = c(
      "Cancer" = "black", 
      "Cardiovascular" = "darkgray", 
      "Injuries" = "gray", 
      "Neuropsychiatric" = "lightgray", 
      "Other Causes" = "dimgray"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Cancer" = 21,  # Circle
      "Cardiovascular" = 22,  # Square
      "Injuries" = 23,  # Diamond
      "Neuropsychiatric" = 24,  # Triangle-up
      "Other Causes" = 25  # Triangle-down
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Cancer" = "solid", 
      "Cardiovascular" = "dashed", 
      "Injuries" = "dotted", 
      "Neuropsychiatric" = "dotdash", 
      "Other Causes" = "longdash"
    )
  ) +
  facet_wrap(~ age_group, labeller = as_labeller(c("1" = "15-29", "2" = "30-44", "3" = "45-59", "4" = "60+"))) +
  labs(
    x = "Year",
    y = "Proportion of total alcohol-attributable deaths by cause",
    linetype = "Category",
    shape = "Category"  # Add shape legend title
  ) +
  theme_minimal()


# TABLES

table_mort_fem <- mortality_results %>%
  # Filter for "Mujer"
  filter(gender == "Mujer") %>% 
  # Group by year and disease (exclude age_group)
  group_by(year, disease) %>%
  # Summarize values across age groups
  summarise(
    mort = sum(round(mort,0), na.rm = TRUE),
    ci = paste0("[", sum(round(ll_mort,0), na.rm = TRUE), "-", sum(round(up_mort,0), na.rm = TRUE), "]"),
    .groups = "drop"
  ) %>%
  # Pivot years into columns
  pivot_wider(
    names_from = year,
    values_from = c(mort, ci)
  ) 


rio::export(table_mort_fem, "Table Mortality Women.xlsx")

table_mort_male <- mortality_results %>%
  # Filter for "Mujer"
  filter(gender == "Hombre") %>% 
  # Group by year and disease (exclude age_group)
  group_by(year, disease) %>%
  # Summarize values across age groups
  summarise(
    mort = sum(round(mort,0), na.rm = TRUE),
    ci = paste0("[", sum(round(ll_mort,0), na.rm = TRUE), "-", sum(round(up_mort,0), na.rm = TRUE), "]"),
    .groups = "drop"
  ) %>%
  # Pivot years into columns
  pivot_wider(
    names_from = year,
    values_from = c(mort, ci)
  ) 


rio::export(table_mort_male, "Table Mortality Male.xlsx")


