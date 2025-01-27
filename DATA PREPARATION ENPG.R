#########################################################
# CLEANING AND MERGING NATIONAL DRUG AND ALCOHOL SURVEY #
#########################################################

# LIBRARIES
library(dplyr)
library(readr)
library(haven)
#############
# ENPG 2008 #
#############

enpg08 <- read_rds("https://github.com/ACC1240138/FONDECYT-REGULAR-/raw/main/rawdata/enpg2008.RDS")
enpg08 <- enpg08 %>% mutate(year = 2008)
data08 <- enpg08 %>% 
  dplyr::select(id, year, region,comuna, exp, sexo, edad,  
         religion = q277, nedu = q279, ecivil = q281, 
         p_originario = q282, ingreso = q284,oh1 = q12, oh2 = q15, oh3 = q16,
         tab1 = q5,tab2 = q8, tab3 = q9, tab4 = q10, mar1 = q33, mar2 = q36,
         coc1 = q79, coc2 = q82, t1 = q133a, t2 = q133b, t3 = q133c, t4= q133d, 
         t5 = q133e, t6 = q133f, t7 = q133g, t8 = q133h, t9 = q133i, t10 = q133j, 
         t1_m = q134, t2_m = q135, t3_m = q136, t4_m = q137, t5_m = q138, t6_m = q139, 
         t7_m = q140, t8_m = q141, t9_m = q142, t10_m = q143,
         audit1 = q18, audit2 = q19, audit3 = q20) %>% 
  mutate(id = factor(id),
         mar1 = factor(mar1, levels = c(1,2), labels = c("Si","No")),
         mar2 = factor(mar2, levels = c (1,2,3), labels = c("<30",">30",">1 año")),
         coc1 = factor(coc1, levels = c(1,2), labels = c("Si","No")),
         coc2 = factor(coc2, levels = c (1,2,3), labels = c("<30",">30",">1 año")),
         tab1 = factor(tab1, levels = c(1,2), labels = c("Si","No")),
         tab2 = factor(tab2, levels = c (1,2,3), labels = c("<30",">30",">1 año")),
         sexo = factor(sexo, levels = c(1,2), labels = c("Hombre","Mujer")),
         oh1 = factor(oh1, levels = c(1,2), labels = c("Si","No")),
         oh2 = factor(oh2, levels = c (1,2,3), labels = c("<30",">30",">1 año")),
         audit1 = factor(audit1, levels = c (0:4), labels = c("Nunca","1 vez al mes o menos"," 2 a 4 veces al mes",
                                                              "2 a 3 veces a la semana","4 o mas veces a la semana")),
         audit2 = factor(audit2, levels = c (0:4), labels = c("0-2","3-4","5-6",
                                                              "7-8", "9 o mas")),
         audit3 = factor(audit3, levels = c (0:4), labels = c("Nunca","menos de 1 vez al mes",
                                                              "mensualmente",
                                                              "semanalmente",
                                                              "todos o casi todos los dias")),
         ecivil = factor(case_when(ecivil <= 2 ~ "soltero",
                                   ecivil == 3 ~ "casado",
                                   ecivil == 4 | ecivil == 5 ~ "separado",
                                   ecivil == 6 | ecivil == 7 ~ "viudo")),
         religion = factor(religion, levels = c(1:4), labels = c("Catolico",
                                                                 "evangelico/protestante",
                                                                 "otra",
                                                                 "ninguna")),
         p_originario = factor(ifelse(p_originario == 9,"no pertenece","pertenece")),
         nedu = factor(case_when(nedu == 1 ~ "basica incompleta",
                                 nedu == 2 ~ "basica completa",
                                 nedu == 3 ~ "media incompleta",
                                 nedu == 4 ~ "media completa",
                                 nedu == 5 | nedu == 7 ~ "superior incompleta",
                                 nedu == 6 | nedu == 8 ~ "superior completa",
                                 TRUE~NA)),
         ingreso = na_if(ingreso, 9),
         ingreso = na_if(ingreso, 0)) %>% 
  rowwise() %>%
  mutate(tranq_vida = if_else(any(c_across(t1:t10) == 1), 1, 0),
         tranq_mes = if_else(any(c_across(t1_m:t10_m) %in% c(1, 2)), 1, 0)) %>%
  ungroup() %>% 
  dplyr::select(-matches("^t[1-9]$|^t10$"), -matches("^t[1-9]_m$|^t10_m$"))

#############
# ENPG 2010 #
#############

enpg10 <- read_rds("https://github.com/ACC1240138/FONDECYT-REGULAR-/raw/main/rawdata/enpg2010.RDS")
data10 <- enpg10 %>% mutate(year = 2010) %>% 
  dplyr::select(id = folio, year,region = pregion, comuna = pcodcom,tab1 = p005, tab2 = p008,
         tab3 = p009, tab4 = p010, mar1 = p033, mar2 = p036, coc1 = p077, coc2 = p080, 
         t1 = p169, t2 = p171, t3 = p173, t4 = p175, t5 = p177,  t6 = p179, t7 = p181,
         t8 = p183, t9 = p185, t1_m = p170, t2_m = p172, t3_m = p174, t4_m = p176,  
         t5_m = p178, t6_m = p180, t7_m = p182,  t8_m = p184,  t9_m = p186, 
         exp = factor_ajustado_com, sexo = psexoent,edad = pedadent, nedu = p288a, 
         religion = p285,ecivil = p282, p_originario = p283, ingreso = p291,
         oh1 = p013, oh2 = p016, oh3 = p017, audit1 = p019, audit2 = p020, audit3 = p021) %>% 
  mutate(id = factor(id),
         sexo = factor(sexo, levels = c(1,2), labels = c("Hombre","Mujer")),
         mar1 = factor(mar1, levels = c(1,2), labels = c("Si","No")),
         mar2 = factor(mar2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
         tab1 = factor(tab1, levels = c(1,2), labels = c("Si","No")),
         tab2 = factor(tab2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
         coc1 = factor(coc1, levels = c(1,2), labels = c("Si","No")),
         coc2 = factor(coc2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
         oh1 = factor(oh1, levels = c(1,2), labels = c("Si","No")),
         oh2 = factor(oh2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
         audit1 = factor(audit1, levels = c (0:4), labels = c("Nunca","1 vez al mes o menos"," 2 a 4 veces al mes",
                                                              "2 a 3 veces a la semana","4 o mas veces a la semana")),
         audit2 = factor(audit2, levels = c (0:4), labels = c("0-2","3-4","5-6",
                                                              "7-8", "9 o mas")),
         audit3 = factor(audit3, levels = c (0:4), labels = c("Nunca","menos de 1 vez al mes",
                                                              "mensualmente",
                                                              "semanalmente",
                                                              "todos o casi todos los dias")),
         ecivil = factor(case_when(ecivil == 0 ~ NA,
                                   ecivil <= 2 ~ "soltero",
                                   ecivil == 3 ~ "casado",
                                   ecivil == 4 | ecivil == 5 ~ "separado",
                                   ecivil == 6 | ecivil == 7 ~ "viudo")),
         religion = factor(religion, levels = c(1:4), labels = c("Catolico",
                                                                 "evangelico/protestante",
                                                                 "otra",
                                                                 "ninguna")),
         p_originario = factor(ifelse(p_originario == 9,"no pertenece","pertenece")),
         nedu = factor(case_when(nedu < 2 | nedu == 16 ~ "basica incompleta",
                                 nedu > 2 & nedu < 5 ~ "basica completa",
                                 nedu >= 5 & nedu <= 8 ~ "media incompleta",
                                 nedu == 4 ~ "media completa",
                                 nedu == 9 | nedu == 11 | nedu == 13 ~ "superior incompleta",
                                 nedu == 10 | nedu == 12 | nedu == 14 | nedu == 15 ~ "superior completa",
                                 TRUE~NA)),
         ingreso = na_if(ingreso, 9)) %>% 
  rowwise() %>%
  mutate(
    tranq_vida = if_else(any(c_across(t1:t9) == 1), 1, 0),
    tranq_mes = if_else(any(c_across(t1_m:t9_m) %in% c(1, 2)), 1, 0)
  ) %>%
  ungroup() %>%
  dplyr::select(-t1:-t9, -t1_m:-t9_m)


#############
# ENPG 2012 #
#############

enpg12 <- read_rds("https://github.com/ACC1240138/FONDECYT-REGULAR-/raw/main/rawdata/enpg2012.RDS")
data12 <- enpg12 %>% mutate(year = 2012) %>% 
  dplyr::select(id = idencuesta, year, region = "región",comuna = "código_comuna",
         exp = PONDERADOR, sexo, edad, nedu1 = p184_1, nedu2 = p185_1, 
         religion = p180, ecivil = p177, p_originario = p179, ingreso = p188,
         tab1 = p2, tab2 = p5, tab3 = p6, tab4 = p7,mar1 = p35, mar2 = p38,
         coc1 = p80, coc2 = p83, t1 = p118_1, t2 = p118_2, t3 = p118_3, t4 = p118_4,
         t5 = p118_5, t6 = p118_6, t7 = p118_7, t8 = p118_8, t9 = p118_9, t1_m = p119_1,
         t2_m = p119_2, t3_m = p119_3, t4_m = p119_4, t5_m = p119_5, t6_m = p119_6,
         t7_m = p119_7, t8_m = p119_8, t9_m = p119_9,
         oh1 = p10, oh2 = p13,oh3 = p14, audit1 = p19, audit2 = p20, audit3 = p21) %>% 
mutate(id = factor(id),
       sexo = factor(sexo, levels = c(1,2), labels = c("Hombre","Mujer")),
       mar1 = factor(mar1, levels = c(1,2), labels = c("Si","No")),
       mar2 = factor(mar2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
       tab1 = factor(tab1, levels = c(1,2), labels = c("Si","No")),
       tab2 = factor(tab2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
       coc1 = factor(coc1, levels = c(1,2), labels = c("Si","No")),
       coc2 = factor(coc2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
       oh1 = factor(oh1, levels = c(1,2), labels = c("Si","No")),
       oh2 = factor(oh2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
       audit1 = factor(audit1, levels = c (0:4), labels = c("Nunca","1 vez al mes o menos"," 2 a 4 veces al mes",
                                                            "2 a 3 veces a la semana","4 o mas veces a la semana")),
       audit2 = factor(audit2, levels = c (0:4), labels = c("0-2","3-4","5-6",
                                                            "7-8", "9 o mas")),
       audit3 = factor(audit3, levels = c (0:4), labels = c("Nunca","menos de 1 vez al mes",
                                                            "mensualmente",
                                                            "semanalmente",
                                                            "todos o casi todos los dias")),
       ecivil = factor(case_when(
                                 ecivil == 1 ~ "soltero",
                                 ecivil == 2 ~ "casado",
                                 ecivil == 4  ~ "viudo",
                                 ecivil == 3 | ecivil == 5 ~ "separado")),
       religion = factor(case_when(religion == 1 ~ "Catolico",
                                   religion == 2 ~ "evangelico/protestante",
                                   religion >= 3 & religion <= 11 ~ "otra",
                                   religion == 12 ~ "ninguna")),
       p_originario = factor(ifelse(p_originario == 10,"no pertenece","pertenece")),
       nedu = factor(case_when(nedu1 < 4 & nedu2 == 2 ~ "basica incompleta",
                               nedu1 >= 2 & nedu1 < 5 & nedu2 == 1~ "basica completa",
                               nedu1 >= 5 & nedu1 <= 8 & nedu2 == 2 ~ "media completa",
                               nedu1 >= 5 & nedu1 <= 8 & nedu2 == 1  ~ "media completa",
                               nedu1 >= 9 & nedu1 <= 10 & nedu2 == 2 ~ "superior incompleta",
                               nedu1 >= 9 & nedu1 <= 10 & nedu2 == 1 ~ "superior completa",
                               nedu1 >= 11 & nedu1 <= 13 ~ "superior completa",
                               TRUE~NA)),
       ingreso = na_if(ingreso, 88),
       ingreso = na_if(ingreso, 99)) %>% 
  dplyr::select(-nedu1, -nedu2) %>% 
rowwise() %>%
  mutate(
    tranq_vida = if_else(any(c_across(t1:t9) == 1), 1, 0),
    tranq_mes = if_else(any(c_across(t1_m:t9_m) %in% c(1, 2)), 1, 0)
  ) %>%
  ungroup() %>%
  dplyr::select(-t1:-t9, -t1_m:-t9_m)

#############
# ENPG 2014 #
#############

enpg14 <- read_rds("https://github.com/ACC1240138/FONDECYT-REGULAR-/raw/main/rawdata/enpg2014.RDS")
data14 <- enpg14 %>% mutate(year = 2014) %>% 
  dplyr::select(id = idencuesta, year, region = Region, comuna = Comuna, exp = RND_F2_MAY_AJUS_com, sexo, edad,
         nedu1 = dp9_1, nedu2 = dp10_1, religion = dp5, ecivil = dp2, 
         p_originario = dp4, ingreso = dp13, oh1,oh2 = oh4, oh3 = oh5, mar1, mar2 = mar4,coc1 , coc2 = coc4, 
         tab1 = st3, tab2 = st6, tab3 = st7, tab4 = st8, t1 = trans1_1, t2 = trans1_2, 
         t3 = trans1_3,t4 = trans1_4, t5 = trans1_5, t6 = trans1_6, t7 = trans1_7,
         t8 = trans1_8, t9 = trans1_9, t1_m = trans3_1, t2_m = trans3_2, t3_m = trans3_3,
         t4_m = trans3_4, t5_m = trans3_5, t6_m = trans3_6, t7_m = trans3_7, t8_m = trans3_8,
         t9_m = trans3_9, audit1 = oh10, audit2 = oh11, audit3 = oh12) %>% 
mutate(id = factor(id),
      sexo = factor(sexo, levels = c(1,2), labels = c("Hombre","Mujer")),
      mar1 = factor(mar1, levels = c(1,2), labels = c("Si","No")),
      mar2 = factor(mar2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
      tab1 = factor(tab1, levels = c(1,2), labels = c("Si","No")),
      tab2 = factor(tab2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
      coc1 = factor(coc1, levels = c(1,2), labels = c("Si","No")),
      coc2 = factor(coc2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
       oh1 = factor(oh1, levels = c(1,2), labels = c("Si","No")),
       oh2 = factor(oh2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
       audit1 = factor(audit1, levels = c (0:4), labels = c("Nunca","1 vez al mes o menos"," 2 a 4 veces al mes",
                                                            "2 a 3 veces a la semana","4 o mas veces a la semana")),
       audit2 = factor(audit2, levels = c (0:4), labels = c("0-2","3-4","5-6",
                                                            "7-8", "9 o mas")),
       audit3 = factor(audit3, levels = c (0:4), labels = c("Nunca","menos de 1 vez al mes",
                                                            "mensualmente",
                                                            "semanalmente",
                                                            "todos o casi todos los dias")),
       ecivil = factor(case_when(
         ecivil == 1 ~ "soltero",
         ecivil == 2 ~ "casado",
         ecivil == 4  ~ "viudo",
         ecivil == 3 | ecivil == 5 ~ "separado")),
       religion = factor(case_when(religion == 1 ~ "Catolico",
                                   religion == 2 ~ "evangelico/protestante",
                                   religion >= 3 & religion <= 11 ~ "otra",
                                   religion == 12 ~ "ninguna")),
       p_originario = factor(ifelse(p_originario == 10,"no pertenece","pertenece")),
       nedu = factor(case_when(nedu1 < 4 & nedu2 == 2 ~ "basica incompleta",
                               nedu1 >= 2 & nedu1 < 5 & nedu2 == 1~ "basica completa",
                               nedu1 >= 5 & nedu1 <= 8 & nedu2 == 2 ~ "media completa",
                               nedu1 >= 5 & nedu1 <= 8 & nedu2 == 1  ~ "media completa",
                               nedu1 >= 9 & nedu1 <= 10 & nedu2 == 2 ~ "superior incompleta",
                               nedu1 >= 9 & nedu1 <= 10 & nedu2 == 1 ~ "superior completa",
                               nedu1 >= 11 & nedu1 <= 13 ~ "superior completa",
                               TRUE~NA)),
       ingreso = na_if(ingreso, 88),
       ingreso = na_if(ingreso, 99)) %>% 
  dplyr::select(-nedu1, -nedu2) %>% 
  rowwise() %>%
  mutate(
    tranq_vida = if_else(any(!is.na(c_across(t1:t9))), 1, 0),
    tranq_mes = if_else(any(c_across(t1_m:t9_m) %in% c(1, 2)), 1, 0)
  ) %>%
  ungroup() %>%
  dplyr::select(-t1:-t9, -t1_m:-t9_m)

#############
# ENPG 2016 #
#############
enpg16 <- read_rds("https://github.com/ACC1240138/FONDECYT-REGULAR-/raw/main/rawdata/enpg2016.RDS")
exp16 <- read_rds("https://github.com/ACC1240138/FONDECYT-REGULAR-/raw/main/rawdata/Expansion16.RDS")

data16 <- enpg16 %>% 
  left_join(exp16, by = "idencuesta") %>% 
  mutate(year = 2016) %>% 
  dplyr::select(id = idencuesta, year, region = región, comuna, exp = Fexp.x, sexo, edad,
         nedu1 = dp_9_a, nedu2 = dp_10_a, religion = dp_5, ecivil = dp_2, 
         p_originario = dp_4, ingreso = dp_13, tab1 = st_3, tab2 = st_6, tab3 = st_7,
         tab4 = st_8, mar1 = mar_1, mar2 = mar_4, coc1 = coc_1, coc2 = coc_4, tranq_vida = trans_1_j,
         t1_m = trans_3_a, t2_m = trans_3_b, t3_m = trans_3_c, t4_m = trans_3_d,
         t5_m = trans_3_e, t6_m = trans_3_f, t7_m = trans_3_g, t8_m = trans_3_h, 
         t9_m = trans_3_i,
         oh1 = oh_1, oh2 = oh_4,oh3 = oh_5, audit1 = oh_14, audit2 = oh_15,
         audit3 = oh_16) %>% 
mutate(id = factor(id),
  sexo = factor(sexo, levels = c(1,2), labels = c("Hombre","Mujer")),
  mar1 = factor(mar1, levels = c(1,2), labels = c("Si","No")),
  mar2 = factor(mar2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
  tab1 = factor(tab1, levels = c(1,2), labels = c("Si","No")),
  tab2 = factor(tab2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
  coc1 = factor(coc1, levels = c(1,2), labels = c("Si","No")),
  coc2 = factor(coc2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
  tranq_vida = ifelse(is.na(tranq_vida), 1,0),
       oh1 = factor(oh1, levels = c(1,2), labels = c("Si","No")),
       oh2 = factor(oh2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
       audit1 = factor(audit1, levels = c (0:4), labels = c("Nunca","1 vez al mes o menos"," 2 a 4 veces al mes",
                                                            "2 a 3 veces a la semana","4 o mas veces a la semana")),
       audit2 = factor(audit2, levels = c (0:4), labels = c("0-2","3-4","5-6",
                                                            "7-8", "9 o mas")),
       audit3 = factor(audit3, levels = c (0:4), labels = c("Nunca","menos de 1 vez al mes",
                                                            "mensualmente",
                                                            "semanalmente",
                                                            "todos o casi todos los dias")),
       ecivil = factor(case_when(
         ecivil == 1 ~ "soltero",
         ecivil == 2 & ecivil == 6 ~ "casado",
         ecivil == 4  ~ "viudo",
         ecivil == 3 | ecivil == 5 ~ "separado")),
       religion = factor(case_when(religion == 1 ~ "Catolico",
                                   religion == 2 ~ "evangelico/protestante",
                                   religion >= 3 & religion <= 11 ~ "otra",
                                   religion == 12 ~ "ninguna")),
       p_originario = factor(ifelse(p_originario == 10,"no pertenece","pertenece")),
       nedu = factor(case_when(nedu1 < 4 & nedu2 == 2 ~ "basica incompleta",
                               nedu1 >= 2 & nedu1 < 5 & nedu2 == 1~ "basica completa",
                               nedu1 >= 5 & nedu1 <= 8 & nedu2 == 2 ~ "media completa",
                               nedu1 >= 5 & nedu1 <= 8 & nedu2 == 1  ~ "media completa",
                               nedu1 >= 9 & nedu1 <= 10 & nedu2 == 2 ~ "superior incompleta",
                               nedu1 >= 9 & nedu1 <= 10 & nedu2 == 1 ~ "superior completa",
                               nedu1 >= 11 & nedu1 <= 13 ~ "superior completa",
                               TRUE~NA)),
       ingreso = na_if(ingreso, 88),
       ingreso = na_if(ingreso, 99)) %>% 
  dplyr::select(-nedu1, -nedu2) %>% 
  rowwise() %>%
  mutate(
    tranq_mes = if_else(any(c_across(t1_m:t9_m) %in% c(1, 2)), 1, 0)
  ) %>%
  ungroup() %>%
  dplyr::select(-t1_m:-t9_m)

#############
# ENPG 2018 #
#############

enpg18 <- read_rds("https://github.com/ACC1240138/FONDECYT-REGULAR-/raw/main/rawdata/enpg2018.RDS")
data18 <- enpg18 %>% 
  mutate(year = 2018) %>% 
  dplyr::select(id = SbjNum, year, region = Region, comuna, exp = Fexp, sexo = S01, edad = S02,
         nedu1 = T_DP_12_1, nedu2 = T_DP_13_1, religion = DP_5, ecivil = DP_2, 
         p_originario = DP_4, ingreso = DP_16, tab1 = ST_2, tab2 = ST_5, tab3 = ST_6,
         tab4 = ST_7, mar1 = MAR_1, mar2 = MAR_4, coc1 = COC_1, coc2 = COC_4, tranq_vida = TRANS_1_O1,
         t1_m = T_TRANS_3_1, t2_m = T_TRANS_3_2, t3_m = T_TRANS_3_3, t4_m = T_TRANS_3_4,
         t5_m = T_TRANS_3_5, t6_m = T_TRANS_3_6, t7_m = T_TRANS_3_7, t8_m = T_TRANS_3_8,
         t9_m = T_TRANS_3_9,oh1 = OH_1, oh2 = OH_4, oh3 = OH_5, audit1 = OH_12, audit2 = OH_13,
         audit3 = OH_14) %>% 
  mutate(id = factor(id),
    sexo = factor(sexo, levels = c(1,2), labels = c("Hombre","Mujer")),
    mar1 = factor(mar1, levels = c(1,2), labels = c("Si","No")),
    mar2 = factor(mar2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
    tab1 = factor(tab1, levels = c(1,2), labels = c("Si","No")),
    tab2 = factor(tab2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
    coc1 = factor(coc1, levels = c(1,2), labels = c("Si","No")),
    coc2 = factor(coc2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
    tranq_vida = ifelse(tranq_vida == 10, 0,1),
         oh1 = factor(oh1, levels = c(1,2), labels = c("Si","No")),
         oh2 = factor(oh2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
         audit1 = factor(audit1, levels = c (0:4), labels = c("Nunca","1 vez al mes o menos"," 2 a 4 veces al mes",
                                                              "2 a 3 veces a la semana","4 o mas veces a la semana")),
         audit2 = factor(audit2, levels = c (0:4), labels = c("0-2","3-4","5-6",
                                                              "7-8", "9 o mas")),
         audit3 = factor(audit3, levels = c (0:4), labels = c("Nunca","menos de 1 vez al mes",
                                                              "mensualmente",
                                                              "semanalmente",
                                                              "todos o casi todos los dias")),
         ecivil = factor(case_when(
           ecivil == 1 ~ "soltero",
           ecivil == 2 & ecivil == 6 ~ "casado",
           ecivil == 4  ~ "viudo",
           ecivil == 3 | ecivil == 5 ~ "separado")),
         religion = factor(case_when(religion == 1 ~ "Catolico",
                                     religion == 2 ~ "evangelico/protestante",
                                     religion >= 3 & religion <= 11 ~ "otra",
                                     religion == 12 ~ "ninguna")),
         p_originario = factor(ifelse(p_originario == 10,"no pertenece","pertenece")),
         nedu = factor(case_when(nedu1 < 4 & nedu2 == 2 ~ "basica incompleta",
                                 nedu1 >= 2 & nedu1 < 5 & nedu2 == 1~ "basica completa",
                                 nedu1 >= 5 & nedu1 <= 8 & nedu2 == 2 ~ "media completa",
                                 nedu1 >= 5 & nedu1 <= 8 & nedu2 == 1  ~ "media completa",
                                 nedu1 >= 9 & nedu1 <= 10 & nedu2 == 2 ~ "superior incompleta",
                                 nedu1 >= 9 & nedu1 <= 10 & nedu2 == 1 ~ "superior completa",
                                 nedu1 >= 11 & nedu1 <= 13 ~ "superior completa",
                                 TRUE~NA)),
         ingreso = na_if(ingreso, 88),
         ingreso = na_if(ingreso, 99)) %>% 
  dplyr::select(-nedu1, -nedu2)%>% 
  rowwise() %>%
  mutate(
    tranq_mes = if_else(any(c_across(t1_m:t9_m) %in% c(1, 2)), 1, 0)
  ) %>%
  ungroup() %>%
  dplyr::select(-t1_m:-t9_m)

#############
# ENPG 2020 #
#############

enpg20 <- read_rds("https://github.com/ACC1240138/FONDECYT-REGULAR-/raw/main/rawdata/enpg2020.RDS")
data20 <- enpg20 %>% 
  mutate(year = 2020) %>% 
  dplyr::select(id = SbjNum, year, region = REGION, comuna = Nom_comuna,exp = FACT_PERS_COMUNA, sexo = S01, edad = S02,
         nedu1 = DP_12, nedu2 = DP_13, religion = DP_5, ecivil = DP_2, 
         p_originario = DP_4, ingreso = DP_16, tab1 = ST_2, tab2 = ST_5, tab3 = ST_6,
         tab4 = ST_7, mar1 = MAR_1, mar2 = MAR_4, coc1 = COC_1, coc2 = COC_4, t1 = TRANS_1_O1,
         t2 = TRANS_1_O2, t3 = TRANS_1_O3, t4 = TRANS_1_O4, t5 = TRANS_1_O5, t6 = TRANS_1_O6,
         t7 = TRANS_1_O7, t8 = TRANS_1_O8, t9 = TRANS_1_O9,
         t1_m = T_TRANS_2_1, t2_m = T_TRANS_2_2, t3_m = T_TRANS_2_3, t4_m = T_TRANS_2_4,
         t5_m = T_TRANS_2_5, t6_m = T_TRANS_2_6, t7_m = T_TRANS_2_7, t8_m = T_TRANS_2_8,
         t9_m = T_TRANS_2_9,oh1 = OH_1, oh2 = OH_4, oh3 = OH_5, audit1 = OH_8, audit2 = OH_9,
         audit3 = OH_10) %>% 
  mutate(id = factor(id),
    sexo = factor(sexo, levels = c(1,2), labels = c("Hombre","Mujer")),
    mar1 = factor(mar1, levels = c(1,2), labels = c("Si","No")),
    mar2 = factor(mar2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
    tab1 = factor(tab1, levels = c(1,2), labels = c("Si","No")),
    tab2 = factor(tab2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
    coc1 = factor(coc1, levels = c(1,2), labels = c("Si","No")),
    coc2 = factor(coc2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
         oh1 = factor(oh1, levels = c(1,2), labels = c("Si","No")),
         oh2 = factor(oh2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
         audit1 = factor(audit1, levels = c (0:4), labels = c("Nunca","1 vez al mes o menos"," 2 a 4 veces al mes",
                                                              "2 a 3 veces a la semana","4 o mas veces a la semana")),
         audit2 = factor(audit2, levels = c (0:4), labels = c("0-2","3-4","5-6",
                                                              "7-8", "9 o mas")),
         audit3 = factor(audit3, levels = c (0:4), labels = c("Nunca","menos de 1 vez al mes",
                                                              "mensualmente",
                                                              "semanalmente",
                                                              "todos o casi todos los dias")),
         ecivil = factor(case_when(
           ecivil == 1 ~ "soltero",
           ecivil == 2 & ecivil == 6 ~ "casado",
           ecivil == 4  ~ "viudo",
           ecivil == 3 | ecivil == 5 ~ "separado")),
         religion = factor(case_when(religion == 1 ~ "Catolico",
                                     religion == 2 ~ "evangelico/protestante",
                                     religion >= 3 & religion <= 11 ~ "otra",
                                     religion == 12 ~ "ninguna")),
         p_originario = factor(ifelse(p_originario == 10,"no pertenece","pertenece")),
         nedu = factor(case_when(nedu1 < 4 & nedu2 == 2 ~ "basica incompleta",
                                 nedu1 >= 2 & nedu1 < 5 & nedu2 == 1~ "basica completa",
                                 nedu1 >= 5 & nedu1 <= 8 & nedu2 == 2 ~ "media completa",
                                 nedu1 >= 5 & nedu1 <= 8 & nedu2 == 1  ~ "media completa",
                                 nedu1 >= 9 & nedu1 <= 10 & nedu2 == 2 ~ "superior incompleta",
                                 nedu1 >= 9 & nedu1 <= 10 & nedu2 == 1 ~ "superior completa",
                                 nedu1 >= 11 & nedu1 <= 13 ~ "superior completa",
                                 TRUE~NA)),
         ingreso = na_if(ingreso, 88),
         ingreso = na_if(ingreso, 99)) %>% 
  dplyr::select(-nedu1, -nedu2) %>% 
  rowwise() %>%
  mutate(
    tranq_vida = if_else(any(c_across(t1:t9) == 1), 1, 0),
    tranq_mes = if_else(any(c_across(t1_m:t9_m) %in% c(1, 2)), 1, 0)
  ) %>%
  ungroup() %>%
  dplyr::select(-t1:-t9, -t1_m:-t9_m)

#############
# ENPG 2022 #
#############

enpg22 <- read_rds("https://github.com/ACC1240138/FONDECYT-REGULAR-/raw/main/rawdata/enpg2022.RDS")
data22 <- enpg22 %>% 
  mutate(year = 2022) %>% 
  dplyr::select(id = FOLIO, year, region = REGION, comuna = COD_COMUNA, exp = FACTOR_EXPANSION, sexo = SEXO, edad = EDAD,
         nedu1 = DP_12, nedu2 = DP_13, religion = DP_5, ecivil = DP_2, 
         p_originario = DP_4, ingreso = DP_16, tab1 = ST_2, tab2 = ST_5, tab3 = ST_6,
         tab4 = ST_7, mar1 = MAR_1, mar2 = MAR_4, coc1 = COC_1, coc2 = COC_4,t1 = TRANS_101,
         t2 = TRANS_102, t3 = TRANS_103, t4 = TRANS_104, t5 = TRANS_105, t6 = TRANS_106,
         t7 = TRANS_107, t8 = TRANS_108, t9 = TRANS_109,
         t1_m = T_TRANS_2_1, t2_m = T_TRANS_2_2, t3_m = T_TRANS_2_3, t4_m = T_TRANS_2_4,
         t5_m = T_TRANS_2_5, t6_m = T_TRANS_2_6, t7_m = T_TRANS_2_7, t8_m = T_TRANS_2_8,
         t9_m = T_TRANS_2_9,
         oh1 = OH_1, oh2 = OH_4, oh3 = OH_5, audit1 = OH_8, audit2 = OH_9,
         audit3 = OH_10) %>% 
  mutate(id = factor(id),
    sexo = factor(sexo, levels = c(1,2), labels = c("Hombre","Mujer")),
    mar1 = factor(mar1, levels = c(1,2), labels = c("Si","No")),
    mar2 = factor(mar2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
    tab1 = factor(tab1, levels = c(1,2), labels = c("Si","No")),
    tab2 = factor(tab2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
    coc1 = factor(coc1, levels = c(1,2), labels = c("Si","No")),
    coc2 = factor(coc2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
         oh1 = factor(oh1, levels = c(1,2), labels = c("Si","No")),
         oh2 = factor(oh2, levels = c (1,2,3), labels = c("30 dias",">30",">1 año")),
         audit1 = factor(audit1, levels = c (0:4), labels = c("Nunca","1 vez al mes o menos"," 2 a 4 veces al mes",
                                                              "2 a 3 veces a la semana","4 o mas veces a la semana")),
         audit2 = factor(audit2, levels = c (0:4), labels = c("0-2","3-4","5-6",
                                                              "7-8", "9 o mas")),
         audit3 = factor(audit3, levels = c (0:4), labels = c("Nunca","menos de 1 vez al mes",
                                                              "mensualmente",
                                                              "semanalmente",
                                                              "todos o casi todos los dias")),
         ecivil = factor(case_when(
           ecivil == 1 ~ "soltero",
           ecivil == 2 & ecivil == 6 ~ "casado",
           ecivil == 4  ~ "viudo",
           ecivil == 3 | ecivil == 5 ~ "separado")),
         religion = factor(case_when(religion == 1 ~ "Catolico",
                                     religion == 2 ~ "evangelico/protestante",
                                     religion >= 3 & religion <= 11 ~ "otra",
                                     religion == 12 ~ "ninguna")),
         p_originario = factor(ifelse(p_originario == 10,"no pertenece","pertenece")),
         nedu = factor(case_when(nedu1 < 4 & nedu2 == 2 ~ "basica incompleta",
                                 nedu1 >= 2 & nedu1 < 5 & nedu2 == 1~ "basica completa",
                                 nedu1 >= 5 & nedu1 <= 8 & nedu2 == 2 ~ "media completa",
                                 nedu1 >= 5 & nedu1 <= 8 & nedu2 == 1  ~ "media completa",
                                 nedu1 >= 9 & nedu1 <= 10 & nedu2 == 2 ~ "superior incompleta",
                                 nedu1 >= 9 & nedu1 <= 10 & nedu2 == 1 ~ "superior completa",
                                 nedu1 >= 11 & nedu1 <= 13 ~ "superior completa",
                                 TRUE~NA)),
         ingreso = na_if(ingreso, 88),
         ingreso = na_if(ingreso, 99)) %>% 
  dplyr::select(-nedu1, -nedu2) %>% 
  rowwise() %>%
  mutate(
    tranq_vida = if_else(any(!is.na(c_across(t1:t9))), 1, 0),
    tranq_mes = if_else(any(c_across(t1_m:t9_m) %in% c(1, 2)), 1, 0)
  ) %>%
  ungroup() %>%
  dplyr::select(-t1:-t9, -t1_m:-t9_m)

#########
# MERGE #
#########

enpg_full <- bind_rows(data08, data10, data12, data14, data16, data18, data20, data22)

write_rds(enpg_full, "ENPG_FULL.RDS", compress = "gz")


########################
# ESTIMATE ALCOHOL USE #
########################

rm(list = ls())
gc()
enpg_full <- readRDS("enpg_full.RDS")
data <- enpg_full %>% 
  filter(edad >=15) %>% 
  mutate(oh3 = case_when(oh1 == "No"| oh2 == ">30" | oh2 == ">1 año"~ 0,
                         TRUE ~ oh3),
    prom_tragos = case_when(oh1 == "No" | oh2 == ">30" | oh2 == ">1 año"  ~ 0,
                                audit2 == "0-2"~ 1,
                                 audit2 == "3-4"~ 3.5,
                                 audit2 == "5-6"~ 5.5,
                                 audit2 == "7-8"~ 7.5,
                                 audit2 == "9 o mas" ~ 9),
         dias_binge = case_when(audit3 == "Nunca" | oh1 == "No"| oh2 == ">30" | oh2 == ">1 año" ~ 0,
                                audit3 == "menos de 1 vez al mes" ~ 0.5,
                                audit3 == "mensualmente" ~ 1,
                                audit3 == "semanalmente" ~ 4,
                                audit3 == "todos o casi todos los dias" ~ 20),
    dias_binge = ifelse(prom_tragos > 5.5, 0, dias_binge),
         diasalchab = oh3-dias_binge,
         diasalchab = ifelse(diasalchab < 0, 0, diasalchab), 
         volalchab = diasalchab*prom_tragos,
        volbinge = dias_binge*6, 
    voltotal = (volbinge + volalchab) * 13,
    voltotMS = (volbinge + volalchab) * 16,
    voltotdia = voltotal/30,
    voltotMINSAL = voltotMS/30,
    catohaj = case_when(
      sexo == "Mujer" & voltotdia == 0 ~ 0,
      sexo == "Mujer" & voltotdia > 0 & voltotdia <= 19.99 ~ 1, 
      sexo == "Mujer" & voltotdia >= 20 & voltotdia <= 39.99 ~ 2,
      sexo == "Mujer" & voltotdia >= 40 & voltotdia <= 100 ~ 3,
      sexo == "Mujer" & voltotdia > 100 ~ 4,
      sexo == "Hombre" & voltotdia == 0 ~ 0,
      sexo == "Hombre" & voltotdia > 0 & voltotdia <= 39.99 ~ 1,
      sexo == "Hombre" & voltotdia >= 40 & voltotdia <= 59.99 ~ 2,
      sexo == "Hombre" & voltotdia >= 60 & voltotdia <= 100 ~ 3,
      sexo == "Hombre" & voltotMINSAL >100 ~ 4,
      TRUE ~ NA_real_),
      catohaj = factor(catohaj, levels = 0:4, 
                              labels = c("Abstinentes", "Categoría 1", "Categoría 2", "Categoría 3", "Categoría 4")),
      catohMS = case_when(
      sexo == "Mujer" & voltotMINSAL == 0 ~ 0,
      sexo == "Mujer" & voltotMINSAL > 0 & voltotMINSAL <= 19.99 ~ 1,
      sexo == "Mujer" & voltotMINSAL >= 20 & voltotMINSAL <= 39.99 ~ 2,
      sexo == "Mujer" & voltotMINSAL >= 40 & voltotMINSAL <= 1000 ~ 3,
      sexo == "Hombre" & voltotMINSAL == 0 ~ 0,
      sexo == "Hombre" & voltotMINSAL > 0 & voltotMINSAL <= 39.99 ~ 1,
      sexo == "Hombre" & voltotMINSAL >= 40 & voltotMINSAL <= 59.99 ~ 2,
      sexo == "Hombre" & voltotMINSAL >= 60 & voltotMINSAL <= 100 ~ 3,
      sexo == "Hombre" & voltotMINSAL > 100 ~ 4,
      TRUE ~ NA_real_),
      catohMS = factor(catohMS, levels = 0:4, 
                       labels = c("Abstinentes", "Categoría 1", "Categoría 2", "Categoría 3", "Categoría 4")),
    volCH = (voltotdia*365),
    volCHMS = (voltotMINSAL*365)) %>% 
  filter(oh3 <=30)

rm(enpg_full)

# APC OMS
total_volCH <- data %>% 
  group_by(year) %>% 
  filter(!is.na(volCH)) %>% 
  summarise(pop = sum(exp),
            pc_totalvolCH = sum(volCH*exp)/pop) 

conversion <- function(x,vol){
  vol_oms = x*0.8
  oms=round((vol_oms*0.789)*1000,2)
  round(oms/vol,2)
}

# 2008 = 7.8
total_volCH[1,3]/(7.8*0.789*1000)

conversion(7.8,total_volCH[1,3])
# 5.4

# 2010 = 7.8
conversion(7.8,total_volCH[2,3])
# 5.18

# 2012 = 7.8
conversion(7.8,total_volCH[3,3])
# 5.26

# 2014 = 7.8
conversion(7.8,total_volCH[4,3])
# 5.37

# 2016 = 6.7
conversion(6.7,total_volCH[5,3])
# 4.13

# 2018 = 6.7
conversion(6.7,total_volCH[6,3])
# 2.52

# 2020 = 7.5
conversion(7.5,total_volCH[7,3])
# 4.83

# 2022 = 7.5 (repetido)
conversion(7.5,total_volCH[8,3])
# 5.53

data <- data %>% 
mutate(volaj = case_when(year == 2008 ~ volCH*5.4,
                         year == 2010 ~ volCH*5.18,
                         year == 2012 ~ volCH*5.26,
                         year == 2014 ~ volCH*5.37,
                         year == 2016 ~ volCH*4.13,
                         year == 2018 ~ volCH*2.52,
                         year == 2020 ~ volCH*4.83,
                         year == 2022 ~ volCH*5.53),
       volajohdia = volaj/365,
       cvolaj = case_when(
  sexo == "Mujer" & volajohdia == 0 ~ 0,
  sexo == "Mujer" & volajohdia > 0 & volajohdia <= 19.99 ~ 1,
  sexo == "Mujer" & volajohdia >= 20 & volajohdia <= 39.99 ~ 2,
  sexo == "Mujer" & volajohdia >= 40 & volajohdia <= 100 ~ 3,
  sexo == "Mujer" & volajohdia > 100 ~ 4,
  sexo == "Hombre" & volajohdia == 0 ~ 0,
  sexo == "Hombre" & volajohdia > 0 & volajohdia <= 39.99 ~ 1,
  sexo == "Hombre" & volajohdia >= 40 & volajohdia <= 59.99 ~ 2,
  sexo == "Hombre" & volajohdia >= 60 & volajohdia <= 100 ~ 3,
  sexo == "Hombre" & volajohdia > 100 ~ 4,
  TRUE ~ NA_real_),
cvolaj = factor(cvolaj, levels = 0:4, 
                labels = c("Abstinentes", "Categoría 1", "Categoría 2", "Categoría 3", "Categoría 4"))) 

write_rds(data, "ENPG_FULL.rds")

# CALCULO DE FACTOR POR AÑO
# REVISAR LA OMS, EL CONSUMO PC EN LITROS SE MULTIPLICA POR 
# LA DENSIDAD DE ALCOHOL, PARA TRANSFORMAR A MASA
# QUE ES 0.789 (PARA TRANSFORMAR A KG)
# Y LUEGO MULTIPLICAR POR 1000 PARA PASAR A GRAMOS
# SE ESTIMA EL TOTAL Y SE DIVIDE POR EL TOTAL DE LA ENCUESTA
# CON LOS FACTORES DE EXPANSION INCLUIDAS
# TOTALIDAD DE ALCOHOL REGISTRADO EN LA ENCUESTA
# ESO DA EL PER CAPITA
# ESE PER CAPITA SERÁ MENOR AL REPORTADO EN LA OMS.
# LUEGO LO QUE HAY QUE HACER ES ENCONTRAR EL FACTOR
# QUE MULTIPLICADO DA EL CONSUMO DE ALCOHOL REPORTADO EN LA OMS
# crear la categoría 4 DE CONSUMO DE ALCOHOL.

# AJUSTAR EL INGRESO SEGUN IPC



