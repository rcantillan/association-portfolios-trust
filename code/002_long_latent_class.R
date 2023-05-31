

# librerías
pacman::p_load(tidyverse,jtools,polycor,ggplot2,ggstatsplot,ggcorrplot,broom,survey,
               kableExtra,scales,panelr,sjPlot,sjlabelled,sjmisc,stargazer,skimr,texreg,
               igraph, signnet, ggraph, extrafont, forcats, xtable, Hmisc, psych, psy,
               nFactors, GPArotation, psychTools, here, LMest)


# load data
load(here::here("data/ELSOC_Long.RData"))

# seleccionar variables 
a=elsoc_long_2016_2022;rm(elsoc_long_2016_2022)
a_full<-a%>%
  filter(ola==1 | ola==3 | ola==6) %>%
  mutate_at(vars(matches("r01")), ~ifelse(. > 2, 1, 0)) %>% 
  mutate_at(vars(matches("c12")), ~ifelse(. < 2, 0, 1)) %>% 
  select(idencuesta, ola, ponderador02, m0_sexo, m0_edad, m01, 
         c12_01, c12_02, c12_03, c12_04, c12_05, c12_06, c12_07, c12_08, 
         r01_01, r01_02, r01_03, r01_04, r01_05, r01_06, r01_07, r01_08, r01_09, 
         r01_10, r01_11, r01_12, r01_13)



