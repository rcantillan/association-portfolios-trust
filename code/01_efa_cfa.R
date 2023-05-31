
# librerías
pacman::p_load(tidyverse,jtools,polycor,ggplot2,ggstatsplot,ggcorrplot,broom,survey,
               kableExtra,scales,panelr,sjPlot,sjlabelled,sjmisc,stargazer,skimr,texreg,
               igraph, signnet, ggraph, extrafont, forcats, xtable, Hmisc, psych, psy,
               nFactors, GPArotation, psychTools)

# load data
load("~/Documents/doctorado-UC/sem/ELSOC_Long.RData")

elsoc_long_2016_2022$c12_09

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

a_na<-a%>%
  filter(ola==1 | ola==3 | ola==6) %>%
  mutate(across(r01_01:r01_13, ~case_when(.x== -888 ~ NA_real_, .x== -999 ~ NA_real_, .x > 1 ~ 1, TRUE ~ 0)))%>%
  mutate(across(c12_01:c12_09, ~case_when(.x== -888 ~ NA_real_, .x== -999 ~ NA_real_, .x < 2 ~ 0, TRUE ~ 1)))%>%
  select(idencuesta, ola, ponderador02, m0_sexo, m0_edad, m01, 
         c12_01, c12_02, c12_03, c12_04, c12_05, c12_06, c12_07, c12_08, c12_09, 
         r01_01, r01_02, r01_03, r01_04, r01_05, r01_06, r01_07, r01_08, r01_09, 
         r01_10, r01_11, r01_12, r01_13)



# Descriptivos `Membresías` 





# EFA ola 1 
## trasnformar a factor  
asoc_ola1<-a_full %>% mutate_at(vars(matches("c12")), as.factor)%>%
  filter(ola==1)%>%
  select(c12_01, c12_02, c12_03, c12_04, c12_05, c12_06, c12_07, c12_08)

## correlation matrix
het_asoc_ola1 <- hetcor(asoc_ola1)$cor
ggcorrplot(het_asoc_ola1, hc.order = F, type = "lower", outline.col = "white")
KMO(het_asoc_ola1) 

## factor analysis
scree.plot(het_asoc_ola1)   
fa.parallel(het_asoc_ola1, n.obs=300)  

fa <- fa(r = het_asoc_ola1, nfactors = 3, n.obs = nrow(asoc_ola1), rotate = "varimax")
fa2latex(fa, caption = "Análisis Factorial Exploratorio Ola 1", heading = " ")

# Ver resultados del modelo que mejor ajusta
fa.diagram(fa, main = "Factores sugeridos")


## alternativa
#fa.1 <- factanal(covmat = het_asoc_ola1, factors = 3, rotation = "varimax")
#fa.1

# EFA ola 3
## trasnformar a factor  
asoc_ola3<-a_full %>% mutate_at(vars(matches("c12")), as.factor)%>%
  filter(ola==3)%>%
  select(c12_01, c12_02, c12_03, c12_04, c12_05, c12_06, c12_07, c12_08)

## correlation matrix
het_asoc_ola3 <- hetcor(asoc_ola3)$cor
ggcorrplot(het_asoc_ola3, hc.order = F, type = "lower", outline.col = "white")
KMO(het_asoc_ola3) 

## factor analysis
scree.plot(het_asoc_ola3)   
fa.parallel(het_asoc_ola3, n.obs=300)  

fa2 <- fa(r = het_asoc_ola3, nfactors = 3, n.obs = nrow(asoc_ola3), rotate = "varimax")
fa2latex(fa2, caption = "Análisis Factorial Exploratorio Ola 3", heading = " ")

# Ver resultados del modelo que mejor ajusta
fa.diagram(fa2, main = "Factores sugeridos")

# EFA ola 6
## trasnformar a factor  
asoc_ola6<-a_full %>% mutate_at(vars(matches("c12")), as.factor)%>%
  filter(ola==6)%>%
  select(c12_01, c12_02, c12_03, c12_04, c12_05, c12_06, c12_07, c12_08)

## correlation matrix
het_asoc_ola6 <- hetcor(asoc_ola6)$cor
ggcorrplot(het_asoc_ola6, hc.order = F, type = "lower", outline.col = "white")
KMO(het_asoc_ola6) 

## factor analysis
scree.plot(het_asoc_ola6)   
fa.parallel(het_asoc_ola6, n.obs=300)  

fa3 <- fa(r = het_asoc_ola6, nfactors = 3, n.obs = nrow(asoc_ola6), rotate = "varimax")
fa2latex(fa3, caption = "Análisis Factorial Exploratorio Ola 6", heading = " ")

# Ver resultados del modelo que mejor ajusta
fa.diagram(fa3, main = "Factores sugeridos")



# Ver cargas de items
EFA_model$loadings
head(EFA_model$scores)

# evaluae estadísticos relativos 
EFA_model $BIC
EFA_model2$BIC
EFA_model3$BIC
EFA_model4$BIC


# CFA 



# latent class analysis 






























