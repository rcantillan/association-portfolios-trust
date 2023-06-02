

# librerías
pacman::p_load(tidyverse,jtools,polycor,ggplot2,ggstatsplot,ggcorrplot,broom,survey,
               kableExtra,scales,panelr,sjPlot,sjlabelled,sjmisc,stargazer,skimr,texreg,
               igraph, signnet, ggraph, extrafont, forcats, xtable, Hmisc, psych, psy,
               nFactors, GPArotation, psychTools, here, LMest, tidyr)

options(knitr.kable.NA = '')
# load data
load(here::here("data/ELSOC_Long.RData"))

# seleccionar variables 
a=elsoc_long_2016_2022;rm(elsoc_long_2016_2022)
a_full<-a%>%
  filter(ola==1 | ola==3 | ola==6) %>%
  mutate(ola=case_when(ola==1~1,ola==3~2,ola==6~3))%>%
  mutate_at(vars(matches("r01")), ~ifelse(. > 2, 1, 0)) %>% 
  mutate_at(vars(matches("c12")), ~ifelse(. < 2, 0, 1)) %>% 
  dplyr::mutate (conf_gral   = case_when(c02       == 1       ~ 1,  
                                         c02       == 2       ~ 0, 
                                         c02       == 3       ~ 1)) %>%
  dplyr::mutate (extranjero  = case_when(m45     %in% 1       ~ 0,
                                         m45     %in% 2 :8    ~ 1)) %>% 
  dplyr::mutate (mujer       = case_when(m0_sexo   == 1       ~ 0,
                                         m0_sexo   == 2       ~ 1)) %>% 
  dplyr::mutate (edad        = case_when(m0_edad %in% 18:24   ~ "18-24",
                                         m0_edad %in% 25:34   ~ "25-34",
                                         m0_edad %in% 35:44   ~ "35-44",
                                         m0_edad %in% 45:54   ~ "45-54",
                                         m0_edad %in% 55:64   ~ "55-64", 
                                         m0_edad %in% 65:88   ~ "65-")) %>%
  dplyr::mutate (nivel_educ  = case_when(m01     %in% 1 :3    ~ "básica",
                                         m01     %in% 4 :5    ~ "media",
                                         m01     %in% 6 :7    ~ "técnica",
                                         m01     %in% 8 :10   ~ "univers")) %>%
  select(idencuesta, ola, ponderador02, mujer, edad, nivel_educ, 
         c12_01, c12_02, c12_03, c12_04, c12_05, c12_06, c12_07, c12_08)

# set NA
a_full[a_full=="-999"] <- NA
a_full[a_full=="-888"] <- NA

# delete NA in covariable
a_full<-a_full %>% drop_na(c("mujer", "edad", "nivel_educ"))

# set panel 
a_full<- panel_data(a_full, id = idencuesta, wave = ola) %>%  
  complete_data(min.waves = "all") %>%
  as.data.frame()
 

# Analysis Latent Markov 
## Modelo 1
modelo1 <-  lmest(responsesFormula = c12_01+c12_02+c12_03+c12_04+c12_05+c12_06+c12_07+c12_08~NULL,
                  latentFormula =~ mujer+edad+nivel_educ,
                  index = c("idencuesta","ola"),
                  output = TRUE,
                  out_se = TRUE,
                  paramLatent = "multilogit",
                  data = a_full,
                  k = 1:5,
                  start = 1,
                  #modBasic = 1,
                  modManifest="LM",
                  seed = 123)

plot(modelo1,what="modSel") # de acuerdo con BIC sugiere modelo con 3 clases 

## Modelo selección de modelos 1:5 estados latentes. 
mod_sel <- lmestSearch(responsesFormula = c12_01+c12_02+c12_03+c12_04+c12_05+c12_06+c12_07+c12_08~NULL,
                   latentFormula =~ mujer+edad+nivel_educ, #mujer+edad+nivel_educ,
                   index = c("idencuesta","ola"), 
                   data = a_full,
                   output = TRUE,
                   out_se = TRUE,
                   version ="categorical",
                   paramLatent = "multilogit",
                   k = 1:5,
                   #modBasic = 1, 
                   seed = 123)

summary(mod_sel) # de acuerdo con el BIC se selecciona modelo de 3 clases

# Crear tabla de ajuste promedio. 
lk<-mod_sel[["lkv"]];lk<-as.data.frame(lk)
aic<-mod_sel[["Aic"]];aic<-as.data.frame(aic)
bic<-mod_sel[["Bic"]];bic<-as.data.frame(bic)
np<-c(8,46,104,102,200);np<-as.data.frame(np)
states<-c(1,2,3,4,5);states<-as.data.frame(states)
t1<-as.data.frame(c(states,lk,np,aic,bic));rm(states,lk,np,aic,bic)
print(xtable(t1, type = "latex"), file = "output/fit1.tex")

# Modelo de mejor ajuste k=3 y tabla de regresión multinomial
modelo3 <- mod_sel$out.single[[3]]
summary(modelo3)

## Resultados de regresión multinomial para describir grupos.
Be<-as.data.frame(modelo3$Be)
seBe<-as.data.frame(modelo3$seBe)
z <- modelo3$Be/modelo3$seBe
p <- (1 - pnorm(abs(z), 0, 1))*2 # two-tailed z test

options(scipen=999)
print(xtable(Be, type = "latex"), file = "output/multinom_coeff.tex")
print(xtable(p, type = "latex"), file = "output/multinom_pvaluestex")

# Probabilidades iniciales 
plot(modelo3, what="CondProb")

## Graficar modelo 3
LMmodelo3 <- reshape2::melt(modelo3$Psi, level=1)
LMmodelo3 = LMmodelo3 %>%
  dplyr::mutate (clase       = case_when(state      == 1 ~ "Class 1\n Closed (36%)",
                                         state      == 2 ~ "Class 2\n Broker (10%)",
                                         state      == 3 ~ "Class 3\n Apathetic (54%)")) %>%
  dplyr::mutate (category    = case_when(category   == 0 ~ "no (not member)",
                                         category   == 1 ~ "yes (member)")) 

LMmodelo3$item <- plyr::mapvalues(LMmodelo3$item, 
                              c('1','2','3','4','5','6','7','8'),
                              c("neighborhood","religious","political","union","professional","charity","sports","student"))
level_order <- c("neighborhood","religious","charity","political","union","professional","sports","student") 


## plot 
p1 <- ggplot(LMmodelo3,aes(x = factor(item, level = level_order), y = value, fill = category))
p1 <- p1 + geom_bar(stat = "identity", position = "stack") 
p1 <- p1 + facet_grid(clase ~ .) 
p1 <- p1 + scale_fill_manual(values = c("#28282B", "#3F00FF")) 
p1 <- p1 + labs(x = "",y="", fill ="") + theme(text = element_text(size=15))
p1 <- p1 + theme(axis.ticks.y=element_blank(),
                   legend.position = "left",
                   panel.grid.major.y=element_blank(),
                   plot.title = element_text(hjust = 0.5, size = 8),
                   axis.title = element_text(size=10),
                   axis.text.x = element_text(size = 11,angle = 45, hjust = 0.9),
                   axis.text.y = element_text(size = 9))
p1 <- p1 + guides(fill = guide_legend(reverse=F))
print(p1)

# Distribución marginal
plot(modelo3, what="marginal")

# Transiciones
plot(modelo3,what="transitions")







