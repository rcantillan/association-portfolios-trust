
# librerías
pacman::p_load(tidyverse,jtools,polycor,ggplot2,ggstatsplot,ggcorrplot,broom,survey,
               kableExtra,scales,panelr,sjPlot,sjlabelled,sjmisc,stargazer,skimr,texreg,
               igraph, signnet, ggraph, extrafont, forcats, xtable, Hmisc, psych, psy,
               nFactors, GPArotation, psychTools, here, semTools, influence.SEM, lavaan,
               semPlot, car, stringr, semTable)


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


# EFA confianza social generalizada 
# load data
load(here::here("data/ELSOC_Long.RData"))

# seleccionar variables 
a=elsoc_long_2016_2022;rm(elsoc_long_2016_2022)
a_full<-a%>%
  filter(ola==1 | ola==3 | ola==6) %>%
  dplyr::mutate (c02 = case_when(c02 == 1 ~ 1, c02 == 2 ~ 0, c02 == 3 ~ 1)) %>% 
  dplyr::mutate (c03 = case_when(c03 == 1 ~ 1, c03 == 2 ~ 0, c03 == 3 ~ 1)) %>%
  dplyr::mutate (c04 = case_when(c04 == 1 ~ 0, c04 == 2 ~ 1, c04 == 3 ~ 1)) %>%
  select(idencuesta, ola, ponderador02, m0_sexo, m0_edad, m01, 
         c02,c03,c04)

a_full[a_full=="-999"] <- NA
a_full[a_full=="-888"] <- NA
a_full=na.omit(a_full)

## plot items
a_plot_ola1<-a_full %>% 
  filter(ola==1)%>%
  select(c02,c03,c04) %>%
  mutate_at(vars(matches("c")), ~ifelse(. == 1, "si", "no"))%>%
  mutate_if(is.character, as.factor)

set_theme(legend.pos = "top")
plot_likert(a_plot_ola1,
            geom.colors = c("#3F00FF","#28282B"),
            axis.labels=c("Confianza gral.", "Altruismo gral.", "La gente trata de ser justa"),
            reverse.colors=F,
            #cat.neutral=NULL,
            values="sum.outside",
            show.prc.sign=F)

# EFA ola 1 
## trasnformar a factor  
conf_ola1<-a_full %>% mutate_at(vars(matches("c")), as.factor)%>%
  filter(ola==1)%>%
  select(c02,c03,c04)

## correlation matrix
het_asoc_ola1 <- hetcor(conf_ola1)$cor
ggcorrplot(het_asoc_ola1, hc.order = F, type = "lower", outline.col = "white",lab = TRUE)
KMO(het_asoc_ola1) 

## factor analysis
scree.plot(het_asoc_ola1)   
fa.parallel(het_asoc_ola1, n.obs=300)  

fa <- fa(r = het_asoc_ola1, nfactors = 1, n.obs = nrow(asoc_ola1), rotate = "varimax")
summary(fa)
fa2latex(fa, caption = "Análisis Factorial Exploratorio Ola 1", heading = " ")
xtable(unclass(fa2$loadings))

# Ver resultados del modelo que mejor ajusta
fa.diagram(fa, main = "Factores sugeridos")


# EFA ola 3
a_plot_ola3<-a_full %>% 
  filter(ola==3)%>%
  select(c02,c03,c04) %>%
  mutate_at(vars(matches("c")), ~ifelse(. == 1, "si", "no"))%>%
  mutate_if(is.character, as.factor)

set_theme(legend.pos = "top")
plot_likert(a_plot_ola3,
            geom.colors = c("#3F00FF","#28282B"),
            axis.labels=c("Confianza gral.", "Altruismo gral.", "La gente trata de ser justa"),
            reverse.colors=F,
            #cat.neutral=NULL,
            values="sum.outside",
            show.prc.sign=F)

## trasnformar a factor  
conf_ola3<-a_full %>% mutate_at(vars(matches("c")), as.factor)%>%
  filter(ola==3)%>%
  select(c02,c03,c04)

## correlation matrix
het_asoc_ola3 <- hetcor(conf_ola3)$cor
ggcorrplot(het_asoc_ola3, hc.order = F, type = "lower", outline.col = "white",lab = TRUE)
KMO(het_asoc_ola3) 

## factor analysis
scree.plot(het_asoc_ola3)   
fa.parallel(het_asoc_ola3, n.obs=300)  

fa2 <- fa(r = het_asoc_ola3, nfactors = 1, n.obs = nrow(asoc_ola3), rotate = "varimax")
fa2latex(fa2, caption = "Análisis Factorial Exploratorio Ola 3", heading = " ")
xtable(unclass(fa2$loadings))


# Ver resultados del modelo que mejor ajusta
fa.diagram(fa2, main = "Factores sugeridos")


# EFA ola 6
a_plot_ola6<-a_full %>% 
  filter(ola==6)%>%
  select(c02,c03,c04) %>%
  mutate_at(vars(matches("c")), ~ifelse(. == 1, "si", "no"))%>%
  mutate_if(is.character, as.factor)

#set_theme(legend.pos = "top")
plot_likert(a_plot_ola6,
            geom.colors = c("#3F00FF","#28282B"),
            axis.labels=c("Confianza gral.", "Altruismo gral.", "La gente trata de ser justa"),
            reverse.colors=F,
            #cat.neutral=NULL,
            values="sum.outside",
            show.prc.sign=F)

## trasnformar a factor  
conf_ola6<-a_full %>% mutate_at(vars(matches("c12")), as.factor)%>%
  filter(ola==6)%>%
  select(c02,c03,c04)

## correlation matrix
het_asoc_ola6 <- hetcor(conf_ola6)$cor
ggcorrplot(het_asoc_ola6, hc.order = F, type = "lower", outline.col = "white",lab = TRUE)
KMO(het_asoc_ola6) 

## factor analysis
scree.plot(het_asoc_ola6)   
fa.parallel(het_asoc_ola6, n.obs=300)  

fa3 <- fa(r = het_asoc_ola6, nfactors = 1, n.obs = nrow(asoc_ola6), rotate = "varimax")
fa
fa2latex(fa3, caption = "Análisis Factorial Exploratorio Ola 6", heading = " ")
xtable(unclass(fa3$loadings))

# Ver resultados del modelo que mejor ajusta
fa.diagram(fa3, main = "Factores sugeridos")


# CFA (confianza generalizada) ------------------------------------------------
## Ajustar modelo de ecuaciones estructurales SEM

m1 <- "
  # Factoriales
  CONF         =~ c02 + w*c03 + w*c04
  "
#CONF         =~ c02 + w*c03 + w*c04

#Evaluar ajuste del modelo ola 1 --------------------------
conf_ola1<-a_full %>% filter(ola==1)%>% select(ponderador02,c02,c03,c04)

m1_fit <- sem(m1, data = conf_ola1, std.lv=T,  
               ordered = c("c02", "c03", "c04"), 
               sampling.weights = "ponderador02", 
               estimator = "ULSMV")

#Evaluar

summary(m1_fit, 
        fit.measure=TRUE, 
        standardized=TRUE, 
        rsquare = TRUE)

vlabs <- c(c02 = "Confianza gral." , c03 = "Altruismo gral." , c04 = "La gente trata de ser justa")
m3_fit_table <- semTable(m1_fit, columns = c("estse", "p"), 
                         paramSets = c("loadings"), 
                         fits = c("chisq" , "rmsea"),  
                         file =  "output/m3_fit_table",
                         varLabels = vlabs,
                         type = "latex", 
                         caption = "CFA generalized trust ola 1", 
                         label = "tab : conf1", 
                         longtable = TRUE,
                         table.float = TRUE)

#Evaluar ajuste del modelo ola 3 --------------------------
conf_ola3<-a_full %>% filter(ola==3)%>% select(ponderador02,c02,c03,c04)

m2_fit <- sem(m1, data = conf_ola3, std.lv=T,  
              ordered = c("c02", "c03", "c04"), 
              sampling.weights = "ponderador02", 
              estimator = "ULSMV")

#Evaluar

summary(m2_fit, 
        fit.measure=TRUE, 
        standardized=TRUE, 
        rsquare = TRUE)

m3_fit_table <- semTable(m2_fit, columns = c("estse", "p"), 
                         paramSets = c("loadings"), 
                         fits = c("chisq" , "rmsea"),  
                         file =  "output/m3_fit_table",
                         varLabels = vlabs,
                         type = "latex", 
                         caption = "CFA generalized trust ola 3", 
                         label = "tab : conf2", 
                         longtable = TRUE,
                         table.float = TRUE)


#Evaluar ajuste del modelo ola 6 --------------------------
conf_ola6<-a_full %>% filter(ola==6)%>% select(ponderador02,c02,c03,c04)
m3_fit <- sem(m1, data = conf_ola6, std.lv=T,  
              ordered = c("c02", "c03", "c04"), 
              sampling.weights = "ponderador02", 
              estimator = "ULSMV")

#Evaluar
summary(m3_fit, 
        fit.measure=TRUE, 
        standardized=TRUE, 
        rsquare = TRUE)


parameterestimates(m3_fit, standardized=TRUE) ##CIs para parámetros
fitted(m3_fit) ##Tabla de covarianzas 
residuals(m3_fit) ##Residuos
fitmeasures(m1_fit) ##Indices de ajuste
modificationindices(m3_fit, sort. = TRUE) ##Indices de modifación 


m3_fit_table <- semTable(m3_fit, columns = c("estse", "p"), 
                   paramSets = c("loadings"), 
                   fits = c("chisq" , "rmsea"),  
                   file =  "output/m3_fit_table",
                   varLabels = vlabs,
                   type = "latex", 
                   caption = "CFA generalized trust", 
                   label = "tab : conf3", 
                   longtable = TRUE,
                   table.float = TRUE)








