

# librerías
pacman::p_load(tidyverse,jtools,polycor,ggplot2,ggstatsplot,ggcorrplot,broom,survey,
               kableExtra,scales,panelr,sjPlot,sjlabelled,sjmisc,stargazer,skimr,texreg,
               igraph, signnet, ggraph, extrafont, forcats, xtable, Hmisc, psych, psy,
               nFactors, GPArotation, psychTools, here, LMest, tidyr, tigerstats)

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
  dplyr::mutate (conf_gral   = case_when(c02 == 1 ~ 1,  
                                         c02 == 2 ~ 0, 
                                         c02 == 3 ~ 1)) %>%
  dplyr::mutate (extranjero  = case_when(m45 %in% 1 ~ 0,
                                         m45 %in% 2 :8 ~ 1)) %>% 
  dplyr::mutate (mujer       = case_when(m0_sexo == 1 ~ 0,
                                         m0_sexo == 2 ~ 1)) %>% 
  dplyr::mutate (edad        = case_when(m0_edad %in% 18:24 ~ "18_24",
                                         m0_edad %in% 25:34 ~ "25_34",
                                         m0_edad %in% 35:44 ~ "35_44",
                                         m0_edad %in% 45:54 ~ "45_54",
                                         m0_edad %in% 55:64 ~ "55_64", 
                                         m0_edad %in% 65:88 ~ "65")) %>%
  dplyr::mutate (nivel_educ  = case_when(m01 %in% 1 :3 ~ "básica",
                                         m01 %in% 4 :5 ~ "media",
                                         m01 %in% 6 :7 ~ "técnica",
                                         m01 %in% 8 :10 ~ "univers")) %>%
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
 

# descriptivos
## sexo
t1<-rowPerc(xtabs(~c12_01+mujer,data=a_full))
t2<-rowPerc(xtabs(~c12_02+mujer,data=a_full))
t3<-rowPerc(xtabs(~c12_03+mujer,data=a_full))
t4<-rowPerc(xtabs(~c12_04+mujer,data=a_full))
t5<-rowPerc(xtabs(~c12_05+mujer,data=a_full))
t6<-rowPerc(xtabs(~c12_06+mujer,data=a_full))
t7<-rowPerc(xtabs(~c12_07+mujer,data=a_full))
t8<-rowPerc(xtabs(~c12_08+mujer,data=a_full))
#t9<-rowPerc(xtabs(~c12_09+mujer,data=a_full))

a<-rbind(t1,t2,t3,t4,t5,t6,t7,t8)
a<-as.data.frame(a)
a$part<-c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1)
a1<-a[a$part==1,]
a1$tipo<-c("JJVV","Religiosa","Deportiva","Caridad","Partido","Profesional","Sindicato","AAEE")
a1<-a1%>%dplyr::select("tipo","0","1")
#a1

## nivel educ. 
t1<-rowPerc(xtabs(~c12_01+nivel_educ,data=a_full))
t2<-rowPerc(xtabs(~c12_02+nivel_educ,data=a_full))
t3<-rowPerc(xtabs(~c12_03+nivel_educ,data=a_full))
t4<-rowPerc(xtabs(~c12_04+nivel_educ,data=a_full))
t5<-rowPerc(xtabs(~c12_05+nivel_educ,data=a_full))
t6<-rowPerc(xtabs(~c12_06+nivel_educ,data=a_full))
t7<-rowPerc(xtabs(~c12_07+nivel_educ,data=a_full))
t8<-rowPerc(xtabs(~c12_08+nivel_educ,data=a_full))
#t9<-rowPerc(xtabs(~+nivel_educ,data=a_full))

a<-rbind(t1,t2,t3,t4,t5,t6,t7,t8)
a<-as.data.frame(a)
a$part<-c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1)
a2<-a[a$part==1,]
a2<-a2%>% mutate("no"= básica + media + técnica)%>% mutate("si"= univers)%>%
   dplyr::select("no", "si")

## edad
t1<-rowPerc(xtabs(~c12_01+edad,data=a_full))
t2<-rowPerc(xtabs(~c12_02+edad,data=a_full))
t3<-rowPerc(xtabs(~c12_03+edad,data=a_full))
t4<-rowPerc(xtabs(~c12_04+edad,data=a_full))
t5<-rowPerc(xtabs(~c12_05+edad,data=a_full))
t6<-rowPerc(xtabs(~c12_06+edad,data=a_full))
t7<-rowPerc(xtabs(~c12_07+edad,data=a_full))
t8<-rowPerc(xtabs(~c12_08+edad,data=a_full))
#t9<-rowPerc(xtabs(~c12_09+edad,data=asociaciones_2016d))

a<-rbind(t1,t2,t3,t4,t5,t6,t7,t8)
a<-as.data.frame(a)
a$part<-c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1)
a3<-a[a$part==1,]

a3<-a3 %>%rowwise() %>% 
  mutate(joven = sum(c_across("18_24":"25_34"), na.rm = T))%>%
  mutate(adulto = sum(c_across("35_44":"45_54"), na.rm = T))%>%
  mutate(mayor = sum(c_across("55_64":"65"), na.rm = T)) %>%
  dplyr::select(joven, adulto, mayor)

tabla1<-cbind(a1,a2,a3)
tabla1<-as.data.frame(tabla1)
tabla1<-tibble::rowid_to_column(tabla1, "id");tabla1$id<-NULL
tabla1<-tabla1 %>% remove_rownames %>% column_to_rownames(var="tipo")
colnames(tabla1)<-c("hombre","mujer","no","si","joven","adulto","mayor")

## tabla
kable(tabla1,
      "latex",
      booktabs = T,
      align = c("r"),
      col.names = c("hombre","mujer","no","si","joven","adulto","mayor")) %>%
  kable_classic("hover", full_width = F)%>%
  add_header_above(c(" " = 1,"Sexo" = 2, "Educ. universitaria" = 2, "Edad"=3))


## plot items
a_plot<-a_full %>% 
  select(c12_01:c12_08) %>%
  mutate_at(vars(matches("c12")), ~ifelse(. == 1, "Miembro", "No miembro"))%>%
  mutate_if(is.character, as.factor)
 
set_theme(legend.pos = "top")
 plot_likert(a_plot,
             geom.colors = c("#3F00FF","#28282B"),
            axis.labels=c("JJVV","Religiosa","Deportiva","Caridad","Partido","Profesional","Sindicato","AAEE"),
            reverse.colors=F,
            #cat.neutral=NULL,
            values="sum.outside",
            show.prc.sign=F)


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

## Estimar errores estandar con bootstrap paramétrico **(no estimar)**
mboot <- bootstrap(modelo3, n = 581, B = 2, seed = 172)

## Tablas
options(scipen=999)
print(xtable(Be, type = "latex"), file = "output/multinom_coeff.tex")
print(xtable(p, type = "latex"), file = "output/multinom_pvaluestex")

# Probabilidades iniciales 
plot(modelo3, what="CondProb")

## Graficar modelo 3
LMmodelo3 <- reshape2::melt(modelo3$Psi, level=1)
LMmodelo3 = LMmodelo3 %>%
  dplyr::mutate (clase = case_when(state == 1 ~ "Class 1\n Closed (36%)",
                                   state == 2 ~ "Class 2\n Broker (10%)",
                                   state == 3 ~ "Class 3\n Apathetic (54%)")) %>%
  dplyr::mutate (category = case_when(category == 0 ~ "no (not member)",
                                      category == 1 ~ "yes (member)")) 

LMmodelo3$item <- plyr::mapvalues(LMmodelo3$item, 
                              c('1','2','3','4','5','6','7','8'),
                              c("neighborhood","religious","political","union","professional","charity","sports","student"))
level_order <- c("neighborhood","religious","charity","political","union","professional","sports","student") 


## plot 
ggplot(LMmodelo3,aes(x = factor(item, level = level_order), y = value, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(clase ~ .) +
  scale_fill_manual(values = c("#28282B", "#3F00FF")) + 
  labs(x = "",y="", fill ="") + theme(text = element_text(size=15)) +
  theme(axis.ticks.y=element_blank(),
                   #legend.position = "top",
                   panel.grid.major.y=element_blank(),
                   plot.title = element_text(hjust = 0.5, size = 8),
                   axis.title = element_text(size=10),
                   axis.text.x = element_text(size = 11,angle = 45, hjust = 0.9),
                   axis.text.y = element_text(size = 9)) +
  guides(fill = guide_legend(reverse=F)) 


# Distribución marginal
plot(modelo3, what="marginal")

# Transiciones
plot(modelo3,what="transitions") # figura 3. 

# Pivot wider and cbind
d1<-modelo3[["data"]]%>%select(idencuesta, ola, ponderador02, mujer, edad, nivel_educ)
d1<- panel_data(d1, id = idencuesta, wave = ola) 
d1_wide<-d1<-panelr::widen_panel(d1, separator = "_")
mseq<-modelo3[["Ul"]] #matrix containing the predicted sequence of latent states by the local decoding method
mseq<-as_tibble(mseq); colnames(mseq)<-c("e1", "e2", "e3")
d1_wide<-cbind(d1_wide,mseq)
save(d1_wide, file = "data/d1_wide.RData")


