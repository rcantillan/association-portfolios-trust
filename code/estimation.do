************************************************************************************************
* Subject		: estimate logit model and panel model
* By 			: Gustavo Ahumada
* Last update	: 05/July/2023	 
************************************************************************************************


clear all 
* global directory
global data "/Volumes/GoogleDrive/My Drive/Reserach_in_progress/social_turst_membership/sem_trabajo/data"


* open dataset
use "$data/ELSOC_Long.dta", clear
* merge data
merge 1:1 idencuesta ola  using "$data/data.long.dta"
* keep match
keep if _merge == 3
drop _merge

* replace values to missing
mvdecode _all, mv(-888)
mvdecode _all, mv(-999)
mvdecode _all, mv(-666)

* generate confianza generalizada
gen conf_gen = 1 if c02 == 1 | c02 == 3
replace conf_gen = 0 if c02 == 2

* generate confianza en vecinos
gen conf_vec = 1 if t01 > 2
replace conf_vec = 0 if t01 < 3

*save 
save "$data/data.dta", replace

* filter
keep if muestra == 1
keep if tipo_atricion == 1
keep if ola == 1 | ola == 3 | ola == 6

* declarate panel data
xtset idencuesta ola

***************
* Estimatations
***************

/* probit model */

* model1
probit conf_gen i.clase i.m0_sexo m0_edad i.ola i.m01 ///
i.m02 i.m36 [pweight=ponderador_long_panel], vce(cluster comuna_cod)
* average marginal effects
margins, dydx(*) post
estimates store model1
estadd local Controls "Yes"
estadd local TE "Yes"

* mode2
probit conf_vec i.clase i.m0_sexo m0_edad i.ola i.m01 ///
i.m02 i.m36 [pweight=ponderador_long_panel], vce(cluster comuna_cod)
* average marginal effects
margins, dydx(*) post
estimates store model2
estadd local Controls "Yes"
estadd local TE "Yes"


/* rando effects probit model */

* model3
xtprobit conf_gen i.clase i.m0_sexo m0_edad i.ola i.m01 ///
i.m02 i.m36, vce(cluster idencuesta)
* average marginal effects
margins, dydx(*) post
estimates store model3
estadd local Controls "Yes"
estadd local TE "Yes"


* model4
xtprobit conf_vec i.clase i.m0_sexo m0_edad i.ola i.m01 ///
i.m02 i.m36, vce(cluster idencuesta)
* average marginal effects
margins, dydx(*) post
estimates store model4
estadd local Controls "Yes"
estadd local TE "Yes"


esttab model1 model2 model3 model4 using "$data/table5.tex", replace b(3) se(3) nonumbers ///
keep(2.clase 3.clase) ///
coeflabels(2.clase "Clase 2 (broker)" 3.clase "Clase 3 (apathetic)") ///
mgroups("Probit model" "Random-effects probit model", pattern(1 0 1 0) ///
prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) alignment(D{.}{.}{-1}) ///
mtitles("Gen. trust" "NBHD trust" "Gen. trust" "NBHD trust") ///
star(* 0.10 ** 0.05 *** 0.01) ///
scalar("TE Time Effects" "Controls Controls") ///
title("Marginal effects of associative behavior on trust \label{tab:table5}")



xtprobit conf_gen L.conf_gen i.m0_sexo m0_edad i.ola i.m01 ///
i.m02 i.m36, re 






