* ==============================================================================
* estimation.do — (Optional) Stata replication (kept for parity with older pipeline)
* ------------------------------------------------------------------------------
* This do-file assumes you have exported data/dt_states.dta from R.
* Paths are relative to the project root.
* ==============================================================================

clear all
set more off

* --- Use relative paths ---
local ROOT = c(pwd)
use "data/dt_states.dta", clear

* Panel setup
xtset id ola

* RE probit for generalized trust
xtprobit trust i.position edad i.woman i.education i.employment i.ola, re vce(cluster id)
margins, dydx(*) post
est store m1

* RE probit for neighborhood trust
xtprobit trust_nh i.position edad i.woman i.education i.employment i.ola, re vce(cluster id)
margins, dydx(*) post
est store m2

* Export table
cap which esttab
if _rc {
    di as error "Package estout not installed. Run: ssc install estout"
}
esttab m1 m2 using "output/trust_models_stata.tex", replace se star(* 0.10 ** 0.05 *** 0.01)
