/*=============================================================================
  03 - Patient characteristics and waiting times
  -----------------------------------------------------------------------------
  Builds the analysis variables (age band, labelled factors, the analytic
  subset) used by every later step, then produces a characteristics table:
  for each level of age, sex, comorbidity, ethnicity and performance status it
  reports the patient count and the mean and SD of each of the three waiting
  time measures.

  Inputs  (in $work):  colon_cohort_analytic.dta
  Outputs (in $work):  colon_analysis.dta            analysis-ready cohort
          (in $out):   colon_characteristics.txt     the table
=============================================================================*/

clear all
set more off

use "$work/colon_cohort_analytic.dta", clear

*------------------------------------------------------------------------------
* Analysis variables
*------------------------------------------------------------------------------
* age band, reference 70-74 (matches the real analysis grouping)
egen age_group = cut(agediag), at(0,50,55,60,65,70,75,80,85,120) icode
label define age_group 0 "<50" 1 "50-54" 2 "55-59" 3 "60-64" 4 "65-69" ///
                       5 "70-74" 6 "75-79" 7 "80-84" 8 "85+"
label values age_group age_group

* comorbidity: numeric Charlson score (0,1,2,3+)
label define rcs 0 "0" 1 "1" 2 "2" 3 "3+"
capture label values rcs_ch_score rcs

* performance status
label define ps 0 "0" 1 "1" 2 "2" 3 "3" 4 "4"
capture label values perf_status ps

* sex
label define sexlab 1 "Male" 2 "Female"
capture confirm string variable sex
if _rc {
    capture label values sex sexlab
}

* route of referral and change of treating trust used later
encode route_combined, gen(route)
label variable route "Route to diagnosis"
label variable change_trust "Diagnosed and treated at different trusts"

* restrict to the analytic subset: reliable dates and a valid in-window wait
keep if inlist(diff_cwt_cr_treat_cat,1,2)
keep if wt_dx_to_dtt > 0 & wt_dx_to_dtt <= 180
keep if wt_dtt_to_tx >= 0 & wt_dx_to_tx <= 180

compress
save "$work/colon_analysis.dta", replace

*------------------------------------------------------------------------------
* Characteristics table: N and mean (SD) of each waiting time by covariate
*------------------------------------------------------------------------------
log using "$out/colon_characteristics.txt", text replace

local wtvars wt_dx_to_dtt wt_dtt_to_tx wt_dx_to_tx
local covars age_group sex rcs_ch_score ethnicity_group_broad perf_status

display _n "Overall"
tabstat `wtvars', stats(n mean sd) columns(statistics)

foreach c of local covars {
    display _n "By `c'"
    tabstat `wtvars', by(`c') stats(n mean sd) nototal columns(statistics)
}

log close
display "Step 03 done: characteristics written to $out/colon_characteristics.txt"
