/*=============================================================================
  04 - Patient characteristics and waiting times
  -----------------------------------------------------------------------------
  Builds the analysis variables (age band, labelled factors, the analytic
  subset) used by every later step, then produces a characteristics table:
  for each level of age, sex, comorbidity, ethnicity and performance status it
  reports the patient count and the mean and SD of each of the three waiting
  time measures.

  Output is written two ways: a human-readable log, and a tidy long-format CSV
  (one row per measure x covariate x level) for reading straight into R/ggplot.

  Inputs  (in $work):  colon_cohort_analytic.dta
  Outputs (in $work):  colon_analysis.dta             analysis-ready cohort
          (in $out):   colon_characteristics.txt      readable log
                       colon_characteristics.csv      tidy table
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

local wtvars wt_dx_to_dtt wt_dtt_to_tx wt_dx_to_tx
local covars age_group sex rcs_ch_score ethnicity_group_broad perf_status

*------------------------------------------------------------------------------
* Readable log
*------------------------------------------------------------------------------
log using "$out/colon_characteristics.txt", text replace

display _n "Overall"
tabstat `wtvars', stats(n mean sd) columns(statistics)

foreach c of local covars {
    display _n "By `c'"
    tabstat `wtvars', by(`c') stats(n mean sd) nototal columns(statistics)
}

log close

*------------------------------------------------------------------------------
* Tidy long-format table for plotting
*  one row per: measure (waiting time) x covariate x level
*------------------------------------------------------------------------------
tempname pf
tempfile results
postfile `pf' str16 measure str24 covariate double level_id str80 level_name ///
    double n double mean double sd using `results', replace

* overall rows (one per waiting-time measure)
foreach wt of local wtvars {
    quietly summarize `wt'
    post `pf' ("`wt'") ("Overall") (.) ("All patients") (r(N)) (r(mean)) (r(sd))
}

* rows for each covariate level
foreach c of local covars {

    * detect whether the covariate is a string (different comparison/label logic)
    capture confirm string variable `c'
    local is_string = (_rc == 0)

    quietly levelsof `c', local(lvls)
    foreach l of local lvls {

        if `is_string' {
            local cond `c' == "`l'"
            local lname "`l'"
            local lid = .
        }
        else {
            local cond `c' == `l'
            local lname : label (`c') `l'
            local lid = `l'
        }

        foreach wt of local wtvars {
            quietly summarize `wt' if `cond'
            post `pf' ("`wt'") ("`c'") (`lid') ("`lname'") (r(N)) (r(mean)) (r(sd))
        }
    }
}
postclose `pf'

* tidy up and export
use `results', clear
replace mean = round(mean, 0.01)
replace sd   = round(sd,   0.01)
order measure covariate level_id level_name n mean sd
export delimited using "$out/colon_characteristics.csv", replace quote

display "Step 04 done: log + tidy CSV written to $out."
