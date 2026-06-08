/*=============================================================================
  02 - Exclusion flow chart
  -----------------------------------------------------------------------------
  Renders the patient-exclusion funnel recorded during the linkage, and adds
  the downstream cohort-definition exclusions (stage 4 / metastases, emergency
  or DCO route, waiting time outside the analysis window). The numbers can be
  read straight into a flow-chart diagram.

  Inputs  (in $work):
    colon_funnel.dta              linkage funnel from step 01
    colon_cohort_analytic.dta     post-linkage cohort
  Output  (in $out):
    colon_flowchart.txt           text summary of the funnel and the drops
=============================================================================*/

clear all
set more off

log using "$out/colon_flowchart.txt", text replace

*------------------------------------------------------------------------------
* Part 1: linkage funnel (from step 01)
*------------------------------------------------------------------------------
display _n "Linkage funnel (registry+HES -> CWT)"
use "$work/colon_funnel.dta", clear

* drop at each step = previous N minus current N
gen long dropped = n_patients[_n-1] - n_patients
replace dropped = . in 1
list step n_patients dropped, noobs sep(0) abbrev(40)

*------------------------------------------------------------------------------
* Part 2: cohort-definition exclusions on the linked cohort
*------------------------------------------------------------------------------
use "$work/colon_cohort_analytic.dta", clear

* a small helper that reports N before and after a drop condition
capture program drop dropstep
program define dropstep
    args label cond
    quietly count
    local before = r(N)
    quietly drop if `cond'
    quietly count
    local after = r(N)
    display as text %-44s "`label'" ///
        as result "  kept " %7.0f `after' "   (dropped " %6.0f (`before'-`after') ")"
end

display _n "Cohort-definition exclusions"
quietly count
display as text %-44s "linked cohort entering definition" ///
    as result "  kept " %7.0f r(N)

* metastatic / stage 4 disease (synthetic cohort is stage 1-3, so these are 0)
dropstep "stage 4 disease"            "stage == \"4\""
dropstep "any recorded metastases"    "any_mets == 1"
dropstep "M1 at staging"              "tnm_m == 1"
dropstep "M1 pre-treatment"           "pretreat_m == 1"

* non-elective presentation
dropstep "emergency presentation"     "emergency == 1"
dropstep "death certificate only"     "dco == 1"

* unreliable waiting time (date disagreement beyond 5 days, kept as cat 1/2 only)
dropstep "CWT/HES dates disagree"     "!inlist(diff_cwt_cr_treat_cat,1,2)"

* analysis window on the primary outcome
dropstep "wait outside (0,180] days"  "!(wt_dx_to_dtt > 0 & wt_dx_to_dtt <= 180)"

quietly count
display _n as text "Final analytic sample: " as result r(N)

log close
display "Step 02 done: flow chart written to $out/colon_flowchart.txt"
