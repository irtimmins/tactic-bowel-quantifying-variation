/*=============================================================================
  01 - Cancer waiting times linkage (the 1b merge)
  -----------------------------------------------------------------------------
  Starts from the registry+HES backbone (Table A) and merges the raw CWT
  records (Table B), reproducing the real linkage so the analytic cohort and
  the exclusion funnel can be built downstream.

  Inputs  (in $syn):
    colon_ncras_hes_synthetic.dta   one row per patient, NCRAS registry + the
                                    linked HES resection (the LEFT side)
    colon_cwt_records_synthetic.dta raw CWT records, several rows per patient,
                                    dates held as dd/mm/YYYY strings (RIGHT side)
  Outputs (in $work):
    colon_cohort_analytic.dta       one row per patient, post-linkage
    colon_funnel.dta                step-by-step counts for the flow chart

  ICON -> rapid name mapping, if you later switch extracts:
    pseudo_patientid  <- patient_pseudo_id
    diagmdy           <- date_of_diagnosis
    tx_date           <- date_of_surgery
    diag_trust        <- trust_code_diag
    tx_trust          <- trust_code_surg
    diag_hosp         <- site_code_of_diagnosis
    wt_dx_to_dtt      <- time_diag_dtt
    wt_dtt_to_tx      <- time_dectotreat_treat
    wt_dx_to_tx       <- diag_treat
=============================================================================*/

clear all
set more off

* linkage constants (match the R generator's merge_const)
local diff_max   = 5                      // max days between HES and CWT treat dates
local mod_keep   "01 23 24"               // CWT modalities kept (surgery)
local mod2324    = date("2020-06-01","YMD")

* hold the funnel counts as we go
tempname fn
postfile `fn' str60 step long n_patients long n_rows using "$work/colon_funnel.dta", replace

*------------------------------------------------------------------------------
* Table A: registry + HES backbone
*------------------------------------------------------------------------------
use "$syn/colon_ncras_hes_synthetic.dta", clear

* keep patients with a recorded surgery date (the resection episode)
drop if missing(tx_date)
count
post `fn' ("NCRAS+HES patients with surgery") (r(N)) (r(N))

* one row per patient at this stage
isid pseudo_patientid

tempfile backbone
save `backbone'

*------------------------------------------------------------------------------
* Table B: raw CWT records, parse the string dates and keep surgery modalities
*------------------------------------------------------------------------------
use "$syn/colon_cwt_records_synthetic.dta", clear

gen double dtt_date    = date(treat_period_start, "DMY")
gen double cwt_tx_date = date(treat_start, "DMY")
gen double mdt_date2   = date(mdt_date, "DMY")
drop mdt_date
rename mdt_date2 mdt_date
format dtt_date cwt_tx_date mdt_date %td

* keep surgical modalities; 23/24 only valid from mid-2020 onwards
gen byte keep_mod = 0
foreach m of local mod_keep {
    replace keep_mod = 1 if modality == "`m'"
}
replace keep_mod = 0 if inlist(modality,"23","24") & cwt_tx_date < `mod2324'
keep if keep_mod == 1
drop keep_mod

tempfile cwt
save `cwt'

*------------------------------------------------------------------------------
* Merge backbone to CWT (1:m: a patient may have several CWT records)
*------------------------------------------------------------------------------
use `backbone', clear
merge 1:m pseudo_patientid using `cwt'

* _merge==1 : surgery patient with no CWT record  -> missing_CWT
* _merge==2 : CWT record with no surgery backbone  -> not part of the cohort
gen byte missing_CWT = (_merge == 1)
drop if _merge == 2
drop _merge

* records with at least one linkable CWT treatment date
gen double tx_date_diff = abs(tx_date - cwt_tx_date)
preserve
    keep if !missing(tx_date_diff)
    count
    local np = r(N)
    bysort pseudo_patientid: keep if _n == 1
    count
    post `fn' ("with a linkable CWT record") (r(N)) (`np')
restore

*------------------------------------------------------------------------------
* Treatment-date agreement between HES and CWT
*------------------------------------------------------------------------------
* signed difference, mirrors Stata diff_cwt_cr_treat = date_of_surgery - treat_start
gen long diff_cwt_cr_treat = tx_date - cwt_tx_date

gen byte diff_cwt_cr_treat_cat = .
replace diff_cwt_cr_treat_cat = 1 if diff_cwt_cr_treat == 0
replace diff_cwt_cr_treat_cat = 2 if inrange(diff_cwt_cr_treat,  1,  4)
replace diff_cwt_cr_treat_cat = 2 if inrange(diff_cwt_cr_treat, -4, -1)
replace diff_cwt_cr_treat_cat = 3 if inrange(diff_cwt_cr_treat,  6, 29)
replace diff_cwt_cr_treat_cat = 3 if inrange(diff_cwt_cr_treat,-29, -6)
replace diff_cwt_cr_treat_cat = 4 if diff_cwt_cr_treat >  30 & !missing(diff_cwt_cr_treat)
replace diff_cwt_cr_treat_cat = 4 if diff_cwt_cr_treat < -30
replace diff_cwt_cr_treat_cat = 5 if missing(cwt_tx_date)
replace diff_cwt_cr_treat_cat = 5 if missing(diff_cwt_cr_treat_cat)

label define diffcat 1 "exact" 2 "within 5d" 3 "5-30d" 4 ">30d" 5 "missing/edge"
label values diff_cwt_cr_treat_cat diffcat

* keep records where the two sources agree to within the tolerance
keep if !missing(tx_date_diff) & tx_date_diff <= `diff_max'
preserve
    bysort pseudo_patientid: keep if _n == 1
    count
    post `fn' ("treatment dates agree (<= `diff_max'd)") (r(N)) (.)
restore

*------------------------------------------------------------------------------
* Waiting-time outcomes (days)
*------------------------------------------------------------------------------
gen double wt_dx_to_dtt = dtt_date - diagmdy   // diagnosis  -> decision to treat
gen double wt_dtt_to_tx = tx_date  - dtt_date  // decision to treat -> surgery
gen double wt_dx_to_tx  = tx_date  - diagmdy   // diagnosis  -> surgery

label variable wt_dx_to_dtt "Days diagnosis to decision-to-treat"
label variable wt_dtt_to_tx "Days decision-to-treat to surgery"
label variable wt_dx_to_tx  "Days diagnosis to surgery"

* tumour site agreement between registry and CWT (used to break ties)
gen byte site_match = (substr(site_icd10,1,3) == substr(sitestr,1,3))

* non-negative total wait
keep if wt_dx_to_tx >= 0
preserve
    bysort pseudo_patientid: keep if _n == 1
    count
    post `fn' ("non-negative total wait") (r(N)) (.)
restore

*------------------------------------------------------------------------------
* Reduce to one record per patient: prefer a site match, then the closest
* treatment-date agreement, then the earliest CWT dates
*------------------------------------------------------------------------------
gsort pseudo_patientid -site_match tx_date_diff dtt_date cwt_tx_date
bysort pseudo_patientid: keep if _n == 1
isid pseudo_patientid

count
post `fn' ("deduplicated to one row per patient") (r(N)) (r(N))
postclose `fn'

* recompute missing_CWT cleanly: every retained patient has a CWT record here
replace missing_CWT = 0

compress
save "$work/colon_cohort_analytic.dta", replace

display "Step 01 done: analytic cohort saved with " _N " patients."
