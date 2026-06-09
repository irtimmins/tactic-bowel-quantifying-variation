/*=============================================================================
  07 - Provider-level characteristics and waiting times
  -----------------------------------------------------------------------------
  Imports the synthetic hospital-characteristics workbook, builds the site-level
  covariates, merges them onto the analytic cohort by hospital of diagnosis, and
  reports mean waiting time (with its standard error and the number of patients
  and hospitals) by each site characteristic.

  This mirrors the real provider-level script; the synthetic workbook is built
  by the R helper Create_synthetic_provider_level_excel_colon.R and holds only
  the columns used here. The competitor status (net patient flow) is merged in
  if its file is present, otherwise that covariate is skipped.

  Inputs:
    $provxlsx                                   hospital characteristics (xlsx)
    $syn/colon_competitor_status_synthetic.dta  optional: winner/loser status
    $work/colon_analysis.dta                    analytic cohort from step 03
  Output (in $out):
    colon_provider_level.txt
=============================================================================*/

clear all
set more off

*------------------------------------------------------------------------------
* Provider characteristics workbook -> one row per site code
*------------------------------------------------------------------------------
import excel "$provxlsx", sheet("Sheet1") firstrow clear

drop if missing(Trust_Name)
drop if missing(Hospital_site_code)

* exclude flagged trust colours (data-quality exclusions in the real file)
drop if inlist(Trust_Name_colour, "Light Red", "Pink Red", "Orange")

* if a site appears more than once, prefer the row carrying surgery information
capture confirm string variable Bowel_ca_surgery
if !_rc destring Bowel_ca_surgery, replace force
bysort Hospital_site_code: gen byte n_rows = _N
drop if n_rows > 1 & missing(Bowel_ca_surgery)
drop n_rows

rename Hospital_site_code diag_hosp

* binary site characteristics
capture confirm string variable Comprehensive_centre
if !_rc destring Comprehensive_centre, replace force
capture confirm string variable Teaching_hospitals
if !_rc destring Teaching_hospitals, replace force

gen byte comprehensive = (Comprehensive_centre == 1)
label define comp 0 "Non-comprehensive" 1 "Comprehensive"
label values comprehensive comp

gen byte teaching = (Teaching_hospitals == 1)
label define teach 0 "Non-teaching" 1 "Teaching"
label values teaching teach

* CQC rating as an ordered score
gen byte cqc_rating = .
replace cqc_rating = 1 if Latest_Rating == "Inadequate"
replace cqc_rating = 2 if Latest_Rating == "Requires Improvement"
replace cqc_rating = 3 if Latest_Rating == "Good"
replace cqc_rating = 4 if Latest_Rating == "Outstanding"
label define cqc 1 "Inadequate" 2 "Requires Improvement" 3 "Good" 4 "Outstanding"
label values cqc_rating cqc

* staff engagement and morale as quintiles
capture confirm string variable Staff_engagement
if !_rc destring Staff_engagement, replace force
capture confirm string variable Moral
if !_rc destring Moral, replace force
xtile staff_eng_cat = Staff_engagement, nquantiles(5)
xtile moral_cat     = Moral,            nquantiles(5)
label define quint 1 "Q1 (lowest)" 2 "Q2" 3 "Q3" 4 "Q4" 5 "Q5 (highest)"
label values staff_eng_cat quint
label values moral_cat     quint

* bed occupancy: high vs normal at the 95% threshold
capture confirm string variable mean
if !_rc destring mean, replace force
gen byte bed_occ_cat = .
replace bed_occ_cat = 0 if mean <  0.95 & !missing(mean)
replace bed_occ_cat = 1 if mean >= 0.95 & !missing(mean)
label define bedocc 0 "Normal (<95%)" 1 "High (>=95%)"
label values bed_occ_cat bedocc

keep diag_hosp comprehensive teaching cqc_rating staff_eng_cat moral_cat bed_occ_cat
tempfile provchar
save `provchar'

*------------------------------------------------------------------------------
* Optional: net patient flow status, for the diagnosing and treating hospital
*------------------------------------------------------------------------------
local have_diag = 0
local have_tx   = 0
capture confirm file "$work/colon_competitor_status_diag.dta"
if !_rc local have_diag = 1
capture confirm file "$work/colon_competitor_status_tx.dta"
if !_rc local have_tx = 1

*------------------------------------------------------------------------------
* Merge characteristics onto the analytic cohort
*------------------------------------------------------------------------------
use "$work/colon_analysis.dta", clear

merge m:1 diag_hosp using `provchar', keep(master match) nogenerate

* diagnosing-hospital flow status keys on diag_hosp
if `have_diag' {
    merge m:1 diag_hosp using "$work/colon_competitor_status_diag.dta", ///
        keepusing(competitor_status_diag) keep(master match) nogenerate
}
* treating-hospital flow status keys on tx_hosp
if `have_tx' {
    merge m:1 tx_hosp using "$work/colon_competitor_status_tx.dta", ///
        keepusing(competitor_status_tx) keep(master match) nogenerate
}

*------------------------------------------------------------------------------
* Mean waiting time by each site characteristic
*  reports mean, SE, patient N and hospital N per group, as a readable log and
*  a tidy long CSV (one row per measure x characteristic x level)
*------------------------------------------------------------------------------
local wtvars wt_dx_to_dtt wt_dtt_to_tx wt_dx_to_tx
local covars comprehensive teaching cqc_rating staff_eng_cat moral_cat bed_occ_cat
if `have_diag' local covars `covars' competitor_status_diag
if `have_tx'   local covars `covars' competitor_status_tx

* tidy output postfile
tempname pf
tempfile results
postfile `pf' str16 measure str24 covariate double level_id str48 level_name ///
    double n_patients double mean double se double n_hospitals ///
    using `results', replace

log using "$out/colon_provider_level.txt", text replace

foreach y of local wtvars {
    foreach c of local covars {

        * distinct-hospital count keys on the role the covariate refers to
        local hospvar diag_hosp
        if "`c'" == "competitor_status_tx" local hospvar tx_hosp

        display _n "Outcome: `y'  |  characteristic: `c'"
        tabstat `y', by(`c') stats(mean semean n) nototal

        display "Hospitals per group (`hospvar'):"
        preserve
            bysort `c' `hospvar': keep if _n == 1
            gen byte one = 1
            tabstat one, by(`c') stats(n) nototal
        restore

        * tidy rows: loop over the levels of this characteristic
        quietly levelsof `c', local(lvls)
        foreach l of local lvls {
            local lname : label (`c') `l'
            if "`lname'" == "" local lname "`l'"

            quietly summarize `y' if `c' == `l'
            local m  = r(mean)
            local se = r(sd) / sqrt(r(N))
            local np = r(N)

            * distinct hospitals contributing to this level
            quietly bysort `hospvar': gen byte _firsth = (_n == 1) if `c' == `l'
            quietly count if _firsth == 1
            local nh = r(N)
            drop _firsth

            post `pf' ("`y'") ("`c'") (`l') ("`lname'") (`np') (`m') (`se') (`nh')
        }
    }
}

log close
postclose `pf'

* export the tidy table
use `results', clear
replace mean = round(mean, 0.01)
replace se   = round(se,   0.01)
order measure covariate level_id level_name n_patients mean se n_hospitals
export delimited using "$out/colon_provider_level.csv", replace quote

display "Step done: provider-level log + tidy CSV written to $out."
