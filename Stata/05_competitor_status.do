/*=============================================================================
  05 - Competitor status: net patient flow from the distance matrix
  -----------------------------------------------------------------------------
  Classifies each hospital as a Winner (net importer of patients), Loser (net
  exporter) or Insignificant, based on whether the observed patient flow
  differs from what would be expected if patients attended their nearest site.

  The classification is produced twice, against the same nearest-site
  benchmark: once for the diagnosing hospital and once for the treating
  hospital, since a patient can be a flow at either point on the pathway.

  Logic:
    1. For each LSOA, the nearest valid hospital is found from the pairwise
       drive-time matrix (minimum drive time to any valid site).
    2. A patient is "core" if their hospital (diagnosis or treatment) is their
       nearest site. If not, they are a leaver from the nearest site and an
       arriver at the actual site.
    3. Net gain = arrivers - leavers per site. A Poisson z-score
       (net / sqrt(arrivers + leavers)) of magnitude >= 1.96 marks a Winner or
       Loser; otherwise Insignificant.

  For the synthetic test data, lsoa11_code is not in the analytic cohort, so
  the patient-lsoa lookup written by the R helper is merged in first. On real
  data lsoa11_code is already present and that block does nothing.

  Inputs  (in $syn):
    colon_pairwise_distance_matrix_synthetic.dta  lsoa11_code x sitecode x time
    colon_valid_sites_synthetic.dta               valid site codes
    colon_patient_lsoa_synthetic.dta              patient -> lsoa (synthetic only)
    $work/colon_analysis.dta                      analytic cohort
  Outputs (in $work):
    colon_competitor_status_diag.dta   diag_hosp + competitor_status_diag
    colon_competitor_status_tx.dta     tx_hosp   + competitor_status_tx
=============================================================================*/

clear all
set more off

local sig_cut = 1.96          // |z| threshold for winner/loser classification

*------------------------------------------------------------------------------
* Nearest valid site per LSOA from the drive-time matrix (shared benchmark)
*------------------------------------------------------------------------------
use "$syn/colon_pairwise_distance_matrix_synthetic.dta", clear

* keep only valid sites; the valid-sites file keys on diag_hosp, so match on that
rename sitecode diag_hosp
merge m:1 diag_hosp using "$syn/colon_valid_sites_synthetic.dta", keep(match) nogenerate
drop valid
rename diag_hosp sitecode

* nearest site = shortest drive time within each LSOA
bysort lsoa11_code (total_drive_time): keep if _n == 1
rename sitecode nearest_site
keep lsoa11_code nearest_site

tempfile nearest
save `nearest'

*------------------------------------------------------------------------------
* Cohort with residence area, both hospital roles, year and nearest site
*------------------------------------------------------------------------------
use "$work/colon_analysis.dta", clear

* --- synthetic only: merge in the patient-lsoa lookup -----------------------
capture confirm file "$syn/colon_patient_lsoa_synthetic.dta"
if !_rc {
    merge 1:1 pseudo_patientid using "$syn/colon_patient_lsoa_synthetic.dta", ///
        keep(master match) nogenerate
}
* --- end synthetic block -----------------------------------------------------

drop if missing(lsoa11_code)
keep pseudo_patientid lsoa11_code diag_hosp tx_hosp ydiag
merge m:1 lsoa11_code using `nearest', keep(master match) nogenerate

tempfile base
save `base'

*------------------------------------------------------------------------------
* Compute the flow classification for each hospital role
*------------------------------------------------------------------------------
foreach h in diag_hosp tx_hosp {

    * short tag for variable and file names
    if "`h'" == "diag_hosp" local tag diag
    else                    local tag tx

    use `base', clear

    *--- volume threshold on this hospital role: >= 10 patients per year
    *    (>= 5 in the most recent partial year), all years must pass -----------
    bysort `h' ydiag: gen long n_cell = _N
    quietly summarize ydiag
    local max_year = r(max)
    gen int threshold = 10
    replace threshold = 5 if ydiag == `max_year'
    gen byte meets = (n_cell >= threshold)
    bysort `h': egen byte all_ok = min(meets)
    keep if all_ok == 1
    drop n_cell threshold meets all_ok

    *--- core vs non-core against the nearest site ---------------------------
    gen byte core = (`h' == nearest_site) if !missing(nearest_site)

    *--- leavers (by nearest site) and arrivers (by actual site) -------------
    preserve
        keep if core == 0
        rename nearest_site site
        contract site, freq(n_leavers)
        tempfile leavers
        save `leavers'
    restore
    preserve
        keep if core == 0
        rename `h' site
        contract site, freq(n_arrivers)
        tempfile arrivers
        save `arrivers'
    restore

    * all distinct sites in this role, including those with no flow
    contract `h', freq(n_total)
    rename `h' site
    merge 1:1 site using `leavers',  nogenerate
    merge 1:1 site using `arrivers', nogenerate
    replace n_leavers  = 0 if missing(n_leavers)
    replace n_arrivers = 0 if missing(n_arrivers)

    *--- net flow, Poisson z-score and classification ------------------------
    gen long   n_net   = n_arrivers - n_leavers
    gen long   n_flow  = n_arrivers + n_leavers
    gen double z_score = n_net / sqrt(n_flow) if n_flow > 0
    replace    z_score = 0 if n_flow == 0

    gen byte competitor_status_`tag' = 3
    replace  competitor_status_`tag' = 1 if z_score >=  `sig_cut' & !missing(z_score)
    replace  competitor_status_`tag' = 2 if z_score <= -`sig_cut'

    label define compstat 1 "Winner" 2 "Loser" 3 "Insignificant diff.", replace
    label values competitor_status_`tag' compstat
    label variable competitor_status_`tag' ///
        "Net patient flow classification (`tag' hospital)"

    rename site `h'
    rename (n_net n_arrivers n_leavers n_flow z_score) ///
           (n_net_`tag' n_arrivers_`tag' n_leavers_`tag' n_flow_`tag' z_`tag')

    *--- QC ------------------------------------------------------------------
    display _n "Flow classification for `tag' hospital:"
    tab competitor_status_`tag', missing
    quietly summarize n_leavers_`tag'
    local tl = r(sum)
    quietly summarize n_arrivers_`tag'
    display "Total leavers " `tl' " vs arrivers " r(sum) " (should match)"

    keep `h' competitor_status_`tag' n_net_`tag' n_arrivers_`tag' ///
         n_leavers_`tag' n_flow_`tag' z_`tag'
    save "$work/colon_competitor_status_`tag'.dta", replace
    display "Saved competitor status for " _N " `tag' hospitals."
}

display "Step 05 done: diagnosing and treating hospital flow status saved."
