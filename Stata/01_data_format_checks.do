/*=============================================================================
  01 - Check input data formats
  -----------------------------------------------------------------------------
  Pre-flight check that the two pipeline inputs are in the shape the rest of
  the scripts expect. Reports every problem it finds (it does not stop at the
  first), then prints a final verdict. Run this before 01.

  Checks:
    - both files exist
    - all variables the pipeline reads are present, with the right storage type
    - the patient id is a string in both files and unique in the backbone
    - the registry date fields are numeric Stata dates in a plausible range
    - the CWT date strings parse as dd/mm/YYYY
    - the CWT modality codes are strings holding the expected values
    - the merge key overlaps between the two files

  Inputs (in $syn):
    colon_ncras_hes_synthetic.dta     Table A, registry + HES backbone
    colon_cwt_records_synthetic.dta   Table B, raw CWT records
=============================================================================*/

clear all
set more off

* set $syn if running this file on its own
if "$syn" == "" global syn "Data/synthetic"

global n_err  = 0
global n_warn = 0

* helper: confirm a variable exists and (optionally) has the expected type
capture program drop chkvar
program define chkvar
    args fname vname vtype
    capture confirm variable `vname'
    if _rc {
        display as error "  FAIL [`fname']: required variable not found - `vname'"
        global n_err = $n_err + 1
        exit
    }
    if "`vtype'" == "string" {
        capture confirm string variable `vname'
        if _rc {
            display as error "  FAIL [`fname']: `vname' should be string but is numeric"
            global n_err = $n_err + 1
        }
    }
    if "`vtype'" == "numeric" {
        capture confirm numeric variable `vname'
        if _rc {
            display as error "  FAIL [`fname']: `vname' should be numeric but is string"
            global n_err = $n_err + 1
        }
    }
    display "  ok [`fname']: `vname'"
end

*==============================================================================
* Table A: registry + HES backbone
*==============================================================================
display _n "{hline 70}"
display "Checking Table A: colon_ncras_hes_synthetic.dta"
display "{hline 70}"

capture confirm file "$syn/colon_ncras_hes_synthetic.dta"
if _rc {
    display as error "  FAIL: file not found in $syn"
    global n_err = $n_err + 1
}
else {
    use "$syn/colon_ncras_hes_synthetic.dta", clear

    * merge key and dates
    chkvar "A" pseudo_patientid string
    chkvar "A" diagmdy          numeric
    chkvar "A" tx_date          numeric

    * variables used downstream (analysis and exclusions)
    chkvar "A" ydiag                            numeric
    chkvar "A" agediag                          numeric
    chkvar "A" sex                              any
    chkvar "A" stage                            string
    chkvar "A" route_combined                   string
    chkvar "A" change_trust                     numeric
    chkvar "A" rcs_ch_score                     numeric
    chkvar "A" cci_group                        any
    chkvar "A" perf_status                      numeric
    chkvar "A" ethnicity_group_broad            string
    chkvar "A" NHSE_reversed_imd_quintile_lsoas numeric
    chkvar "A" diag_trust                       string
    chkvar "A" diag_hosp                        string
    chkvar "A" tx_trust                         string
    chkvar "A" tx_hosp                          string
    chkvar "A" sitestr                          string
    chkvar "A" emergency                        numeric
    chkvar "A" dco                              numeric
    chkvar "A" any_mets                         numeric
    chkvar "A" tnm_m                            numeric
    chkvar "A" pretreat_m                       numeric

    * one row per patient
    capture isid pseudo_patientid
    if _rc {
        display as error "  FAIL [A]: pseudo_patientid is not unique (expected one row per patient)"
        global n_err = $n_err + 1
    }

    * registry dates should be real dates in a sensible range
    capture confirm numeric variable diagmdy
    if !_rc {
        quietly summarize diagmdy
        local dmin = r(min)
        local dmax = r(max)
        display "  diagmdy range: " %td `dmin' " to " %td `dmax'
        if `dmin' < td(01jan2000) | `dmax' > td(31dec2035) {
            display as error "  WARN [A]: diagmdy outside 2000-2035 - is it a true Stata date?"
            global n_warn = $n_warn + 1
        }
    }
    capture confirm numeric variable tx_date
    if !_rc {
        quietly summarize tx_date
        display "  tx_date range: " %td r(min) " to " %td r(max)
    }
}

*==============================================================================
* Table B: raw CWT records
*==============================================================================
display _n "{hline 70}"
display "Checking Table B: colon_cwt_records_synthetic.dta"
display "{hline 70}"

capture confirm file "$syn/colon_cwt_records_synthetic.dta"
if _rc {
    display as error "  FAIL: file not found in $syn"
    global n_err = $n_err + 1
}
else {
    use "$syn/colon_cwt_records_synthetic.dta", clear

    chkvar "B" pseudo_patientid   string
    chkvar "B" treat_period_start string
    chkvar "B" treat_start        string
    chkvar "B" mdt_date           string
    chkvar "B" modality           string
    chkvar "B" site_icd10         string
    chkvar "B" org_treat_start    string

    * CWT date strings must parse as dd/mm/YYYY
    foreach d in treat_period_start treat_start mdt_date {
        capture confirm string variable `d'
        if !_rc {
            quietly count if !missing(`d')
            local n_nonmiss = r(N)
            quietly gen double _chk = date(`d', "DMY")
            quietly count if missing(_chk) & !missing(`d')
            local n_bad = r(N)
            drop _chk
            if `n_bad' > 0 {
                display as error "  FAIL [B]: `n_bad' of `n_nonmiss' `d' values do not parse as dd/mm/YYYY"
                global n_err = $n_err + 1
            }
            else {
                display "  `d': all non-missing values parse as dd/mm/YYYY"
            }
        }
    }

    * modality should hold the surgical codes the merge keeps
    capture confirm string variable modality
    if !_rc {
        display "  modality values (counts):"
        tab modality

        local found_keep = 0
        foreach m in 01 23 24 {
            quietly count if modality == "`m'"
            if r(N) > 0 local found_keep = 1
        }
        if !`found_keep' {
            display as error "  WARN [B]: none of the kept modalities (01, 23, 24) found - check coding"
            global n_warn = $n_warn + 1
        }
    }
}

*==============================================================================
* Cross-file: merge key compatibility and overlap
*==============================================================================
display _n "{hline 70}"
display "Checking the merge key across the two files"
display "{hline 70}"

capture confirm file "$syn/colon_ncras_hes_synthetic.dta"
local haveA = (_rc == 0)
capture confirm file "$syn/colon_cwt_records_synthetic.dta"
local haveB = (_rc == 0)

if `haveA' & `haveB' {
    use pseudo_patientid using "$syn/colon_ncras_hes_synthetic.dta", clear
    quietly count
    local nA = r(N)
    tempfile akeys
    quietly save `akeys'

    use pseudo_patientid using "$syn/colon_cwt_records_synthetic.dta", clear
    quietly bysort pseudo_patientid: keep if _n == 1
    quietly count
    local nB = r(N)

    merge 1:1 pseudo_patientid using `akeys'
    quietly count if _merge == 3
    local n_match = r(N)
    quietly count if _merge == 2
    local n_aonly = r(N)

    display "  Table A patients: `nA'"
    display "  Table B distinct patients: `nB'"
    display "  patients in both files: `n_match'"
    display "  Table A patients with no CWT record: `n_aonly'"

    if `n_match' == 0 {
        display as error "  FAIL: no patient ids overlap - the merge would return nothing"
        global n_err = $n_err + 1
    }
}
else {
    display as error "  skipped (one or both files missing)"
}

*==============================================================================
* Summary
*==============================================================================
display _n "{hline 70}"
if $n_err == 0 & $n_warn == 0 {
    display "Input check passed: no problems found. Safe to run 01."
}
else {
    display "Input check finished with " $n_err " error(s) and " $n_warn " warning(s)."
    if $n_err > 0 display as error "Fix the errors above before running 01."
}
display "{hline 70}"
