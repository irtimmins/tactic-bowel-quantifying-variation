/*=============================================================================
  06 - Trust-level effects (fixed and random)
  -----------------------------------------------------------------------------
  Estimates between-trust variation in waiting time. For each outcome, within
  each comorbidity stratum, four models are fitted:
    (a) fixed effects, unadjusted        trust dummies only
    (b) fixed effects, risk-adjusted     + age (+ comorbidity overall)
    (c) random effects, unadjusted       trust random intercept
    (d) random effects, risk-adjusted    + age (+ comorbidity overall)
  The random-effects models give empirical Bayes (shrunken) trust estimates.

  The analysis runs three times: overall, comorbidity 0/1, and comorbidity 2+.
  Within the stratified models, risk adjustment is for age only.

  Outputs as a readable log, a tidy per-trust CSV (all four estimates, every
  outcome and stratum, for caterpillar plots) and a model-summary CSV (ICC and
  the fixed-effects joint tests).

  Input  (in $work):  colon_analysis.dta
  Outputs (in $out):  colon_trust_effects.txt
                      colon_trust_estimates.csv      per-trust, tidy
                      colon_trust_model_summary.csv  ICC + FE joint tests
          (in $work): colon_trust_eb_<outcome>_<stratum>.dta
=============================================================================*/

clear all
set more off

local wtvars wt_dx_to_dtt wt_dtt_to_tx wt_dx_to_tx

* model-summary postfile (long format)
tempname psum
tempfile sum_file
postfile `psum' str16 measure str10 stratum str10 model str10 statistic ///
    double value using `sum_file', replace

log using "$out/colon_trust_effects.txt", text replace

* three strata: overall, comorbidity 0/1, comorbidity 2+
foreach stratum in overall cci01 cci2plus {

    use "$work/colon_analysis.dta", clear

    if "`stratum'" == "cci01" {
        keep if inlist(rcs_ch_score, 0, 1)
        local strat_label "comorbidity 0/1"
        local adj_covars  "ib5.age_group"
    }
    else if "`stratum'" == "cci2plus" {
        keep if rcs_ch_score >= 2 & !missing(rcs_ch_score)
        local strat_label "comorbidity 2+"
        local adj_covars  "ib5.age_group"
    }
    else {
        local strat_label "overall"
        local adj_covars  "ib5.age_group i.rcs_ch_score"
    }

    * encode trust after filtering so unused levels are dropped
    encode diag_trust, gen(trust_id)

    quietly count
    display _n(3) "====================================================="
    display        "Stratum: `strat_label'  (N = " r(N) ")"
    display        "====================================================="

    foreach y of local wtvars {

        display _n(2) "---- Outcome: `y'  |  stratum: `strat_label' ----"

        *----------------------------------------------------------------------
        * (a) fixed effects model, unadjusted (trust dummies only)
        *----------------------------------------------------------------------
        display _n "(a) Fixed effects model, unadjusted"
        regress `y' ib1.trust_id
        testparm i.trust_id
        post `psum' ("`y'") ("`stratum'") ("FE_unadj") ("F") (r(F))
        post `psum' ("`y'") ("`stratum'") ("FE_unadj") ("p") (r(p))

        quietly gen double fe_unadj = 0     // reference trust = 0
        quietly levelsof trust_id, local(tlevs)
        foreach t of local tlevs {
            capture quietly replace fe_unadj = _b[`t'.trust_id] if trust_id == `t'
        }

        *----------------------------------------------------------------------
        * (b) fixed effects model, risk-adjusted
        *----------------------------------------------------------------------
        display _n "(b) Fixed effects model, risk-adjusted (`adj_covars')"
        regress `y' ib1.trust_id `adj_covars'
        testparm i.trust_id
        post `psum' ("`y'") ("`stratum'") ("FE_adj") ("F") (r(F))
        post `psum' ("`y'") ("`stratum'") ("FE_adj") ("p") (r(p))

        quietly gen double fe_adj = 0
        foreach t of local tlevs {
            capture quietly replace fe_adj = _b[`t'.trust_id] if trust_id == `t'
        }

        *----------------------------------------------------------------------
        * (c) random effects model, unadjusted
        *----------------------------------------------------------------------
        display _n "(c) Random effects model, unadjusted"
        mixed `y' || trust_id:, reml
        estat icc
        local icc_u = .
        capture local icc_u = r(icc2)
        post `psum' ("`y'") ("`stratum'") ("RE_unadj") ("icc") (`icc_u')

        predict eb_unadj, reffects
        predict se_unadj, reses

        *----------------------------------------------------------------------
        * (d) random effects model, risk-adjusted
        *----------------------------------------------------------------------
        display _n "(d) Random effects model, risk-adjusted (`adj_covars')"
        mixed `y' `adj_covars' || trust_id:, reml
        estat icc
        local icc_a = .
        capture local icc_a = r(icc2)
        post `psum' ("`y'") ("`stratum'") ("RE_adj") ("icc") (`icc_a')

        predict eb_adj, reffects
        predict se_adj, reses

        *----------------------------------------------------------------------
        * one row per trust: all four estimates, tagged with measure + stratum
        *----------------------------------------------------------------------
        preserve
            bysort trust_id: keep if _n == 1
            keep trust_id diag_trust fe_unadj fe_adj eb_unadj se_unadj eb_adj se_adj
            rename (eb_adj se_adj) (eb_riskadj se_riskadj)
            gen str16 measure = "`y'"
            gen str10 stratum = "`stratum'"
            label variable fe_unadj   "Fixed effect, unadjusted (days)"
            label variable fe_adj     "Fixed effect, risk-adjusted (days)"
            label variable eb_unadj   "EB trust effect, unadjusted (days)"
            label variable eb_riskadj "EB trust effect, risk-adjusted (days)"
            gsort -eb_unadj
            save "$work/colon_trust_eb_`y'_`stratum'.dta", replace

            display _n "FE vs FE adjusted (correlation):"
            correlate fe_unadj fe_adj
            display _n "RE unadjusted vs RE adjusted (correlation):"
            correlate eb_unadj eb_riskadj
            display _n "FE vs RE unadjusted (shrinkage check):"
            correlate fe_unadj eb_unadj
        restore

        drop fe_unadj fe_adj eb_unadj se_unadj eb_adj se_adj
    }
}

log close
postclose `psum'

*------------------------------------------------------------------------------
* Combine all per-trust files into one tidy CSV
*------------------------------------------------------------------------------
clear
local first = 1
foreach stratum in overall cci01 cci2plus {
    foreach y of local wtvars {
        if `first' {
            use "$work/colon_trust_eb_`y'_`stratum'.dta", clear
            local first = 0
        }
        else {
            append using "$work/colon_trust_eb_`y'_`stratum'.dta"
        }
    }
}
foreach v in fe_unadj fe_adj eb_unadj se_unadj eb_riskadj se_riskadj {
    replace `v' = round(`v', 0.001)
}
order measure stratum trust_id diag_trust fe_unadj fe_adj eb_unadj se_unadj eb_riskadj se_riskadj
export delimited using "$out/colon_trust_estimates.csv", replace quote

*------------------------------------------------------------------------------
* Model-summary CSV (ICC and fixed-effects joint tests)
*------------------------------------------------------------------------------
use `sum_file', clear
replace value = round(value, 0.0001)
order measure stratum model statistic value
export delimited using "$out/colon_trust_model_summary.csv", replace quote

display "Step done: trust effects log + tidy CSVs written to $out."
