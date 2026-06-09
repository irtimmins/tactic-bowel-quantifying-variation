/*=============================================================================
  05 - Trust-level random effects (empirical Bayes)
  -----------------------------------------------------------------------------
  Estimates between-trust variation in waiting time using a random-intercept
  model on trust of diagnosis, and extracts the empirical Bayes (shrunken)
  trust effects. This is done twice for each outcome:
    (a) unadjusted          - trust intercepts only
    (b) risk-adjusted       - adjusting for age band and comorbidity only
  The shrunken effects and their standard errors are saved so the before/after
  comparison and any caterpillar plot can be drawn.

  Input  (in $work):  colon_analysis.dta
  Outputs (in $out):  colon_trust_effects.txt          model output + ICCs
          (in $work): colon_trust_eb_<outcome>.dta     EB effects per trust
=============================================================================*/

clear all
set more off

use "$work/colon_analysis.dta", clear

* numeric trust identifier for the random effect
encode diag_trust, gen(trust_id)

log using "$out/colon_trust_effects.txt", text replace

local wtvars wt_dx_to_dtt wt_dtt_to_tx wt_dx_to_tx

foreach y of local wtvars {

    display _n(2) "==== Outcome: `y' ===="

    *--------------------------------------------------------------------------
    * (a) unadjusted random-intercept model
    *--------------------------------------------------------------------------
    display _n "Unadjusted model"
    mixed `y' || trust_id:, reml
    estat icc

    * empirical Bayes (shrunken) trust effects and their standard errors
    predict eb_raw, reffects
    predict se_raw, reses

    *--------------------------------------------------------------------------
    * (b) risk-adjusted model (age band and comorbidity only)
    *--------------------------------------------------------------------------
    display _n "Risk-adjusted model (age and comorbidity)"
    mixed `y' ib5.age_group i.rcs_ch_score || trust_id:, reml
    estat icc

    predict eb_adj, reffects
    predict se_adj, reses

    *--------------------------------------------------------------------------
    * one row per trust: shrunken effect before and after adjustment
    *--------------------------------------------------------------------------
    preserve
        bysort trust_id: keep if _n == 1
        keep trust_id diag_trust eb_raw se_raw eb_adj se_adj
        rename (eb_raw se_raw eb_adj se_adj) ///
               (eb_unadj se_unadj eb_riskadj se_riskadj)
        label variable eb_unadj   "EB trust effect, unadjusted (days)"
        label variable eb_riskadj "EB trust effect, risk-adjusted (days)"
        gsort -eb_unadj
        save "$work/colon_trust_eb_`y'.dta", replace

        display _n "Correlation of trust effects before vs after adjustment (`y'):"
        correlate eb_unadj eb_riskadj
    restore

    drop eb_raw se_raw eb_adj se_adj
}

log close
display "Step 06 done: trust effects written to $out and $work."
