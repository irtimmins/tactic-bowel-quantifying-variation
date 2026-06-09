/*=============================================================================
  04 - Interaction: change of trust x route to diagnosis
  -----------------------------------------------------------------------------
  Tests whether the effect of being diagnosed and treated at different trusts
  on waiting time depends on the route to diagnosis. Fits a linear model with
  the two-way interaction for each waiting-time outcome, reports the joint test
  of the interaction, and the adjusted means for every change-trust x route
  combination.

  Output is written as a readable log and two tidy CSVs: the adjusted cell
  means (for plotting the interaction) and the joint interaction test per
  outcome.

  Input  (in $work):  colon_analysis.dta
  Outputs (in $out):  colon_interaction.txt              readable log
                      colon_interaction_means.csv        adjusted means by cell
                      colon_interaction_tests.csv        joint interaction tests
=============================================================================*/

clear all
set more off

use "$work/colon_analysis.dta", clear

local wtvars wt_dx_to_dtt wt_dtt_to_tx wt_dx_to_tx

* postfiles for the tidy outputs
tempname pmeans ptests
tempfile means_file tests_file
postfile `pmeans' str16 measure double change_trust str40 change_trust_lbl ///
    double route str40 route_lbl double adj_mean double se double ll double ul ///
    using `means_file', replace
postfile `ptests' str16 measure double F double df_num double df_den double p ///
    using `tests_file', replace

log using "$out/colon_interaction.txt", text replace

foreach y of local wtvars {

    display _n(2) "==== Outcome: `y' ===="

    * main effects plus the change-trust x route interaction,
    * lightly adjusted for age and comorbidity
    regress `y' i.change_trust##i.route ib5.age_group i.rcs_ch_score

    * joint test that the interaction adds nothing
    display _n "Joint test of the change_trust x route interaction:"
    testparm i.change_trust#i.route
    post `ptests' ("`y'") (r(F)) (r(df)) (r(df_r)) (r(p))

    * adjusted mean wait for each combination
    display _n "Adjusted means by change_trust and route:"
    margins change_trust#route

    * effect of changing trust within each route (for the log)
    display _n "Effect of changing trust, within each route:"
    margins route, dydx(change_trust)

    *--------------------------------------------------------------------------
    * tidy adjusted means: one predictive margin per change_trust x route cell
    *--------------------------------------------------------------------------
    quietly levelsof change_trust, local(ct_levels)
    quietly levelsof route,        local(rt_levels)

    foreach ct of local ct_levels {
        local ct_lbl : label (change_trust) `ct'
        foreach rt of local rt_levels {
            local rt_lbl : label (route) `rt'

            quietly margins, at(change_trust = `ct' route = `rt')
            matrix b = r(b)
            matrix V = r(V)
            local est = b[1,1]
            local se  = sqrt(V[1,1])
            local ll  = `est' - 1.96 * `se'
            local ul  = `est' + 1.96 * `se'

            post `pmeans' ("`y'") (`ct') ("`ct_lbl'") (`rt') ("`rt_lbl'") ///
                (`est') (`se') (`ll') (`ul')
        }
    }
}

log close
postclose `pmeans'
postclose `ptests'

* export the tidy tables
use `means_file', clear
foreach v in adj_mean se ll ul {
    replace `v' = round(`v', 0.01)
}
order measure change_trust change_trust_lbl route route_lbl adj_mean se ll ul
export delimited using "$out/colon_interaction_means.csv", replace quote

use `tests_file', clear
replace F = round(F, 0.001)
replace p = round(p, 0.0001)
order measure F df_num df_den p
export delimited using "$out/colon_interaction_tests.csv", replace quote

display "Step 04 done: log + tidy CSVs written to $out."
