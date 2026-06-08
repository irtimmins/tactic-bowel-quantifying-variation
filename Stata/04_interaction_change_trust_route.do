/*=============================================================================
  04 - Interaction: change of trust x route to diagnosis
  -----------------------------------------------------------------------------
  Tests whether the effect of being diagnosed and treated at different trusts
  on waiting time depends on the route to diagnosis. Fits a linear model with
  the two-way interaction for each waiting-time outcome, reports the joint test
  of the interaction, and the adjusted means for every change-trust x route
  combination.

  Input  (in $work):  colon_analysis.dta
  Output (in $out):   colon_interaction.txt
=============================================================================*/

clear all
set more off

use "$work/colon_analysis.dta", clear

log using "$out/colon_interaction.txt", text replace

local wtvars wt_dx_to_dtt wt_dtt_to_tx wt_dx_to_tx

foreach y of local wtvars {

    display _n(2) "==== Outcome: `y' ===="

    * main effects plus the change-trust x route interaction,
    * lightly adjusted for age and comorbidity
    regress `y' i.change_trust##i.route ib5.age_group i.rcs_ch_score

    * joint test that the interaction adds nothing
    display _n "Joint test of the change_trust x route interaction:"
    testparm i.change_trust#i.route

    * adjusted mean wait for each combination
    display _n "Adjusted means by change_trust and route:"
    margins change_trust#route

    * difference made by changing trust, within each route
    display _n "Effect of changing trust, within each route:"
    margins route, dydx(change_trust)
}

log close
display "Step 04 done: interaction results written to $out/colon_interaction.txt"
