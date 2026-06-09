/*=============================================================================
  Colon waiting times - master script
  -----------------------------------------------------------------------------
  Runs the full pipeline on the synthetic test data. Each numbered step is a
  self-contained do-file that reads what it needs and writes what the next step
  expects, so they can also be run one at a time while developing.

  Variable names follow the ICON convention (the older comparable data). The
  newer "rapid" data uses different names; see the mapping note in 01 if you
  later want to run the same pipeline on the rapid extract.

  Edit the three paths below, then run this file.
=============================================================================*/

clear all
set more off
version 16

/* folders ------------------------------------------------------------------ */
* root path - update this if the project moves
global root "D:/Projects/#2045_ICON_TACTIC/Project1_interim_bowel/tactic-bowel-quantifying-variation"
cd "$root"
* where the synthetic .dta files from the R generator live
global syn   "Data/synthetic"
* where to put derived datasets and the analytic cohort
global work  "Data/temp"
* where to put tables, figures and logs
global out   "Output"
* the provider characteristics workbook (synthetic, from the R helper)
global provxlsx "Data/synthetic/NHSHospitals_services_colon_synthetic.xlsx"

cap mkdir "$work"
cap mkdir "$out"

/* run the pipeline --------------------------------------------------------- */
do "Stata/01_merge_cwt.do"
do "Stata/02_flowchart_exclusions.do"
do "Stata/03_patient_characteristics.do"
do "Stata/04_interaction_change_trust_route.do"
do "Stata/05_trust_random_effects.do"
do "Stata/06_trust_random_effects.do"
do "Stata/07_provider_level.do"

display "Pipeline complete."
