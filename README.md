# Code for reproducing simulations and analyses

## Library files

The NIMBLE R package is necessary to run these scripts, and is freely available from r-nimble.org.
The script sim-cred-bands.R contains methods for computing quantile-based simultaneous credible bands from posterior samples, and simulation-functions.R contains helper methods used in the simulation study script.

## Simulations

To reproduce the simulation study, run simulation.R.

## Analyses

The paper includes two analyses run from separate files.
The "standard" analysis, which tests, in sequence, the overall treatment effect, treatment-covariate interactions, and subgroups defined by significant interactions, is run from the following files, respectively:
   1. example-standard-overall.R
   2. example-standard-interaction.R
   3. example-standard-subset.R
The main analysis is run via example.R, and the interaction prior specification may be switched among shrinkage (standard analysis), flat (prior sensitivity analysis), and spike-and-slab (variable selection analysis) by (un)commenting the so-labeled code blocks in the model specification.
