# BayesianSIHVRModel
This repository provides code to fit the Bayesian SHIVR Model for estimating the burden of COVID-19 on local healthcare systems described in 'A Bayesian Susceptible-Infectious-Hospitalized-Ventilated-Recovered Model to Predict Demand for COVID-19 Inpatient Care in a Large Healthcare System' Self et al., 2020.

MCMC Code.R is an Rscript file containing the code which implements the Markov chain Monte Carlo algorithm for fitting the Bayesian SHIVR model described in 'A Bayesian Susceptible-Infectious-Hospitalized-Ventilated-Recovered Model to Predict Demand for COVID-19 Inpatient Care in a Large Healthcare System', by Self et. al.

Example Data.csv contains simulated data formatted appropriately for the MCMC Code.R script.

Example Population Data.csv contains the population of each of the geographic areas referenced in Example Data.csv.

For the code to run properly, the area names in Example Data.csv and Example Population Data.csv much match the names given in the 'county.list' vector defined on line 29 of MCMC Code.R

For further details, see 'A Bayesian Susceptible-Infectious-Hospitalized-Ventilated-Recovered Model to Predict Demand for COVID-19 Inpatient Care in a Large Healthcare System', by Self et. al.

Some users have had trouble installing the package on R version 4.1.  The source code for all functions in the package is provided in the Source_Code file.
