READ ME
##################################################################################################################
Source Code and data for the manuscript "A non-parametric proportional risk model to assess a treatment effect in time-to-event data", by L. Ameis, O. Kuss, A. Hoyer and K. Moellenhoff.
##################################################################################################################

For questions or comments please contact K. Moellenhoff at the University of Cologne: 
kathrin.moellenhoff@uni-koeln.de

This code has been written on R version 4.3.0. 

For each file (NPPR_base.R, NPPR_simulation.R, NPPR_simulation_results.R, NPPR_case_study.R, NPPR_plot_generation.R and NPPR_master.R) please start with setting the working directory to the same location. Precisely, set it to the storage location of the file NPPR_base.R. Do so e.g. by using the function setwd() and inserting the address of the respective folder (e.g. setwd("C:/Users/example/storage_location")). We prepared the function for this (line 66 in the file NPPR_master.R and line 36 in all others). Here, substitute "Please define" for the correct address and remove the hashtag at the beginning of the line. 

Functions for applying the NPPR estimator as well as all functions used for the simulation study and its analysis are provided in the file NPPR_base.R. This must be run first in any case.

To reproduce the simulation study run the file NPPR_simulation.R stored in the subfolder 'simulation' to generate the simulation study and run the estimation. To generate the results run the file NPPR_simulation_results.R stored in the same subfolder. Tables are automatically stored in the subfolder 'tables'.

WARNING: The calculation of the confidence intervals for the NPPR estimator is very time consuming. The time can be reduced by changing the argument 'L' from 500 to a smaller number. For details we refer to the beginning of section 'Coverage of the NPPR, PPR and LL models' in file NPPR_simulation.R.
Also, due to the parallelization results of re-runs might differ slightly. 
The estimation of all simulations also takes time due to the large number of different scenarios. For a detailed instruction see the explanation in the file NPPR_simulation.R starting in line 737.

To reproduce the case study run the file NPPR_case_study.R stored in the subfolder 'case_study'.

To reproduce the figures presented in the manuscript please run all other files first. Afterwards run the file NPPR_plot_generation.R stored in the subfile 'plot_generation'. The figures are automatically stored in the subfolder 'plot' as pdfs. All figures are also stored as eps in the subfolder ‘plot_eps’.

The folder 'data' contains the both DapaHF datasets used in the manuscript. The dataset with its primary outcome, that is used as base for the simulation study, is stored with the name dapa_hf_primary_outcome.sas7bdat. The data set with its secondary outcome, that is analysed as a case study, is stored with the name DAPA-HF_Death_from_any_cause.csv.

The file NPPR_master.R contains the master script and can be used to run all other files without manually opening them.

For the .RData files including the older simulations we refer to Zenodo: https://zenodo.org/records/11103494
