READ ME
##################################################################################################################
Source Code and data for the manuscript "A non-parametric proportional risk model to assess a treatment effect in time-to-event data", by L. Ameis, O. Kuss, A. Hoyer and K. Moellenhoff.
##################################################################################################################

For questions or comments pleas contact K. Moellenhoff at the Univeristy of Cologne: 
kathrin.moellenhoff@uni-koeln.de

This code as it is has been written on R version 4.3.0. 
Of note, the results regarding the NPPR estimator for PR and PH data as well as the PPR estimator for PR data (Tables 4, 5 for Bias and MSE, Tables C, D for Coverage and Tables G, H for numerical robustness) were generated on R version 4.2.2. Therefore, when run on R version 4.3.0, the results may differ slightly. The tables as saved in the respectiv folder present the data as listed in the manuscript.

Functions for applying the NPPR estimator as well as all functions used for the simulation study and its analysis are provided in the file NPPR_base.r. This has to be run first in any case.

To reproduce the simulation study run the file NPPR_simulation.r stored in the subfolder 'simulation'. The results are stored in the subfolder 'tables'.
WARNING: The calculation of the confidence intervals for the NPPR estimator is very time consuming. The time can be reduced by changing the argument 'L' from 500 to a smaller number. This is commented on in the file as well.
Also, due to the parallelization results of re-runs might differ slightly.

To reproduce the case study run the file NPPR_case_study.r stored in the subfolder 'case_study'.

To reproduce the figure presented in the manuscript please run all other files first. Afterwards run the file NPPR_plot_generation.r stored in the subfile 'plot_generation'. The figures are stored in the subfile 'plot'.

The folder 'data' contains the both DapaHF datasets used in the manuscript. The dataset with its primary outcome, that is used as base for the simulation study, is stored with the name dapa_hf_primary_outcome.sas7bdat. The data set with its secondary outcome, that is analysed as a case study, is stored with the name DAPA-HF_Death_from_any_cause.csv.

Please always start with setting the working directory to the folder containing the NPPR_base.r file in the provided line.
