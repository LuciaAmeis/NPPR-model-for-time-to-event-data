#########################################################################
#                                                                       #
#    A non-parametric proportional risk model to assess a treatment     #
#                    effect in time-to-event data                       #
#                                                                       #
#         L. Ameis, O. Kuss, A. Hoyer and K. Moellenhoff                #
#                                                                       #
#                           MASTER                                      #
#                                                                       #
#########################################################################
#This master file summarizes the results of the analysis
#Please run all other files beforehand (see Read Me for instructions)
####################
####Session Info####
####################
#session Info of the device the code was run on
#sessionInfo() 
#R version 4.3.0 (2023-04-21 ucrt)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 19045)

#Matrix products: default

#locale:
#[1] LC_COLLATE=German_Germany.utf8  LC_CTYPE=German_Germany.utf8   
#LC_MONETARY=German_Germany.utf8 LC_NUMERIC=C                   
#[5] LC_TIME=German_Germany.utf8    

#time zone: Europe/Berlin
#tzcode source: internal

#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base     

#loaded via a namespace (and not attached):
#[1] compiler_4.3.0    parallel_4.3.0    tools_4.3.0       
#rstudioapi_0.15.0 codetools_0.2-19  doParallel_1.0.17 iterators_1.0.14  
#foreach_1.5.2 


#Of note: Some of the results were produces using R version 4.2
#Please refer to the Read Me file for details

###############
####Library####
###############
library("Rlab")
library("rlist")
library("survival")
library("haven")
library("sas7bdat")
library("dplyr")
library("parallel")
library("foreach")
library("doParallel")
library("doSNOW")
library("numDeriv")
library("ggplot2")
library("flexsurv")

source ("https://git.io/fjinW")     #CDF/PDF of the EU distribution
#Provided DOI: 10.5433/1679-0375.2019v40n2p107
#Dec. 2019

#set WD
#setwd("Please define")

#######################################
##### Part 1: Intermediate results ####
#######################################
#The NPPR_base.R file contains all the functions used for analysis as well as
#the preparation of the simulations. This also includes the analysis of 
#the DapaHF data set with the studies primary outcome used for the simulations. 
source("NPPR_base.R")
#Please run this first in any case.

#The NPPR_case_study.R file contains the case study: The DapaHF data set with 
#the secondary outcome 'death of all causes'.
#See Section 4 of the manuscript.
source("case_study/NPPR_case_study.R")


#The NPPR_simulation.R file contains the estimation part of the simulation
#study. 

#WARNING: As is, the code contained in NPPR_simulation.R will run a long time!!!

#We recommend to skip this subfile and just run the files below. In this
#case the provided .RData file will be used for the following subfiles.
#The .RData can be manually loaded using:
#load("./simulation/NPPR_simulation.RData")

#A description how to reduce the run time was added to both the file at the 
#specific points in the code and the Read Me. Please refer to the code
#for specific instructions.
source("simulation/NPPR_simulation.R")


####################################
#### Part 2 - Results and Plots ####
####################################
#The NPPR_simulation_results.R file contains the collected results of the 
#simulation study including the Tables presented in section 3 of the manuscript
#and the supplementary material.

#If the instructions how to reduce the run time of the file NPPR_simulation.R were
#followed, some changes have to be made to this file as well. Please follow the 
#instructions provided at every beginning of a section.
source("simulation/NPPR_simulation_results.R")
#Please run the NPPR_simulation.R file beforehand.

#The NPPR_plot_generation.R file contains the code for all plots presented in the
#manuscript and the supplementary material.
#The plots are generated as pdfs.
source("plot_generation/NPPR_plot_generation.R")
#Please run all other files beforehand.