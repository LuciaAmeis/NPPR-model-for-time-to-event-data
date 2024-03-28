#########################################################################
#                                                                       #
#    A non-parametric proportional risk model to assess a treatment     #
#                    effect in time-to-event data                       #
#                                                                       #
#         L. Ameis, O. Kuss, A. Hoyer and K. Moellenhoff                #
#                                                                       #
#                             CASE STUDY                                #
#                                                                       #
#########################################################################


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


setwd("Please define")

#load workspace
load("NPPR_base.RData")


#######################################################################
#### Analysis of the DapaHF dataset death of all causes (Section 4)####
#######################################################################
DapaHF_ac <- read.csv("./data/DAPA-HF_Death_from_any_cause.csv")
#create binary variable for the treatment group
DapaHF_ac$Treat_bin<-ifelse(DapaHF_ac$Treatment=="Placebo",0,1)
View(DapaHF_ac)
#col_time is 3, col_cens is 2 and col_target is 4

#Basic information
nrow(DapaHF_ac)
nrow(subset(DapaHF_ac,DapaHF_ac$Treat_bin==1))
nrow(subset(DapaHF_ac,DapaHF_ac$Treat_bin==1&DapaHF_ac$Event==1))
nrow(subset(DapaHF_ac,DapaHF_ac$Treat_bin==0))
nrow(subset(DapaHF_ac,DapaHF_ac$Treat_bin==0&DapaHF_ac$Event==1))


#calculate the HR to validate digitized data
summary(coxph(Surv(SurvivalTimeMonths,Event)~Treat_bin,data=DapaHF_ac))


#estimating beta and RR
estimate_ac<-B_easy_fast(DapaHF_ac,col=c(col_time=3,col_cens=2),col_target=4)
exp(-estimate_ac)

#confidence interval for beta and RR
ci_fast_ac<-CI_fast(DapaHF_ac,col=c(col_time=3,col_cens=2),col_target=4,L=500,alpha=0.05)
exp(-ci_fast_ac)


#nnt at time point 24 (end of study)
nnt_dapahfac_24<-NNT_fast(DapaHF_ac,col=c(col_time=3,col_cens=2),col_target=4,24)

#CI for nnt at time point 24 (end of study)
ci_nnt_dapahfac_24<-CI_nnt_fast(DapaHF_ac,col=c(col_time=3,col_cens=2),col_target=4,L=500,alpha=0.05,24)


#####Save Worspace####
save.image(file = "./case_study/NPPR_case_study.RData")
