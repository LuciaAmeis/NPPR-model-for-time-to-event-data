#########################################################################
#                                                                       #
#    A non-parametric proportional risk model to assess a treatment     #
#                    effect in time-to-event data                       #
#                                                                       #
#         L. Ameis, O. Kuss, A. Hoyer and K. Moellenhoff                #
#                                                                       #
#                     SIMULATION  RESULTS                               #
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


#set WD
#setwd("Please define")


#load workspace
load("NPPR_base.RData")
load("./simulation/NPPR_simulation_subset.RData")

############################
####Numerical robustness####
############################
#WARNING: If only a subset of simulation scenarios was regarded in the file
#NPPR_simulaion.R some adjustment have to be made to the code below. First 
#remove the hash tags from the code in lines 52 to 125 and run it. Afterwards,
#run the tables as is. The cases that were not regarded will be filled with NAs.

#PR cases
#NPPR
#NA_PR_NPPR_filled<-rep(NA,45) #since there are 45 PR scenarios (with censoring)
                              #we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  NA_PR_NPPR_filled[j]<-NA_PR_NPPR[i] #order the numerical robustness to the index
#}
#NA_PR_NPPR<-NA_PR_NPPR_filled #change name to use following code

#PPR
#NA_PR_PPR_filled<-rep(NA,45) #since there are 45 PR scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  NA_PR_PPR_filled[j]<-NA_PR_PPR[i] #order the numerical robustness to the index
#}
#NA_PR_PPR<-NA_PR_PPR_filled #change name to use following code

#Cox
#NA_PR_Cox_filled<-rep(NA,45) #since there are 45 PR scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  NA_PR_Cox_filled[j]<-NA_PR_Cox[i] #order the numerical robustness to the index
#}
#NA_PR_Cox<-NA_PR_Cox_filled #change name to use following code

#LL
#NA_PR_LL_filled<-rep(NA,45) #since there are 45 PR scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  NA_PR_LL_filled[j]<-NA_PR_LL[i] #order the numerical robustness to the index
#}
#NA_PR_LL<-NA_PR_LL_filled #change name to use following code


#PH cases
#NPPR
#NA_PH_NPPR_filled<-rep(NA,45) #since there are 45 PH scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  NA_PH_NPPR_filled[j]<-NA_PH_NPPR[i] #order the numerical robustness to the index
#}
#NA_PH_NPPR<-NA_PH_NPPR_filled #change name to use following code

#Cox
#NA_PH_Cox_filled<-rep(NA,45) #since there are 45 PH scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  NA_PH_Cox_filled[j]<-NA_PH_Cox[i] #order the numerical robustness to the index
#}
#NA_PH_Cox<-NA_PH_Cox_filled #change name to use following code


#PO cases
#NPPR
#NA_PO_NPPR_filled<-rep(NA,45) #since there are 45 PO scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  NA_PO_NPPR_filled[j]<-NA_PO_NPPR[i] #order the numerical robustness to the index
#}
#NA_PO_NPPR<-NA_PO_NPPR_filled #change name to use following code

#LL
#NA_PO_LL_filled<-rep(NA,45) #since there are 45 PO scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  NA_PO_LL_filled[j]<-NA_PO_LL[i] #order the numerical robustness to the index
#}
#NA_PO_LL<-NA_PO_LL_filled #change name to use following code

#With this we can produce the Tables as is, with NA for all cases, that were
#not regarded.

#Numerical robustness
#Create Tables
#Of note, the Tables in the paper have another order and need to be reordered
#Table G Numerical Robustness for PR data
TableG<-data.frame(
  "Effect"=c(rep(0,9),rep(0.5,9),rep(0.25,9),rep(-0.25,9),rep(-0.5,9)),
  "Censoring (%)"=rep(c(rep(30,3),rep(50,3),rep(70,3)),5),
  "Participants"=rep(c(500,100,50),15),
  "NPPR"=c(1000-NA_PR_NPPR[1:18],1000-NA_PR_NPPR[28:36],1000-NA_PR_NPPR[37:45],1000-NA_PR_NPPR[19:27]),
  "PPR"=c(NA_PR_PPR[1:18],NA_PR_PPR[28:36],NA_PR_PPR[37:45],NA_PR_PPR[19:27]),
  "Cox"=c(1000-NA_PR_Cox[1:18],1000-NA_PR_Cox[28:36],1000-NA_PR_Cox[37:45],1000-NA_PR_Cox[19:27]),
  "LL"=c(NA_PR_LL[1:18],NA_PR_LL[28:36],NA_PR_LL[37:45],NA_PR_LL[19:27])
)
write.csv2(TableG,file="./simulation/tables/TableG.csv")

#Table H Numerical Robustness for PH data
TableH<-data.frame(
  "Effect"=c(rep(0,9),rep(0.5,9),rep(0.25,9),rep(-0.25,9),rep(-0.5,9)),
  "Censoring (%)"=rep(c(rep(30,3),rep(50,3),rep(70,3)),5),
  "Participants"=rep(c(500,100,50),15),
  "NPPR"=c(1000-NA_PH_NPPR[1:18],1000-NA_PH_NPPR[28:36],1000-NA_PH_NPPR[37:45],1000-NA_PH_NPPR[19:27]),
  "Cox"=c(1000-NA_PH_Cox[1:18],1000-NA_PH_Cox[28:36],1000-NA_PH_Cox[37:45],1000-NA_PH_Cox[19:27])
)
write.csv2(TableH,file="./simulation/tables/TableH.csv")

#Table I Numerical Robustness for PO data
TableI<-data.frame(
  "Effect"=c(rep(0,9),rep(0.5,9),rep(0.25,9),rep(-0.25,9),rep(-0.5,9)),
  "Censoring (%)"=rep(c(rep(30,3),rep(50,3),rep(70,3)),5),
  "Participants"=rep(c(500,100,50),15),
  "NPPR"=c(1000-NA_PO_NPPR[1:18],1000-NA_PO_NPPR[28:36],1000-NA_PO_NPPR[37:45],1000-NA_PO_NPPR[19:27]),
  "LL"=c(NA_PO_LL[1:18],NA_PO_LL[28:36],NA_PO_LL[37:45],NA_PO_LL[19:27])
)
write.csv2(TableI,file="./simulation/tables/TableI.csv")

#WARNING:
#As above, the nocens cases have to be prepared if only a subset of the simulation 
#scenarios was regarded. Please remove the hash tags and run the code from line 172
#to 236

#PR cases
#NPPR
#NA_PR_NPPR_nocens_filled<-rep(NA,15) #since there are 15 PR scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  NA_PR_NPPR_nocens_filled[j]<-NA_PR_NPPR_nocens[i] #order the numerical robustness to the index
#}
#NA_PR_NPPR_nocens<-NA_PR_NPPR_nocens_filled #change name to use following code

#Cox
#NA_PR_Cox_nocens_filled<-rep(NA,15) #since there are 15 PR scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  NA_PR_Cox_nocens_filled[j]<-NA_PR_Cox_nocens[i] #order the numerical robustness to the index
#}
#NA_PR_Cox_nocens<-NA_PR_Cox_nocens_filled #change name to use following code

#LL
#NA_PR_LL_nocens_filled<-rep(NA,15) #since there are 15 PR scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  NA_PR_LL_nocens_filled[j]<-NA_PR_LL_nocens[i] #order the numerical robustness to the index
#}
#NA_PR_LL_nocens<-NA_PR_LL_nocens_filled #change name to use following code


#PH cases
#NPPR
#NA_PH_NPPR_nocens_filled<-rep(NA,15) #since there are 15 PH scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  NA_PH_NPPR_nocens_filled[j]<-NA_PH_NPPR_nocens[i] #order the numerical robustness to the index
#}
#NA_PH_NPPR_nocens<-NA_PH_NPPR_nocens_filled #change name to use following code

#Cox
#NA_PH_Cox_nocens_filled<-rep(NA,15) #since there are 15 PH scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  NA_PH_Cox_nocens_filled[j]<-NA_PH_Cox_nocens[i] #order the numerical robustness to the index
#}
#NA_PH_Cox_nocens<-NA_PH_Cox_nocens_filled #change name to use following code


#PO cases
#NPPR
#NA_PO_NPPR_nocens_filled<-rep(NA,15) #since there are 15 PO scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  NA_PO_NPPR_nocens_filled[j]<-NA_PO_NPPR_nocens[i] #order the numerical robustness to the index
#}
#NA_PO_NPPR_nocens<-NA_PO_NPPR_nocens_filled #change name to use following code

#LL
#NA_PO_LL_nocens_filled<-rep(NA,15) #since there are 15 PO scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  NA_PO_LL_nocens_filled[j]<-NA_PO_LL_nocens[i] #order the numerical robustness to the index
#}
#NA_PO_LL_nocens<-NA_PO_LL_nocens_filled #change name to use following code

#With this we can produce the Tables as is, with NA for all cases, that were
#not regarded.


#Table J Numerical Robustness for no censoring
TableJ<-data.frame(
  "Underlying Model"=c(rep("PPR",15),rep("Weibull PH",15),rep("LL PO",15)),
  "Effect"=rep(c(rep(0,3),rep(0.5,3),rep(0.25,3),rep(-0.25,3),rep(-0.5,3)),3),
  "Participants"=rep(c(500,100,50),15),
  "NPPR"=c(1000-NA_PR_NPPR_nocens[1:6],1000-NA_PR_NPPR_nocens[10:12],1000-NA_PR_NPPR_nocens[13:15],1000-NA_PR_NPPR_nocens[7:9],1000-NA_PH_NPPR_nocens[1:6],1000-NA_PH_NPPR_nocens[10:12],1000-NA_PH_NPPR_nocens[13:15],1000-NA_PH_NPPR_nocens[7:9],1000-NA_PO_NPPR_nocens[1:6],1000-NA_PO_NPPR_nocens[10:12],1000-NA_PO_NPPR_nocens[13:15],1000-NA_PO_NPPR_nocens[7:9]),
  "Cox"=c(1000-NA_PR_Cox_nocens[1:6],1000-NA_PR_Cox_nocens[10:12],1000-NA_PR_Cox_nocens[13:15],1000-NA_PR_Cox_nocens[7:9],1000-NA_PH_Cox_nocens[1:6],1000-NA_PH_Cox_nocens[10:12],1000-NA_PH_Cox_nocens[13:15],1000-NA_PH_Cox_nocens[7:9],rep("-",15)),
  "LL"=c(NA_PR_LL_nocens[1:6],NA_PR_LL_nocens[10:12],NA_PR_LL_nocens[13:15],NA_PR_LL_nocens[7:9],rep("-",15),NA_PO_LL_nocens[1:6],NA_PO_LL_nocens[10:12],NA_PO_LL_nocens[13:15],NA_PO_LL_nocens[7:9])
)
write.csv2(TableJ,file="./simulation/tables/TableJ.csv")

############################
####Output (Section 3.3)####
############################
#PR case
#NPPR
#List to store results
Out_PR_NPPR<-list()
for(i in 1:length(Est_PR_NPPR_nna)){
  b1<-true.beta[i]          #extract true beta
  Est<-Est_PR_NPPR_nna[[i]] #extract estimates
  Out_PR_NPPR[[i]]<-Sim_out_log(b_est=Est,b=b1) #data.frame including mean, Bias and MSE of the given case
}

#PPR
Out_PR_PPR<-list()
for(i in 1:length(Est_PR_PPR_nna)){
  b1<-true.beta[i]
  Est<-Est_PR_PPR_nna[[i]]
  Out_PR_PPR[[i]]<-Sim_EU_out_log(ests=Est,b=b1) #data.frame including mean, bias and MSE of the given case
}

#Cox's model
Out_PR_Cox<-list()
for(i in 1:length(Est_PR_Cox_nna)){
  b1<-true.beta[i]
  Est<-Est_PR_Cox_nna[[i]]
  Out_PR_Cox[[i]]<-Sim_PH_out(est=Est,b=b1) #data.frame including mean, bias, MSE and coverage of the given case
}

#LL
Out_PR_LL<-list()
for(i in 1:length(Est_PR_LL_nna)){
  b1<-true.beta[i]
  Est<-Est_PR_LL_nna[[i]]
  Out_PR_LL[[i]]<-Sim_PO_out_log(ests=Est,b=b1) #data.frame including mean, bias and MSE of the given case
}


#PH case
#NPPR
Out_PH_NPPR<-list()
for(i in 1:length(Est_PH_NPPR_nna)){
  b1<-true.beta[i]
  Est<-Est_PH_NPPR_nna[[i]]
  Out_PH_NPPR[[i]]<-Sim_out_log(b_est=Est,b=b1) #data.frame including mean, Bias and MSE of the given case
}

#Cox's model
Out_PH_Cox<-list()
for(i in 1:length(Est_PH_Cox_nna)){
  b1<-true.beta[i]
  Est<-Est_PH_Cox_nna[[i]]
  Out_PH_Cox[[i]]<-Sim_PH_out(est=Est,b=b1) #data.frame including mean, bias, MSE and coverage of the given case
}


#PO case
#NPPR
Out_PO_NPPR<-list()
for(i in 1:length(Est_PO_NPPR_nna)){
  b1<-true.beta[i]
  Est<-Est_PO_NPPR_nna[[i]]
  Out_PO_NPPR[[i]]<-Sim_out_log(b_est=Est,b=b1) #data.frame including mean, Bias and MSE of the given case
}

#LL
Out_PO_LL<-list()
for(i in 1:length(Est_PO_LL_nna)){
  b1<-true.beta[i]
  Est<-Est_PR_LL_nna[[i]]
  Out_PO_LL[[i]]<-Sim_PO_out_log(ests=Est,b=b1) #data.frame including mean, bias and MSE of the given case
}


#No censoring
#PR case
#NPPR
Out_PR_NPPR_nocens<-list()
for(i in 1:length(Est_PR_NPPR_nocens_nna)){
  b1<-true.beta_nocens[i]
  Est<-Est_PR_NPPR_nocens_nna[[i]]
  Out_PR_NPPR_nocens[[i]]<-Sim_out_log(b_est=Est,b=b1) #data.frame including mean, Bias and MSE of the given case
}

#Cox's model
Out_PR_Cox_nocens<-list()
for(i in 1:length(Est_PR_Cox_nocens_nna)){
  b1<-true.beta_nocens[i]
  Est<-Est_PR_Cox_nocens_nna[[i]]
  Out_PR_Cox_nocens[[i]]<-Sim_PH_out(est=Est,b=b1) #data.frame including mean, bias, MSE and coverage of the given case
}

#LL
Out_PR_LL_nocens<-list()
for(i in 1:length(Est_PR_LL_nocens_nna)){
  b1<-true.beta_nocens[i]
  Est<-Est_PR_LL_nocens_nna[[i]]
  Out_PR_LL_nocens[[i]]<-Sim_PO_out_log(ests=Est,b=b1) #data.frame including mean, bias and MSE of the given case
}

#PH case
#NPPR
Out_PH_NPPR_nocens<-list()
for(i in 1:length(Est_PH_NPPR_nocens_nna)){
  b1<-true.beta_nocens[i]
  Est<-Est_PH_NPPR_nocens_nna[[i]]
  Out_PH_NPPR_nocens[[i]]<-Sim_out_log(b_est=Est,b=b1) #data.frame including mean, Bias and MSE of the given case
}

#Cox's model
Out_PH_Cox_nocens<-list()
for(i in 1:length(Est_PH_Cox_nocens_nna)){
  b1<-true.beta_nocens[i]
  Est<-Est_PH_Cox_nocens_nna[[i]]
  Out_PH_Cox_nocens[[i]]<-Sim_PH_out(est=Est,b=b1) #data.frame including mean, bias, MSE and coverage of the given case
}


#PO case
#NPPR
Out_PO_NPPR_nocens<-list()
for(i in 1:length(Est_PO_NPPR_nocens_nna)){
  b1<-true.beta_nocens[i]
  Est<-Est_PO_NPPR_nocens_nna[[i]]
  Out_PO_NPPR_nocens[[i]]<-Sim_out_log(b_est=Est,b=b1) #data.frame including mean, Bias and MSE of the given case
}

#LL
Out_PO_LL_nocens<-list()
for(i in 1:length(Est_PO_LL_nocens_nna)){
  b1<-true.beta_nocens[i]
  Est<-Est_PR_LL_nocens_nna[[i]]
  Out_PO_LL_nocens[[i]]<-Sim_PO_out_log(ests=Est,b=b1) #data.frame including mean, bias and MSE of the given case
}

#After collecting the results, we will now reproduce the Tables as represented
#first we create data frames for each model each including the PR, PH and PO scenario
#respectively
#vector of censoring

#Next we combine the list to one large data frame 
#PR case
#NPPR model
Results_PR_NPPR<-data.frame() #empty data frame to store results
for(i in 1:length(Out_PR_NPPR)){
  data<-Out_PR_NPPR[[i]]      #extract output
  Results_PR_NPPR<-rbind(Results_PR_NPPR,data) #binde to data frame
}


#PPR model
Results_PR_PPR<-data.frame()
for(i in 1:length(Out_PR_PPR)){
  data<-Out_PR_PPR[[i]]
  Results_PR_PPR<-rbind(Results_PR_PPR,data)
}


#Cox's model
Results_PR_Cox<-data.frame()
for(i in 1:length(Out_PR_Cox)){
  data<-Out_PR_Cox[[i]]
  Results_PR_Cox<-rbind(Results_PR_Cox,data)
}

#LL model
Results_PR_LL<-data.frame()
for(i in 1:length(Out_PR_LL)){
  data<-Out_PR_LL[[i]]
  Results_PR_LL<-rbind(Results_PR_LL,data)
}


#PH case
#NPPR model
Results_PH_NPPR<-data.frame()
for(i in 1:length(Out_PH_NPPR)){
  data<-Out_PH_NPPR[[i]]
  Results_PH_NPPR<-rbind(Results_PH_NPPR,data)
}


#Cox's model
Results_PH_Cox<-data.frame()
for(i in 1:length(Out_PH_Cox)){
  data<-Out_PH_Cox[[i]]
  Results_PH_Cox<-rbind(Results_PH_Cox,data)
}


#PO case
#NPPR model
Results_PO_NPPR<-data.frame()
for(i in 1:length(Out_PO_NPPR)){
  data<-Out_PO_NPPR[[i]]
  Results_PO_NPPR<-rbind(Results_PO_NPPR,data)
}


#LL model
Results_PO_LL<-data.frame()
for(i in 1:length(Out_PO_LL)){
  data<-Out_PO_LL[[i]]
  Results_PO_LL<-rbind(Results_PO_LL,data)
}


#No censoring
#PR case
#NPPR model
Results_PR_NPPR_nocens<-data.frame()
for(i in 1:length(Out_PR_NPPR_nocens)){
  data<-Out_PR_NPPR_nocens[[i]]
  Results_PR_NPPR_nocens<-rbind(Results_PR_NPPR_nocens,data)
}


#Cox's model
Results_PR_Cox_nocens<-data.frame()
for(i in 1:length(Out_PR_Cox_nocens)){
  data<-Out_PR_Cox_nocens[[i]]
  Results_PR_Cox_nocens<-rbind(Results_PR_Cox_nocens,data)
}


#LL model
Results_PR_LL_nocens<-data.frame()
for(i in 1:length(Out_PR_LL_nocens)){
  data<-Out_PR_LL_nocens[[i]]
  Results_PR_LL_nocens<-rbind(Results_PR_LL_nocens,data)
}


#PH case
#NPPR model
Results_PH_NPPR_nocens<-data.frame()
for(i in 1:length(Out_PH_NPPR_nocens)){
  data<-Out_PH_NPPR_nocens[[i]]
  Results_PH_NPPR_nocens<-rbind(Results_PH_NPPR_nocens,data)
}


#Cox's model
Results_PH_Cox_nocens<-data.frame()
for(i in 1:length(Out_PH_Cox_nocens)){
  data<-Out_PH_Cox_nocens[[i]]
  Results_PH_Cox_nocens<-rbind(Results_PH_Cox_nocens,data)
}


#PO case
#NPPR model
Results_PO_NPPR_nocens<-data.frame()
for(i in 1:length(Out_PO_NPPR_nocens)){
  data<-Out_PO_NPPR_nocens[[i]]
  Results_PO_NPPR_nocens<-rbind(Results_PO_NPPR_nocens,data)
}


#LL model
Results_PO_LL_nocens<-data.frame()
for(i in 1:length(Out_PO_LL_nocens)){
  data<-Out_PO_LL_nocens[[i]]
  Results_PO_LL_nocens<-rbind(Results_PO_LL_nocens,data)
}


#now we have to ad the coverage for the NPPR, PPR and LL models 
#PR case
#NPPR model
Results_PR_NPPR$Coverage<-CI_PR_NPPR

#PPR model
Results_PR_PPR$Coverage<-CI_PR_PPR

#LL model
Results_PR_LL$Coverage<-CI_PR_LL


#PH case
#NPPR model
Results_PH_NPPR$Coverage<-CI_PH_NPPR


#PO case
#NPPR model
Results_PO_NPPR$Coverage<-CI_PO_NPPR

#LL model
Results_PO_LL$Coverage<-CI_PO_LL


#No censoring
#PR case
#NPPR model
Results_PR_NPPR_nocens$Coverage<-CI_PR_NPPR_nocens

#LL model
Results_PR_LL_nocens$Coverage<-CI_PR_LL_nocens


#PH case
#NPPR model
Results_PH_NPPR_nocens$Coverage<-CI_PH_NPPR_nocens


#PO case
#NPPR model
Results_PO_NPPR_nocens$Coverage<-CI_PO_NPPR_nocens

#LL model
Results_PO_LL_nocens$Coverage<-CI_PO_LL_nocens


#The respective results are collected next

#WARNING: If only a subset of scenarios was regarded at this point some changes 
#have to be made to the code in order to produce the following tables. Precisely,
#we will fill the result data frames with NAs for the indices of the simulation
#scenarios, that were not regarded. Please remove the hash tags and run the code
#from line to

#PR cases
#NPPR
#Results_PR_NPPR_filled<-data.frame(
#"Variable"=rep(1,45),
#"Mean"=rep(NA,45),
#"Bias"=rep(NA,45),
#"MSE"=rep(NA,45),
#"Coverage"=rep(NA,45)
#) #since there are 15 PR scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  Results_PR_NPPR_filled[j,]<-Results_PR_NPPR[i,] #order the numerical robustness to the index
#}
#Results_PR_NPPR<-Results_PR_NPPR_filled #change name to use following code

#PPR
#Results_PR_PPR_filled<-data.frame(
#"Variable"=rep(1,45),
#"Mean"=rep(NA,45),
#"Bias"=rep(NA,45),
#"MSE"=rep(NA,45),
#"Coverage"=rep(NA,45)
#) #since there are 15 PR scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  Results_PR_PPR_filled[j,]<-Results_PR_PPR[i,] #order the numerical robustness to the index
#}
#Results_PR_PPR<-Results_PR_PPR_filled #change name to use following code

#Cox
#Results_PR_Cox_filled<-data.frame(
#"Variable"=rep(1,45),
#"Mean"=rep(NA,45),
#"Bias"=rep(NA,45),
#"MSE"=rep(NA,45),
#"Coverage"=rep(NA,45)
#) #since there are 15 PR scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  Results_PR_Cox_filled[j,]<-Results_PR_Cox[i,] #order the numerical robustness to the index
#}
#Results_PR_Cox<-Results_PR_Cox_filled #change name to use following code


#LL
#Results_PR_LL_filled<-data.frame(
#"Variable"=rep(1,45),
#"Mean"=rep(NA,45),
#"Bias"=rep(NA,45),
#"MSE"=rep(NA,45),
#"Coverage"=rep(NA,45)
#) #since there are 15 PR scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  Results_PR_LL_filled[j,]<-Results_PR_LL[i,] #order the numerical robustness to the index
#}
#Results_PR_LL<-Results_PR_LL_filled #change name to use following code


#PH cases
#NPPR
#Results_PH_NPPR_filled<-data.frame(
#"Variable"=rep(1,45),
#"Mean"=rep(NA,45),
#"Bias"=rep(NA,45),
#"MSE"=rep(NA,45),
#"Coverage"=rep(NA,45)
#) #since there are 15 PH scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  Results_PH_NPPR_filled[j,]<-Results_PH_NPPR[i,] #order the numerical robustness to the index
#}
#Results_PH_NPPR<-Results_PH_NPPR_filled #change name to use following code

#Cox
#Results_PH_Cox_filled<-data.frame(
#"Variable"=rep(1,45),
#"Mean"=rep(NA,45),
#"Bias"=rep(NA,45),
#"MSE"=rep(NA,45),
#"Coverage"=rep(NA,45)
#) #since there are 15 PH scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  Results_PH_Cox_filled[j,]<-Results_PH_Cox[i,] #order the numerical robustness to the index
#}
#Results_PH_Cox<-Results_PH_Cox_filled #change name to use following code


#PO cases
#NPPR
#Results_PO_NPPR_filled<-data.frame(
#"Variable"=rep(1,45),
#"Mean"=rep(NA,45),
#"Bias"=rep(NA,45),
#"MSE"=rep(NA,45),
#"Coverage"=rep(NA,45)
#) #since there are 15 PO scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  Results_PO_NPPR_filled[j,]<-Results_PO_NPPR[i,] #order the numerical robustness to the index
#}
#Results_PO_NPPR<-Results_PO_NPPR_filled #change name to use following code

#LL
#Results_PO_LL_filled<-data.frame(
#"Variable"=rep(1,45),
#"Mean"=rep(NA,45),
#"Bias"=rep(NA,45),
#"MSE"=rep(NA,45),
#"Coverage"=rep(NA,45)
#) #since there are 15 PO scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  Results_PO_LL_filled[j,]<-Results_PO_LL[i,] #order the numerical robustness to the index
#}
#Results_PO_LL<-Results_PO_LL_filled #change name to use following code

#Now we can run the tables.

#Of note, again we need to rearrange the order to fit the paper
#PR case
#Bias Table 4
#MSE Table 4
Table4<-data.frame(
  "Effect"=c(rep(0,9),rep(0.5,9),rep(0.25,9),rep(-0.25,9),rep(-0.5,9)),
  "Censoring (%)"=rep(c(rep(30,3),rep(50,3),rep(70,3)),5),
  "Participants"=rep(c(500,100,50),15),
  "Bias_NPPR"=c(Results_PR_NPPR$Bias[1:18],Results_PR_NPPR$Bias[28:36],Results_PR_NPPR$Bias[37:45],Results_PR_NPPR$Bias[19:27]),
  "Bias_PPR"=c(Results_PR_PPR$Bias[1:18],Results_PR_PPR$Bias[28:36],Results_PR_PPR$Bias[37:45],Results_PR_PPR$Bias[19:27]),
  "Bias_Cox"=c(Results_PR_Cox$Bias[1:18],Results_PR_Cox$Bias[28:36],Results_PR_Cox$Bias[37:45],Results_PR_Cox$Bias[19:27]),
  "Bias_LL"=c(Results_PR_LL$Bias[1:18],Results_PR_LL$Bias[28:36],Results_PR_LL$Bias[37:45],Results_PR_LL$Bias[19:27]),
  "MSE_NPPR"=c(Results_PR_NPPR$MSE[1:18],Results_PR_NPPR$MSE[28:36],Results_PR_NPPR$MSE[37:45],Results_PR_NPPR$MSE[19:27]),
  "MSE_PPR"=c(Results_PR_PPR$MSE[1:18],Results_PR_PPR$MSE[28:36],Results_PR_PPR$MSE[37:45],Results_PR_PPR$MSE[19:27]),
  "MSE_Cox"=c(Results_PR_Cox$MSE[1:18],Results_PR_Cox$MSE[28:36],Results_PR_Cox$MSE[37:45],Results_PR_Cox$MSE[19:27]),
  "MSE_LL"=c(Results_PR_LL$MSE[1:18],Results_PR_LL$MSE[28:36],Results_PR_LL$MSE[37:45],Results_PR_LL$MSE[19:27])
)
write.csv2(Table4,file="./simulation/tables/Table4.csv")

#Coverage 
#Table C
TableC<-data.frame(
  "Effect"=c(rep(0,9),rep(0.5,9),rep(0.25,9),rep(-0.25,9),rep(-0.5,9)),
  "Censoring (%)"=rep(c(rep(30,3),rep(50,3),rep(70,3)),5),
  "Participants"=rep(c(500,100,50),15),
  "Coverage_NPPR"=c(Results_PR_NPPR$Coverage[1:18],Results_PR_NPPR$Coverage[28:36],Results_PR_NPPR$Coverage[37:45],Results_PR_NPPR$Coverage[19:27]),
  "Coverage_PPR"=c(Results_PR_PPR$Coverage[1:18],Results_PR_PPR$Coverage[28:36],Results_PR_PPR$Coverage[37:45],Results_PR_PPR$Coverage[19:27]),
  "Coverage_Cox"=c(Results_PR_Cox$Coverage[1:18],Results_PR_Cox$Coverage[28:36],Results_PR_Cox$Coverage[37:45],Results_PR_Cox$Coverage[19:27]),
  "Coverage_LL"=c(Results_PR_LL$Coverage[1:18],Results_PR_LL$Coverage[28:36],Results_PR_LL$Coverage[37:45],Results_PR_LL$Coverage[19:27])
)
write.csv2(TableC,file="./simulation/tables/TableC.csv")



#PH case 
#Bias Table 5
#MSE Table 5
Table5<-data.frame(
  "Effect"=c(rep(0,9),rep(0.5,9),rep(0.25,9),rep(-0.25,9),rep(-0.5,9)),
  "Censoring (%)"=rep(c(rep(30,3),rep(50,3),rep(70,3)),5),
  "Participants"=rep(c(500,100,50),15),
  "Bias_NPPR"=c(Results_PH_NPPR$Bias[1:18],Results_PH_NPPR$Bias[28:36],Results_PH_NPPR$Bias[37:45],Results_PH_NPPR$Bias[19:27]),
  "Bias_Cox"=c(Results_PH_Cox$Bias[1:18],Results_PH_Cox$Bias[28:36],Results_PH_Cox$Bias[37:45],Results_PH_Cox$Bias[19:27]),
  "MSE_NPPR"=c(Results_PH_NPPR$MSE[1:18],Results_PH_NPPR$MSE[28:36],Results_PH_NPPR$MSE[37:45],Results_PH_NPPR$MSE[19:27]),
  "MSE_Cox"=c(Results_PH_Cox$MSE[1:18],Results_PH_Cox$MSE[28:36],Results_PH_Cox$MSE[37:45],Results_PH_Cox$MSE[19:27])
)
write.csv2(Table5,file="./simulation/tables/Table5.csv")

#Coverage 
#Table D
TableD<-data.frame(
  "Effect"=c(rep(0,9),rep(0.5,9),rep(0.25,9),rep(-0.25,9),rep(-0.5,9)),
  "Censoring (%)"=rep(c(rep(30,3),rep(50,3),rep(70,3)),5),
  "Participants"=rep(c(500,100,50),15),
  "Coverage_NPPR"=c(Results_PH_NPPR$Coverage[1:18],Results_PH_NPPR$Coverage[28:36],Results_PH_NPPR$Coverage[37:45],Results_PH_NPPR$Coverage[19:27]),
  "Coverage_Cox"=c(Results_PH_Cox$Coverage[1:18],Results_PH_Cox$Coverage[28:36],Results_PH_Cox$Coverage[37:45],Results_PH_Cox$Coverage[19:27])
)
write.csv2(TableD,file="./simulation/tables/TableD.csv")


#PO case
#Bias Table 6
#MSE Table 6
Table6<-data.frame(
  "Effect"=c(rep(0,9),rep(0.5,9),rep(0.25,9),rep(-0.25,9),rep(-0.5,9)),
  "Censoring (%)"=rep(c(rep(30,3),rep(50,3),rep(70,3)),5),
  "Participants"=rep(c(500,100,50),15),
  "Bias_NPPR"=c(Results_PO_NPPR$Bias[1:18],Results_PO_NPPR$Bias[28:36],Results_PO_NPPR$Bias[37:45],Results_PO_NPPR$Bias[19:27]),
  "Bias_LL"=c(Results_PO_LL$Bias[1:18],Results_PO_LL$Bias[28:36],Results_PO_LL$Bias[37:45],Results_PO_LL$Bias[19:27]),
  "MSE_NPPR"=c(Results_PO_NPPR$MSE[1:18],Results_PO_NPPR$MSE[28:36],Results_PO_NPPR$MSE[37:45],Results_PO_NPPR$MSE[19:27]),
  "MSE_LL"=c(Results_PO_LL$MSE[1:18],Results_PO_LL$MSE[28:36],Results_PO_LL$MSE[37:45],Results_PO_LL$MSE[19:27])
)
write.csv2(Table6,file="./simulation/tables/Table6.csv")

#Coverage 
#Table E
TableE<-data.frame(
  "Effect"=c(rep(0,9),rep(0.5,9),rep(0.25,9),rep(-0.25,9),rep(-0.5,9)),
  "Censoring (%)"=rep(c(rep(30,3),rep(50,3),rep(70,3)),5),
  "Participants"=rep(c(500,100,50),15),
  "Coverage_NPPR"=c(Results_PO_NPPR$Coverage[1:18],Results_PO_NPPR$Coverage[28:36],Results_PO_NPPR$Coverage[37:45],Results_PO_NPPR$Coverage[19:27]),
  "Coverage_LL"=c(Results_PO_LL$Coverage[1:18],Results_PO_LL$Coverage[28:36],Results_PO_LL$Coverage[37:45],Results_PO_LL$Coverage[19:27])
)
write.csv2(TableE,file="./simulation/tables/TableE.csv")


#No censoring
#WARNING: If only a subset of scenarios was regarded, we need to change the code
#as above. Please remove the hash tags and run the code from line to and change the 
#example vector 1:9 to the vector with indices of interest

#PR cases
#NPPR
#Results_PR_NPPR_nocens_filled<-data.frame(
#"Variable"=rep(1,15),
#"Mean"=rep(NA,15),
#"Bias"=rep(NA,15),
#"MSE"=rep(NA,15),
#"Coverage"=rep(NA,15)
#) #since there are 15 PR scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  Results_PR_NPPR_nocens_filled[j,]<-Results_PR_NPPR_nocens[i,] #order the numerical robustness to the index
#}
#Results_PR_NPPR_nocens<-Results_PR_NPPR_nocens_filled #change name to use following code

#Cox
#Results_PR_Cox_nocens_filled<-data.frame(
#"Variable"=rep(1,15),
#"Mean"=rep(NA,15),
#"Bias"=rep(NA,15),
#"MSE"=rep(NA,15),
#"Coverage"=rep(NA,15)
#) #since there are 15 PR scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  Results_PR_Cox_nocens_filled[j,]<-Results_PR_Cox_nocens[i,] #order the numerical robustness to the index
#}
#Results_PR_Cox_nocens<-Results_PR_Cox_nocens_filled #change name to use following code


#LL
#Results_PR_LL_nocens_filled<-data.frame(
#"Variable"=rep(1,15),
#"Mean"=rep(NA,15),
#"Bias"=rep(NA,15),
#"MSE"=rep(NA,15),
#"Coverage"=rep(NA,15)
#) #since there are 15 PR scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  Results_PR_LL_nocens_filled[j,]<-Results_PR_LL_nocens[i,] #order the numerical robustness to the index
#}
#Results_PR_LL_nocens<-Results_PR_LL_nocens_filled #change name to use following code


#PH cases
#NPPR
#Results_PH_NPPR_nocens_filled<-data.frame(
#"Variable"=rep(1,15),
#"Mean"=rep(NA,15),
#"Bias"=rep(NA,15),
#"MSE"=rep(NA,15),
#"Coverage"=rep(NA,15)
#) #since there are 15 PH scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  Results_PH_NPPR_nocens_filled[j,]<-Results_PH_NPPR_nocens[i,] #order the numerical robustness to the index
#}
#Results_PH_NPPR_nocens<-Results_PH_NPPR_nocens_filled #change name to use following code

#Cox
#Results_PH_Cox_nocens_filled<-data.frame(
#"Variable"=rep(1,15),
#"Mean"=rep(NA,15),
#"Bias"=rep(NA,15),
#"MSE"=rep(NA,15),
#"Coverage"=rep(NA,15)
#) #since there are 15 PH scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  Results_PH_Cox_nocens_filled[j,]<-Results_PH_Cox_nocens[i,] #order the numerical robustness to the index
#}
#Results_PH_Cox_nocens<-Results_PH_Cox_nocens_filled #change name to use following code


#PO cases
#NPPR
#Results_PO_NPPR_nocens_filled<-data.frame(
#"Variable"=rep(1,15),
#"Mean"=rep(NA,15),
#"Bias"=rep(NA,15),
#"MSE"=rep(NA,15),
#"Coverage"=rep(NA,15)
#) #since there are 15 PO scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  Results_PO_NPPR_nocens_filled[j,]<-Results_PO_NPPR_nocens[i,] #order the numerical robustness to the index
#}
#Results_PO_NPPR_nocens<-Results_PO_NPPR_nocens_filled #change name to use following code

#LL
#Results_PO_LL_nocens_filled<-data.frame(
#"Variable"=rep(1,15),
#"Mean"=rep(NA,15),
#"Bias"=rep(NA,15),
#"MSE"=rep(NA,15),
#"Coverage"=rep(NA,15)
#) #since there are 15 PO scenarios (with censoring)
#we will fill up with NAs
#for(i in 1:length(1:9)){ #please change 1:9 to vector with indices of interest
#  j<-c(1:9)[i]
#  Results_PO_LL_nocens_filled[j,]<-Results_PO_LL_nocens[i,] #order the numerical robustness to the index
#}
#Results_PO_LL_nocens<-Results_PO_LL_nocens_filled #change name to use following code

#Now we can run the tables.
#Bias Table B
#MSE Table B
TableB<-data.frame(
  "Underlying Model"=c(rep("PPR",15),rep("Weibull PH",15),rep("LL PO",15)),
  "Effect"=rep(c(rep(0,3),rep(0.5,3),rep(0.25,3),rep(-0.25,3),rep(-0.5,3)),3),
  "Participants"=rep(c(500,100,50),15),
  "Bias_NPPR"=c(Results_PR_NPPR_nocens$Bias[1:6],Results_PR_NPPR_nocens$Bias[10:12],Results_PR_NPPR_nocens$Bias[13:15],Results_PR_NPPR_nocens$Bias[7:9],Results_PH_NPPR_nocens$Bias[1:6],Results_PH_NPPR_nocens$Bias[10:12],Results_PH_NPPR_nocens$Bias[13:15],Results_PH_NPPR_nocens$Bias[7:9],Results_PO_NPPR_nocens$Bias[1:6],Results_PO_NPPR_nocens$Bias[10:12],Results_PO_NPPR_nocens$Bias[13:15],Results_PO_NPPR_nocens$Bias[7:9]),
  "Bias_Cox"=c(Results_PR_Cox_nocens$Bias[1:6],Results_PR_Cox_nocens$Bias[10:12],Results_PR_Cox_nocens$Bias[13:15],Results_PR_Cox_nocens$Bias[7:9],Results_PH_Cox_nocens$Bias[1:6],Results_PH_Cox_nocens$Bias[10:12],Results_PH_Cox_nocens$Bias[13:15],Results_PH_Cox_nocens$Bias[7:9],rep("-",15)),
  "Bias_LL"=c(Results_PR_LL_nocens$Bias[1:6],Results_PR_LL_nocens$Bias[10:12],Results_PR_LL_nocens$Bias[13:15],Results_PR_LL_nocens$Bias[7:9],rep("-",15),Results_PO_LL_nocens$Bias[1:6],Results_PO_LL_nocens$Bias[10:12],Results_PO_LL_nocens$Bias[13:15],Results_PO_LL_nocens$Bias[7:9]),
  "MSE_NPPR"=c(Results_PR_NPPR_nocens$MSE[1:6],Results_PR_NPPR_nocens$MSE[10:12],Results_PR_NPPR_nocens$MSE[13:15],Results_PR_NPPR_nocens$MSE[7:9],Results_PH_NPPR_nocens$MSE[1:6],Results_PH_NPPR_nocens$MSE[10:12],Results_PH_NPPR_nocens$MSE[13:15],Results_PH_NPPR_nocens$MSE[7:9],Results_PO_NPPR_nocens$MSE[1:6],Results_PO_NPPR_nocens$MSE[10:12],Results_PO_NPPR_nocens$MSE[13:15],Results_PO_NPPR_nocens$MSE[7:9]),
  "MSE_Cox"=c(Results_PR_Cox_nocens$MSE[1:6],Results_PR_Cox_nocens$MSE[10:12],Results_PR_Cox_nocens$MSE[13:15],Results_PR_Cox_nocens$MSE[7:9],Results_PH_Cox_nocens$MSE[1:6],Results_PH_Cox_nocens$MSE[10:12],Results_PH_Cox_nocens$MSE[13:15],Results_PH_Cox_nocens$MSE[7:9],rep("-",15)),
  "MSE_LL"=c(Results_PR_LL_nocens$MSE[1:6],Results_PR_LL_nocens$MSE[10:12],Results_PR_LL_nocens$MSE[13:15],Results_PR_LL_nocens$MSE[7:9],rep("-",15),Results_PO_LL_nocens$MSE[1:6],Results_PO_LL_nocens$MSE[10:12],Results_PO_LL_nocens$MSE[13:15],Results_PO_LL_nocens$MSE[7:9])
)
write.csv2(TableB,file="./simulation/tables/TableB.csv")


#Coverage Table F
TableF<-data.frame(
  "Underlying Model"=c(rep("PPR",15),rep("Weibull PH",15),rep("LL PO",15)),
  "Effect"=rep(c(rep(0,3),rep(0.5,3),rep(0.25,3),rep(-0.25,3),rep(-0.5,3)),3),
  "Participants"=rep(c(500,100,50),15),
  "Coverage_NPPR"=c(Results_PR_NPPR_nocens$Coverage[1:6],Results_PR_NPPR_nocens$Coverage[10:12],Results_PR_NPPR_nocens$Coverage[13:15],Results_PR_NPPR_nocens$Coverage[7:9],Results_PH_NPPR_nocens$Coverage[1:6],Results_PH_NPPR_nocens$Coverage[10:12],Results_PH_NPPR_nocens$Coverage[13:15],Results_PH_NPPR_nocens$Coverage[7:9],Results_PO_NPPR_nocens$Coverage[1:6],Results_PO_NPPR_nocens$Coverage[10:12],Results_PO_NPPR_nocens$Coverage[13:15],Results_PO_NPPR_nocens$Coverage[7:9]),
  "Coverage_Cox"=c(Results_PR_Cox_nocens$Coverage[1:6],Results_PR_Cox_nocens$Coverage[10:12],Results_PR_Cox_nocens$Coverage[13:15],Results_PR_Cox_nocens$Coverage[7:9],Results_PH_Cox_nocens$Coverage[1:6],Results_PH_Cox_nocens$Coverage[10:12],Results_PH_Cox_nocens$Coverage[13:15],Results_PH_Cox_nocens$Coverage[7:9],rep("-",15)),
  "Coverage_LL"=c(round(Results_PR_LL_nocens$Coverage,2)[1:6],round(Results_PR_LL_nocens$Coverage,2)[10:12],round(Results_PR_LL_nocens$Coverage,2)[13:15],round(Results_PR_LL_nocens$Coverage,2)[7:9],rep("-",15),round(Results_PO_LL_nocens$Coverage,2)[1:6],round(Results_PO_LL_nocens$Coverage,2)[10:12],round(Results_PO_LL_nocens$Coverage,2)[13:15],round(Results_PO_LL_nocens$Coverage,2)[7:9])
)
write.csv2(TableF,file="./simulation/tables/TableF.csv")

#####Save Worspace####
save.image(file = "./simulation/NPPR_simulation_results.RData")
