#########################################################################
#                                                                       #
#    A non-parametric proportional risk model to assess a treatment     #
#                    effect in time-to-event data                       #
#                                                                       #
#         L. Ameis, O. Kuss, A. Hoyer and K. Moellenhoff                #
#                                                                       #
#                           SIMULATION                                  #
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

#seed
set.seed(60)

#######################
####Parallelization####
#######################
#Preperation of the parallelization
n.cores<-parallel::detectCores()-1   #number of cores
my.cluster<-parallel::makeCluster(   #specify cluster
  n.cores,
  type="PSOCK"
)
###########################################
####Simulation of the data (Section 3) ####
###########################################
m1<-1000 #number of simulations per scenario
participants<-c(500,100,50)#for 500, 100 and 50 participants respectively


#first we will manually test for the parameters of the UD that result in the desired
#censoring rate 
#the parameters are listed in Table A (Supplementary Material)
#afterwards we will simulate studies for 500, 100 and 50 participants
#The respective scale parameter of the simulated treatment groups calculated from the 
#desired true are listed in Table 1


#PR data
#b=0
b1<-0

#30%
#censoring parameters
rateC_EU_DapaHF_0_censSC30<-c(0,175)  
#test censoring rate
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_0_censSC30,prob=0.5,b=b1)          
#to store simulations
Sim_PR_b0_30cens<-list()
#simualate
for(i in 1:3){
  Sim_PR_b0_30cens[[i]]<-Sim_EUSC(m=m1,N=participants[i],param=ML_EU$par,rateC =rateC_EU_DapaHF_0_censSC30,prob=0.5,b=b1)
}

#50%
rateC_EU_DapaHF_0_censSC50<-c(0,100)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_0_censSC50,prob=0.5,b=b1)
Sim_PR_b0_50cens<-list()
for(i in 1:3){
  Sim_PR_b0_50cens[[i]]<-Sim_EUSC(m=m1,N=participants[i],param=ML_EU$par,rateC =rateC_EU_DapaHF_0_censSC50,prob=0.5,b=b1)
}

#70%
rateC_EU_DapaHF_0_censSC70<-c(0,60)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_0_censSC70,prob=0.5,b=b1)
Sim_PR_b0_70cens<-list()
for(i in 1:3){
  Sim_PR_b0_70cens[[i]]<-Sim_EUSC(m=m1,N=participants[i],param=ML_EU$par,rateC =rateC_EU_DapaHF_0_censSC70,prob=0.5,b=b1)
}


#b=0.5
b1<-0.5

#30%
rateC_EU_DapaHF_.5_censSC30<-c(0,250)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_.5_censSC30,prob=0.5,b=b1)
Sim_PR_b0.5_30cens<-list()
for(i in 1:3){
  Sim_PR_b0.5_30cens[[i]]<-Sim_EUSC(m=m1,N=participants[i],param=ML_EU$par,rateC =rateC_EU_DapaHF_.5_censSC30,prob=0.5,b=b1)
}

#50%
rateC_EU_DapaHF_.5_censSC50<-c(0,140)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_.5_censSC50,prob=0.5,b=b1)
Sim_PR_b0.5_50cens<-list()
for(i in 1:3){
  Sim_PR_b0.5_50cens[[i]]<-Sim_EUSC(m=m1,N=participants[i],param=ML_EU$par,rateC =rateC_EU_DapaHF_.5_censSC50,prob=0.5,b=b1)
}

#70%
rateC_EU_DapaHF_.5_censSC70<-c(0,75)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_.5_censSC70,prob=0.5,b=b1)
Sim_PR_b0.5_70cens<-list()
for(i in 1:3){
  Sim_PR_b0.5_70cens[[i]]<-Sim_EUSC(m=m1,N=participants[i],param=ML_EU$par,rateC =rateC_EU_DapaHF_.5_censSC70,prob=0.5,b=b1)
}


#b=-0.5
b1<--0.5

#30%
rateC_EU_DapaHF_n.5_censSC30<-c(0,130)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_n.5_censSC30,prob=0.5,b=b1)
Sim_PR_bn0.5_30cens<-list()
for(i in 1:3){
  Sim_PR_bn0.5_30cens[[i]]<-Sim_EUSC(m=m1,N=participants[i],param=ML_EU$par,rateC =rateC_EU_DapaHF_n.5_censSC30,prob=0.5,b=b1)
}

#50%
rateC_EU_DapaHF_n.5_censSC50<-c(0,80)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_n.5_censSC50,prob=0.5,b=b1)
Sim_PR_bn0.5_50cens<-list()
for(i in 1:3){
  Sim_PR_bn0.5_50cens[[i]]<-Sim_EUSC(m=m1,N=participants[i],param=ML_EU$par,rateC =rateC_EU_DapaHF_n.5_censSC50,prob=0.5,b=b1)
}

#70
rateC_EU_DapaHF_n.5_censSC70<-c(0,40)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_n.5_censSC70,prob=0.5,b=b1)
Sim_PR_bn0.5_70cens<-list()
for(i in 1:3){
  Sim_PR_bn0.5_70cens[[i]]<-Sim_EUSC(m=m1,N=participants[i],param=ML_EU$par,rateC =rateC_EU_DapaHF_n.5_censSC70,prob=0.5,b=b1)
}


#b=0.25
b1<-0.25 

#30%
rateC_EU_DapaHF_.25_censSC30<-c(0,200)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_.25_censSC30,prob=0.5,b=b1)
Sim_PR_b0.25_30cens<-list()
for(i in 1:3){
  Sim_PR_b0.25_30cens[[i]]<-Sim_EUSC(m=m1,N=participants[i],param=ML_EU$par,rateC =rateC_EU_DapaHF_.25_censSC30,prob=0.5,b=b1)
}

#50%
rateC_EU_DapaHF_.25_censSC50<-c(0,125)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_.25_censSC50,prob=0.5,b=b1)
Sim_PR_b0.25_50cens<-list()
for(i in 1:3){
  Sim_PR_b0.25_50cens[[i]]<-Sim_EUSC(m=m1,N=participants[i],param=ML_EU$par,rateC =rateC_EU_DapaHF_.25_censSC50,prob=0.5,b=b1)
}

#70%
rateC_EU_DapaHF_.25_censSC70<-c(0,70)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_.25_censSC70,prob=0.5,b=b1)
Sim_PR_b0.25_70cens<-list()
for(i in 1:3){
  Sim_PR_b0.25_70cens[[i]]<-Sim_EUSC(m=m1,N=participants[i],param=ML_EU$par,rateC =rateC_EU_DapaHF_.25_censSC70,prob=0.5,b=b1)
}

#b=-0.25
b1<--0.25
#30
rateC_EU_DapaHF_n.25_censSC30<-c(0,150)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_n.25_censSC30,prob=0.5,b=b1)
Sim_PR_bn0.25_30cens<-list()
for(i in 1:3){
  Sim_PR_bn0.25_30cens[[i]]<-Sim_EUSC(m=m1,N=participants[i],param=ML_EU$par,rateC =rateC_EU_DapaHF_n.25_censSC30,prob=0.5,b=b1)
}

#50
rateC_EU_DapaHF_n.25_censSC50<-c(0,90)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_n.25_censSC50,prob=0.5,b=b1)
Sim_PR_bn0.25_50cens<-list()
for(i in 1:3){
  Sim_PR_bn0.25_50cens[[i]]<-Sim_EUSC(m=m1,N=participants[i],param=ML_EU$par,rateC =rateC_EU_DapaHF_n.25_censSC50,prob=0.5,b=b1)
}

#70
rateC_EU_DapaHF_n.25_censSC70<-c(0,50)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_n.25_censSC70,prob=0.5,b=b1)
Sim_PR_bn0.25_70cens<-list()
for(i in 1:3){
  Sim_PR_bn0.25_70cens[[i]]<-Sim_EUSC(m=m1,N=participants[i],param=ML_EU$par,rateC =rateC_EU_DapaHF_n.25_censSC70,prob=0.5,b=b1)
}


#Collect all simulations
Sim_PR<-c(Sim_PR_b0_30cens,Sim_PR_b0_50cens,Sim_PR_b0_70cens,
          Sim_PR_b0.5_30cens,Sim_PR_b0.5_50cens,Sim_PR_b0.5_70cens,
          Sim_PR_bn0.5_30cens,Sim_PR_bn0.5_50cens,Sim_PR_bn0.5_70cens,
          Sim_PR_b0.25_30cens,Sim_PR_b0.25_50cens,Sim_PR_b0.25_70cens,
          Sim_PR_bn0.25_30cens,Sim_PR_bn0.25_50cens,Sim_PR_bn0.25_70cens
)

#remove the redundant simulations
rm(Sim_PR_b0_30cens,Sim_PR_b0_50cens,Sim_PR_b0_70cens,
   Sim_PR_b0.5_30cens,Sim_PR_b0.5_50cens,Sim_PR_b0.5_70cens,
   Sim_PR_bn0.5_30cens,Sim_PR_bn0.5_50cens,Sim_PR_bn0.5_70cens,
   Sim_PR_b0.25_30cens,Sim_PR_b0.25_50cens,Sim_PR_b0.25_70cens,
   Sim_PR_bn0.25_30cens,Sim_PR_bn0.25_50cens,Sim_PR_bn0.25_70cens)

#PH
#b=0
b1<-0

#30%
rateC_PH_DapaHF_0_censSC30<-c(0,275)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_0_censSC30,prob=0.5,b=b1)
Sim_PH_b0_30cens<-list()
for(i in 1:3){
  Sim_PH_b0_30cens[[i]]<-Sim_PHSC(m=m1,N=participants[i],param=ML_0$par,rateC =rateC_PH_DapaHF_0_censSC30,prob=0.5,b=b1)
}


#50%
rateC_PH_DapaHF_0_censSC50<-c(0,130)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_0_censSC50,prob=0.5,b=b1)
Sim_PH_b0_50cens<-list()
for(i in 1:3){
  Sim_PH_b0_50cens[[i]]<-Sim_PHSC(m=m1,N=participants[i],param=ML_0$par,rateC =rateC_PH_DapaHF_0_censSC50,prob=0.5,b=b1)
}


#70%
rateC_PH_DapaHF_0_censSC70<-c(0,60)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_0_censSC70,prob=0.5,b=b1)
Sim_PH_b0_70cens<-list()
for(i in 1:3){
  Sim_PH_b0_70cens[[i]]<-Sim_PHSC(m=m1,N=participants[i],param=ML_0$par,rateC =rateC_PH_DapaHF_0_censSC70,prob=0.5,b=b1)
}


#b=0.5
b1<-0.5 

#30%
rateC_PH_DapaHF_.5_censSC30<-c(0,350)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_.5_censSC30,prob=0.5,b=b1)
Sim_PH_b0.5_30cens<-list()
for(i in 1:3){
  Sim_PH_b0.5_30cens[[i]]<-Sim_PHSC(m=m1,N=participants[i],param=ML_0$par,rateC =rateC_PH_DapaHF_.5_censSC30,prob=0.5,b=b1)
}

#50%
rateC_PH_DapaHF_.5_censSC50<-c(0,150)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_.5_censSC50,prob=0.5,b=b1)
Sim_PH_b0.5_50cens<-list()
for(i in 1:3){
  Sim_PH_b0.5_50cens[[i]]<-Sim_PHSC(m=m1,N=participants[i],param=ML_0$par,rateC =rateC_PH_DapaHF_.5_censSC50,prob=0.5,b=b1)
}

#70%
rateC_PH_DapaHF_.5_censSC70<-c(0,75)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_.5_censSC70,prob=0.5,b=b1)
Sim_PH_b0.5_70cens<-list()
for(i in 1:3){
  Sim_PH_b0.5_70cens[[i]]<-Sim_PHSC(m=m1,N=participants[i],param=ML_0$par,rateC =rateC_PH_DapaHF_.5_censSC70,prob=0.5,b=b1)
}


#b=-0.5
b1<--0.5

#30%
rateC_PH_DapaHF_n.5_censSC30<-c(0,220)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_n.5_censSC30,prob=0.5,b=b1)
Sim_PH_bn0.5_30cens<-list()
for(i in 1:3){
  Sim_PH_bn0.5_30cens[[i]]<-Sim_PHSC(m=m1,N=participants[i],param=ML_0$par,rateC =rateC_PH_DapaHF_n.5_censSC30,prob=0.5,b=b1)
}

#50%
rateC_PH_DapaHF_n.5_censSC50<-c(0,110)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_n.5_censSC50,prob=0.5,b=b1)
Sim_PH_bn0.5_50cens<-list()
for(i in 1:3){
  Sim_PH_bn0.5_50cens[[i]]<-Sim_PHSC(m=m1,N=participants[i],param=ML_0$par,rateC =rateC_PH_DapaHF_n.5_censSC50,prob=0.5,b=b1)
}

#70%
rateC_PH_DapaHF_n.5_censSC70<-c(0,45)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_n.5_censSC70,prob=0.5,b=b1)
Sim_PH_bn0.5_70cens<-list()
for(i in 1:3){
  Sim_PH_bn0.5_70cens[[i]]<-Sim_PHSC(m=m1,N=participants[i],param=ML_0$par,rateC =rateC_PH_DapaHF_n.5_censSC70,prob=0.5,b=b1)
}


#b=0.25
b1<-0.25 

#30%
rateC_PH_DapaHF_.25_censSC30<-c(0,320)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_.25_censSC30,prob=0.5,b=b1)
Sim_PH_b0.25_30cens<-list()
for(i in 1:3){
  Sim_PH_b0.25_30cens[[i]]<-Sim_PHSC(m=m1,N=participants[i],param=ML_0$par,rateC =rateC_PH_DapaHF_.25_censSC30,prob=0.5,b=b1)
}

#50%
rateC_PH_DapaHF_.25_censSC50<-c(0,155)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_.25_censSC50,prob=0.5,b=b1)
Sim_PH_b0.25_50cens<-list()
for(i in 1:3){
  Sim_PH_b0.25_50cens[[i]]<-Sim_PHSC(m=m1,N=participants[i],param=ML_0$par,rateC =rateC_PH_DapaHF_.25_censSC50,prob=0.5,b=b1)
}

#70%
rateC_PH_DapaHF_.25_censSC70<-c(0,70)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_.25_censSC70,prob=0.5,b=b1)
Sim_PH_b0.25_70cens<-list()
for(i in 1:3){
  Sim_PH_b0.25_70cens[[i]]<-Sim_PHSC(m=m1,N=participants[i],param=ML_0$par,rateC =rateC_PH_DapaHF_.25_censSC70,prob=0.5,b=b1)
}


#b=-0.25
b1<--0.25

#30%
rateC_PH_DapaHF_n.25_censSC30<-c(0,250)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_n.25_censSC30,prob=0.5,b=b1)
Sim_PH_bn0.25_30cens<-list()
for(i in 1:3){
  Sim_PH_bn0.25_30cens[[i]]<-Sim_PHSC(m=m1,N=participants[i],param=ML_0$par,rateC =rateC_PH_DapaHF_n.25_censSC30,prob=0.5,b=b1)
}

#50%
rateC_PH_DapaHF_n.25_censSC50<-c(0,115)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_n.25_censSC50,prob=0.5,b=b1)
Sim_PH_bn0.25_50cens<-list()
for(i in 1:3){
  Sim_PH_bn0.25_50cens[[i]]<-Sim_PHSC(m=m1,N=participants[i],param=ML_0$par,rateC =rateC_PH_DapaHF_n.25_censSC50,prob=0.5,b=b1)
}

#70%
rateC_PH_DapaHF_n.25_censSC70<-c(0,50)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_n.25_censSC70,prob=0.5,b=b1)
Sim_PH_bn0.25_70cens<-list()
for(i in 1:3){
  Sim_PH_bn0.25_70cens[[i]]<-Sim_PHSC(m=m1,N=participants[i],param=ML_0$par,rateC =rateC_PH_DapaHF_n.25_censSC70,prob=0.5,b=b1)
}

#Collect all simulations
Sim_PH<-c(Sim_PH_b0_30cens,Sim_PH_b0_50cens,Sim_PH_b0_70cens,
          Sim_PH_b0.5_30cens,Sim_PH_b0.5_50cens,Sim_PH_b0.5_70cens,
          Sim_PH_bn0.5_30cens,Sim_PH_bn0.5_50cens,Sim_PH_bn0.5_70cens,
          Sim_PH_b0.25_30cens,Sim_PH_b0.25_50cens,Sim_PH_b0.25_70cens,
          Sim_PH_bn0.25_30cens,Sim_PH_bn0.25_50cens,Sim_PH_bn0.25_70cens
)

#remove redundant simulations
rm(Sim_PH_b0_30cens,Sim_PH_b0_50cens,Sim_PH_b0_70cens,
   Sim_PH_b0.5_30cens,Sim_PH_b0.5_50cens,Sim_PH_b0.5_70cens,
   Sim_PH_bn0.5_30cens,Sim_PH_bn0.5_50cens,Sim_PH_bn0.5_70cens,
   Sim_PH_b0.25_30cens,Sim_PH_b0.25_50cens,Sim_PH_b0.25_70cens,
   Sim_PH_bn0.25_30cens,Sim_PH_bn0.25_50cens,Sim_PH_bn0.25_70cens)

#Set another seed
set.seed(29) 

#PO data
#b=0                                                                                                 #true value
b1<-0

#30%                                                                                         #censoring rate
rateC_PO_DapaHF_0_censSC30<-c(0,500)                                                                 #parameters for censoring
simulPOLL_SC_censrate(N=10000,param=ML_LL$par,rateC=rateC_PO_DapaHF_0_censSC30,prob=0.5,b=b1)          #test censoring rate
Sim_PO_b0_30cens<-list()
for(i in 1:3){
  Sim_PO_b0_30cens[[i]]<-Sim_POSC(m=m1,N=participants[i],param=ML_LL$par,rateC =rateC_PO_DapaHF_0_censSC30,prob=0.5,b=b1)
}

#50%
rateC_PO_DapaHF_0_censSC50<-c(0,175)
simulPOLL_SC_censrate(N=10000,param=ML_LL$par,rateC=rateC_PO_DapaHF_0_censSC50,prob=0.5,b=b1)
Sim_PO_b0_50cens<-list()
for(i in 1:3){
  Sim_PO_b0_50cens[[i]]<-Sim_POSC(m=m1,N=participants[i],param=ML_LL$par,rateC =rateC_PO_DapaHF_0_censSC50,prob=0.5,b=b1)
}

#70%
rateC_PO_DapaHF_0_censSC70<-c(0,60)
simulPOLL_SC_censrate(N=10000,param=ML_LL$par,rateC=rateC_PO_DapaHF_0_censSC70,prob=0.5,b=b1)
Sim_PO_b0_70cens<-list()
for(i in 1:3){
  Sim_PO_b0_70cens[[i]]<-Sim_POSC(m=m1,N=participants[i],param=ML_LL$par,rateC =rateC_PO_DapaHF_0_censSC70,prob=0.5,b=b1)
}


#b=0.5
b1<-0.5 

#30%
rateC_PO_DapaHF_.5_censSC30<-c(0,650)
simulPOLL_SC_censrate(N=10000,param=ML_LL$par,rateC=rateC_PO_DapaHF_.5_censSC30,prob=0.5,b=b1)
Sim_PO_b0.5_30cens<-list()
for(i in 1:3){
  Sim_PO_b0.5_30cens[[i]]<-Sim_POSC(m=m1,N=participants[i],param=ML_LL$par,rateC =rateC_PO_DapaHF_.5_censSC30,prob=0.5,b=b1)
}

#50%
rateC_PO_DapaHF_.5_censSC50<-c(0,225)
simulPOLL_SC_censrate(N=10000,param=ML_LL$par,rateC=rateC_PO_DapaHF_.5_censSC50,prob=0.5,b=b1)
Sim_PO_b0.5_50cens<-list()
for(i in 1:3){
  Sim_PO_b0.5_50cens[[i]]<-Sim_POSC(m=m1,N=participants[i],param=ML_LL$par,rateC =rateC_PO_DapaHF_.5_censSC50,prob=0.5,b=b1)
}

#70%
rateC_PO_DapaHF_.5_censSC70<-c(0,75)
simulPOLL_SC_censrate(N=10000,param=ML_LL$par,rateC=rateC_PO_DapaHF_.5_censSC70,prob=0.5,b=b1)
Sim_PO_b0.5_70cens<-list()
for(i in 1:3){
  Sim_PO_b0.5_70cens[[i]]<-Sim_POSC(m=m1,N=participants[i],param=ML_LL$par,rateC =rateC_PO_DapaHF_.5_censSC70,prob=0.5,b=b1)
}


#b=-0.5
b1<--0.5

#30%
rateC_PO_DapaHF_n.5_censSC30<-c(0,400)
simulPOLL_SC_censrate(N=10000,param=ML_LL$par,rateC=rateC_PO_DapaHF_n.5_censSC30,prob=0.5,b=b1)
Sim_PO_bn0.5_30cens<-list()
for(i in 1:3){
  Sim_PO_bn0.5_30cens[[i]]<-Sim_POSC(m=m1,N=participants[i],param=ML_LL$par,rateC =rateC_PO_DapaHF_n.5_censSC30,prob=0.5,b=b1)
}

#50%
rateC_PO_DapaHF_n.5_censSC50<-c(0,125)
simulPOLL_SC_censrate(N=10000,param=ML_LL$par,rateC=rateC_PO_DapaHF_n.5_censSC50,prob=0.5,b=b1)
Sim_PO_bn0.5_50cens<-list()
for(i in 1:3){
  Sim_PO_bn0.5_50cens[[i]]<-Sim_POSC(m=m1,N=participants[i],param=ML_LL$par,rateC =rateC_PO_DapaHF_n.5_censSC50,prob=0.5,b=b1)
}

#70%
rateC_PO_DapaHF_n.5_censSC70<-c(0,50)
simulPOLL_SC_censrate(N=10000,param=ML_LL$par,rateC=rateC_PO_DapaHF_n.5_censSC70,prob=0.5,b=b1)
Sim_PO_bn0.5_70cens<-list()
for(i in 1:3){
  Sim_PO_bn0.5_70cens[[i]]<-Sim_POSC(m=m1,N=participants[i],param=ML_LL$par,rateC =rateC_PO_DapaHF_n.5_censSC70,prob=0.5,b=b1)
}


#b=0.25
b1<-0.25 

#30%
rateC_PO_DapaHF_.25_censSC30<-c(0,550)
simulPOLL_SC_censrate(N=10000,param=ML_LL$par,rateC=rateC_PO_DapaHF_.25_censSC30,prob=0.5,b=b1)
Sim_PO_b0.25_30cens<-list()
for(i in 1:3){
  Sim_PO_b0.25_30cens[[i]]<-Sim_POSC(m=m1,N=participants[i],param=ML_LL$par,rateC =rateC_PO_DapaHF_.25_censSC30,prob=0.5,b=b1)
}

#50%
rateC_PO_DapaHF_.25_censSC50<-c(0,200)
simulPOLL_SC_censrate(N=10000,param=ML_LL$par,rateC=rateC_PO_DapaHF_.25_censSC50,prob=0.5,b=b1)
Sim_PO_b0.25_50cens<-list()
for(i in 1:3){
  Sim_PO_b0.25_50cens[[i]]<-Sim_POSC(m=m1,N=participants[i],param=ML_LL$par,rateC =rateC_PO_DapaHF_.25_censSC50,prob=0.5,b=b1)
}

#70%
rateC_PO_DapaHF_.25_censSC70<-c(0,75)
simulPOLL_SC_censrate(N=10000,param=ML_LL$par,rateC=rateC_PO_DapaHF_.25_censSC70,prob=0.5,b=b1)
Sim_PO_b0.25_70cens<-list()
for(i in 1:3){
  Sim_PO_b0.25_70cens[[i]]<-Sim_POSC(m=m1,N=participants[i],param=ML_LL$par,rateC =rateC_PO_DapaHF_.25_censSC70,prob=0.5,b=b1)
}


#b=-0.25
b1<--0.25

#30%
rateC_PO_DapaHF_n.25_censSC30<-c(0,400)
simulPOLL_SC_censrate(N=10000,param=ML_LL$par,rateC=rateC_PO_DapaHF_n.25_censSC30,prob=0.5,b=b1)
Sim_PO_bn0.25_30cens<-list()
for(i in 1:3){
  Sim_PO_bn0.25_30cens[[i]]<-Sim_POSC(m=m1,N=participants[i],param=ML_LL$par,rateC =rateC_PO_DapaHF_n.25_censSC30,prob=0.5,b=b1)
}

#50%
rateC_PO_DapaHF_n.25_censSC50<-c(0,150)
simulPOLL_SC_censrate(N=10000,param=ML_LL$par,rateC=rateC_PO_DapaHF_n.25_censSC50,prob=0.5,b=b1)
Sim_PO_bn0.25_50cens<-list()
for(i in 1:3){
  Sim_PO_bn0.25_50cens[[i]]<-Sim_POSC(m=m1,N=participants[i],param=ML_LL$par,rateC =rateC_PO_DapaHF_n.25_censSC50,prob=0.5,b=b1)
}

#70%
rateC_PO_DapaHF_n.25_censSC70<-c(0,55)
simulPOLL_SC_censrate(N=10000,param=ML_LL$par,rateC=rateC_PO_DapaHF_n.25_censSC70,prob=0.5,b=b1)
Sim_PO_bn0.25_70cens<-list()
for(i in 1:3){
  Sim_PO_bn0.25_70cens[[i]]<-Sim_POSC(m=m1,N=participants[i],param=ML_LL$par,rateC =rateC_PO_DapaHF_n.25_censSC70,prob=0.5,b=b1)
}

#Collect simulations
Sim_PO<-c(Sim_PO_b0_30cens,Sim_PO_b0_50cens,Sim_PO_b0_70cens,
          Sim_PO_b0.5_30cens,Sim_PO_b0.5_50cens,Sim_PO_b0.5_70cens,
          Sim_PO_bn0.5_30cens,Sim_PO_bn0.5_50cens,Sim_PO_bn0.5_70cens,
          Sim_PO_b0.25_30cens,Sim_PO_b0.25_50cens,Sim_PO_b0.25_70cens,
          Sim_PO_bn0.25_30cens,Sim_PO_bn0.25_50cens,Sim_PO_bn0.25_70cens
)

#remove redundant simulations
rm(Sim_PO_b0_30cens,Sim_PO_b0_50cens,Sim_PO_b0_70cens,
   Sim_PO_b0.5_30cens,Sim_PO_b0.5_50cens,Sim_PO_b0.5_70cens,
   Sim_PO_bn0.5_30cens,Sim_PO_bn0.5_50cens,Sim_PO_bn0.5_70cens,
   Sim_PO_b0.25_30cens,Sim_PO_b0.25_50cens,Sim_PO_b0.25_70cens,
   Sim_PO_bn0.25_30cens,Sim_PO_bn0.25_50cens,Sim_PO_bn0.25_70cens)

#No Censoring
#set another seed
set.seed(678)

#PR data
#b=0
b1<-0
Sim_PR_b0_nocens<-list()
for(i in 1:3){
  Sim_PR_b0_nocens[[i]]<-Sim_EUnocens(m=m1,N=participants[i],param=ML_EU$par,prob=0.5,b=b1)
}

#b=0.5
b1<-0.5
Sim_PR_b0.5_nocens<-list()
for(i in 1:3){
  Sim_PR_b0.5_nocens[[i]]<-Sim_EUnocens(m=m1,N=participants[i],param=ML_EU$par,prob=0.5,b=b1)
}

#b=-0.5
b1<--0.5
Sim_PR_bn0.5_nocens<-list()
for(i in 1:3){
  Sim_PR_bn0.5_nocens[[i]]<-Sim_EUnocens(m=m1,N=participants[i],param=ML_EU$par,prob=0.5,b=b1)
}

#b=0.25
b1<-0.25
Sim_PR_b0.25_nocens<-list()
for(i in 1:3){
  Sim_PR_b0.25_nocens[[i]]<-Sim_EUnocens(m=m1,N=participants[i],param=ML_EU$par,prob=0.5,b=b1)
}

#b=-0.25
b1<--0.25
Sim_PR_bn0.25_nocens<-list()
for(i in 1:3){
  Sim_PR_bn0.25_nocens[[i]]<-Sim_EUnocens(m=m1,N=participants[i],param=ML_EU$par,prob=0.5,b=b1)
}

#collect Simulations
Sim_PR_nocens<-c(Sim_PR_b0_nocens,
                 Sim_PR_b0.5_nocens,
                 Sim_PR_bn0.5_nocens,
                 Sim_PR_b0.25_nocens,
                 Sim_PR_bn0.25_nocens
)

#remove redundant simulations
rm(Sim_PR_b0_nocens,
   Sim_PR_b0.5_nocens,
   Sim_PR_bn0.5_nocens,
   Sim_PR_b0.25_nocens,
   Sim_PR_bn0.25_nocens)


#PH data
#b=0
b1<-0
Sim_PH_b0_nocens<-list()
for(i in 1:3){
  Sim_PH_b0_nocens[[i]]<-Sim_PHnocens(m=m1,N=participants[i],param=ML_0$par,prob=0.5,b=b1)
}

#b=0.5
b1<-0.5
Sim_PH_b0.5_nocens<-list()
for(i in 1:3){
  Sim_PH_b0.5_nocens[[i]]<-Sim_PHnocens(m=m1,N=participants[i],param=ML_0$par,prob=0.5,b=b1)
}

#b=-0.5
b1<--0.5
Sim_PH_bn0.5_nocens<-list()
for(i in 1:3){
  Sim_PH_bn0.5_nocens[[i]]<-Sim_PHnocens(m=m1,N=participants[i],param=ML_0$par,prob=0.5,b=b1)
}

#b=0.25
b1<-0.25
Sim_PH_b0.25_nocens<-list()
for(i in 1:3){
  Sim_PH_b0.25_nocens[[i]]<-Sim_PHnocens(m=m1,N=participants[i],param=ML_0$par,prob=0.5,b=b1)
}

#b=-0.25
b1<--0.25
Sim_PH_bn0.25_nocens<-list()
for(i in 1:3){
  Sim_PH_bn0.25_nocens[[i]]<-Sim_PHnocens(m=m1,N=participants[i],param=ML_0$par,prob=0.5,b=b1)
}


#collect simulations
Sim_PH_nocens<-c(Sim_PH_b0_nocens,
                 Sim_PH_b0.5_nocens,
                 Sim_PH_bn0.5_nocens,
                 Sim_PH_b0.25_nocens,
                 Sim_PH_bn0.25_nocens
)

#remove redundant simulations
rm(Sim_PH_b0_nocens,
   Sim_PH_b0.5_nocens,
   Sim_PH_bn0.5_nocens,
   Sim_PH_b0.25_nocens,
   Sim_PH_bn0.25_nocens)

#PO data
#b=0
b1<-0
Sim_PO_b0_nocens<-list()
for(i in 1:3){
  Sim_PO_b0_nocens[[i]]<-Sim_POnocens(m=m1,N=participants[i],param=ML_LL$par,prob=0.5,b=b1)
}

#b=0.5
b1<-0.5
Sim_PO_b0.5_nocens<-list()
for(i in 1:3){
  Sim_PO_b0.5_nocens[[i]]<-Sim_POnocens(m=m1,N=participants[i],param=ML_LL$par,prob=0.5,b=b1)
}

#b=-0.5
b1<--0.5
Sim_PO_bn0.5_nocens<-list()
for(i in 1:3){
  Sim_PO_bn0.5_nocens[[i]]<-Sim_POnocens(m=m1,N=participants[i],param=ML_LL$par,prob=0.5,b=b1)
}

#b=0.25
b1<-0.25
Sim_PO_b0.25_nocens<-list()
for(i in 1:3){
  Sim_PO_b0.25_nocens[[i]]<-Sim_POnocens(m=m1,N=participants[i],param=ML_LL$par,prob=0.5,b=b1)
}

#b=-0.25
b1<--0.25
Sim_PO_bn0.25_nocens<-list()
for(i in 1:3){
  Sim_PO_bn0.25_nocens[[i]]<-Sim_POnocens(m=m1,N=participants[i],param=ML_LL$par,prob=0.5,b=b1)
}


#collect simulations
Sim_PO_nocens<-c(Sim_PO_b0_nocens,
                 Sim_PO_b0.5_nocens,
                 Sim_PO_bn0.5_nocens,
                 Sim_PO_b0.25_nocens,
                 Sim_PO_bn0.25_nocens
)

#remove redundant simulations
rm(Sim_PO_b0_nocens,
   Sim_PO_b0.5_nocens,
   Sim_PO_bn0.5_nocens,
   Sim_PO_b0.25_nocens,
   Sim_PO_bn0.25_nocens)

#Storing parameters into tables
#Table 3 (distribution parameters)
#Of note, since there are different parametrizations of the PPR model used for
#the simulation of the data and the evaluate reparemetrization is necessary
Table3<-data.frame("beta"=c(0,0.5,0.25,-0.25,-0.5),
           "RR/HR/OR"=exp(-c(0,0.5,0.25,-0.25,-0.5)),
           "PPR_alpha"=rep(ML_EU$par[1],5),
           "PPR_beta0"=rep(1/ML_EU$par[2],5),
           "PPR_beta1"=exp(c(0,0.5,0.25,-0.25,-0.5))^(-1/ML_EU$par[1])*(1/ML_EU$par[2]),
           "WeibullPH_k"=rep(ML_0$par[1],times=5),
           "WeibullPH_lambda1"=ML_0$par[2]*exp(c(0,0.5,0.25,-0.25,-0.5))^(1/ML_0$par[1]),
           "WeibullPH_lambda0"=rep(ML_0$par[2],times=5),
           "LL_a"=rep(ML_LL$par[1],5),
           "LL_b0"=rep(ML_LL$par[2],5),
           "LL_b1"=ML_LL$par[2]*exp(c(0,0.5,0.25,-0.25,-0.5))^(1/ML_LL$par[1])
)
write.csv2(Table3,file="./simulation/tables/Table3.csv")


#Table A Supplement (parameters UD, censoring)
#PPR model
TableA_a<-data.frame(
  "Number of censored participants (%)"=c(rep(30,5),rep(50,5),rep(70,5)),
  "True underlying effect beta"=rep(c(0,0.5,0.25,-0.25,-0.5),3),
  "ParametersUD_min"=rep(0,15),
  "ParametersUD_max"=c(rateC_EU_DapaHF_0_censSC30[2],rateC_EU_DapaHF_.5_censSC30[2],
                       rateC_EU_DapaHF_.25_censSC30[2],rateC_EU_DapaHF_n.25_censSC30[2],
                       rateC_EU_DapaHF_n.5_censSC30[2],rateC_EU_DapaHF_0_censSC50[2],
                       rateC_EU_DapaHF_.5_censSC50[2],rateC_EU_DapaHF_.25_censSC50[2],
                       rateC_EU_DapaHF_n.25_censSC50[2],rateC_EU_DapaHF_n.5_censSC50[2],
                       rateC_EU_DapaHF_0_censSC70[2],rateC_EU_DapaHF_.5_censSC70[2],
                       rateC_EU_DapaHF_.25_censSC70[2],rateC_EU_DapaHF_n.25_censSC70[2],
                       rateC_EU_DapaHF_n.5_censSC70[2])
)
write.csv2(TableA_a,file="./simulation/tables/TableA_a.csv")

#WeibullPH model
TableA_b<-data.frame(
  "Number of censored participants (%)"=c(rep(30,5),rep(50,5),rep(70,5)),
  "True underlying effect beta"=rep(c(0,0.5,0.25,-0.25,-0.5),3),
  "ParametersUD_min"=rep(0,15),
  "ParametersUD_max"=c(rateC_PH_DapaHF_0_censSC30[2],rateC_PH_DapaHF_.5_censSC30[2],
                       rateC_PH_DapaHF_.25_censSC30[2],rateC_PH_DapaHF_n.25_censSC30[2],
                       rateC_PH_DapaHF_n.5_censSC30[2],rateC_PH_DapaHF_0_censSC50[2],
                       rateC_PH_DapaHF_.5_censSC50[2],rateC_PH_DapaHF_.25_censSC50[2],
                       rateC_PH_DapaHF_n.25_censSC50[2],rateC_PH_DapaHF_n.5_censSC50[2],
                       rateC_PH_DapaHF_0_censSC70[2],rateC_PH_DapaHF_.5_censSC70[2],
                       rateC_PH_DapaHF_.25_censSC70[2],rateC_PH_DapaHF_n.25_censSC70[2],
                       rateC_PH_DapaHF_n.5_censSC70[2])
)
write.csv2(TableA_b,file="./simulation/tables/TableA_b.csv")

#LL PO model
TableA_c<-data.frame(
  "Number of censored participants (%)"=c(rep(30,5),rep(50,5),rep(70,5)),
  "True underlying effect beta"=rep(c(0,0.5,0.25,-0.25,-0.5),3),
  "ParametersUD_min"=rep(0,15),
  "ParametersUD_max"=c(rateC_PO_DapaHF_0_censSC30[2],rateC_PO_DapaHF_.5_censSC30[2],
                       rateC_PO_DapaHF_.25_censSC30[2],rateC_PO_DapaHF_n.25_censSC30[2],
                       rateC_PO_DapaHF_n.5_censSC30[2],rateC_PO_DapaHF_0_censSC50[2],
                       rateC_PO_DapaHF_.5_censSC50[2],rateC_PO_DapaHF_.25_censSC50[2],
                       rateC_PO_DapaHF_n.25_censSC50[2],rateC_PO_DapaHF_n.5_censSC50[2],
                       rateC_PO_DapaHF_0_censSC70[2],rateC_PO_DapaHF_.5_censSC70[2],
                       rateC_PO_DapaHF_.25_censSC70[2],rateC_PO_DapaHF_n.25_censSC70[2],
                       rateC_PO_DapaHF_n.5_censSC70[2])
)
write.csv2(TableA_c,file="./simulation/tables/TableA_c.csv")

#############################################################################
####Application of the different models to the simulated data (Section 3)####
#############################################################################
#WARNING: Running all estimations will take some hours! The run time can be
#reduced by running only a subset of estimations
#In order to run a subset please change the index of the for loop below to the
#vector containing the indices of the studies of interest.
#The studies appear in the following order
#True value     Participants        Censoring       INDEX
#0                  500                 30            1
#0                  100                 30            2
#0                   50                 30            3
#0                  500                 50            4
#0                  100                 50            5
#0                   50                 50            6
#0                  500                 70            7
#0                  100                 70            8
#0                   50                 70            9
#0.5                500                 30            10
#0.5                100                 30            11
#0.5                 50                 30            12
#0.5                500                 50            13
#0.5                100                 50            14
#0.5                 50                 50            15
#0.5                500                 70            16
#0.5                100                 70            17
#0.5                 50                 70            18
#-0.5               500                 30            19
#-0.5               100                 30            20
#-0.5                50                 30            21
#-0.5               500                 50            22
#-0.5               100                 50            23
#-0.5                50                 50            24
#-0.5               500                 70            25
#-0.5               100                 70            26
#-0.5                50                 70            27
#0.25               500                 30            28
#0.25               100                 30            29
#0.25                50                 30            30
#0.25               500                 50            31
#0.25               100                 50            32
#0.25                50                 50            33
#0.25               500                 70            34
#0.25               100                 70            35
#0.25                50                 70            36
#-0.25              500                 30            37
#-0.25              100                 30            38
#-0.25               50                 30            39
#-0.25              500                 50            40
#-0.25              100                 50            41
#-0.25               50                 50            42
#-0.25              500                 70            43
#-0.25              100                 70            44
#-0.25               50                 70            45
#E.g. if only the PR cases with true value 0 evaluated using the NPPR estimator
#are of interest, change:
#Est_PR_NPPR<-list()
#for(i in 1:length(1:9)){
# x<-c(1:9)[i]
#  Sim<-Sim_PR[[x]]
#  Est_PR_NPPR[[i]]<-Sim_est(sims=Sim)
#}
#For instructions how to proceed with a subset of estimations see comments at the
#beginning of each following chapter
#OF NOTE: To reduce the runtime chose studies with fewer patients. The larger 
#the number of simulated patients, the longer the run time.

#PR data
#NPPR model
#empty list to fill with estimation of PR data using the NPPR model
Est_PR_NPPR<-list()
for(i in 1:length(Sim_PR)){#for every simulation generated
  Sim<-Sim_PR[[i]]#extract simulation
  Est_PR_NPPR[[i]]<-Sim_est(sims=Sim)#estimate and save
}

#PPR model
#fist we will have to set the starting values for optim
#we will use the true value, to test the NPPR model against the optimal used 
#competitor
#the true underlying beta as used for the data generation are
true.beta<-c(rep(0,times=9), #500,100,50 participants per 30%,50%,70% censoring each
             rep(0.5,times=9),
             rep(-0.5,times=9),
             rep(0.25,times=9),
             rep(-0.25,times=9)
)
start_PPR<-list()#empty list for start values
for(i in 1:length(true.beta)){#for every simulation
  b1<-true.beta[i]#extract true beta
  start_PPR[[i]]<-c(ML_EU$par[1],ML_EU$par[2],exp(-b1)^(-1/ML_EU$par[1])*ML_EU$par[2])    
  #save start values
  #the first two are
  #as estimated
  #from the DapaHF
  #control group
  #the third is
  #the scale parameter
  #of the treatment group
  #calculated from the 
  #true beta
}

#empty list to fill with estimation of PR data using the PPR model
Est_PR_PPR<-list()
for(i in 1:length(Sim_PR)){#for every simulation generated
  Sim<-Sim_PR[[i]]#extract simulation
  staASV<-start_PPR[[i]]#extract start value
  Est_PR_PPR[[i]]<-Sim_EU_ASV_Est(sims=Sim,start=staASV)#estimate and save
}

#Cox's model
Est_PR_Cox<-list()
for(i in 1:length(Sim_PR)){#for every simulation generated
  Sim<-Sim_PR[[i]]#extract simulation
  Est_PR_Cox[[i]]<-Sim_PH_est(sims=Sim)#estimate and save
}

#LL model
#we will use the true values of the similar PO case as true starting values
start_LL<-list()#empty list for start values
for(i in 1:length(true.beta)){#for every simulation
  b1<-true.beta[i]#extract true beta
  start_LL[[i]]<-c(ML_LL$par[1],ML_LL$par[2],exp(-b1)^(-1/ML_LL$par[1])*ML_LL$par[2])     
  #save start values
  #the first two are
  #as estimated
  #from the DapaHF
  #control group
  #the third is
  #the scale parameter
  #of the treatment group
  #calculated from the 
  #true beta
}

Est_PR_LL<-list()
for(i in 1:length(Sim_PR)){#for every simulation generated
  Sim<-Sim_PR[[i]]#extract simulation
  staV<-start_LL[[i]]#extract start value
  Est_PR_LL[[i]]<-Sim_PO_Est(sims=Sim,start = staV)#estimate and save
}


#PH data
#NPPR model
#empty list to fill with estimation of PR data using the NPPR model
Est_PH_NPPR<-list()
for(i in 1:length(Sim_PH)){#for every simulation generated
  Sim<-Sim_PH[[i]]#extract simulation
  Est_PH_NPPR[[i]]<-Sim_est(sims=Sim)#estimate and save
}

#Cox's model
Est_PH_Cox<-list()
for(i in 1:length(Sim_PH)){#for every simulation generated
  Sim<-Sim_PH[[i]]#extract simulation
  Est_PH_Cox[[i]]<-Sim_PH_est(sims=Sim)#estimate and save
}



#PO data
Est_PO_NPPR<-list()
for(i in 1:length(Sim_PO)){#for every simulation generated
  Sim<-Sim_PO[[i]]#extract simulation
  Est_PO_NPPR[[i]]<-Sim_est(sims=Sim)#estimate and save
}


#LL model
Est_PO_LL<-list()
for(i in 1:length(Sim_PO)){#for every simulation generated
  Sim<-Sim_PO[[i]]#extract simulation
  staV<-start_LL[[i]]#extract start value
  Est_PO_LL[[i]]<-Sim_PO_Est(sims=Sim,start = staV)#estimate and save
}


#No censoring
#WARNING: Here, the table of the simulation scenarios is slightly different.
#The studies appear in the following order
#True value     Participants         INDEX
#0                  500                  1
#0                  100                  2
#0                   50                  3
#0.5                500                  4
#0.5                100                  5
#0.5                 50                  6
#-0.5               500                  7
#-0.5               100                  8
#-0.5                50                  9
#0.25               500                 10
#0.25               100                 11
#0.25                50                 12
#-0.25              500                 13
#-0.25              100                 14
#-0.25               50                 15

#PR data
#NPPR model
#empty list to fill with estimation of PR data using the NPPR model
Est_PR_NPPR_nocens<-list()
for(i in 1:length(Sim_PR_nocens)){#for every simulation generated
  Sim<-Sim_PR_nocens[[i]]#extract simulation
  Est_PR_NPPR_nocens[[i]]<-Sim_est(sims=Sim)#estimate and save
}

#Cox's model
Est_PR_Cox_nocens<-list()
for(i in 1:length(Sim_PR_nocens)){#for every simulation generated
  Sim<-Sim_PR_nocens[[i]]#extract simulation
  Est_PR_Cox_nocens[[i]]<-Sim_PH_est(sims=Sim)#estimate and save
}


#We have to fit the start value list to the number of the nocens simulations
start_LL_nocens<-start_LL[c(1,2,3,10,11,12,19,20,21,28,29,30,37,38,39)]
#LL model
Est_PR_LL_nocens<-list()
for(i in 1:length(Sim_PR_nocens)){#for every simulation generated
  Sim<-Sim_PR_nocens[[i]]#extract simulation
  staV<-start_LL_nocens[[i]]#extract start value
  Est_PR_LL_nocens[[i]]<-Sim_PO_Est(sims=Sim,start = staV)#estimate and save
}


#PH data
#NPPR model
#empty list to fill with estimation of PR data using the NPPR model
Est_PH_NPPR_nocens<-list()
for(i in 1:length(Sim_PH_nocens)){#for every simulation generated
  Sim<-Sim_PH_nocens[[i]]#extract simulation
  Est_PH_NPPR_nocens[[i]]<-Sim_est(sims=Sim)#estimate and save
}


#Cox's model
Est_PH_Cox_nocens<-list()
for(i in 1:length(Sim_PH_nocens)){#for every simulation generated
  Sim<-Sim_PH_nocens[[i]]#extract simulation
  Est_PH_Cox_nocens[[i]]<-Sim_PH_est(sims=Sim)#estimate and save
}


#PO data
#NPPR model
#empty list to fill with estimation of PR data using the NPPR model
Est_PO_NPPR_nocens<-list()
for(i in 1:length(Sim_PO_nocens)){#for every simulation generated
  Sim<-Sim_PO_nocens[[i]]#extract simulation
  Est_PO_NPPR_nocens[[i]]<-Sim_est(sims=Sim)#estimate and save
}


#LL model
Est_PO_LL_nocens<-list()
for(i in 1:length(Sim_PO_nocens)){#for every simulation generated
  Sim<-Sim_PO_nocens[[i]]#extract simulation
  staV<-start_LL_nocens[[i]]#extract start value
  Est_PO_LL_nocens[[i]]<-Sim_PO_Est(sims=Sim,start = staV)#estimate and save
}



############################################################################
####Pre-processing of the estimation results (Section 3.3.4, Tables G-J)####
############################################################################
#WARNING: If only a subset of the simulation studies was analysed in the previous
#section, some of the code has to be adjusted in the following.
#Precisely, go to lines 1018, 1022, 1098, 1102, 1146 and 1150 and add the vector
#containing the indices of the simulations of interest as restriction.
#E.g. to only include the cases with a true effect of 0 
#Sim_PR_NPPR_nna<-Sim_PR[1:9]
#Est_PR_NPPR_nna<-Est_PR_NPPR
#for(i in 1:length(which(NA_PR_NPPR>0))){#for all cases
#  x<-which(NA_PR_NPPR>0)[i] #extract one case
#  Sim_PR_NPPR_nna[[x]]<-Sim_PR[1:9][[x]][!is.na(Est_PR_NPPR[[x]])]  
#  Est_PR_NPPR_nna[[x]]<-Est_PR_NPPR[[x]][!is.na(Est_PR_NPPR[[x]])]  
#}
#If which(NA_PR_NPPR>0) is empty is suffices to run only the first two lines.
#In this case please comment out the for loop.

#PR data
#NPPR model
#first we will check for NAs
NA_PR_NPPR<-numeric(length(Est_PR_NPPR))
for(i in 1:length(Est_PR_NPPR)){
  est<-Est_PR_NPPR[[i]]
  NA_PR_NPPR[i]<-sum(is.na(est))
}
#1000-result is collected in Table G under numerical robustness  and 
#discussed in section 3.3.4

#exclude NA cases
#List of all estimations and simulations in which we will exclude the NA cases
Sim_PR_NPPR_nna<-Sim_PR
Est_PR_NPPR_nna<-Est_PR_NPPR
for(i in 1:length(which(NA_PR_NPPR>0))){#for all cases
  x<-which(NA_PR_NPPR>0)[i] #extract one case
  Sim_PR_NPPR_nna[[x]]<-Sim_PR[[x]][!is.na(Est_PR_NPPR[[x]])]  
  Est_PR_NPPR_nna[[x]]<-Est_PR_NPPR[[x]][!is.na(Est_PR_NPPR[[x]])]  
}

#we will check the ranges of the estimations for all cases
range_PR_NPPR<-list()
for(i in 1:length(Est_PR_NPPR_nna)){
  est<-Est_PR_NPPR_nna[[i]]
  range_PR_NPPR[[i]]<-range(est)
}
#all ranges are in between the cutoff values at -3 and 3


#PPR model
#first we will test for the number of estimations with numerical problems
NA_PR_PPR<-numeric(length(Est_PR_PPR))
for(i in 1:length(Est_PR_PPR)){
  est<-Est_PR_PPR[[i]]
  NA_PR_PPR[i]<-Sim_EU_nna_test_log(est) 
}
#The results are listed in Table G under numerical robustness and discussed in section 3.3.4

#exclude the problematic cases and transform estimated parameters to -log(RR)
Est_PR_PPR_nna<-Est_PR_PPR
for(i in 1:length(Est_PR_PPR)){
  Est_PR_PPR_nna[[i]]<-Sim_EU_nna_log(Est_PR_PPR[[i]]) 
}


#Cox's model
#check for estimations indicating numerical problems
NA_PR_Cox<-numeric(length(Est_PR_Cox))
for(i in 1:length(Est_PR_Cox)){
  est<-Est_PR_Cox[[i]] #one simulation case
  NA_PR_Cox[i]<-length(which(estcox(est)>=(3)))
}
#1000-results are listed in Table G under numerical robustness and discussed in section 3.3.4

#exclude those with numerical problems
Est_PR_Cox_nna<-Est_PR_Cox
for(i in 1:length(which(NA_PR_Cox>0))){
  x<-which(NA_PR_Cox>0)[i]
  est<-Est_PR_Cox[[x]]
  Est_PR_Cox_nna[[x]]<-est[-which(estcox(est)>=(3))]
}


#LL model
#first we will test for the number of estimations with numerical problems
NA_PR_LL<-numeric(length(Est_PR_LL))
for(i in 1:length(Est_PR_LL)){
  est<-Est_PR_LL[[i]]
  NA_PR_LL[i]<-Sim_PO_nna_test_log(est) 
}
#The results are listed in Table G under numerical robustness and discussed in section 3.3.4

#exclude the problematic cases and transform estimated parameters to -log(RR)
Est_PR_LL_nna<-Est_PR_LL
for(i in 1:length(Est_PR_LL)){
  Est_PR_LL_nna[[i]]<-Sim_PO_nna_log(Est_PR_LL[[i]]) 
}



#PH data
#NPPR model
#first we will check for NAs
NA_PH_NPPR<-numeric(length(Est_PH_NPPR))
for(i in 1:length(Est_PH_NPPR)){
  est<-Est_PH_NPPR[[i]]
  NA_PH_NPPR[i]<-sum(is.na(est))
}
#1000-result is collected in Table H under numerical robustness
#here we will exclude the NA cases

#List of all estimations and simulations in which we will exclude the NA cases
Sim_PH_NPPR_nna<-Sim_PH
Est_PH_NPPR_nna<-Est_PH_NPPR
for(i in 1:length(which(NA_PH_NPPR>0))){#for all cases
  x<-which(NA_PH_NPPR>0)[i] #extract one case
  Sim_PH_NPPR_nna[[x]]<-Sim_PH[[x]][!is.na(Est_PH_NPPR[[x]])]  
  Est_PH_NPPR_nna[[x]]<-Est_PH_NPPR[[x]][!is.na(Est_PH_NPPR[[x]])]  
}

#we will check the ranges of the estimations for all cases
range_PH_NPPR<-list()
for(i in 1:length(Est_PH_NPPR_nna)){
  est<-Est_PH_NPPR_nna[[i]]
  range_PH_NPPR[[i]]<-range(est)
}
#all ranges are in between the cutoff values at -3 and 3


#Cox's model
#check for estimations indicating numerical problems
NA_PH_Cox<-numeric(length(Est_PH_Cox))
for(i in 1:length(Est_PH_Cox)){
  est<-Est_PH_Cox[[i]] #one simulation case
  NA_PH_Cox[i]<-length(which(estcox(est)>=(3)))
}
#1000-results are listed in Table H under numerical robustness and discussed in section 3.3.4

#exclude those with numerical problems
Est_PH_Cox_nna<-Est_PH_Cox
for(i in 1:length(which(NA_PH_Cox>0))){
  x<-which(NA_PH_Cox>0)[i]
  est<-Est_PH_Cox[[x]]
  Est_PH_Cox_nna[[x]]<-est[-which(estcox(est)>=(3))]
}



#PO data
#NPPR model
#first we will check for NAs
NA_PO_NPPR<-numeric(length(Est_PO_NPPR))
for(i in 1:length(Est_PO_NPPR)){
  est<-Est_PO_NPPR[[i]]
  NA_PO_NPPR[i]<-sum(is.na(est))
}
#1000-result is collected in Table I under numerical robustness
#here we will exclude the NA cases

#List of all estimations and simulations in which we will exclude the NA cases
Sim_PO_NPPR_nna<-Sim_PO
Est_PO_NPPR_nna<-Est_PO_NPPR
for(i in 1:length(which(NA_PO_NPPR>0))){#for all cases
  x<-which(NA_PO_NPPR>0)[i] #extract one case
  Sim_PO_NPPR_nna[[x]]<-Sim_PO[[x]][!is.na(Est_PO_NPPR[[x]])]  
  Est_PO_NPPR_nna[[x]]<-Est_PO_NPPR[[x]][!is.na(Est_PO_NPPR[[x]])]  
}

#we will check the ranges of the estimations for all cases
range_PO_NPPR<-list()
for(i in 1:length(Est_PO_NPPR_nna)){
  est<-Est_PO_NPPR_nna[[i]]
  range_PO_NPPR[[i]]<-range(est)
}
#all ranges are in between the cutoff values at -3 and 3


#LL model
#first we will test for the number of estimations with numerical problems
NA_PO_LL<-numeric(length(Est_PO_LL))
for(i in 1:length(Est_PO_LL)){
  est<-Est_PO_LL[[i]]
  NA_PO_LL[i]<-Sim_PO_nna_test_log(est) 
}
#The results are listed in Table I under numerical robustness and discussed in section 3.3.4

#exclude the problematic cases and transform estimated parameters to -log(RR)
Est_PO_LL_nna<-Est_PO_LL
for(i in 1:length(Est_PO_LL)){
  Est_PO_LL_nna[[i]]<-Sim_PO_nna_log(Est_PO_LL[[i]]) 
}


#No censoring
#Please note the incides for the no censoring cases

#PR data
#NPPR model
#first we will check for NAs
NA_PR_NPPR_nocens<-numeric(length(Est_PR_NPPR_nocens))
for(i in 1:length(Est_PR_NPPR_nocens)){
  est<-Est_PR_NPPR_nocens[[i]]
  NA_PR_NPPR_nocens[i]<-sum(is.na(est))
}
#1000-result is collected in Table J under numerical robustness  and discussed in section 3.3.4

#exclude NA cases
#List of all estimations and simulations in which we will exclude the NA cases
Sim_PR_NPPR_nocens_nna<-Sim_PR_nocens
Est_PR_NPPR_nocens_nna<-Est_PR_NPPR_nocens
#Since no case produces NAs we can leave it as is

#we will check the ranges of the estimations for all cases
range_PR_NPPR_nocens<-list()
for(i in 1:length(Est_PR_NPPR_nocens_nna)){
  est<-Est_PR_NPPR_nocens_nna[[i]]
  range_PR_NPPR_nocens[[i]]<-range(est)
}
#all ranges are in between the cutoff values at -3 and 3



#Cox's model
#check for estimations indicating numerical problems
NA_PR_Cox_nocens<-numeric(length(Est_PR_Cox_nocens))
for(i in 1:length(Est_PR_Cox_nocens)){
  est<-Est_PR_Cox_nocens[[i]] #one simulation case
  NA_PR_Cox_nocens[i]<-length(which(estcox(est)>=(3)))
}
#1000-results are listed in Table J under numerical robustness and discussed in section 3.3.4

#exclude those with numerical problems
Est_PR_Cox_nocens_nna<-Est_PR_Cox_nocens
#Since no case produces NAs we can leave it as is

#LL model
#first we will test for the number of estimations with numerical problems
NA_PR_LL_nocens<-numeric(length(Est_PR_LL_nocens))
for(i in 1:length(Est_PR_LL_nocens)){
  est<-Est_PR_LL_nocens[[i]]
  NA_PR_LL_nocens[i]<-Sim_PO_nna_test_log(est) 
}
#The results are listed in Table J under numerical robustness and discussed in section 3.3.4

#exclude the problematic cases and transform estimated parameters to -log(RR)
Est_PR_LL_nocens_nna<-Est_PR_LL_nocens
for(i in 1:length(Est_PR_LL_nocens)){
  Est_PR_LL_nocens_nna[[i]]<-Sim_PO_nna_log(Est_PR_LL_nocens[[i]]) 
}



#PH data
#NPPR model
#first we will check for NAs
NA_PH_NPPR_nocens<-numeric(length(Est_PH_NPPR_nocens))
for(i in 1:length(Est_PH_NPPR_nocens)){
  est<-Est_PH_NPPR_nocens[[i]]
  NA_PH_NPPR_nocens[i]<-sum(is.na(est))
}
#1000-result is collected in Table J under numerical robustness
#here we will exclude the NA cases

#List of all estimations and simulations in which we will exclude the NA cases
Sim_PH_NPPR_nocens_nna<-Sim_PH_nocens
Est_PH_NPPR_nocens_nna<-Est_PH_NPPR_nocens
#Since no case produced any NAs we can leave it as is

#we will check the ranges of the estimations for all cases
range_PH_NPPR_nocens<-list()
for(i in 1:length(Est_PH_NPPR_nocens_nna)){
  est<-Est_PH_NPPR_nocens_nna[[i]]
  range_PH_NPPR_nocens[[i]]<-range(est)
}
#all ranges are in between the cutoff values at -3 and 3


#Cox's model
#check for estimations indicating numerical problems
NA_PH_Cox_nocens<-numeric(length(Est_PH_Cox_nocens))
for(i in 1:length(Est_PH_Cox_nocens)){
  est<-Est_PH_Cox_nocens[[i]] #one simulation case
  NA_PH_Cox_nocens[i]<-length(which(estcox(est)>=(3)))
}
#1000-results are listed in Table J under numerical robustness and discussed in section 3.3.4

#exclude those with numerical problems
Est_PH_Cox_nocens_nna<-Est_PH_Cox_nocens
#SInce no case produces any NAs we can leave it as is


#PO data
#NPPR model
#first we will check for NAs
NA_PO_NPPR_nocens<-numeric(length(Est_PO_NPPR_nocens))
for(i in 1:length(Est_PO_NPPR_nocens)){
  est<-Est_PO_NPPR_nocens[[i]]
  NA_PO_NPPR_nocens[i]<-sum(is.na(est))
}
#1000-result is collected in Table J under numerical robustness
#here we will exclude the NA cases

#List of all estimations and simulations in which we will exclude the NA cases
Sim_PO_NPPR_nocens_nna<-Sim_PO_nocens
Est_PO_NPPR_nocens_nna<-Est_PO_NPPR_nocens
#Since no case produced any NAs we can leave it as is

#we will check the ranges of the estimations for all cases
range_PO_NPPR_nocens<-list()
for(i in 1:length(Est_PO_NPPR_nocens_nna)){
  est<-Est_PO_NPPR_nocens_nna[[i]]
  range_PO_NPPR_nocens[[i]]<-range(est)
}
#all ranges are in between the cutoff values at -3 and 3


#LL model
#first we will test for the number of estimations with numerical problems
NA_PO_LL_nocens<-numeric(length(Est_PO_LL_nocens))
for(i in 1:length(Est_PO_LL_nocens)){
  est<-Est_PO_LL_nocens[[i]]
  NA_PO_LL_nocens[i]<-Sim_PO_nna_test_log(est) 
}
#The results are listed in Table J under numerical robustness and discussed in section 3.3.4

#exclude the problematic cases and transform estimated parameters to -log(RR)
Est_PO_LL_nocens_nna<-Est_PO_LL_nocens
for(i in 1:length(Est_PO_LL_nocens)){
  Est_PO_LL_nocens_nna[[i]]<-Sim_PO_nna_log(Est_PO_LL_nocens[[i]]) 
}


##############################################################################
####Coverage of the NPPR, PPR and LL models (Section 3.3.3., Tables (C-F))####
##############################################################################
#We will now calculate the Coverage of the NPPR, PPR and LL models for all different scenarios
#the coverage of Cox's model will be included later, as the respective confidence intervals
#were already estimated in the estimation step

#WARNING:
#The computation using the NPPR estimator is 
#time-consuming.
#This depends mostly on the number of participants in the simulation study and "L"
#the number of repeats in the bootstrap (and also the number of cores used for 
# the palatalization)
#Setting the number of participants to a smaller value and reducing "L" in the function
#below 
#Of note, this can change the results drastically if L is set very small.

#WARNING 2:
#If only a subset of the simulation scenarios was regarded in the previous 
#subsections, the vector true.beta has to be adjusted. Change the subscription
#1:9 in the following example to the vector of indices of the simulation studies
#of interest and remove the hash tag
#true.beta<-true.beta[1:9]
#1:9 again represents the example of only regarding the cases with a true effect 
#of 0
#The same has to be done below to the vector number.df (see line 1357)

#PR case
#NPPR model
CI_PR_NPPR<-numeric(length(Sim_PR_NPPR_nna))
for(i in 1:length(Sim_PR_NPPR_nna)){#for every simulated case
  b1<-true.beta[i]#true value
  Sim<-Sim_PR_NPPR_nna[[i]]
  CI_PR_NPPR[i]<-Sim_CI_fast(sim=Sim,b=b1,L=500,alpha=0.05)    #coverage for the given case
}

#PPR model
#here we need the number of participants for the df for each scenario
number.df<-rep(c(500,100,50),times=15)
#number.df<-rep(c(500,100,50),times=15)[1:9] #if only the cases with true effect 0
#were regarded. If another subset was considered change 1:9 to the fitting vector
#of indices
CI_PR_PPR<-numeric(length(Est_PR_PPR))
for(i in 1:length(Est_PR_PPR)){#for every simulated case
  b1<-true.beta[i]#true value
  df1<-number.df[i]#number of participants
  Est<-Est_PR_PPR[[i]]
  CI_PR_PPR[i]<-Sim_EU_CI(ests=Sim_EU_nna_CI(Est),b=b1,df=df1)    #coverage for the given case
}


#LL model
CI_PR_LL<-numeric(length(Est_PR_LL))
for(i in 1:length(Est_PR_LL)){#for every simulated case
  b1<-true.beta[i]#true value
  df1<-number.df[i]#number of participants
  Est<-Est_PR_LL[[i]]
  CI_PR_LL[i]<-Sim_PO_CI(ests=Sim_PO_nna_CI(Est),b=b1,df=df1)    #coverage for the given case
}


#PH case
#NPPR model
CI_PH_NPPR<-numeric(length(Sim_PH_NPPR_nna))
for(i in 1:length(Sim_PH_NPPR_nna)){#for every simulated case
  b1<-true.beta[i]#true value
  Sim<-Sim_PH_NPPR_nna[[i]]
  CI_PH_NPPR[i]<-Sim_CI_fast(sim=Sim,b=b1,L=500,alpha=0.05)    #coverage for the given case
}


#PO case
#NPPR model
CI_PO_NPPR<-numeric(length(Sim_PO_NPPR_nna))
for(i in 1:length(Sim_PO_NPPR_nna)){#for every simulated case
  b1<-true.beta[i]#true value
  Sim<-Sim_PO_NPPR_nna[[i]]
  CI_PO_NPPR[i]<-Sim_CI_fast(sim=Sim,b=b1,L=500,alpha=0.05)    #coverage for the given case
}


#LL model
CI_PO_LL<-numeric(length(Est_PO_LL))
for(i in 1:length(Est_PO_LL)){#for every simulated case
  b1<-true.beta[i]#true value
  df1<-number.df[i]#number of participants
  Est<-Est_PO_LL[[i]]
  CI_PO_LL[i]<-Sim_PO_CI(ests=Sim_PO_nna_CI(Est),b=b1,df=df1)    #coverage for the given case
}


#No censoring
#We have to fit the true beta list to the number of the nocens simulations
true.beta_nocens<-true.beta[c(1,2,3,10,11,12,19,20,21,28,29,30,37,38,39)]
#Again please add a subscription with the indices of the simulations of interes
#e.g. (replace 1:3 with fitting vector)
#true.beta_nocens<-true.beta_nocens[1:3]

#PR case
#NPPR model
CI_PR_NPPR_nocens<-numeric(length(Sim_PR_NPPR_nocens_nna))
for(i in 1:length(Sim_PR_NPPR_nocens_nna)){#for every simulated case
  b1<-true.beta_nocens[i]#true value
  Sim<-Sim_PR_NPPR_nocens_nna[[i]]
  CI_PR_NPPR_nocens[i]<-Sim_CI_fast(sim=Sim,b=b1,L=500,alpha=0.05)    #coverage for the given case
}


#LL model
#here we need the number of participants for the df for each scenario
number.df_nocens<-rep(c(500,100,50),times=5)
#number.df_nocens<-number.df_nocens[1:3] #if only the cases with true effect 0 were
#regarded. If another subset was considered please change 1:3 to the vector of indices
CI_PR_LL_nocens<-numeric(length(Est_PR_LL_nocens))
for(i in 1:length(Est_PR_LL_nocens)){#for every simulated case
  b1<-true.beta_nocens[i]#true value
  df1<-number.df_nocens[i]#number of participants
  Est<-Est_PR_LL_nocens[[i]]
  CI_PR_LL_nocens[i]<-Sim_PO_CI(ests=Sim_PO_nna_CI(Est),b=b1,df=df1)    #coverage for the given case
}



#PH case
#NPPR model
CI_PH_NPPR_nocens<-numeric(length(Sim_PH_NPPR_nocens_nna))
for(i in 1:length(Sim_PH_NPPR_nocens_nna)){#for every simulated case
  b1<-true.beta_nocens[i]#true value
  Sim<-Sim_PH_NPPR_nocens_nna[[i]]
  CI_PH_NPPR_nocens[i]<-Sim_CI_fast(sim=Sim,b=b1,L=500,alpha=0.05)    #coverage for the given case
}



#PO case
#NPPR model
CI_PO_NPPR_nocens<-numeric(length(Sim_PO_NPPR_nocens_nna))
for(i in 1:length(Sim_PO_NPPR_nocens_nna)){#for every simulated case
  b1<-true.beta_nocens[i]#true value
  Sim<-Sim_PO_NPPR_nocens_nna[[i]]
  CI_PO_NPPR_nocens[i]<-Sim_CI_fast(sim=Sim,b=b1,L=500,alpha=0.05)    #coverage for the given case
}


#LL model
CI_PO_LL_nocens<-numeric(length(Est_PO_LL_nocens))
for(i in 1:length(Est_PO_LL_nocens)){#for every simulated case
  b1<-true.beta_nocens[i]#true value
  df1<-number.df_nocens[i]#number of participants
  Est<-Est_PO_LL_nocens[[i]]
  CI_PO_LL_nocens[i]<-Sim_PO_CI(ests=Sim_PO_nna_CI(Est),b=b1,df=df1)    #coverage for the given case
}


#######################
####Old simulations####
#######################
#As stated in the Read Me, the results of the NPPR and PPR estimators in case
#of PR data and of the NPPR estimator in case of PH data were produces using an 
#earlier R version. The code was the same, but the seed appears to produce different
#simulations. Load the following and re-run the part of the code above to produce 
#the results.
#load("./old_simulation/Sim_PR_old.RData")
#Sim_PR<-Sim_PR_old
#load("./old_simulation/Sim_PH_old.RData")
#Sim_PH<-Sim_PH_old
#####Save Worspace####
save.image(file = "./simulation/NPPR_simulation.RData")
save(true.beta,true.beta_nocens,rateC_EU_DapaHF_0_censSC30,
     NA_PR_NPPR,NA_PR_PPR,NA_PR_Cox,NA_PR_LL,
     NA_PH_NPPR,NA_PH_Cox,NA_PO_NPPR,NA_PO_LL,
     
     NA_PR_NPPR_nocens,NA_PR_Cox_nocens,NA_PR_LL_nocens,
     NA_PH_NPPR_nocens,NA_PH_Cox_nocens,
     NA_PO_NPPR_nocens,NA_PO_LL_nocens,
     
     Est_PR_NPPR_nna,Est_PR_PPR_nna,Est_PR_Cox_nna,Est_PR_LL_nna,
     Est_PH_NPPR_nna,Est_PH_Cox_nna,
     Est_PO_NPPR_nna,Est_PO_LL_nna,
     
     Est_PR_NPPR_nocens_nna,Est_PR_Cox_nocens_nna,Est_PR_LL_nocens_nna,
     Est_PH_NPPR_nocens_nna,Est_PH_Cox_nocens_nna,
     Est_PO_NPPR_nocens_nna,Est_PO_LL_nocens_nna,
     
     CI_PR_NPPR,CI_PR_PPR,CI_PR_LL,
     CI_PH_NPPR,CI_PO_NPPR,CI_PO_LL,
     
     CI_PR_NPPR_nocens,CI_PR_LL_nocens,
     CI_PH_NPPR_nocens,CI_PO_NPPR_nocens,CI_PO_LL_nocens,
     
     file="./simulation/NPPR_simulation_subset.RData")
