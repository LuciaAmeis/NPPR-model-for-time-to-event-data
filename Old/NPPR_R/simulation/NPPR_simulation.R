#########################################################################
#                                                                       #
#    A non-parametric proportional risk model to assess a treatment     #
#                    effect in time-to-event data                       #
#                                                                       #
#         L. Ameis, O. Kuss, A. Hoyer and K. Moellenhoff                #
#                                                                       #
#                           SIMULATION                                    #
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
setwd("Please define")

#load workspace
load("NPPR_base.RData")

#seed
set.seed(60)

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

#Table I Numerical Robustness for PH data
TableI<-data.frame(
  "Effect"=c(rep(0,9),rep(0.5,9),rep(0.25,9),rep(-0.25,9),rep(-0.5,9)),
  "Censoring (%)"=rep(c(rep(30,3),rep(50,3),rep(70,3)),5),
  "Participants"=rep(c(500,100,50),15),
  "NPPR"=c(1000-NA_PO_NPPR[1:18],1000-NA_PO_NPPR[28:36],1000-NA_PO_NPPR[37:45],1000-NA_PO_NPPR[19:27]),
  "LL"=c(NA_PO_LL[1:18],NA_PO_LL[28:36],NA_PO_LL[37:45],NA_PO_LL[19:27])
)
write.csv2(TableI,file="./simulation/tables/TableI.csv")


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
  data<-Out_PR_NPPR[[i]]
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
  data<-Out_PR_NPPR_nocens[[i]]
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

#Of note, again we need to rearange the order to fit the paper
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
save.image(file = "./simulation/NPPR_simulation.RData")
