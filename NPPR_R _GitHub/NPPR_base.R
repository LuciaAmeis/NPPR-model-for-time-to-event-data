#########################################################################
#                                                                       #
#    A non-parametric proportional risk model to assess a treatment     #
#                    effect in time-to-event data                       #
#                                                                       #
#         L. Ameis, O. Kuss, A. Hoyer and K. Moellenhoff                #
#                                                                       #
#                             BASICS                                    #
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

#######################
####Parallelization####
#######################
#Preperation of the parallelization
n.cores<-parallel::detectCores()-1   #number of cores
my.cluster<-parallel::makeCluster(   #specify cluster
  n.cores,
  type="PSOCK"
)


#####################
#### Preparation ####
#####################

#Survival curve estimated with the Kaplan-Meier (KM) estimator
KM_surv<-function(t,daten,col){
  #t is time
  #daten is the data frame
  #col is a vector consisting of named values corresponding to the columns with 
  #information about survival time (col_time) and censoring (col_cens) in the data frame
  #e.g. col=c(col_time=1,col_cens=3)
  
  #Extract the values corresponding to the columns time and censoring
  x1=as.numeric(col["col_time"])
  x2=as.numeric(col["col_cens"])
  
  #KM-estimator included in the package survival
  km_fit <- survfit(Surv(daten[,x1], daten[,x2])~ 1, data=data.frame(daten))
  
  #Extract the vectors consisting of the time and survival probability 
  km_fit_time<-km_fit$time
  km_fit_surv<-km_fit$surv
  
  #length of the time vector
  m_time<-length(km_fit$time)
  
  #The KM-estimator is a step function
  #We define the probabilities of survival accordingly
  #For convenience we fix the values for t larger then max(k_fit_time)
  #at 0. Those values will not be needed for further analysis
  
  if(t<min(km_fit_time)){
    surv_prob<-1                     #1 before any event
  } else{ if(t>max(km_fit_time)){
    surv_prob<-0                     #set to 0 after all events (see above)
  }else{
    for(i in 1:m_time){              #Steps as estimated using survfit
      y1=km_fit_time[i]              #from event
      y2=km_fit_time[i+1]            #to event   
      y3=km_fit_surv[i]              #estimated survival probability
      ymax=max(km_fit_time)          #last event time
      y4=km_fit_surv[m_time]         #estimated probability at last event
      if(i==m_time){
        if(t==ymax){
          surv_prob<-y4              #assign last probability to max event  
        }} else{
          if(y1<=t& t<y2){
            surv_prob<-y3            #assign probability to each interval 
          }                          #between two events
          
        }}}}
                                     #return probability of the time point 
  return(surv_prob)                  #passed to function 
}                                      

#Pointwise standard error of the KM-estimation approximated by Greenwood's formula
KM_se<-function(t,daten,col){
  #t is time
  #daten is the data frame
  #col is a vector consisting of named values corresponding to the columns with 
  #information about survival time (col_time) and censoring (col_cens) in the data frame
  
  #Extract the values corresponding to the columns time and censoring
  x1=as.numeric(col["col_time"])
  x2=as.numeric(col["col_cens"])
  
  #KM-Estimator included in the package survival
  km_fit <- survfit(Surv(daten[,x1], daten[,x2])~ 1, data=data.frame(daten))
  
  #Extract the vectors consisting of the time and standard error (Greenwood's formula)
  km_fit_time<-km_fit$time
  km_fit_se<-km_fit$std.err
  
  #length of the time vector
  m_time<-length(km_fit$time)
  
  #The KM-estimator is a step function
  #We define the SE of survival at any time point accordingly
  #For convenience we fix the values for t larger then max(k_fit_time) to 1. Those will
  #not factor into further analysis
  
  if(t<min(km_fit_time)){
    surv_se<-1
  } else{ if(t>max(km_fit_time)){
    surv_se<-1
  }else{
    for(i in 1:m_time){                #Steps as estimated using survfit
      y1=km_fit_time[i]                #from event
      y2=km_fit_time[i+1]              #to event 
      y3=km_fit_se[i]                  #approximated SE
      ymax=max(km_fit_time)            #last event time
      y4=km_fit_se[m_time]             #SE at last event
      if(i==m_time){
        if(t==ymax){
          surv_se<-y4                  #assign last SE to max event 
        }} else{
          if(y1<=t& t<y2){             #assign SE to each interval 
            surv_se<-y3                #between two events
          } 
        }}}}
                                       #return SE of the time point  
  return(surv_se)                      #passed to function 
}

#Estimated CDFs using the KM estimator
F_est<-function(t,daten,col){
  #t is time
  #daten is the data frame
  #col is a vector consisting of named values corresponding to the columns with 
  #information about survival time (col_time) and censoring (col_cens) in the data frame
  #e.g. col=c(col_time=1,col_cens=3)
  1-KM_surv(t,daten,col)       #CDF=1-survival function
}

#Extracting the set \tilde{T}
T_jump<-function(events1,events0){
  #events1 set of events observed in group 1 (e.g. treatment)
  #events2 set of events observed in group 2 (e.g. control)

  t_min<-max(min(events1),min(events0))#lower bound
  t_max<-min(max(events1),max(events0))#upper bound
  
  events<-c(events1,events0)#all events
  
  #events between lower and upper bound
  events_bound<-subset(events,events>=t_min&events<=t_max)
  
  #return events between bounds
  return(events_bound)
}

#Estimator
Beta_est<-function(t_vector,Fn1,Fn0,W1,W0){
  #t_vector is the vector of event time points (later extracted with T_jump)
  #Fn1 and Fn0 are functions (KM-estimated CDFs)
  #W1 and W0 are the weights (KM-estimated SE)
  
  w<-numeric(length(t_vector))     #to store weights for each event
  d<-numeric(length(t_vector))     #to store RR for each event
  for(i in 1:length(t_vector)){
    t<-t_vector[i]                                #for every time point
    Wi<-(1/Fn1(t)^2)*W1(t)^2+(1/Fn0(t)^2)*W0(t)^2 #calculate the weight at time point
    d[i]<-(1/Wi)*log(Fn1(t)/Fn0(t))               #weighted log of estimated RR at time point
    
    w[i]<-1/Wi                                    #save the weights
  }
  W<-sum(w)                                       #Sum of weights                                   
  D<-sum(d)                                       #Sum of weighted log of estimated RR
  beta_est=-(((1/W)*D))                           #Combine to weighted mean
  
  #return estimated beta
  return(beta_est)              
}

#Estimator function (combined)
B_easy<-function(daten,col,col_target){
  #daten is the data frame
  #col is as above
  #col_target is the number corresponding to the column with the variable of 
  #interest
  
  #Extract the numbers corresponding to the columns survival time and censoring
  x1=as.numeric(col["col_time"])
  x2=as.numeric(col["col_cens"])
  #and target
  x3=col_target
  
  #create data frames for the different strata 
  group_test<-subset(daten,daten[,x3]==1)
  group_control<-subset(daten,daten[,x3]==0)
  
  #create the corresponding event vectors
  events_test<-subset(daten[,x1],daten[,x2]==1&daten[,x3]==1)
  events_control<-subset(daten[,x1],daten[,x2]==1&daten[,x3]==0)
  
  #create event time vector using T_jump
  t_jump<-T_jump(events1  = events_test,events0  = events_control)
  
  #create the cdf-estimates corresponding to the strata
  F_est_test<-function(t){
    F_est(t,group_test,col)
  }
  F_est_control<-function(t){
    F_est(t,group_control,col)
  }
  
  #weights
  W_test<-function(t){
    KM_se(t,daten=group_test,col=col)
  }
  W_control<-function(t){
    KM_se(t,daten=group_control,col=col)
  }
  
  #Use the estimator defined above
  if(length(t_jump)>0){
    beta_est<-Beta_est(t_vector=t_jump,Fn1=F_est_test,Fn0=F_est_control,W1=W_test,W0=W_control)
  } else{beta_est=NA #Exclude cases, if \tilde{T} is empty
  }
  
  #return estimated beta
  out<-c(beta_target=beta_est)
  return(out)
}

#Estimator function (fast)
#The same function as above without nesting of functions, therefore faster
B_easy_fast<-function(daten,col,col_target){
  #daten is the data frame
  #col is as above
  #col_target is the number corresponding to the column with the variable of 
  #interest
  
  #we extract the numbers corresponding to the columns survival time and censoring
  x1=as.numeric(col["col_time"])
  x2=as.numeric(col["col_cens"])
  #and target
  x3=col_target
  
  #stratify 
  group_test<-subset(daten,daten[,x3]==1)
  group_control<-subset(daten,daten[,x3]==0)
  
  #create the event time vectors
  events_test<-subset(daten[,x1],daten[,x2]==1&daten[,x3]==1)
  events_control<-subset(daten[,x1],daten[,x2]==1&daten[,x3]==0)
  
  t_min<-max(min(events_test),min(events_control))#lower bound
  t_max<-min(max(events_test),max(events_control))#upper bound
  
  events<-c(events_test,events_control)#all events
  
  #Event time vectors
  events_bound<-subset(events,events>=t_min&events<=t_max)#all events in bounds
  
  #CDF and weight function
  #test
  #We use the KM-Estimator included in the package survival
  km_fit_test<- survfit(Surv(group_test[,x1], group_test[,x2])~ 1, data=data.frame(group_test))
  
  #We extract the vectors consisting of the time and survival probability 
  km_fit_time_test<-km_fit_test$time
  km_fit_surv_test<-km_fit_test$surv
  km_fit_se_test<-km_fit_test$std.err
  
  #length of the time vector
  m_time_test<-length(km_fit_test$time)
  
  #control
  #We use the KM-Estimator included in the package survival
  km_fit_control<- survfit(Surv(group_control[,x1], group_control[,x2])~ 1, data=data.frame(group_control))
  
  #We extract the vectors consisting of the time and survival probability 
  km_fit_time_control<-km_fit_control$time
  km_fit_surv_control<-km_fit_control$surv
  km_fit_se_control<-km_fit_control$std.err
  
  #length of the time vector
  m_time_control<-length(km_fit_control$time)
  
  #estimating CDFs and with that beta as above
  Beta_est<-function(t){                   
    if(t<min(km_fit_time_test)){
      surv_prob_test<-1
      surv_se_test<-1
    } else{ if(t>max(km_fit_time_test)){
      surv_prob_test<-0
      surv_se_test<-1
    }else{
      for(i in 1:m_time_test){
        y1=km_fit_time_test[i]
        y2=km_fit_time_test[i+1]
        y3=km_fit_surv_test[i]
        y4=km_fit_se_test[i]
        ymax=max(km_fit_time_test)
        y5=km_fit_surv_test[m_time_test]
        y6=km_fit_se_test[m_time_test]
        if(i==m_time_test){
          if(t==ymax){
            surv_prob_test<-y5
            surv_se_test<-y6
          }} else{
            if(y1<=t& t<y2){
              surv_prob_test<-y3
              surv_se_test<-y4
            } 
          }}}}
    
    if(t<min(km_fit_time_control)){
      surv_prob_control<-1
      surv_se_control<-1
    } else{ if(t>max(km_fit_time_control)){
      surv_prob_control<-0
      surv_se_control<-1
    }else{
      for(i in 1:m_time_control){
        y1=km_fit_time_control[i]
        y2=km_fit_time_control[i+1]
        y3=km_fit_surv_control[i]
        y4=km_fit_se_control[i]
        ymax=max(km_fit_time_control)
        y5=km_fit_surv_control[m_time_control]
        y6=km_fit_se_control[m_time_control]
        if(i==m_time_control){
          if(t==ymax){
            surv_prob_control<-y5
            surv_se_control<-y6
          }} else{
            if(y1<=t& t<y2){
              surv_prob_control<-y3
              surv_se_control<-y4
            } 
          }}}}
    
    W<-(1/(1-surv_prob_test)^2)*surv_se_test^2+(1/(1-surv_prob_control)^2)*surv_se_control^2
    B<-(1/W)*log((1-surv_prob_test)/(1-surv_prob_control))
    
    return(c(B,(1/W)))
  }
  
  
  #use estimator for events
  if(length(events_bound)>0){
    #Use the estimator defined above
    beta_est_matrix<-sapply(X=events_bound, FUN=Beta_est, simplify = TRUE, USE.NAMES = TRUE)
    B_est<-sum(beta_est_matrix[1,])
    W_est<-sum(beta_est_matrix[2,])
    beta_est<--(1/W_est)*B_est
  } else{beta_est=NA}
   
  #return estimated beta 
  return(beta_est)
}

#Confidence interval for beta
CI_fast<-function(daten,col,col_target,L,alpha){
  #daten is the data set
  #L is the number of samples we draw for the boostrap
  #col and col_target as above
  #alpha is the confidence level
  
  #extracting numbers of columns
  x1=as.numeric(col["col_time"])
  x2=as.numeric(col["col_cens"])
  x3=col_target
  
  #fitting the estimator function to the situation of the data set
  B_easy_fast_sample<-function(daten){
    B_easy_fast(daten,col=c(col_time=x1,col_cens=x2),col_target=x3)
  }
  
  #parallelization
  #(Of note, windows implementation results in problems with set.seed. If implemented
  #this might cause slightly different results, as "sample" draws randomly)
  doParallel::registerDoParallel(cl=my.cluster)
  y<-foreach(1:L, .combine="rbind", .export=c("B_easy_fast"),.packages="survival") %dopar% {
    samp<-sample(1:nrow(daten),replace =T)               #draw sample with replacements
    x<-as.numeric(B_easy_fast_sample(daten[samp,]))      #estimate beta for sample
    x
  }
  stopImplicitCluster()
  
  #alpha for the lower/upper quantile
  y1<-alpha/2                                         
  y2<-1-alpha/2
  
  #excluding NA cases
  y_nna<-y[!is.na(y)]
  
  #calculate the empirical quantile
  Q<-quantile(y_nna,c(y1,y2))
  
  #return the cofindence limits (quantiles)
  return(Q)
}

#Number needed to treat
NNT_fast<-function(daten,col,col_target,tp){
  #daten is the data frame
  #col is as above
  #col_target is the number corresponding to the column with the variable of 
  #interest
  #tp time point at which NNT is to be evaluated at
  
  #extract the numbers corresponding to the columns survival time and censoring
  x1=as.numeric(col["col_time"])
  x2=as.numeric(col["col_cens"])
  #and target
  x3=col_target
  
  #create data frames for the different strata 
  group_test<-subset(daten,daten[,x3]==1)
  group_control<-subset(daten,daten[,x3]==0)
  
  #create the corresponding event vectors
  events_test<-subset(daten[,x1],daten[,x2]==1&daten[,x3]==1)
  events_control<-subset(daten[,x1],daten[,x2]==1&daten[,x3]==0)
  
  t_min<-max(min(events_test),min(events_control))#lower bound
  t_max<-min(max(events_test),max(events_control))#upper bound
  
  events<-c(events_test,events_control)#all events
  
  #event time points of interest (\tilde{T})
  events_bound<-subset(events,events>=t_min&events<=t_max)#
  
  #CDF and weight function
  #treatment
  #KM-estimator included in the package survival as above
  km_fit_test<- survfit(Surv(group_test[,x1], group_test[,x2])~ 1, data=data.frame(group_test))
  
  #Extract the vectors consisting of the time and survival probability 
  km_fit_time_test<-km_fit_test$time
  km_fit_surv_test<-km_fit_test$surv
  km_fit_se_test<-km_fit_test$std.err
  
  #length of the time vector
  m_time_test<-length(km_fit_test$time)
  
  #control
  #KM-estimator included in the package survival
  km_fit_control<- survfit(Surv(group_control[,x1], group_control[,x2])~ 1, data=data.frame(group_control))
  
  #extract the vectors consisting of the time and survival probability/SE 
  km_fit_time_control<-km_fit_control$time
  km_fit_surv_control<-km_fit_control$surv
  km_fit_se_control<-km_fit_control$std.err
  
  #length of the time vector
  m_time_control<-length(km_fit_control$time)
  
  #KM-estimated CDF of control
  F0_est<-function(t){
    if(t<min(km_fit_time_control)){
      surv_prob_control<-1
    } else{ if(t>max(km_fit_time_control)){
      surv_prob_control<-0
    }else{
      for(i in 1:m_time_control){
        y1=km_fit_time_control[i]
        y2=km_fit_time_control[i+1]
        y3=km_fit_surv_control[i]
        ymax=max(km_fit_time_control)
        y5=km_fit_surv_control[m_time_control]
        if(i==m_time_control){
          if(t==ymax){
            surv_prob_control<-y5
          }} else{
            if(y1<=t& t<y2){
              surv_prob_control<-y3
            }} } 
    }
      
    }
    return(1-surv_prob_control)
  }
  
  #estimator for beta as above
  Beta_est<-function(t){
    if(t<min(km_fit_time_test)){
      surv_prob_test<-1
      surv_se_test<-1
    } else{ if(t>max(km_fit_time_test)){
      surv_prob_test<-0
      surv_se_test<-1
    }else{
      for(i in 1:m_time_test){
        y1=km_fit_time_test[i]
        y2=km_fit_time_test[i+1]
        y3=km_fit_surv_test[i]
        y4=km_fit_se_test[i]
        ymax=max(km_fit_time_test)
        y5=km_fit_surv_test[m_time_test]
        y6=km_fit_se_test[m_time_test]
        if(i==m_time_test){
          if(t==ymax){
            surv_prob_test<-y5
            surv_se_test<-y6
          }} else{
            if(y1<=t& t<y2){
              surv_prob_test<-y3
              surv_se_test<-y4
            } 
          }}}}
    if(t<min(km_fit_time_control)){
      surv_prob_control<-1
      surv_se_control<-1
    } else{ if(t>max(km_fit_time_control)){
      surv_prob_control<-0
      surv_se_control<-1
    }else{
      for(i in 1:m_time_control){
        y1=km_fit_time_control[i]
        y2=km_fit_time_control[i+1]
        y3=km_fit_surv_control[i]
        y4=km_fit_se_control[i]
        ymax=max(km_fit_time_control)
        y5=km_fit_surv_control[m_time_control]
        y6=km_fit_se_control[m_time_control]
        if(i==m_time_control){
          if(t==ymax){
            surv_prob_control<-y5
            surv_se_control<-y6
          }} else{
            if(y1<=t& t<y2){
              surv_prob_control<-y3
              surv_se_control<-y4
            } 
          }}}}
    
    W<-(1/(1-surv_prob_test)^2)*surv_se_test^2+(1/(1-surv_prob_control)^2)*surv_se_control^2
    B<-(1/W)*log((1-surv_prob_test)/(1-surv_prob_control))
    
    return(c(B,(1/W)))
  }
  
  #estimator evaluated using the bound events
  if(length(events_bound)>0){
    #Use the estimator defined above
    beta_est_matrix<-sapply(X=events_bound, FUN=Beta_est, simplify = TRUE, USE.NAMES = TRUE)
    B_est<-sum(beta_est_matrix[1,])
    W_est<-sum(beta_est_matrix[2,])
    beta_est<--(1/W_est)*B_est
  } else{beta_est=NA}
  
  #CDF at time point tp
  F0_nnt<-F0_est(tp) 
  #Risk difference
  rd<-(1-exp(-beta_est))*F0_nnt
  #nnt
  nnt<-1/rd
  
  #return nnt at time point
  return(nnt)
}


#Confidence interval for nnt at time point tp
CI_nnt_fast<-function(daten,col,col_target,L,alpha,tp){
  #daten is the data set
  #L is the number of samples we draw
  #col and col_target as above
  #tp time point
  #alpha confidence level
  
  #extract columns as above
  x1=as.numeric(col["col_time"])
  x2=as.numeric(col["col_cens"])
  x3=col_target
  
  #fit NNT to sample
  NNT_fast_sample<-function(daten){
    NNT_fast(daten,col=c(col_time=x1,col_cens=x2),col_target=x3,tp=tp)
  }
  
  #again as seeds do cause problems with paralleization on Windows, results can slightly
  #differ
  doParallel::registerDoParallel(cl=my.cluster)
  y<-foreach(1:L, .combine="rbind", .export=c("NNT_fast"),.packages="survival") %dopar% {
    samp<-sample(1:nrow(daten),replace =T)        #simulate sample
    x<-as.numeric(NNT_fast_sample(daten[samp,]))  #estimate for every sample
    x
  }
  stopImplicitCluster()
  
  #quantile
  y1<-alpha/2
  y2<-1-alpha/2
  
  #exclude NAs
  y_nna<-y[!is.na(y)]
  
  #empirical quantile
  Q<-quantile(y_nna,c(y1,y2))
  
  #return limits of the confidence interval (quantiles)
  return(Q)
}


##################################################################################
#### Analysis of the DapaHF dataset primary outcome (Section 3.1.1 Supplement)####
##################################################################################
#Analysis of the data set regarding the primary outcome
#This data set is used as base for the simulation study
#Results of analysis using the NPPR estimator are recorded in Section 3.1.1.
#of the supplementary material

#load data set 
DapaHF <- read_sas("./data/dapa_hf_primary_outcome.sas7bdat")
DapaHF<-as.data.frame(DapaHF)
View(DapaHF)
#col_time is 3, col_cens is 1 and col_target is 2

#estimating beta and RR
estimate<-B_easy_fast(DapaHF,col=c(col_time=3,col_cens=1),col_target=2)
exp(-estimate)

#confidence interval for beta and RR
ci_fast<-CI_fast(DapaHF,col=c(col_time=3,col_cens=1),col_target=2,L=500,alpha=0.05)
exp(-ci_fast)

#nnt at time point 10
nnt_dapahf_10<-NNT_fast(DapaHF,col=c(col_time=3,col_cens=1),col_target=2,10)

#CI for nnt at time point 10
ci_nnt_dapahf_10<-CI_nnt_fast(DapaHF,col=c(col_time=3,col_cens=1),col_target=2,L=500,alpha=0.05,10)

#########################################################
####Preparation of the simulation study (Section 3.1)####
#########################################################
#As preparation for the simulation study we will fit the exponentiated uniform
#distribution, the Weibull distribution and the log-logistic distribution the the
#control group of the dataset
#The estimated parameters are listed in Table 1 under the respective shape parameters
#and the scale parameters of the respective control groups

#load data set if the analysis of the DapaHF dataset was not run before this
#DapaHF <- read_sas("./DapaHF/dapa_hf_primary_outcome.sas7bdat")
#DapaHF<-as.data.frame(DapaHF)
#View(DapaHF)

#Preparation of the PR data using the EU distribution
#likelihood
likeEU<-function(param, daten,col){
  #param shape and form parameters
  #the parametrization of the scale parameter in the used is slightly different 
  #than in the github functions
  #theta=1/(t1-t0)=1/t1, whereby t0 is always 0 here. 
  #daten data set
  #col=c(col_time,col_cens)
  #col_time column with survival times
  #col_cens column with information about censoring as above
  
  alpha<-as.numeric(param[1])  #shape
  t1<-as.numeric(param[2])     #scale
  
  #extract column information
  x1<-as.numeric(col["col_time"])
  x2<-as.numeric(col["col_cens"])
  
  #split data set into censored and not censored
  daten_event<-subset(daten,daten[,x2]==1)
  daten_cens<-subset(daten,daten[,x2]==0)
  
  #time vectors
  time_event<-daten_event[,x1]
  time_cens<-daten_cens[,x1]
  
  #number of events/censored observations
  n1<-nrow(daten_event)
  n2<-nrow(daten_cens)
  
  #distribution function
  EU_event<-function(t){
    dEU(t,alpha=alpha,a=0,b=t1)
  }
  
  #evaluate distribution function at event time points
  eu_event<-lapply(time_event,EU_event)
  eu_event_corr<-as.numeric(eu_event)*50 #due to the number of observations
                                         #the product would be to small for R
                                         #to work with, therefore we correct
                                         #this does not influence the estimated values
  c1<-prod(as.numeric(eu_event_corr))    #product part of the likelihood include events
  
  
  #CDF
  EU_cens<-function(t){
    (1-pEU(t,alpha=alpha,a=0,b=t1))
  }
  
  #eveluate CDF at censored time points
  eu_cens<-lapply(time_cens,EU_cens)
  c2<-prod(as.numeric(eu_cens))         #product part of the likelihood include censored 
                                        #time points
  
  #product combining results above to likelihood
  like<-c1*c2
  
  #return likelihood
  return(like)
}

#negative log likelihood
negloglikeEU<-function(param, daten,col){
  #param shape and form parameters
  #the parametrization of the scale parameter in the used is slightly different 
  #than in the github functions
  #theta=1/(t1-t0)=1/t1, whereby t0 is always 0 here. 
  #daten data set
  #col=c(col_time,col_cens)
  #col_time column with survival times
  #col_cens column with information about censoring as above
  
  -log(likeEU(param, daten,col))  #evaluate negative log likelihood
}


#Fitting a EU distribution to the control group of the data set
#starting values were determined by plotting  different EU distributions 
#against the KM-estimated CDF of the control group
ML_EU<-optim(par=c(0.8,130),negloglikeEU,daten=subset(DapaHF,DapaHF$Arm==0),col=c(col_time=3,col_cens=1),method="Nelder-Mead")

#Preparation of the PH data using the Weibull distribution
#likelihood
like_Weib<-function(param, daten,col){
  #param parameters to be estimated
  #daten data set
  #col=c(col_time,col_cens)
  #col_time column with survival times
  #col_cens column with information about censoring as above
  
  a1<-as.numeric(param[1]) #shape
  a2<-as.numeric(param[2]) #scale
  
  #extract column information
  x1<-as.numeric(col["col_time"])
  x2<-as.numeric(col["col_cens"])
  
  #split data set into censored and event
  daten_event<-subset(daten,daten[,x2]==1)
  daten_cens<-subset(daten,daten[,x2]==0)
  
  #time vectors
  time_event<-daten_event[,x1]
  time_cens<-daten_cens[,x1]
  
  #number of events/censored observations
  n1<-nrow(daten_event)
  n2<-nrow(daten_cens)
  
  #events in likelihood
  #distribution function
  Weib_event<-function(t){
    dweibull(t,shape=a1,scale=a2,log=F)
  }
  
  #evaluate distribution function at event time points
  weib_event<-lapply(time_event,Weib_event)
  weib_event_corr<-as.numeric(weib_event)*50  #due to the large number of observations
                                              #the product would be to small for R
                                              #to work with, therefore we correct
                                              #this does not influence the estimated values
  c1<-prod(as.numeric(weib_event_corr))       #product part of the likelihood include events
  
  #censored in likelihood
  #CDF
  Weib_cens<-function(t){
    (1-pweibull(t,shape=a1,scale=a2,lower.tail=T,log.p=F))
  }
  
  #eveluate CDF at censored time points
  weib_cens<-lapply(time_cens,Weib_cens)
  weib_cens_corr<-as.numeric(weib_cens)*1.45 #due to the number of observations
                                             #the product would be to small for R
                                             #to work with, therefore we correct
                                             #this does not influence the estimated values
  c2<-prod(as.numeric(weib_cens_corr))       #product part of the likelihood include censored 
                                             #time points
  
  #product combining results above to likelihood
  like<-c1*c2
  
  #return likelihood
  return(like)
}

#negative log likelihood
negloglike_Weib<-function(param, daten,col){
  #param parameters to be estimated
  #daten data set
  #col=c(col_time,col_cens)
  #col_time column with survival times
  #col_cens column with information about censoring as above
  -log(like_Weib(param,daten,col))  #evaluate negative log likelihood
}

#Fitting a Weibull distribution to the control group of the data set
#starting values were determined by plotting  different Weibull distributions 
#against the KM-estimated CDF of the control group
ML_0<-optim(par=c(1.2,70),fn=negloglike_Weib,daten=subset(DapaHF,DapaHF$Arm==0),col=c(col_time=3,col_cens=1),lower = c(0, 0))


#likelihood
like_LL<-function(param, daten,col){
  #param parameters to be estimated
  #daten data set
  #col=c(col_time=,col_cens=)
  #col_time column with survival times
  #col_cens column with information about censoring as above
  
  a1<-as.numeric(param[1]) #shape
  a2<-as.numeric(param[2]) #scale
  
  
  #extract column information
  x1<-as.numeric(col["col_time"])
  x2<-as.numeric(col["col_cens"])
  
  #split data set into censored and event
  daten_event<-subset(daten,daten[,x2]==1)
  daten_cens<-subset(daten,daten[,x2]==0)
  
  #time vectors
  time_event<-daten_event[,x1]
  time_cens<-daten_cens[,x1]
  
  #number of events/censored observations
  n1<-nrow(daten_event)
  n2<-nrow(daten_cens)
  
  #events in likelihood
  #distribution
  LL_event<-function(t){
    dllogis(t,shape=a1,scale=a2,log=F)
  }
  
  #eveluate distribution function at event time points
  LL_event<-lapply(time_event,LL_event)
  LL_event_corr<-as.numeric(LL_event)*50 #due to the large number of observations
                                         #the product would be to small for R
                                         #to work with, therefore we correct
                                         #this does not influence the estimated values
  c1<-prod(as.numeric(LL_event_corr))    #product part of the likelihood include events
  
  #censored in likelihood
  #CDF
  LL_cens<-function(t){
    (1-pllogis(t,shape=a1,scale=a2,lower.tail=T,log.p=F))
  }
  
  #eveluate CDF at censored time point
  LL_cens<-lapply(time_cens,LL_cens)
  c2<-prod(as.numeric(LL_cens))          #product part of the likelihood include censored 
                                         #time points
  
  #product combining results above to likelihood
  like<-c1*c2
  
  #return likelihood  
  return(like)
}

#negative log likelihood
negloglike_LL<-function(param, daten,col){
  #param parameters to be estimated
  #daten data set
  #col=c(col_time=,col_cens=)
  #col_time column with survival times
  #col_cens column with information about censoring as above
  -log(like_LL(param,daten,col))  #evaluate negative log likelihood
}

#Preparation of the PH data using the Weibull distribution
#Fitting a LL distribution to the control group of the data set
#starting values were determined by plotting  different LL distributions 
#against the KM-estimated CDF of the control group
ML_LL<-optim(par=c(1,70),fn=negloglike_LL,daten=subset(DapaHF,DapaHF$Arm==0),col=c(col_time=3,col_cens=1))



#Next we will construct the function later used to create the simulated data

#PR case
#simulation of one study
simulEU_SC<-function(N,param, rateC,prob, b){
  #N is the number of participants
  #param are shape and scale for the EU distribution 
  #rateC parameters for censoring
  #prob is a vector of probabilities for the Bernoulli distributions, it determines 
  #the group distribution treatment:control
  #b is a vector with the real values of the parameters we wish to estimate later
  
  a2<-as.numeric(param[1]) #shape
  a1<-as.numeric(param[2]) #scale
  
  #if we do not state otherwise we assume that patients are assignment to the 
  #groups with equal probability
  if(missing(prob)){prob=0.5}
  
  #group assignment
  expressions<-rbern(N, prob=prob)
  
  #to store list sub data.frames, one for each of the two groups
  datalist<-list()
  
  #vector of possible group assignments
  expr_vector<-c(1,0)
  for(i in 1:2){#for every possible group
    #indicating the group     
    pers_expr<-expr_vector[i]
    
    #only regard the subset assigned to the current group
    single_expr<-subset(expressions,expressions==pers_expr)
    
    #number of patient in this group
    Nn<-length(single_expr)
    
    #calculating the scale parameter according to the group
    b1<-exp(-b*pers_expr)^(-1/a2)*a1
    
    #draw survival times 
    survival<-rEU(n=Nn,alpha=a2,a=0,b = b1)
    
    #add to the sub data set
    daten_single_expr<-data.frame(expressions=single_expr,survival=survival)
    
    #add to list of data frames
    datalist[[i]]<-daten_single_expr
  }
  
  #combining the data sets
  sim_all<-datalist[[1]]
  sim_all<-rbind(sim_all,datalist[[2]])
  
  #censoring
  #extracting the parameters for the UD 
  minC<-rateC[1]   #always 0
  maxC<-rateC[2]
  
  #draw censoring times
  cens<-runif(n=N,min=minC,max=maxC) 
  
  #add to data frame
  sim_all$cens<-cens
  
  #calculate time and status from survival and censoring times
  time <- pmin(sim_all$survival, sim_all$cens)
  #determine status
  status <- as.numeric(sim_all$survival <= sim_all$cens) 
  
  #extracting the important information
  sim<-data.frame(
    id=1:N,                #simulated patient ID
    time=time,             #time of event/censoring
    status=status,         #status
    expr=sim_all[1]        #group assignment
  )
  
  #data set with only the information we need
  return(sim)
}

#to test out correct parameters for the UD we use the following function
simulEU_SC_censrate<-function(N,param, rateC,prob, b){
  #parameters as above
  #proceed as above until censoring times are drawn
  
  a2<-as.numeric(param[1]) 
  a1<-as.numeric(param[2]) 
  
  if(missing(prob)){prob=0.5}
  
  expressions<-rbern(N, prob=prob)
  
  datalist<-list()
  
  expr_vector<-c(1,0)
  for(i in 1:2){
    pers_expr<-expr_vector[i]
    
    single_expr<-subset(expressions,expressions==pers_expr)
    
    Nn<-length(single_expr)
    
    b1<-exp(-b*pers_expr)^(-1/a2)*a1
    
    survival<-rEU(n=Nn,alpha=a2,a=0,b = b1)
    
    daten_single_expr<-data.frame(expressions=single_expr,survival=survival)
    
    datalist[[i]]<-daten_single_expr
  }
  
  sim_all<-datalist[[1]]
  sim_all<-rbind(sim_all,datalist[[2]])
  
  minC<-rateC[1]
  maxC<-rateC[2]
  cens<-runif(n=N,min=minC,max=maxC) 
  
  sim_all$cens<-cens
  
  #now we can calculate the percentage of censoring
  Cens<-as.numeric(sim_all$cens<sim_all$survival)
  CensRate<-sum(Cens)/N
  
  #returne percentage
  return(CensRate)
}

#multiple simulations
Sim_EUSC<-function(m,N,param,rateC,prob,b){
  #m number of simulations
  #N,param,rateC,prob,b for the single simulation
  
  #empty list for the simulated data sets
  sims<-list()
  
  #create m simulation studies
  for(i in 1:m){
    daten<-simulEU_SC(N=N,param = param,rateC=rateC,prob=prob,b=b)
    sims[[i]]<-daten
  }
  
  #return list of simulation studies
  return(sims)
}


#No censoring
#simulation of one study
simulEU_nocens<-function(N,param,prob, b){
  #N is the number of participants
  #param are shape and scale for the EU distribution 
  #prob is a vector of probabilities for the Bernoulli distributions, it determines 
  #the group distribution treatment:control
  #b is a vector with the real values of the parameters we wish to estimate later
  
  a2<-as.numeric(param[1]) #shape
  a1<-as.numeric(param[2]) #scale
  
  #if we do not state otherwise we assume that patients are assignment to the 
  #groups with equal probability
  if(missing(prob)){prob=0.5}
  
  #group assignment
  expressions<-rbern(N, prob=prob)
  
  #list to store sub data.frames, one for each of the two groups
  datalist<-list()
  
  #vector of possible group assignments
  expr_vector<-c(1,0)
  for(i in 1:2){#for every possible group
    #indicating the group     
    pers_expr<-expr_vector[i]
    
    #only regard the subset assigned to the current group
    single_expr<-subset(expressions,expressions==pers_expr)
    
    #number of patient in this group
    Nn<-length(single_expr)
    
    #calculating the scale parameter according to the group
    b1<-exp(-b*pers_expr)^(-1/a2)*a1
    
    #draw survival times 
    survival<-rEU(n=Nn,alpha=a2,a=0,b = b1)
    
    #add to the sub data set
    daten_single_expr<-data.frame(expressions=single_expr,survival=survival)
    
    #add to list of data frames
    datalist[[i]]<-daten_single_expr
  }
  
  #combining the data sets
  sim_all<-datalist[[1]]
  sim_all<-rbind(sim_all,datalist[[2]])
  
  
  #calculate time and status from survival 
  time <- sim_all$survival
  status <- rep(1,times=N) #all events since there is no censoring
  
  #extracting the important information
  sim<-data.frame(
    id=1:N,         #simulated patient ID
    time=time,      #event time
    status=status,  #censoring stats (allways event)
    expr=sim_all[1] #group assignment
  )
  
  #data set with only the information we need
  return(sim)
}

#multiple simulations
Sim_EUnocens<-function(m,N,param,prob,b){
  #m number of simulations
  #N,lambda,gamma,prob,b for the single simulation
  
  #empty list for the simulated data sets
  sims<-list()
  
  #create m simulation studies
  for(i in 1:m){
    daten<-simulEU_nocens(N=N,param = param,prob=prob,b=b)
    sims[[i]]<-daten
  }
  
  #return list of simulation studies
  return(sims)
}



#PH case
#simulation of one study
simulPHWeib_SC<-function(N,param, rateC,prob, b){
  #N is the number of participants
  #param are shape and scale for the Weibull
  #rateC are the parameters used for censoring
  #prob is a vector of probabilities for the Bernoulli distributions 
  #that determines group assignments
  #b is a vector with the real values of the parameters we wish to estimate later
  
  a1<-as.numeric(param[1]) #shape
  a2<-as.numeric(param[2]) #scale
  
  #if we do not state otherwise we assume that patients are assigned to the 
  #groups with equal probability
  if(missing(prob)){prob=0.5}
  
  #draw group assignment
  expressions<-rbern(N, prob=prob)
  
  #to store created sub data.frames for each of the two groups
  datalist<-list()
  
  #vector with both possible groups
  expr_vector<-c(1,0)
  for(i in 1:2){#for every possible expressions
    #group we regard now
    pers_expr<-expr_vector[i]
    
    #subset assigned to this group
    single_expr<-subset(expressions,expressions==pers_expr)
    
    #number of patients assigned to this group
    Nn<-length(single_expr)
    
    #draw survival times according to group assignment
    survival<-rweibull(n=Nn,shape=a1,scale = a2*exp(b*pers_expr)^(1/a1))
    
    #add to the subdataset
    daten_single_expr<-data.frame(expressions=single_expr,survival=survival)
    
    #add to list of data frames
    datalist[[i]]<-daten_single_expr
  }
  
  #combine subdatasets
  sim_all<-datalist[[1]]
  sim_all<-rbind(sim_all,datalist[[2]])
  
  #censoring
  #extracting parameters of UD
  minC<-rateC[1] #always 0
  maxC<-rateC[2]
  
  #draw censoring times
  cens<-runif(n=N,min=minC,max=maxC) 
  
  #add to data frame
  sim_all$cens<-cens
  
  #combine survival times and censoring times to time and status
  time <- pmin(sim_all$survival, sim_all$cens)            
  status <- as.numeric(sim_all$survival <= sim_all$cens) 
  
  #extracting the important information
  sim<-data.frame(
    id=1:N,         #simulated patient ID
    time=time,      #event/censoring time 
    status=status,  #censoring status
    expr=sim_all[1] #group assignment
  )
  
  #data set with only the information we need
  return(sim)
}

#to test out correct parameters for the UD we use the following function
simulPHWeib_SC_censrate<-function(N,param, rateC,prob, b){
  #same as above until censoring times are drawn
  
  a1<-as.numeric(param[1]) 
  a2<-as.numeric(param[2]) 
  
  if(missing(prob)){prob=0.5}
  
  expressions<-rbern(N, prob=prob)
  
  datalist<-list()
  
  expr_vector<-c(1,0)
  for(i in 1:2){
    pers_expr<-expr_vector[i]
    single_expr<-subset(expressions,expressions==pers_expr)
    
    Nn<-length(single_expr)
    
    survival<-rweibull(n=Nn,shape=a1,scale = a2*exp(b*pers_expr)^(1/a1))
    
    daten_single_expr<-data.frame(expressions=single_expr,survival=survival)
    
    datalist[[i]]<-daten_single_expr
  }
  
  sim_all<-datalist[[1]]
  sim_all<-rbind(sim_all,datalist[[2]])
  
  minC<-rateC[1]
  maxC<-rateC[2]
  cens<-runif(n=N,min=minC,max=maxC) 
  
  sim_all$cens<-cens
  
  #now we can calculate the percentage of censoring
  Cens<-as.numeric(sim_all$cens<sim_all$survival)
  CensRate<-sum(Cens)/N
  
  #return percentage
  return(CensRate)
}

#multiple simulations
Sim_PHSC<-function(m,N,param,rateC,prob,b){
  #m number of simulation studies
  #N,param, rateC,prob,b for the single simulation
  
  #simulation
  sims<-list()#empty list for the simulated data sets
  for(i in 1:m){
    #simulate m data sets
    daten<-simulPHWeib_SC(N=N,param = param,rateC=rateC,prob=prob,b=b)
    sims[[i]]<-daten #store
  }
  
  #return list of simulated data sets
  return(sims)
}


#No censoring
#simulation of one study
simulPHWeib_nocens<-function(N,param,prob, b){
  #N is the number of participants
  #param are shape and scale for the Weibull
  #prob is a vector of probabilities for the Bernoulli distributions 
  #that determines group assignments
  #b is a vector with the real values of the parameters we wish to estimate later
  
  a1<-as.numeric(param[1]) #shape
  a2<-as.numeric(param[2]) #scale
  
  #if we do not state otherwise we assume that patients are assigned to the 
  #groups with equal probability
  if(missing(prob)){prob=0.5}
  
  #draw group assignment
  expressions<-rbern(N, prob=prob)
  
  #to store sub data.frames for each of the two groups
  datalist<-list()
  
  #vector with both possible groups
  expr_vector<-c(1,0)
  for(i in 1:2){#for every possible expressions
    #group we regard now
    pers_expr<-expr_vector[i]
    
    #subset assigned to this group
    single_expr<-subset(expressions,expressions==pers_expr)
    
    #number assigned to this group
    Nn<-length(single_expr)
    
    #draw survival times according to group assignment
    survival<-rweibull(n=Nn,shape=a1,scale = a2*exp(b*pers_expr)^(1/a1))
    
    #add to the subdataset
    daten_single_expr<-data.frame(expressions=single_expr,survival=survival)
    
    #add to list of data frames
    datalist[[i]]<-daten_single_expr
  }
  
  #combine subdatasets
  sim_all<-datalist[[1]]
  sim_all<-rbind(sim_all,datalist[[2]])
  
  #combine survival times and censoring times to time and status
  time <- sim_all$survival
  status <- rep(1,times=N) #all events
  
  #extracting the important information
  sim<-data.frame(
    id=1:N,         #simulated patient ID
    time=time,      #event time
    status=status,  #censoring status (allways event)
    expr=sim_all[1] #group assignment
  )
  
  #data set with only the information we need
  return(sim)
}

#multiple simulations
Sim_PHnocens<-function(m,N,param,prob,b){
  #m number of simulation studies
  #N,lambda,gamma,prob,b for the single simulation
  
  #simulation
  sims<-list()#empty list for the simulated data sets
  for(i in 1:m){
    #simulate m data sets
    daten<-simulPHWeib_nocens(N=N,param = param,prob=prob,b=b)
    sims[[i]]<-daten
  }
  
  #return list of simulated data sets
  return(sims)
}



#PO case
#simulation of one study
simulPOLL_SC<-function(N,param, rateC,prob, b){
  #N is the number of participants
  #param are shape and scale for the LL
  #rateC are the parameters used for censoring
  #prob is a vector of probabilities for the Bernoulli distributions that determines group assignments
  #b is a vector with the real values of the parameters we wish to estimate later
  
  a1<-as.numeric(param[1]) #shape
  a2<-as.numeric(param[2]) #scale
  
  #if we do not state otherwise we assume that patients are assigned to the 
  #groups with equal probability
  if(missing(prob)){prob=0.5}
  
  #draw group assignment
  expressions<-rbern(N, prob=prob)
  
  #to store sub data.frames for each of the two groups
  datalist<-list()
  
  #vector with both possible groups
  expr_vector<-c(1,0)
  for(i in 1:2){#for every possible expressions
    #group we regard now
    pers_expr<-expr_vector[i]
    
    #subset assigned to this group
    single_expr<-subset(expressions,expressions==pers_expr)
    
    #number assigned to this group
    Nn<-length(single_expr)
    
    #draw survival times according to group assignment
    survival<-rllogis(n=Nn,shape=a1,scale = a2*exp(b*pers_expr)^(1/a1))
    
    #add to the subdataset
    daten_single_expr<-data.frame(expressions=single_expr,survival=survival)
    
    #add to list of data frames
    datalist[[i]]<-daten_single_expr
  }
  
  #combine subdatasets
  sim_all<-datalist[[1]]
  sim_all<-rbind(sim_all,datalist[[2]])
  
  #censoring
  #extracting parameters of UD
  minC<-rateC[1]     #always 0
  maxC<-rateC[2]
  
  #draw censoring times
  cens<-runif(n=N,min=minC,max=maxC) 
  
  #add to data frame
  sim_all$cens<-cens
  
  #combine survival times and censoring times to time and status
  time <- pmin(sim_all$survival, sim_all$cens)
  status <- as.numeric(sim_all$survival <= sim_all$cens) 
  
  #extracting the important information
  sim<-data.frame(
    id=1:N,         #simulated patient ID
    time=time,      #event or censoring time
    status=status,  #censoring status
    expr=sim_all[1] #group assignment
  )
  
  #data set with only the information we need
  return(sim)
}

#to test out correct parameters for the UD we use the following function
simulPOLL_SC_censrate<-function(N,param, rateC,prob, b){
  #same as above until censoring times are drawn
  
  a1<-as.numeric(param[1]) 
  a2<-as.numeric(param[2]) 
  
  if(missing(prob)){prob=0.5}
  
  expressions<-rbern(N, prob=prob)
  
  datalist<-list()
  
  expr_vector<-c(1,0)
  for(i in 1:2){
    pers_expr<-expr_vector[i]
    single_expr<-subset(expressions,expressions==pers_expr)
    
    Nn<-length(single_expr)
    
    survival<-rllogis(n=Nn,shape=a1,scale = a2*exp(b*pers_expr)^(1/a1))
    
    daten_single_expr<-data.frame(expressions=single_expr,survival=survival)
    
    datalist[[i]]<-daten_single_expr
  }
  
  sim_all<-datalist[[1]]
  sim_all<-rbind(sim_all,datalist[[2]])
  
  minC<-rateC[1]
  maxC<-rateC[2]
  cens<-runif(n=N,min=minC,max=maxC) 
  
  sim_all$cens<-cens
  
  #now we can calculate the percentage of censoring
  Cens<-as.numeric(sim_all$cens<sim_all$survival)
  CensRate<-sum(Cens)/N
  
  #retunr percentage
  return(CensRate)
}

#multiple simulations
Sim_POSC<-function(m,N,param,rateC,prob,b){
  #m number of simulation studies
  #N,param, rateC, prob, b for the single simulation
  
  #simulation
  sims<-list()#empty list for the simulated data sets
  for(i in 1:m){
    #simulate m data sets
    daten<-simulPOLL_SC(N=N,param = param,rateC=rateC,prob=prob,b=b)
    sims[[i]]<-daten
  }
  
  #return list of simulated data sets
  return(sims)
}


#No censoring
#simulation of one study
simulPOLL_nocens<-function(N,param, prob, b){
  #N is the number of participants
  #param are shape and scale for the LL
  #prob is a vector of probabilities for the Bernoulli distributions that determines group assignments
  #b is a vector with the real values of the parameters we wish to estimate later
  
  a1<-as.numeric(param[1]) #shape
  a2<-as.numeric(param[2]) #scale
  
  #if we do not state otherwise we assume that patients are assigned to the 
  #groups with equal probability
  if(missing(prob)){prob=0.5}
  
  #draw group assignment
  expressions<-rbern(N, prob=prob)
  
  #create of sub data.frames for each of the two groups
  datalist<-list()
  
  #vector with both possible groups
  expr_vector<-c(1,0)
  for(i in 1:2){#for every possible expressions
    #group we regard now
    pers_expr<-expr_vector[i]
    
    #subset assigned to this group
    single_expr<-subset(expressions,expressions==pers_expr)
    
    #number assigned to this group
    Nn<-length(single_expr)
    
    #draw survival times according to group assignment
    survival<-rllogis(n=Nn,shape=a1,scale = a2*exp(b*pers_expr)^(1/a1))
    
    #add to the subdataset
    daten_single_expr<-data.frame(expressions=single_expr,survival=survival)
    
    #add to list of data frames
    datalist[[i]]<-daten_single_expr
  }
  
  #combine subdatasets
  sim_all<-datalist[[1]]
  sim_all<-rbind(sim_all,datalist[[2]])
  
  #combine survival times and censoring times to time and status
  time <- sim_all$survival
  status <- rep(1,times=N) #all events
  
  #extracting the important information
  sim<-data.frame(
    id=1:N,         #simulated patient ID
    time=time,      #event time
    status=status,  #censoring status (all events)
    expr=sim_all[1] #group assignment
  )
  
  #data set with only the information we need
  return(sim)
}

#multiple simulations
Sim_POnocens<-function(m,N,param,prob,b){
  #m number of simulation studies
  #N,param, prob, b for the single simulation
  
  #simulation
  sims<-list()#empty list for the simulated data sets
  for(i in 1:m){
    #simulate m data sets
    daten<-simulPOLL_nocens(N=N,param = param,prob=prob,b=b)
    sims[[i]]<-daten
  }
  
  #return list of simulated data sets
  return(sims)
}




#Next we will prepare the analysis of the simulated data 

#NPPR model
#Estimation
Sim_est<-function(sims){
  #sims list of simulated studies
  
  #number of sims
  m<-length(sims)
  
  #estimation
  #vector to store results
  b_est<-numeric(m) 
  for(i in 1:m){
    #estimate for every simulation study
    est<-B_easy_fast(daten=sims[[i]],col=c(col_time=2,col_cens=3),col_target=4)
    b_est[i]<-as.numeric(est)   #store results
  }
  
  #return vector of all estimated beta
  return(b_est)
}

#Coverage 
Sim_CI_fast<-function(sim,b,L,alpha){
  #sim list of simulated studies
  #b,L, alpha for CI_fast
  
  #number of simulations
  m<-length(sim)
  
  #list to store CIs
  lim<-list()
  for(i in 1:m){
    #Estimate CI for every simulated study
    limits<-CI_fast(sim[[i]],col=c(col_time=2,col_cens=3),col_target=4,L=L,alpha=alpha)
    lim[[length(lim)+1]]<-limits #store
  }
  
  #coverage
  #vector to store information about coverage
  cov<-numeric(m)
  for(i in 1:m){
    #extract case
    limit<-lim[[i]]
    
    #extract limits
    lim_low<-as.numeric(limit[1])
    lim_up<-as.numeric(limit[2])
    
    #Check for every CI if b is covered
    cov[i]<-ifelse(b>=lim_low & b<=lim_up,1,0)
  }
  
  #percentage of coverage
  coverage<-(sum(cov)/m)*100
  
  #return percentage
  return(coverage)
}

#Output
Sim_out_log<-function(b,b_est){
  #b true value
  #b_est vector of estimated beta
  
  #mean
  Mean<-mean(b_est)
  
  #Bias
  Bias<-(mean(b_est)-b)
  
  #MSE
  MSE<-(var(b_est)+Bias^2)
  
  #return
  out<-data.frame(
    variable=1,
    Mean=Mean,
    Bias=Bias,
    MSE=MSE
  )
  
  #return results
  return(out)
}



#PPR model
#constructing the likelihood for the two group situation
#density function
d_EU<-function(t,alpha,theta){
  #alpha shape parameter
  #theta scale parameter
  #t time
  
  if(alpha<=0){
    return(NA)        #per definition, alpha is positive, otherwise the CDF would not
                      #be monotonically increasing
                      #to avoid problems using optim in the two group situation, 
                      #we manually exclude these cases 
  } else{
    p<-alpha*((theta)^alpha)*t^(alpha-1)    #parametrization
    return(p)}
}

#cdf
p_EU<-function(t,alpha,theta){
  #alpha shape parameter
  #theta scale parameter
  #t time
  
  
  if(alpha<=0){
    return(NA)          #see above
  } else{
    p<-(theta*t)^alpha  #parametrization
    return(p)
  }
}

#likelihood
like_EU_RCT<-function(param, daten, col_target, col){
  #param shape and scale parameters
  #daten data set
  #col=c(col_time,col_cens)
  #col_time column with survival times
  #col_cens column with information about censoring as above
  
  #shape parameter for both groups
  alpha<-as.numeric(param[1])
  
  #extract column
  x1<-col_target
  x2<-as.numeric(col["col_time"])
  x3<-as.numeric(col["col_cens"])
  
  #stratify data
  datenstrat<-list(subset(daten,daten[,x1]==0),subset(daten,daten[,x1]==1))
  
  #vector to store interim results
  likevector<-c(NA,NA,NA,NA)
  for(i in 1:2){                            #for both groups
    subdaten<-datenstrat[[i]]
    
    #extracting parameter according to group
    t1<-as.numeric(param[1+i])
    t2<-1/t1
    
    #split data set into censored and event
    daten_event<-subset(subdaten,subdaten[,x3]==1)
    daten_cens<-subset(subdaten,subdaten[,x3]==0)
    
    #time vectors
    time_event<-daten_event[,x2]
    time_cens<-daten_cens[,x2]
    
    #number of events, censored observations
    n1<-nrow(daten_event)
    n2<-nrow(daten_cens)
    
    #events
    #distribution function
    EU_event<-function(t){
      d_EU(t,alpha=alpha,theta=t2)
    }
    
    #eveluate at event time points
    eu_event<-lapply(time_event,EU_event)
    eu_event_corr<-as.numeric(eu_event)*100 #due to the large number of observations
                                            #the product would be to small for R
                                            #to work with, therefore we correct
                                            #this does not influence the estimated values
    c1<-prod(as.numeric(eu_event_corr))     #product part of likelihood for events
    
    #censored
    #CDF
    EU_cens<-function(t){
      (1-p_EU(t,alpha=alpha,theta=t2))
    }
    
    #evaluate at censored time points
    eu_cens<-lapply(time_cens,EU_cens)
    c2<-prod(as.numeric(eu_cens))         #product part of likelihood for censored times
    
    likevector[2*i-1]<-c1                 #store both interim results
    likevector[2*i]<-c2
  }
  
  #product regarding all cases
  like<-prod(likevector)
  
  #return likelihood
  return(like)
}

#(negative) log likelihood
loglike_EU_RCT<-function(param, daten, col_target, col){
  #param shape and scale parameters
  #daten data set
  #col=c(col_time,col_cens)
  #col_time column with survival times
  #col_cens column with information about censoring as above
  
  log(like_EU_RCT(param, daten, col_target, col)) #evaluate log likelihood
}

negloglike_EU_RCT<-function(param, daten, col_target, col){
  #param shape and scale parameters
  #daten data set
  #col=c(col_time,col_cens)
  #col_time column with survival times
  #col_cens column with information about censoring as above
  
  -log(like_EU_RCT(param, daten, col_target, col)) #evaluate negative log likelihood
}

#estimation
Sim_EU_ASV_Est<-function(sims,start){
  #sims list of simulated studies
  #start start parameters for optim
  
  #number of simulations
  m<-length(sims)
  
  #to store results
  EU<-list()
  
  for(i in 1:m){#for every simulation
    sim<-sims[[i]]#extract sim
    
    #estimate the parameters
    est<-optim(par=start,fn=negloglike_EU_RCT,daten=sim,col=c(col_time=2,col_cens=3),col_target=4,hessian=F)
    
    #estimate the Hessmatrix
    hess<-optimHess(par=start,fn=loglike_EU_RCT,daten=sim,col=c(col_time=2,col_cens=3),col_target=4)
    
    #store results
    EU[[2*i-1]]<-est
    EU[[2*i]]<-hess
  }
  
  #return list with estimated parameters and Hessmatrix
  return(EU)
}

#extracting estimations with numerical problems, before further analysis, 
#keeping the Hessmatrix
Sim_EU_nna_CI<-function(ests){
  #ests list with output from optim and Hessmatrix
  #number of simulations  
   m<-length(ests)/2
  
  #vector of estimated RR
  RR<-rep(NA,length=m)
  for(i in 1:m){
    est<-ests[[2*i-1]]
    alpha<-est$par[1]
    t0<-est$par[2]
    t1<-est$par[3]
    RR[i]<-(t0/t1)^alpha
  }
  
  #exclusion of those with numerical problems
  nna<-which(-log(RR)>-3&-log(RR)<3)
  nna_vector<-numeric(2*length(nna))
  
  for(i in 1:length(nna)){         
    nna_vector[2*i-1]<-2*nna[i]-1   #exclude which result
    nna_vector[2*i]<-2*nna[i]       #exclude which Hessmatrix 
  }
  
  #exclude identifier results
  ests_nna<-ests[nna_vector]
  
  #return only estimated optim output/Hessmatrix without numerical problems
  return(ests_nna)
}

#Coverage
Sim_EU_CI<-function(ests,b,df){
  #ests output of estimation including Hessematrix
  #b true value
  #df degree of freedom
  
  #number of simulation runs
  m<-length(ests)/2
  
  #gradient function
  grad_h<-function(param){
    alpha<-param[1]
    t0<-param[2]
    t1<-param[3]
    
    h<-c(NA,NA,NA)
    h[1]<-log(t0/t1)*(t0/t1)^alpha
    h[2]<-alpha*((1/t1)^alpha)*t0^(alpha-1)
    h[3]<-(t0^alpha)*(-alpha)*(t1^(-alpha-1))
    
    return(h)
  }
  
  #quantiles
  crit <- qt(1-0.05/2, df = df-1)
  
  #vectors to store information
  cov<-rep(NA,times=m)
  grad<-list()
  Fish<-list()
  s2<-rep(NA,times=m)
  
  for( i in 1:m ){
    #Hessian negative as we minimized the neg log likelihood
    est<-ests[[2*i-1]]
    grad_h_est<-grad_h(est$par)
    hess<--ests[[2*i]]
    
    #estimator of covariance matrix
    invFisher <- solve(hess)
    
    #RR
    alpha<-est$par[1]
    t0<-est$par[2]
    t1<-est$par[3]
    RR<-(t0/t1)^alpha
    
    #estimated variance
    s2_est<-grad_h_est%*%(invFisher%*%grad_h_est)
    if(s2_est<0){
      cov[i]<-NA
    }else{
      s_est<-sqrt(s2_est)
      
      lower<-RR-crit*s_est
      upper<-RR+crit*s_est
      
      s2[i]<-s2_est
      grad[[i]]<-grad_h_est
      Fish[[i]]<-invFisher
      
      cov[i]<-ifelse(exp(-b)>=lower & exp(-b)<=upper,1,0)
    }}
  
  #exclude NAs
  cov_nna<-cov[!is.na(cov)]
  
  #calculate percentage of coverage
  C<-sum(cov_nna)
  coverage<-(C*100)/length(cov_nna)
  
  #return percentage
  return(coverage)
}

#test how many estimations have numerical problems, these numbers are later collected 
#in the table "numerical robustness"
Sim_EU_nna_test_log<-function(ests){
  #ests output optim and Hessmatrix
  
  #number of simulations
  m<-length(ests)/2
  
  #RR
  RR<-rep(NA,length=m)
  for(i in 1:m){
    est<-ests[[2*i-1]]
    alpha<-est$par[1]
    t0<-est$par[2]
    t1<-est$par[3]
    RR[i]<-(t0/t1)^alpha
  }
  
  #-log(RR)
  nlogRR<--log(RR)
  
  #exclude those with numerical problems
  nlogRR_nna<-subset(nlogRR,nlogRR>-3&nlogRR<3)
  
  #number without numerical problems
  return(length(nlogRR_nna))
}

#excluding estimations with numerical problems, transforming ests to RR
Sim_EU_nna_log<-function(ests){
  #ests output optim and Hessmatrix
  
  #number of simulations
  m<-length(ests)/2
  
  #RR
  RR<-rep(NA,length=m)
  for(i in 1:m){
    est<-ests[[2*i-1]]
    alpha<-est$par[1]
    t0<-est$par[2]
    t1<-est$par[3]
    RR[i]<-(t0/t1)^alpha
  }
  
  #excluding estimations with numerical Problems
  RR_naa<-subset(RR,-log(RR)>-3&-log(RR)<3)
  
  #return only RRs, excluding those with numerical problems
  return(RR_naa)
}

#Output
Sim_EU_out_log<-function(ests,b){
  #ests estimated RRs
  #b true value
  
  #The values of interest are on the -log(RR) scale
  #and are transformed accordingly 
  
  #mean
  Mean<-mean(-log(ests))
  
  #Bias
  Bias<-(mean(-log(ests))-b)
  
  #MSE
  MSE<-(var(-log(ests))+Bias^2)
  
  #return
  out<-data.frame(
    variable=1,
    Mean=Mean,
    Bias=Bias,
    MSE=MSE
  )
  
  return(out)
}



#Cox's model
#estimation
Sim_PH_est<-function(sims){
  #sims list of simulations
  
  #number of simulations
  m<-length(sims)
  
  #estimate PH cox
  cox<-list()
  for(i in 1:m){
    sim<-sims[[i]]  #extract simulation
    #estimate
    cox[[i]]<-summary(coxph(Surv(sim$time,sim$status)~expressions,data=sim))
  }
  #return estimated information
  return(cox)
}

#check for simulations with high estimations indicating numerical problems
maxvalue.ph<-function(coxlist){
  #coxlist output from the estimation
  
  #vector to store information
  b_est<-numeric(length(coxlist))
  for(i in 1:length(coxlist)){
    #extract estimated information
    cox<-coxlist[[i]]
    #extract coefficients (-log(HR))
    b_est[i]<--cox$coefficients[1]
  }
  #maximal value
  out<-max(abs(b_est))
  #return maximal value
  return(out)
}

#later used to exclude high estimations
estcox<-function(cox){
  #cox results from estimation
  
  #vector to store results
  b_est<-numeric(length(cox))
  for(i in 1:length(cox)){
    c<-cox[[i]] #extract estimations
    b_est[i]<-abs(-c$coefficients[1]) #extract coefficients (-log(HR))
  }
  #return absolute -log(HR)
  return(b_est)
}

#Output
Sim_PH_out<-function(est,b){
  #est output from estimation
  #b true value
  
  #number of simulations
  m<-length(est)

  #estimates
  #vector to store estimates
  b_est<-numeric(m)
  for(i in 1:m){
    cox<-est[[i]] #extract estimation
    b_est[i]<--cox$coefficients[1] #extract coefficients (-log(HR))
  }
  
  
  #coverage
  #vector to store information
  cov<-numeric(m)
  for(i in 1:m){
    cox<-est[[i]] #extract estimation
    #extract CI
    lim_low<-cox$conf.int[3] 
    lim_up<-cox$conf.int[4]
    #check if true value is covered
    cov[i]<-ifelse(exp(-b)<=lim_up&exp(-b)>=lim_low,1,0)
  }
  
  #calculate percentage
  coverage<-(sum(cov)/m)*100
  
  
  #mean
  Mean<-mean(b_est)
  
  #Bias
  Bias<-(mean(b_est)-b)
  
  
  #MSE
  MSE<-(var(b_est)+Bias^2)
  
  
  #return results
  out<-data.frame(
    variable=1,
    Mean=Mean,
    Bias=Bias,
    MSE=MSE,
    Coverage=coverage
  )
  return(out)
}




#LL model
#likelihood for the two group situation
like_PO_RCT<-function(param, daten, col_target, col){
  #param shape and scale parameters
  #daten data set
  #col=c(col_times,col_cens)
  #col_time column with survival times
  #col_cens column with information about censoring as above
  
  #shape parameter for both groups
  alpha<-as.numeric(param[1])
  
  #extract column
  x1<-col_target
  x2<-as.numeric(col["col_time"])
  x3<-as.numeric(col["col_cens"])
  
  #stratify data
  datenstrat<-list(subset(daten,daten[,x1]==0),subset(daten,daten[,x1]==1))
  
  #vector to store interim results
  likevector<-c(NA,NA,NA,NA)
  for(i in 1:2){                            #for both groups
    subdaten<-datenstrat[[i]]
    
    #extracting parameter according to group
    t1<-as.numeric(param[1+i])
    
    #split data set into censored and event
    daten_event<-subset(subdaten,subdaten[,x3]==1)
    daten_cens<-subset(subdaten,subdaten[,x3]==0)
    
    #time vectors
    time_event<-daten_event[,x2]
    time_cens<-daten_cens[,x2]
    
    #number of events, censored observations
    n1<-nrow(daten_event)
    n2<-nrow(daten_cens)
    
    #events
    #distribution function
    PO_event<-function(t){
      dllogis(t,shape=alpha,scale=t1)
    }
    
    #evaluate at events
    po_event<-lapply(time_event,PO_event)
    po_event_corr<-as.numeric(po_event)*500 #correcting since otherwise the product
                                            #will get to small for R to evaluate
    c1<-prod(as.numeric(po_event_corr))     #product for event time points
    
    #censored
    #CDFs
    PO_cens<-function(t){
      (1- pllogis(t,shape=alpha,scale=t1))
    }
    
    #evaluate at censoring times
    po_cens<-lapply(time_cens,PO_cens)
    c2<-prod(as.numeric(po_cens))           #product for censoring times
    
    #store results
    likevector[2*i-1]<-c1
    likevector[2*i]<-c2
  }
  
  #product regarding all cases
  like<-prod(likevector)
  
  #return likelihood
  return(like)
}

#negative log likelihood
negloglike_PO_RCT<-function(param, daten, col_target, col){
  #param shape and scale parameters
  #daten data set
  #col=c(col_times,col_cens)
  #col_time column with survival times
  #col_cens column with information about censoring as above
  
  -log(like_PO_RCT(param, daten, col_target, col))  #evaluate negative log likelihood
}

#estimation
Sim_PO_Est<-function(sims,start){
  #sims list of simulated studies
  #start start parameters for optim
  
  #number of simulations
  m<-length(sims)
  
  #vector to store results
  PO<-list()
  
  for(i in 1:m){#for every simulation
    sim<-sims[[i]] #extract simulations
    
    #estimate the parameters
    est<-optim(par=start,fn=negloglike_PO_RCT,daten=sim,col=c(col_time=2,col_cens=3),col_target=4,hessian=F)
    
    #estimate the Hessmatrix
    hess<-optimHess(par=start,fn=negloglike_PO_RCT,daten=sim,col=c(col_time=2,col_cens=3),col_target=4)
    
    #store results
    PO[[2*i-1]]<-est
    PO[[2*i]]<-hess
  }
  
  #return list with estimated parameters and Hessmatrix
  return(PO)
}

#test how many estimations have numerical problems, these numbers are later collected 
#in the table "numerical robustness"
Sim_PO_nna_test_log<-function(ests){
  #ests output optim and Hessmatrix
  
  #number of simulations
  m<-length(ests)/2
  
  #OR
  ODD<-rep(NA,length=m)
  for(i in 1:m){
    est<-ests[[2*i-1]]
    alpha<-est$par[1]
    t0<-est$par[2]
    t1<-est$par[3]
    ODD[i]<-(t0/t1)^alpha
  }
  
  #-log(OR)
  nlogODD<--log(ODD)
  
  #exclude those with numerical problems
  nlogODD_nna<-subset(nlogODD,nlogODD>-3&nlogODD<3)
  
  #number without numerical problems
  return(length(nlogODD_nna))
}

#excluding estimations with numerical problems, transforming ests to -log(PO)
Sim_PO_nna_log<-function(ests){
  #ests output optim and Hessmatrix
  
  #number of sumulations
  m<-length(ests)/2
  
  #OR
  ODD<-rep(NA,length=m)
  for(i in 1:m){
    est<-ests[[2*i-1]]
    alpha<-est$par[1]
    t0<-est$par[2]
    t1<-est$par[3]
    ODD[i]<-(t0/t1)^alpha
  }
  
  #excluding estimations with numerical Problems
  ODD_naa<-subset(ODD,-log(ODD)>-3&-log(ODD)<3)
  
  #return only ORs, excluding those with numerical problems
  return(ODD_naa)
}

#extracting estimations with numerical problems, before further analysis, 
#keeping the Hessmatrix
Sim_PO_nna_CI<-function(ests){
  #ests list with output from optim and Hessmatrix
  
  #number of simulations
  m<-length(ests)/2
  
  #vector of estimated OR
  ODD<-rep(NA,length=m)
  for(i in 1:m){
    est<-ests[[2*i-1]]
    alpha<-est$par[1]
    t0<-est$par[2]
    t1<-est$par[3]
    ODD[i]<-(t0/t1)^alpha
  }
  
  #exclusion of those with numerical problems
  nna<-which(-log(ODD)>-3&-log(ODD)<3)
  nna_vector<-numeric(2*length(nna))
  
  for(i in 1:length(nna)){         
    nna_vector[2*i-1]<-2*nna[i]-1   #exclude which result
    nna_vector[2*i]<-2*nna[i]       #exclude which Hessmatrix 
  }
  
  #exclude those identifiert with numerical problems
  ests_nna<-ests[nna_vector]
  
  #return only estimated optim output/Hessmatrix without numerical problems
  return(ests_nna)
}

#Coverage
Sim_PO_CI<-function(ests,b,df){
  #ests output of estimation including Hessematrix
  #b true value
  #df degree of freedom
  
  #number of simulations
  m<-length(ests)/2
  
  #gradient function
  grad_h<-function(param){
    alpha<-param[1]
    t0<-param[2]
    t1<-param[3]
    
    h<-c(NA,NA,NA)
    h[1]<-log(t0/t1)*(t0/t1)^alpha
    h[2]<-alpha*((1/t1)^alpha)*t0^(alpha-1)
    h[3]<-(t0^alpha)*(-alpha)*(t1^(-alpha-1))
    
    return(h)
  }
  
  #quantils
  crit <- qt(1-0.05/2, df = df-1)
  
  #vectors to store information
  cov<-rep(NA,times=m)
  grad<-list()
  Fish<-list()
  s2<-rep(NA,times=m)
  for( i in 1:m ){
    #Hessian negative as we minimized the neg log likelihood
    est<-ests[[2*i-1]]
    grad_h_est<-grad_h(est$par)
    hess<-ests[[2*i]]
    
    #estimator of covariance matrix
    invFisher <- solve(hess)
    
    #OR
    alpha<-est$par[1]
    t0<-est$par[2]
    t1<-est$par[3]
    ODD<-(t0/t1)^alpha
    
    #estimated variance
    s2_est<-grad_h_est%*%(invFisher%*%grad_h_est)
    if(s2_est<0){
      cov[i]<-NA
    }else{
      s_est<-sqrt(s2_est)
      
      lower<-ODD-crit*s_est
      upper<-ODD+crit*s_est
      
      s2[i]<-s2_est
      grad[[i]]<-grad_h_est
      Fish[[i]]<-invFisher
      cov[i]<-ifelse(exp(-b)>=lower & exp(-b)<=upper,1,0)
    }}
  
  #exclude NAs
  cov_nna<-cov[!is.na(cov)]
  
  #calculate percentage of coverage
  C<-sum(cov_nna)
  coverage<-(C*100)/length(cov_nna)
  
  #retunr percentage
  return(coverage)
}

#Out
Sim_PO_out_log<-function(ests,b){
  #ests estimated ODD
  #b true value
  
  #The values of interest are on the -log(ODD) scale
  #and are transformed accordingly 
  
  #mean
  Mean<-mean(-log(ests))
  
  #Bias
  Bias<-(mean(-log(ests))-b)
  
  #MSE
  MSE<-(var(-log(ests))+Bias^2)
  
  #return results
  out<-data.frame(
    variable=1,
    Mean=Mean,
    Bias=Bias,
    MSE=MSE
  )
  return(out)
}


#####Save Worspace####
save.image(file = "./NPPR_base.RData")