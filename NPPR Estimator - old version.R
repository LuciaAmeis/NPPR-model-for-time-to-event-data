####Packages####
library("Rlab")
library("rlist")
library("survival")
library("haven")
library("sas7bdat")
source ("https://git.io/fjinW")
library("dplyr")
library("parallel")
library("foreach")
library("doParallel")
library("doSNOW")
library("numDeriv")
library("ggplot2")




####Parallelization####
n.cores<-parallel::detectCores()-1
my.cluster<-parallel::makeCluster(
  n.cores,
  type="PSOCK"
)




####Seed####
set.seed(60)




####Stepwise construction of the estimator####
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
    surv_prob<-1
  } else{ if(t>max(km_fit_time)){
    surv_prob<-0
  }else{
    for(i in 1:m_time){
      y1=km_fit_time[i]
      y2=km_fit_time[i+1]
      y3=km_fit_surv[i]
      ymax=max(km_fit_time)
      y4=km_fit_surv[m_time]
      if(i==m_time){
        if(t==ymax){
          surv_prob<-y4
        }} else{
          if(y1<=t& t<y2){
            surv_prob<-y3
          } 
          
        }}}}
  
  return(surv_prob)
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
    for(i in 1:m_time){
      y1=km_fit_time[i]
      y2=km_fit_time[i+1]
      y3=km_fit_se[i]
      ymax=max(km_fit_time)
      y4=km_fit_se[m_time]
      if(i==m_time){
        if(t==ymax){
          surv_se<-y4
        }} else{
          if(y1<=t& t<y2){
            surv_se<-y3
          } 
        }}}}
  
  return(surv_se)
}

#Estimated CDFs using the KM estimator
F_est<-function(t,daten,col){
  1-KM_surv(t,daten,col)
}

#Extracting the set \tilde{T}
T_jump<-function(events1,events0){
  #function of two event vectors, test an control
  
  t_min<-max(min(events1),min(events0))#lower bound
  t_max<-min(max(events1),max(events0))#upper bound
  
  events<-c(events1,events0)#all events
  
  #events between lower and upper bound
  events_bound<-subset(events,events>=t_min&events<=t_max)#all events in bounds
  
  #return
  return(events_bound)
}

#Estimator
Beta_est<-function(t_vector,Fn1,Fn0,W1,W0){
  #t_vector is the vector of event time points (later extracted with T_jump)
  #Fn1 and Fn0 are functions (KM-estimated CDFs)
  #W1 and W0 are the weights (KM-estimated SE)
  
  w<-numeric(length(t_vector))
  d<-numeric(length(t_vector))
  for(i in 1:length(t_vector)){
    t<-t_vector[i]                                    #for every time point
    Wi<-(1/Fn1(t)^2)*W1(t)^2+(1/Fn0(t)^2)*W0(t)^2     #calculate the weight at time point
    d[i]<-(1/Wi)*log(Fn1(t)/Fn0(t))                   #weighted log of estimated RR at time point

        w[i]<-1/Wi                                    #save the weights
  }
  W<-sum(w)                                           #Sum of weights                                   
  D<-sum(d)                                           #Sum of weighted log of estimated RR
  beta_est=-(((1/W)*D)) #take the weighted mean       #Combine 
  
  #return
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
  
  #return
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
  
  
  #estimator
  if(length(events_bound)>0){
    #Use the estimator defined above
    beta_est_matrix<-sapply(X=events_bound, FUN=Beta_est, simplify = TRUE, USE.NAMES = TRUE)
    B_est<-sum(beta_est_matrix[1,])
    W_est<-sum(beta_est_matrix[2,])
    beta_est<--(1/W_est)*B_est
  } else{beta_est=NA}
  
  return(beta_est)
}

#Confidence interval for beta
CI_fast<-function(daten,col,col_target,L,alpha){
  #daten is the data set
  #L is the number of samples we draw
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
  
  return(Q)
}

#Number needed to treat
NNT_fast<-function(daten,col,col_target,tp){
  #daten is the data frame
  #col is as above
  #col_target is the number corresponding to the column with the variable of 
               #interest
  #tp time
  
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
  events_bound<-subset(events,events>=t_min&events<=t_max)#all events in bounds
  
  #CDF and weight function
  #treatment
  #KM-estimator included in the package survival
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
  
  #estimator
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
    samp<-sample(1:nrow(daten),replace =T)
    #estimate for every sample
    x<-as.numeric(NNT_fast_sample(daten[samp,]))
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
  
  return(Q)
}




####Analysis of the DAPA-HF data set####
DapaHF <- read_sas("Please refer to repository")
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
NNT_fast(DapaHF,col=c(col_time=3,col_cens=1),col_target=2,10)

#CI for nnt at time point 10
CI_nnt_fast(DapaHF,col=c(col_time=3,col_cens=1),col_target=2,L=500,alpha=0.05,10)





