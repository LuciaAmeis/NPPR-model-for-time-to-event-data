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




####Preparation: Setting of the simulation studies####
#Weibull
#likelihood
like_Weib<-function(param, daten,col){
  #param parameters to be estimated
  #daten data set
  #col_time column with survival times
  #col_cens column with information about censoring as above
  
  a1<-as.numeric(param[1]) #shape
  a2<-as.numeric(param[2]) #scale
  
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
  Weib_event<-function(t){
    dweibull(t,shape=a1,scale=a2,log=F)
  }
  
  weib_event<-lapply(time_event,Weib_event)
  weib_event_corr<-as.numeric(weib_event)*50  #due to the large number of observations
                                              #the product would be to small for R
                                              #to work with, therefore we correct
                                              #this does not influence the estimated values
  c1<-prod(as.numeric(weib_event_corr))
  
  #censored in likelihood
  Weib_cens<-function(t){
    (1-pweibull(t,shape=a1,scale=a2,lower.tail=T,log.p=F))
  }
  
  weib_cens<-lapply(time_cens,Weib_cens)
  weib_cens_corr<-as.numeric(weib_cens)*1.45 #due to the number of observations
                                             #the product would be to small for R
                                             #to work with, therefore we correct
                                             #this does not influence the estimated values
  c2<-prod(as.numeric(weib_cens_corr))
  
  #products
  like<-c1*c2
  
  return(like)
}

#negative log likelihood
negloglike_Weib<-function(param, daten,col){
  -log(like_Weib(param,daten,col))
}

#Fitting a Weibull distribution to the control group of the data set
#starting values were determined by plotting  different Weibull distributions 
#against the KM-estimated CDF of the control group
ML_0<-optim(par=c(1.2,70),fn=negloglike_Weib,daten=DapaHF_0,col=c(col_time=3,col_cens=1),lower = c(0, 0))



#PPR
#likelihood
likeEU<-function(param, daten,col){
  #param shape and form parameters
  #alpha is the shape
  #the parametrization of the scale parameter in the used is slightly different than in the paper
  #theta=1/(t1-t0)=1/t1, whereby t0 is always 0 here. 
  #daten data set
  #col_time column with survival times
  #col_cens column with information about censoring as above
  
  alpha<-as.numeric(param[1])
  t1<-as.numeric(param[2])
  
  x1<-as.numeric(col["col_time"])
  x2<-as.numeric(col["col_cens"])
  
  #split data set into censored and not censord
  daten_event<-subset(daten,daten[,x2]==1)
  daten_cens<-subset(daten,daten[,x2]==0)
  
  #time vectors
  time_event<-daten_event[,x1]
  time_cens<-daten_cens[,x1]
  
  #number of events/censored observations
  n1<-nrow(daten_event)
  n2<-nrow(daten_cens)
  
  #event data
  EU_event<-function(t){
    dEU(t,alpha=alpha,a=0,b=t1)
  }
  
  eu_event<-lapply(time_event,EU_event)
  eu_event_corr<-as.numeric(eu_event)*50 #due to the number of observations
                                         #the product would be to small for R
                                         #to work with, therefore we correct
                                         #this does not influence the estimated values
  c1<-prod(as.numeric(eu_event_corr))
  
  #censored
  EU_cens<-function(t){
    (1-pEU(t,alpha=alpha,a=0,b=t1))
  }
  
  eu_cens<-lapply(time_cens,EU_cens)
  c2<-prod(as.numeric(eu_cens))
  
  #product
  like<-c1*c2
  
  return(like)
}

#negative log likelihood
negloglikeEU<-function(param, daten,col){
  -log(likeEU(param, daten,col))
}

#Fitting a EU distribution to the control group of the data set
#starting values were determined by plotting  different EU distributions 
#against the KM-estimated CDF of the control group
ML_EU<-optim(par=c(0.8,130),negloglikeEU,daten=DapaHF_0,col=c(col_time=3,col_cens=1),method="Nelder-Mead")




####Preparation: Simulation####
#EU
#simulation of one study
simulEU_SC<-function(N,param, rateC,prob, b){
  #N is the number of participants
  #param are shape and scale for the EU distribution 
  #rateC parameters for censoring
  #prob is a vector of probabilities for the Bernoulli distributions, it determines 
         #the group distribution treatment:control
  #b is a vector with the real values of the parameters we wish to estimate later
  
  a2<-as.numeric(param[1]) #alpha
  a1<-as.numeric(param[2]) #theta
  
  #if we do not state otherwise we assume that patients are assignment to the 
  #groups with equal probability
  if(missing(prob)){prob=0.5}
  
  #group assignment
  expressions<-rbern(N, prob=prob)
  
  #list sub data.frames, one for each of the two groups
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
  minC<-rateC[1]
  maxC<-rateC[2]
  
  #draw censoring times
  cens<-runif(n=N,min=minC,max=maxC) 
  
  #add to data frame
  sim_all$cens<-cens
  
  #calculate time and status from survival and censoring times
  time <- pmin(sim_all$survival, sim_all$cens)
  status <- as.numeric(sim_all$survival <= sim_all$cens) #uncensored (experiences an event)
  
  #extracting the important information
  sim<-data.frame(
    id=1:N,
    time=time,
    status=status,
    expr=sim_all[1]
  )
  
  #data set with only the information we need
  return(sim)
}

#to test out correct parameters for the UD we use the following function
simulEU_SC_censrate<-function(N,param, rateC,prob, b){
  #parameters as above
  #proceed as above until censoring times are drawn
  
  a2<-as.numeric(param[1]) #alpha
  a1<-as.numeric(param[2]) #theta
  
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
  
  return(CensRate)
}

#multiple simulations
Sim_EUSC<-function(m,N,param,rateC,prob,b){
  #m number of simulations
  #N,lambda,gamma,rateC,prob,b for the single simulation
  
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



#Weibull
#simulation of one study
simulPHWeib_SC<-function(N,param, rateC,prob, b){#N is the number of participants
  #param are shape and scale for the Weibull
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
    survival<-rweibull(n=Nn,shape=a1,scale = a2*exp(b*pers_expr))
    
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
  minC<-rateC[1]
  maxC<-rateC[2]
  
  #draw censoring times
  cens<-runif(n=N,min=minC,max=maxC) 
  
  #add to data frame
  sim_all$cens<-cens
  
  #combine survival times and censoring times to time and status
  time <- pmin(sim_all$survival, sim_all$cens)
  status <- as.numeric(sim_all$survival <= sim_all$cens) #uncensored (experiences an event)
  
  #extracting the important information
  sim<-data.frame(
    id=1:N,
    time=time,
    status=status,
    expr=sim_all[1]
  )
  
  #data set with only the information we need
  return(sim)
}

#to test out correct parameters for the UD we use the following function
simulPHWeib_SC_censrate<-function(N,param, rateC,prob, b){
  #same as above until censoring times are drawn
  
  a1<-as.numeric(param[1]) #shape
  a2<-as.numeric(param[2]) #scale
  
  if(missing(prob)){prob=0.5}

  expressions<-rbern(N, prob=prob)
  
  datalist<-list()
  
  expr_vector<-c(1,0)
  for(i in 1:2){
   pers_expr<-expr_vector[i]
   single_expr<-subset(expressions,expressions==pers_expr)
   
   Nn<-length(single_expr)
    
   survival<-rweibull(n=Nn,shape=a1,scale = a2*exp(b*pers_expr))
  
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
  
  return(CensRate)
}

#multiple simulations
Sim_PHSC<-function(m,N,param,rateC,prob,b){
  #m number of simulation studies
  #N,lambda,gamma,rateC,prob,b for the single simulation
  
  #simulation
  sims<-list()#empty list for the simulated data sets
  for(i in 1:m){
    #simulate m data sets
    daten<-simulPHWeib_SC(N=N,param = param,rateC=rateC,prob=prob,b=b)
    sims[[i]]<-daten
  }
  
  #return list of simulated data sets
  return(sims)
}




####Preparation: Analysis####
#NPPR model
#Estimation
Sim_est<-function(sims){
  #sims list of simulated studies
  
  m<-length(sims)
  
  #estimation
  b_est<-numeric(m)
  for(i in 1:m){
    #estimate for every simulation study
    est<-B_easy_fast(daten=sims[[i]],col=c(col_time=2,col_cens=3),col_target=4)
    b_est[i]<-as.numeric(est)
  }
  
  #return vector of all estimated beta
  return(b_est)
}

#Coverage 
Sim_CI_fast<-function(sim,b,L,alpha){
  #sim list of simulated studies
  #b,L, alpha for CI_fast
  
  m<-length(sim)
  
  #CI
  lim<-list()
  for(i in 1:m){
    #Estimate CI for every simulated study
    limits<-CI_fast(sim[[i]],col=c(col_time=2,col_cens=3),col_target=4,L=L,alpha=alpha)
    lim[[length(lim)+1]]<-limits
  }
  
  #coverage
  cov<-numeric(m)
  for(i in 1:m){
    limit<-lim[[i]]
    
    #extract limits
    lim_low<-as.numeric(limit[1])
    lim_up<-as.numeric(limit[2])
    
    #Check for every CI if b is covered
    cov[i]<-ifelse(b>=lim_low & b<=lim_up,1,0)
  }
  
  #percentage of coverage
  coverage<-(sum(cov)/m)*100
  
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
    EU_event<-function(t){
      d_EU(t,alpha=alpha,theta=t2)
    }
    
    eu_event<-lapply(time_event,EU_event)
    eu_event_corr<-as.numeric(eu_event)*100
    c1<-prod(as.numeric(eu_event_corr))
    
    #censored
    EU_cens<-function(t){
      (1-p_EU(t,alpha=alpha,theta=t2))
    }
    
    eu_cens<-lapply(time_cens,EU_cens)
    eu_cens_corr<-as.numeric(eu_cens)
    c2<-prod(as.numeric(eu_cens_corr))
    
    likevector[2*i-1]<-c1
    likevector[2*i]<-c2
  }
  
  #product regarding all cases
  like<-prod(likevector)
  
  return(like)
}

#negative log likelihood
negloglike_EU_RCT<-function(param, daten, col_target, col){
  -log(like_EU_RCT(param, daten, col_target, col))
}

#estimation
Sim_EU_ASV_Est<-function(sims,start){
  #sims list of simulated studies
  #start start parameters for optim
  
  m<-length(sims)
  EU<-list()
  
  for(i in 1:m){#for every simulation
    sim<-sims[[i]]
    
    #estimate the parameters
    est<-optim(par=start,fn=negloglike_EU_RCT,daten=sim,col=c(col_time=2,col_cens=3),col_target=4,hessian=F)
    
    #estimate the Hessmatrix
    hess<-optimHess(par=start,fn=loglike_EU_RCT,daten=sim,col=c(col_time=2,col_cens=3),col_target=4)
    EU[[2*i-1]]<-est
    EU[[2*i]]<-hess
  }
  
  #return list with estimated parameters and Hessmatrix
  return(EU)
}


#extracting estimations with numerical problems, before further analysis, keeping the Hessmatrix
Sim_EU_nna_CI<-function(ests){
  #ests list with output from optim and Hessmatrix
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
    nna_vector[2*i-1]<-2*nna[i]-1   #exclude result
    nna_vector[2*i]<-2*nna[i]       #exclude Hessmatrix 
  }
  
  ests_nna<-ests[nna_vector]
  
  #return only estimated optim output/Hessmatrix without numerical problems
  return(ests_nna)
}

#Coverage
Sim_EU_CI<-function(ests,b,df){
  #ests output of estimation including Hessematrix
  #b true value
  #df degree of freedom
  
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
  
  return(coverage)
}


#pre-processing output
#test how many estimations have numerical problems, these numbers are collected 
#in the table "numerical robustness"
Sim_EU_nna_test_log<-function(ests){
  #ests output optim and Hessmatrix
  
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



#Putting results together
#CI
results_CI<-function(reslist){
  #reslists list of outputs created with CIOut functions above
  
  m<-length(reslist)
  
  #the outputs of the CI functions are combined to one vector
  Coverage<-rep(NA,m)
  for(i in 1:m){
    Coverage[i]<-reslist[[i]]
  }
  
  result<-data.frame(
    "Effect"=c(rep(0,9),rep(0.5,9),rep(0.25,9),rep(-0.25,9),rep(-0.5,9)),
    "Censoring"=c(rep(30,3),rep(50,3),rep(70,3),rep(30,3),rep(50,3),rep(70,3),rep(30,3),rep(50,3),rep(70,3),rep(30,3),rep(50,3),rep(70,3),rep(30,3),rep(50,3),rep(70,3)),
    "Paricipants"=c(rep(c(500,100,50),3),rep(c(500,100,50),3),rep(c(500,100,50),3),rep(c(500,100,50),3),rep(c(500,100,50),3)),
    "Coverage"=Coverage
  )
  
  return(result)
}

#Mean,Bias,MSE
results_log<-function(reslist){
  #reslist output from out functions above
  
  m<-length(reslist)
  
  #combine parameters to one vector each
  Mean<-rep(NA,m)
  Bias<-rep(NA,m)
  MSE<-rep(NA,m)
  for(i in 1:m){
    res<-reslist[[i]]
    Mean[i]<-res$Mean
    Bias[i]<-res$Bias
    MSE[i]<-res$MSE
  }
  
  result<-data.frame(
    "Effect"=c(rep(0,9),rep(0.5,9),rep(0.25,9),rep(-0.25,9),rep(-0.5,9)),
    "Censoring"=c(rep(30,3),rep(50,3),rep(70,3),rep(30,3),rep(50,3),rep(70,3),rep(30,3),rep(50,3),rep(70,3),rep(30,3),rep(50,3),rep(70,3),rep(30,3),rep(50,3),rep(70,3)),
    "Paricipants"=c(rep(c(500,100,50),3),rep(c(500,100,50),3),rep(c(500,100,50),3),rep(c(500,100,50),3),rep(c(500,100,50),3)),
    "Mean"=Mean,
    "Bias"=Bias,
    "MSE"=MSE
  )
  
  return(result)
}




####Simulation####
#PR
m1<-1000
#b=0                                                                                                 #true value
b1<-0
#30                                                                                                  #censoring rate
rateC_EU_DapaHF_0_censSC30<-c(0,175)                                                                 #parameters for censoring
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_0_censSC30,prob=0.5,b=b1)          #test censoring rate
#simulations
Sim_EU_500_DapaHF_0_censSC30<-Sim_EUSC(m=m1,N=500,param=ML_EU$par,rateC =rateC_EU_DapaHF_0_censSC30,prob=0.5,b=b1)
Sim_EU_100_DapaHF_0_censSC30<-Sim_EUSC(m=m1,N=100,param=ML_EU$par,rateC =rateC_EU_DapaHF_0_censSC30,prob=0.5,b=b1)
Sim_EU_50_DapaHF_0_censSC30<-Sim_EUSC(m=m1,N=50,param=ML_EU$par,rateC =rateC_EU_DapaHF_0_censSC30,prob=0.5,b=b1)

#50
rateC_EU_DapaHF_0_censSC50<-c(0,100)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_0_censSC50,prob=0.5,b=b1)
Sim_EU_500_DapaHF_0_censSC50<-Sim_EUSC(m=m1,N=500,param=ML_EU$par,rateC =rateC_EU_DapaHF_0_censSC50,prob=0.5,b=b1)
Sim_EU_100_DapaHF_0_censSC50<-Sim_EUSC(m=m1,N=100,param=ML_EU$par,rateC =rateC_EU_DapaHF_0_censSC50,prob=0.5,b=b1)
Sim_EU_50_DapaHF_0_censSC50<-Sim_EUSC(m=m1,N=50,param=ML_EU$par,rateC =rateC_EU_DapaHF_0_censSC50,prob=0.5,b=b1)

#70
rateC_EU_DapaHF_0_censSC70<-c(0,60)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_0_censSC70,prob=0.5,b=b1)
Sim_EU_500_DapaHF_0_censSC70<-Sim_EUSC(m=m1,N=500,param=ML_EU$par,rateC =rateC_EU_DapaHF_0_censSC70,prob=0.5,b=b1)
Sim_EU_100_DapaHF_0_censSC70<-Sim_EUSC(m=m1,N=100,param=ML_EU$par,rateC =rateC_EU_DapaHF_0_censSC70,prob=0.5,b=b1)
Sim_EU_50_DapaHF_0_censSC70<-Sim_EUSC(m=m1,N=50,param=ML_EU$par,rateC =rateC_EU_DapaHF_0_censSC70,prob=0.5,b=b1)


#b=0.5
b1<-0.5 
#30
rateC_EU_DapaHF_.5_censSC30<-c(0,250)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_.5_censSC30,prob=0.5,b=b1)
Sim_EU_500_DapaHF_.5_censSC30<-Sim_EUSC(m=m1,N=500,param=ML_EU$par,rateC =rateC_EU_DapaHF_.5_censSC30,prob=0.5,b=b1)
Sim_EU_100_DapaHF_.5_censSC30<-Sim_EUSC(m=m1,N=100,param=ML_EU$par,rateC =rateC_EU_DapaHF_.5_censSC30,prob=0.5,b=b1)
Sim_EU_50_DapaHF_.5_censSC30<-Sim_EUSC(m=m1,N=50,param=ML_EU$par,rateC =rateC_EU_DapaHF_.5_censSC30,prob=0.5,b=b1)

#50
rateC_EU_DapaHF_.5_censSC50<-c(0,140)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_.5_censSC50,prob=0.5,b=b1)
Sim_EU_500_DapaHF_.5_censSC50<-Sim_EUSC(m=m1,N=500,param=ML_EU$par,rateC =rateC_EU_DapaHF_.5_censSC50,prob=0.5,b=b1)
Sim_EU_100_DapaHF_.5_censSC50<-Sim_EUSC(m=m1,N=100,param=ML_EU$par,rateC =rateC_EU_DapaHF_.5_censSC50,prob=0.5,b=b1)
Sim_EU_50_DapaHF_.5_censSC50<-Sim_EUSC(m=m1,N=50,param=ML_EU$par,rateC =rateC_EU_DapaHF_.5_censSC50,prob=0.5,b=b1)

#70
rateC_EU_DapaHF_.5_censSC70<-c(0,75)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_.5_censSC70,prob=0.5,b=b1)
Sim_EU_500_DapaHF_.5_censSC70<-Sim_EUSC(m=m1,N=500,param=ML_EU$par,rateC =rateC_EU_DapaHF_.5_censSC70,prob=0.5,b=b1)
Sim_EU_100_DapaHF_.5_censSC70<-Sim_EUSC(m=m1,N=100,param=ML_EU$par,rateC =rateC_EU_DapaHF_.5_censSC70,prob=0.5,b=b1)
Sim_EU_50_DapaHF_.5_censSC70<-Sim_EUSC(m=m1,N=50,param=ML_EU$par,rateC =rateC_EU_DapaHF_.5_censSC70,prob=0.5,b=b1)


#b=-0.5
b1<--0.5
#30
rateC_EU_DapaHF_n.5_censSC30<-c(0,130)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_n.5_censSC30,prob=0.5,b=b1)
Sim_EU_500_DapaHF_n.5_censSC30<-Sim_EUSC(m=m1,N=500,param=ML_EU$par,rateC =rateC_EU_DapaHF_n.5_censSC30,prob=0.5,b=b1)
Sim_EU_100_DapaHF_n.5_censSC30<-Sim_EUSC(m=m1,N=100,param=ML_EU$par,rateC =rateC_EU_DapaHF_n.5_censSC30,prob=0.5,b=b1)
Sim_EU_50_DapaHF_n.5_censSC30<-Sim_EUSC(m=m1,N=50,param=ML_EU$par,rateC =rateC_EU_DapaHF_n.5_censSC30,prob=0.5,b=b1)

#50
rateC_EU_DapaHF_n.5_censSC50<-c(0,80)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_n.5_censSC50,prob=0.5,b=b1)
Sim_EU_500_DapaHF_n.5_censSC50<-Sim_EUSC(m=m1,N=500,param=ML_EU$par,rateC =rateC_EU_DapaHF_n.5_censSC50,prob=0.5,b=b1)
Sim_EU_100_DapaHF_n.5_censSC50<-Sim_EUSC(m=m1,N=100,param=ML_EU$par,rateC =rateC_EU_DapaHF_n.5_censSC50,prob=0.5,b=b1)
Sim_EU_50_DapaHF_n.5_censSC50<-Sim_EUSC(m=m1,N=50,param=ML_EU$par,rateC =rateC_EU_DapaHF_n.5_censSC50,prob=0.5,b=b1)

#70
rateC_EU_DapaHF_n.5_censSC70<-c(0,40)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_n.5_censSC70,prob=0.5,b=b1)
Sim_EU_500_DapaHF_n.5_censSC70<-Sim_EUSC(m=m1,N=500,param=ML_EU$par,rateC =rateC_EU_DapaHF_n.5_censSC70,prob=0.5,b=b1)
Sim_EU_100_DapaHF_n.5_censSC70<-Sim_EUSC(m=m1,N=100,param=ML_EU$par,rateC =rateC_EU_DapaHF_n.5_censSC70,prob=0.5,b=b1)
Sim_EU_50_DapaHF_n.5_censSC70<-Sim_EUSC(m=m1,N=50,param=ML_EU$par,rateC =rateC_EU_DapaHF_n.5_censSC70,prob=0.5,b=b1)


#b=0.25
b1<-0.25 
#30
rateC_EU_DapaHF_.25_censSC30<-c(0,200)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_.25_censSC30,prob=0.5,b=b1)
Sim_EU_500_DapaHF_.25_censSC30<-Sim_EUSC(m=m1,N=500,param=ML_EU$par,rateC =rateC_EU_DapaHF_.25_censSC30,prob=0.5,b=b1)
Sim_EU_100_DapaHF_.25_censSC30<-Sim_EUSC(m=m1,N=100,param=ML_EU$par,rateC =rateC_EU_DapaHF_.25_censSC30,prob=0.5,b=b1)
Sim_EU_50_DapaHF_.25_censSC30<-Sim_EUSC(m=m1,N=50,param=ML_EU$par,rateC =rateC_EU_DapaHF_.25_censSC30,prob=0.5,b=b1)

#50
rateC_EU_DapaHF_.25_censSC50<-c(0,125)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_.25_censSC50,prob=0.5,b=b1)
Sim_EU_500_DapaHF_.25_censSC50<-Sim_EUSC(m=m1,N=500,param=ML_EU$par,rateC =rateC_EU_DapaHF_.25_censSC50,prob=0.5,b=b1)
Sim_EU_100_DapaHF_.25_censSC50<-Sim_EUSC(m=m1,N=100,param=ML_EU$par,rateC =rateC_EU_DapaHF_.25_censSC50,prob=0.5,b=b1)
Sim_EU_50_DapaHF_.25_censSC50<-Sim_EUSC(m=m1,N=50,param=ML_EU$par,rateC =rateC_EU_DapaHF_.25_censSC50,prob=0.5,b=b1)

#70
rateC_EU_DapaHF_.25_censSC70<-c(0,70)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_.25_censSC70,prob=0.5,b=b1)
Sim_EU_500_DapaHF_.25_censSC70<-Sim_EUSC(m=m1,N=500,param=ML_EU$par,rateC =rateC_EU_DapaHF_.25_censSC70,prob=0.5,b=b1)
Sim_EU_100_DapaHF_.25_censSC70<-Sim_EUSC(m=m1,N=100,param=ML_EU$par,rateC =rateC_EU_DapaHF_.25_censSC70,prob=0.5,b=b1)
Sim_EU_50_DapaHF_.25_censSC70<-Sim_EUSC(m=m1,N=50,param=ML_EU$par,rateC =rateC_EU_DapaHF_.25_censSC70,prob=0.5,b=b1)


#b=-0.25
b1<--0.25
#30
rateC_EU_DapaHF_n.25_censSC30<-c(0,150)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_n.25_censSC30,prob=0.5,b=b1)
Sim_EU_500_DapaHF_n.25_censSC30<-Sim_EUSC(m=m1,N=500,param=ML_EU$par,rateC =rateC_EU_DapaHF_n.25_censSC30,prob=0.5,b=b1)
Sim_EU_100_DapaHF_n.25_censSC30<-Sim_EUSC(m=m1,N=100,param=ML_EU$par,rateC =rateC_EU_DapaHF_n.25_censSC30,prob=0.5,b=b1)
Sim_EU_50_DapaHF_n.25_censSC30<-Sim_EUSC(m=m1,N=50,param=ML_EU$par,rateC =rateC_EU_DapaHF_n.25_censSC30,prob=0.5,b=b1)

#50
rateC_EU_DapaHF_n.25_censSC50<-c(0,90)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_n.25_censSC50,prob=0.5,b=b1)
Sim_EU_500_DapaHF_n.25_censSC50<-Sim_EUSC(m=m1,N=500,param=ML_EU$par,rateC =rateC_EU_DapaHF_n.25_censSC50,prob=0.5,b=b1)
Sim_EU_100_DapaHF_n.25_censSC50<-Sim_EUSC(m=m1,N=100,param=ML_EU$par,rateC =rateC_EU_DapaHF_n.25_censSC50,prob=0.5,b=b1)
Sim_EU_50_DapaHF_n.25_censSC50<-Sim_EUSC(m=m1,N=50,param=ML_EU$par,rateC =rateC_EU_DapaHF_n.25_censSC50,prob=0.5,b=b1)

#70
rateC_EU_DapaHF_n.25_censSC70<-c(0,50)
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=rateC_EU_DapaHF_n.25_censSC70,prob=0.5,b=b1)
Sim_EU_500_DapaHF_n.25_censSC70<-Sim_EUSC(m=m1,N=500,param=ML_EU$par,rateC =rateC_EU_DapaHF_n.25_censSC70,prob=0.5,b=b1)
Sim_EU_100_DapaHF_n.25_censSC70<-Sim_EUSC(m=m1,N=100,param=ML_EU$par,rateC =rateC_EU_DapaHF_n.25_censSC70,prob=0.5,b=b1)
Sim_EU_50_DapaHF_n.25_censSC70<-Sim_EUSC(m=m1,N=50,param=ML_EU$par,rateC =rateC_EU_DapaHF_n.25_censSC70,prob=0.5,b=b1)



#PH
m1<-1000
#b=0
b1<-0
#30
rateC_PH_DapaHF_0_censSC30<-c(0,275)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_0_censSC30,prob=0.5,b=b1)
Sim_PH_500_DapaHF_0_censSC30<-Sim_PHSC(m=m1,N=500,param=ML_0$par,rateC =rateC_PH_DapaHF_0_censSC30,prob=0.5,b=b1)
Sim_PH_100_DapaHF_0_censSC30<-Sim_PHSC(m=m1,N=100,param=ML_0$par,rateC =rateC_PH_DapaHF_0_censSC30,prob=0.5,b=b1)
Sim_PH_50_DapaHF_0_censSC30<-Sim_PHSC(m=m1,N=50,param=ML_0$par,rateC =rateC_PH_DapaHF_0_censSC30,prob=0.5,b=b1)

#50
rateC_PH_DapaHF_0_censSC50<-c(0,130)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_0_censSC50,prob=0.5,b=b1)
Sim_PH_500_DapaHF_0_censSC50<-Sim_PHSC(m=m1,N=500,param=ML_0$par,rateC =rateC_PH_DapaHF_0_censSC50,prob=0.5,b=b1)
Sim_PH_100_DapaHF_0_censSC50<-Sim_PHSC(m=m1,N=100,param=ML_0$par,rateC =rateC_PH_DapaHF_0_censSC50,prob=0.5,b=b1)
Sim_PH_50_DapaHF_0_censSC50<-Sim_PHSC(m=m1,N=50,param=ML_0$par,rateC =rateC_PH_DapaHF_0_censSC50,prob=0.5,b=b1)

#70
rateC_PH_DapaHF_0_censSC70<-c(0,60)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_0_censSC70,prob=0.5,b=b1)
Sim_PH_500_DapaHF_0_censSC70<-Sim_PHSC(m=m1,N=500,param=ML_0$par,rateC =rateC_PH_DapaHF_0_censSC70,prob=0.5,b=b1)
Sim_PH_100_DapaHF_0_censSC70<-Sim_PHSC(m=m1,N=100,param=ML_0$par,rateC =rateC_PH_DapaHF_0_censSC70,prob=0.5,b=b1)
Sim_PH_50_DapaHF_0_censSC70<-Sim_PHSC(m=m1,N=50,param=ML_0$par,rateC =rateC_PH_DapaHF_0_censSC70,prob=0.5,b=b1)


#b=0.5
b1<-0.5 
#30
rateC_PH_DapaHF_.5_censSC30<-c(0,350)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_.5_censSC30,prob=0.5,b=b1)
Sim_PH_500_DapaHF_.5_censSC30<-Sim_PHSC(m=m1,N=500,param=ML_0$par,rateC =rateC_PH_DapaHF_.5_censSC30,prob=0.5,b=b1)
Sim_PH_100_DapaHF_.5_censSC30<-Sim_PHSC(m=m1,N=100,param=ML_0$par,rateC =rateC_PH_DapaHF_.5_censSC30,prob=0.5,b=b1)
Sim_PH_50_DapaHF_.5_censSC30<-Sim_PHSC(m=m1,N=50,param=ML_0$par,rateC =rateC_PH_DapaHF_.5_censSC30,prob=0.5,b=b1)

#50
rateC_PH_DapaHF_.5_censSC50<-c(0,150)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_.5_censSC50,prob=0.5,b=b1)
Sim_PH_500_DapaHF_.5_censSC50<-Sim_PHSC(m=m1,N=500,param=ML_0$par,rateC =rateC_PH_DapaHF_.5_censSC50,prob=0.5,b=b1)
Sim_PH_100_DapaHF_.5_censSC50<-Sim_PHSC(m=m1,N=100,param=ML_0$par,rateC =rateC_PH_DapaHF_.5_censSC50,prob=0.5,b=b1)
Sim_PH_50_DapaHF_.5_censSC50<-Sim_PHSC(m=m1,N=50,param=ML_0$par,rateC =rateC_PH_DapaHF_.5_censSC50,prob=0.5,b=b1)

#70
rateC_PH_DapaHF_.5_censSC70<-c(0,75)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_.5_censSC70,prob=0.5,b=b1)
Sim_PH_500_DapaHF_.5_censSC70<-Sim_PHSC(m=m1,N=500,param=ML_0$par,rateC =rateC_PH_DapaHF_.5_censSC70,prob=0.5,b=b1)
Sim_PH_100_DapaHF_.5_censSC70<-Sim_PHSC(m=m1,N=100,param=ML_0$par,rateC =rateC_PH_DapaHF_.5_censSC70,prob=0.5,b=b1)
Sim_PH_50_DapaHF_.5_censSC70<-Sim_PHSC(m=m1,N=50,param=ML_0$par,rateC =rateC_PH_DapaHF_.5_censSC70,prob=0.5,b=b1)


#b=-0.5
b1<--0.5
#30
rateC_PH_DapaHF_n.5_censSC30<-c(0,220)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_n.5_censSC30,prob=0.5,b=b1)
Sim_PH_500_DapaHF_n.5_censSC30<-Sim_PHSC(m=m1,N=500,param=ML_0$par,rateC =rateC_PH_DapaHF_n.5_censSC30,prob=0.5,b=b1)
Sim_PH_100_DapaHF_n.5_censSC30<-Sim_PHSC(m=m1,N=100,param=ML_0$par,rateC =rateC_PH_DapaHF_n.5_censSC30,prob=0.5,b=b1)
Sim_PH_50_DapaHF_n.5_censSC30<-Sim_PHSC(m=m1,N=50,param=ML_0$par,rateC =rateC_PH_DapaHF_n.5_censSC30,prob=0.5,b=b1)

#50
rateC_PH_DapaHF_n.5_censSC50<-c(0,110)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_n.5_censSC50,prob=0.5,b=b1)
Sim_PH_500_DapaHF_n.5_censSC50<-Sim_PHSC(m=m1,N=500,param=ML_0$par,rateC =rateC_PH_DapaHF_n.5_censSC50,prob=0.5,b=b1)
Sim_PH_100_DapaHF_n.5_censSC50<-Sim_PHSC(m=m1,N=100,param=ML_0$par,rateC =rateC_PH_DapaHF_n.5_censSC50,prob=0.5,b=b1)
Sim_PH_50_DapaHF_n.5_censSC50<-Sim_PHSC(m=m1,N=50,param=ML_0$par,rateC =rateC_PH_DapaHF_n.5_censSC50,prob=0.5,b=b1)

#70
rateC_PH_DapaHF_n.5_censSC70<-c(0,45)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_n.5_censSC70,prob=0.5,b=b1)
Sim_PH_500_DapaHF_n.5_censSC70<-Sim_PHSC(m=m1,N=500,param=ML_0$par,rateC =rateC_PH_DapaHF_n.5_censSC70,prob=0.5,b=b1)
Sim_PH_100_DapaHF_n.5_censSC70<-Sim_PHSC(m=m1,N=100,param=ML_0$par,rateC =rateC_PH_DapaHF_n.5_censSC70,prob=0.5,b=b1)
Sim_PH_50_DapaHF_n.5_censSC70<-Sim_PHSC(m=m1,N=50,param=ML_0$par,rateC =rateC_PH_DapaHF_n.5_censSC70,prob=0.5,b=b1)


#b=0.25
b1<-0.25 
#30
rateC_PH_DapaHF_.25_censSC30<-c(0,320)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_.25_censSC30,prob=0.5,b=b1)
Sim_PH_500_DapaHF_.25_censSC30<-Sim_PHSC(m=m1,N=500,param=ML_0$par,rateC =rateC_PH_DapaHF_.25_censSC30,prob=0.5,b=b1)
Sim_PH_100_DapaHF_.25_censSC30<-Sim_PHSC(m=m1,N=100,param=ML_0$par,rateC =rateC_PH_DapaHF_.25_censSC30,prob=0.5,b=b1)
Sim_PH_50_DapaHF_.25_censSC30<-Sim_PHSC(m=m1,N=50,param=ML_0$par,rateC =rateC_PH_DapaHF_.25_censSC30,prob=0.5,b=b1)

#50
rateC_PH_DapaHF_.25_censSC50<-c(0,155)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_.25_censSC50,prob=0.5,b=b1)
Sim_PH_500_DapaHF_.25_censSC50<-Sim_PHSC(m=m1,N=500,param=ML_0$par,rateC =rateC_PH_DapaHF_.25_censSC50,prob=0.5,b=b1)
Sim_PH_100_DapaHF_.25_censSC50<-Sim_PHSC(m=m1,N=100,param=ML_0$par,rateC =rateC_PH_DapaHF_.25_censSC50,prob=0.5,b=b1)
Sim_PH_50_DapaHF_.25_censSC50<-Sim_PHSC(m=m1,N=50,param=ML_0$par,rateC =rateC_PH_DapaHF_.25_censSC50,prob=0.5,b=b1)

#70
rateC_PH_DapaHF_.25_censSC70<-c(0,70)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_.25_censSC70,prob=0.5,b=b1)
Sim_PH_500_DapaHF_.25_censSC70<-Sim_PHSC(m=m1,N=500,param=ML_0$par,rateC =rateC_PH_DapaHF_.25_censSC70,prob=0.5,b=b1)
Sim_PH_100_DapaHF_.25_censSC70<-Sim_PHSC(m=m1,N=100,param=ML_0$par,rateC =rateC_PH_DapaHF_.25_censSC70,prob=0.5,b=b1)
Sim_PH_50_DapaHF_.25_censSC70<-Sim_PHSC(m=m1,N=50,param=ML_0$par,rateC =rateC_PH_DapaHF_.25_censSC70,prob=0.5,b=b1)


#b=-0.25
b1<--0.25
#30
rateC_PH_DapaHF_n.25_censSC30<-c(0,250)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_n.25_censSC30,prob=0.5,b=b1)
Sim_PH_500_DapaHF_n.25_censSC30<-Sim_PHSC(m=m1,N=500,param=ML_0$par,rateC =rateC_PH_DapaHF_n.25_censSC30,prob=0.5,b=b1)
Sim_PH_100_DapaHF_n.25_censSC30<-Sim_PHSC(m=m1,N=100,param=ML_0$par,rateC =rateC_PH_DapaHF_n.25_censSC30,prob=0.5,b=b1)
Sim_PH_50_DapaHF_n.25_censSC30<-Sim_PHSC(m=m1,N=50,param=ML_0$par,rateC =rateC_PH_DapaHF_n.25_censSC30,prob=0.5,b=b1)

#50
rateC_PH_DapaHF_n.25_censSC50<-c(0,115)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_n.25_censSC50,prob=0.5,b=b1)
Sim_PH_500_DapaHF_n.25_censSC50<-Sim_PHSC(m=m1,N=500,param=ML_0$par,rateC =rateC_PH_DapaHF_n.25_censSC50,prob=0.5,b=b1)
Sim_PH_100_DapaHF_n.25_censSC50<-Sim_PHSC(m=m1,N=100,param=ML_0$par,rateC =rateC_PH_DapaHF_n.25_censSC50,prob=0.5,b=b1)
Sim_PH_50_DapaHF_n.25_censSC50<-Sim_PHSC(m=m1,N=50,param=ML_0$par,rateC =rateC_PH_DapaHF_n.25_censSC50,prob=0.5,b=b1)

#70
rateC_PH_DapaHF_n.25_censSC70<-c(0,50)
simulPHWeib_SC_censrate(N=10000,param=ML_0$par,rateC=rateC_PH_DapaHF_n.25_censSC70,prob=0.5,b=b1)
Sim_PH_500_DapaHF_n.25_censSC70<-Sim_PHSC(m=m1,N=500,param=ML_0$par,rateC =rateC_PH_DapaHF_n.25_censSC70,prob=0.5,b=b1)
Sim_PH_100_DapaHF_n.25_censSC70<-Sim_PHSC(m=m1,N=100,param=ML_0$par,rateC =rateC_PH_DapaHF_n.25_censSC70,prob=0.5,b=b1)
Sim_PH_50_DapaHF_n.25_censSC70<-Sim_PHSC(m=m1,N=50,param=ML_0$par,rateC =rateC_PH_DapaHF_n.25_censSC70,prob=0.5,b=b1)



####Estimation: PR data####
#NPPR
#b=0                                                                       #true value
#30                                                                        #censoring rate
Est_EU_500_DapaHF_0_censSC30<-Sim_est(sims=Sim_EU_500_DapaHF_0_censSC30)
Est_EU_100_DapaHF_0_censSC30<-Sim_est(sims=Sim_EU_100_DapaHF_0_censSC30)
Est_EU_50_DapaHF_0_censSC30<-Sim_est(sims=Sim_EU_50_DapaHF_0_censSC30)

#50
Est_EU_500_DapaHF_0_censSC50<-Sim_est(sims=Sim_EU_500_DapaHF_0_censSC50)
Est_EU_100_DapaHF_0_censSC50<-Sim_est(sims=Sim_EU_100_DapaHF_0_censSC50)
Est_EU_50_DapaHF_0_censSC50<-Sim_est(sims=Sim_EU_50_DapaHF_0_censSC50)

#70
Est_EU_500_DapaHF_0_censSC70<-Sim_est(sims=Sim_EU_500_DapaHF_0_censSC70)
Est_EU_100_DapaHF_0_censSC70<-Sim_est(sims=Sim_EU_100_DapaHF_0_censSC70)
Est_EU_50_DapaHF_0_censSC70<-Sim_est(sims=Sim_EU_50_DapaHF_0_censSC70)


#b=0.5
#30
Est_EU_500_DapaHF_.5_censSC30<-Sim_est(sims=Sim_EU_500_DapaHF_.5_censSC30)
Est_EU_100_DapaHF_.5_censSC30<-Sim_est(sims=Sim_EU_100_DapaHF_.5_censSC30)
Est_EU_50_DapaHF_.5_censSC30<-Sim_est(sims=Sim_EU_50_DapaHF_.5_censSC30)

#50
Est_EU_500_DapaHF_.5_censSC50<-Sim_est(sims=Sim_EU_500_DapaHF_.5_censSC50)
Est_EU_100_DapaHF_.5_censSC50<-Sim_est(sims=Sim_EU_100_DapaHF_.5_censSC50)
Est_EU_50_DapaHF_.5_censSC50<-Sim_est(sims=Sim_EU_50_DapaHF_.5_censSC50)

#70
Est_EU_500_DapaHF_.5_censSC70<-Sim_est(sims=Sim_EU_500_DapaHF_.5_censSC70)
Est_EU_100_DapaHF_.5_censSC70<-Sim_est(sims=Sim_EU_100_DapaHF_.5_censSC70)
Est_EU_50_DapaHF_.5_censSC70<-Sim_est(sims=Sim_EU_50_DapaHF_.5_censSC70)


#b=-0.5
#30
Est_EU_500_DapaHF_n.5_censSC30<-Sim_est(sims=Sim_EU_500_DapaHF_n.5_censSC30)
Est_EU_100_DapaHF_n.5_censSC30<-Sim_est(sims=Sim_EU_100_DapaHF_n.5_censSC30)
Est_EU_50_DapaHF_n.5_censSC30<-Sim_est(sims=Sim_EU_50_DapaHF_n.5_censSC30)

#50
Est_EU_500_DapaHF_n.5_censSC50<-Sim_est(sims=Sim_EU_500_DapaHF_n.5_censSC50)
Est_EU_100_DapaHF_n.5_censSC50<-Sim_est(sims=Sim_EU_100_DapaHF_n.5_censSC50)
Est_EU_50_DapaHF_n.5_censSC50<-Sim_est(sims=Sim_EU_50_DapaHF_n.5_censSC50)

#70
Est_EU_500_DapaHF_n.5_censSC70<-Sim_est(sims=Sim_EU_500_DapaHF_n.5_censSC70)
Est_EU_100_DapaHF_n.5_censSC70<-Sim_est(sims=Sim_EU_100_DapaHF_n.5_censSC70)
Est_EU_50_DapaHF_n.5_censSC70<-Sim_est(sims=Sim_EU_50_DapaHF_n.5_censSC70)


#b=0.25
#30
Est_EU_500_DapaHF_.25_censSC30<-Sim_est(sims=Sim_EU_500_DapaHF_.25_censSC30)
Est_EU_100_DapaHF_.25_censSC30<-Sim_est(sims=Sim_EU_100_DapaHF_.25_censSC30)
Est_EU_50_DapaHF_.25_censSC30<-Sim_est(sims=Sim_EU_50_DapaHF_.25_censSC30)

#50
Est_EU_500_DapaHF_.25_censSC50<-Sim_est(sims=Sim_EU_500_DapaHF_.25_censSC50)
Est_EU_100_DapaHF_.25_censSC50<-Sim_est(sims=Sim_EU_100_DapaHF_.25_censSC50)
Est_EU_50_DapaHF_.25_censSC50<-Sim_est(sims=Sim_EU_50_DapaHF_.25_censSC50)

#70
Est_EU_500_DapaHF_.25_censSC70<-Sim_est(sims=Sim_EU_500_DapaHF_.25_censSC70)
Est_EU_100_DapaHF_.25_censSC70<-Sim_est(sims=Sim_EU_100_DapaHF_.25_censSC70)
Est_EU_50_DapaHF_.25_censSC70<-Sim_est(sims=Sim_EU_50_DapaHF_.25_censSC70)


#b=-0.25
#30
Est_EU_500_DapaHF_n.25_censSC30<-Sim_est(sims=Sim_EU_500_DapaHF_n.25_censSC30)
Est_EU_100_DapaHF_n.25_censSC30<-Sim_est(sims=Sim_EU_100_DapaHF_n.25_censSC30)
Est_EU_50_DapaHF_n.25_censSC30<-Sim_est(sims=Sim_EU_50_DapaHF_n.25_censSC30)

#50
Est_EU_500_DapaHF_n.25_censSC50<-Sim_est(sims=Sim_EU_500_DapaHF_n.25_censSC50)
Est_EU_100_DapaHF_n.25_censSC50<-Sim_est(sims=Sim_EU_100_DapaHF_n.25_censSC50)
Est_EU_50_DapaHF_n.25_censSC50<-Sim_est(sims=Sim_EU_50_DapaHF_n.25_censSC50)

#70
Est_EU_500_DapaHF_n.25_censSC70<-Sim_est(sims=Sim_EU_500_DapaHF_n.25_censSC70)
Est_EU_100_DapaHF_n.25_censSC70<-Sim_est(sims=Sim_EU_100_DapaHF_n.25_censSC70)
Est_EU_50_DapaHF_n.25_censSC70<-Sim_est(sims=Sim_EU_50_DapaHF_n.25_censSC70)



#PPR
#b=0                                                                          #true value
b1<-0                                                                         
staASV<-c(ML_EU$par[1],ML_EU$par[2],exp(-b1)^(-1/ML_EU$par[1])*ML_EU$par[2])  #start values are true values
#30                                                                #censoring rate
EU_ASV_Est_EU_500_DapaHF_0_censSC30<-Sim_EU_ASV_Est(sims=Sim_EU_500_DapaHF_0_censSC30,start=staASV)
EU_ASV_Est_EU_100_DapaHF_0_censSC30<-Sim_EU_ASV_Est(sims=Sim_EU_100_DapaHF_0_censSC30,start=staASV)
EU_ASV_Est_EU_50_DapaHF_0_censSC30<-Sim_EU_ASV_Est(sims=Sim_EU_50_DapaHF_0_censSC30,start=staASV)

#50
EU_ASV_Est_EU_500_DapaHF_0_censSC50<-Sim_EU_ASV_Est(sims=Sim_EU_500_DapaHF_0_censSC50,start=staASV)
EU_ASV_Est_EU_100_DapaHF_0_censSC50<-Sim_EU_ASV_Est(sims=Sim_EU_100_DapaHF_0_censSC50,start=staASV)
EU_ASV_Est_EU_50_DapaHF_0_censSC50<-Sim_EU_ASV_Est(sims=Sim_EU_50_DapaHF_0_censSC50,start=staASV)

#70
EU_ASV_Est_EU_500_DapaHF_0_censSC70<-Sim_EU_ASV_Est(sims=Sim_EU_500_DapaHF_0_censSC70,start=staASV)
EU_ASV_Est_EU_100_DapaHF_0_censSC70<-Sim_EU_ASV_Est(sims=Sim_EU_100_DapaHF_0_censSC70,start=staASV)
EU_ASV_Est_EU_50_DapaHF_0_censSC70<-Sim_EU_ASV_Est(sims=Sim_EU_50_DapaHF_0_censSC70,start=staASV)


#b=0.5
b1<-0.5
staASV<-c(ML_EU$par[1],ML_EU$par[2],exp(-b1)^(-1/ML_EU$par[1])*ML_EU$par[2])
#30
EU_ASV_Est_EU_500_DapaHF_.5_censSC30<-Sim_EU_ASV_Est(sims=Sim_EU_500_DapaHF_.5_censSC30,start=staASV)
EU_ASV_Est_EU_100_DapaHF_.5_censSC30<-Sim_EU_ASV_Est(sims=Sim_EU_100_DapaHF_.5_censSC30,start=staASV)
EU_ASV_Est_EU_50_DapaHF_.5_censSC30<-Sim_EU_ASV_Est(sims=Sim_EU_50_DapaHF_.5_censSC30,start=staASV)

#50
EU_ASV_Est_EU_500_DapaHF_.5_censSC50<-Sim_EU_ASV_Est(sims=Sim_EU_500_DapaHF_.5_censSC50,start=staASV)
EU_ASV_Est_EU_100_DapaHF_.5_censSC50<-Sim_EU_ASV_Est(sims=Sim_EU_100_DapaHF_.5_censSC50,start=staASV)
EU_ASV_Est_EU_50_DapaHF_.5_censSC50<-Sim_EU_ASV_Est(sims=Sim_EU_50_DapaHF_.5_censSC50,start=staASV)

#70
EU_ASV_Est_EU_500_DapaHF_.5_censSC70<-Sim_EU_ASV_Est(sims=Sim_EU_500_DapaHF_.5_censSC70,start=staASV)
EU_ASV_Est_EU_100_DapaHF_.5_censSC70<-Sim_EU_ASV_Est(sims=Sim_EU_100_DapaHF_.5_censSC70,start=staASV)
EU_ASV_Est_EU_50_DapaHF_.5_censSC70<-Sim_EU_ASV_Est(sims=Sim_EU_50_DapaHF_.5_censSC70,start=staASV)


#b=-0.5
b1<--0.5
staASV<-c(ML_EU$par[1],ML_EU$par[2],exp(-b1)^(-1/ML_EU$par[1])*ML_EU$par[2])
#30
EU_ASV_Est_EU_500_DapaHF_n.5_censSC30<-Sim_EU_ASV_Est(sims=Sim_EU_500_DapaHF_n.5_censSC30,start=staASV)
EU_ASV_Est_EU_100_DapaHF_n.5_censSC30<-Sim_EU_ASV_Est(sims=Sim_EU_100_DapaHF_n.5_censSC30,start=staASV)
EU_ASV_Est_EU_50_DapaHF_n.5_censSC30<-Sim_EU_ASV_Est(sims=Sim_EU_50_DapaHF_n.5_censSC30,start=staASV)

#50
EU_ASV_Est_EU_500_DapaHF_n.5_censSC50<-Sim_EU_ASV_Est(sims=Sim_EU_500_DapaHF_n.5_censSC50,start=staASV)
EU_ASV_Est_EU_100_DapaHF_n.5_censSC50<-Sim_EU_ASV_Est(sims=Sim_EU_100_DapaHF_n.5_censSC50,start=staASV)
EU_ASV_Est_EU_50_DapaHF_n.5_censSC50<-Sim_EU_ASV_Est(sims=Sim_EU_50_DapaHF_n.5_censSC50,start=staASV)

#70
EU_ASV_Est_EU_500_DapaHF_n.5_censSC70<-Sim_EU_ASV_Est(sims=Sim_EU_500_DapaHF_n.5_censSC70,start=staASV)
EU_ASV_Est_EU_100_DapaHF_n.5_censSC70<-Sim_EU_ASV_Est(sims=Sim_EU_100_DapaHF_n.5_censSC70,start=staASV)
EU_ASV_Est_EU_50_DapaHF_n.5_censSC70<-Sim_EU_ASV_Est(sims=Sim_EU_50_DapaHF_n.5_censSC70,start=staASV)


#b=0.25
b1<-0.25
staASV<-c(ML_EU$par[1],ML_EU$par[2],exp(-b1)^(-1/ML_EU$par[1])*ML_EU$par[2])
#30
EU_ASV_Est_EU_500_DapaHF_.25_censSC30<-Sim_EU_ASV_Est(sims=Sim_EU_500_DapaHF_.25_censSC30,start=staASV)
EU_ASV_Est_EU_100_DapaHF_.25_censSC30<-Sim_EU_ASV_Est(sims=Sim_EU_100_DapaHF_.25_censSC30,start=staASV)
EU_ASV_Est_EU_50_DapaHF_.25_censSC30<-Sim_EU_ASV_Est(sims=Sim_EU_50_DapaHF_.25_censSC30,start=staASV)

#50
EU_ASV_Est_EU_500_DapaHF_.25_censSC50<-Sim_EU_ASV_Est(sims=Sim_EU_500_DapaHF_.25_censSC50,start=staASV)
EU_ASV_Est_EU_100_DapaHF_.25_censSC50<-Sim_EU_ASV_Est(sims=Sim_EU_100_DapaHF_.25_censSC50,start=staASV)
EU_ASV_Est_EU_50_DapaHF_.25_censSC50<-Sim_EU_ASV_Est(sims=Sim_EU_50_DapaHF_.25_censSC50,start=staASV)

#70
EU_ASV_Est_EU_500_DapaHF_.25_censSC70<-Sim_EU_ASV_Est(sims=Sim_EU_500_DapaHF_.25_censSC70,start=staASV)
EU_ASV_Est_EU_100_DapaHF_.25_censSC70<-Sim_EU_ASV_Est(sims=Sim_EU_100_DapaHF_.25_censSC70,start=staASV)
EU_ASV_Est_EU_50_DapaHF_.25_censSC70<-Sim_EU_ASV_Est(sims=Sim_EU_50_DapaHF_.25_censSC70,start=staASV)


#b=-0.25
b1<--0.25
staASV<-c(ML_EU$par[1],ML_EU$par[2],exp(-b1)^(-1/ML_EU$par[1])*ML_EU$par[2])
#30
EU_ASV_Est_EU_500_DapaHF_n.25_censSC30<-Sim_EU_ASV_Est(sims=Sim_EU_500_DapaHF_n.25_censSC30,start=staASV)
EU_ASV_Est_EU_100_DapaHF_n.25_censSC30<-Sim_EU_ASV_Est(sims=Sim_EU_100_DapaHF_n.25_censSC30,start=staASV)
EU_ASV_Est_EU_50_DapaHF_n.25_censSC30<-Sim_EU_ASV_Est(sims=Sim_EU_50_DapaHF_n.25_censSC30,start=staASV)

#50
EU_ASV_Est_EU_500_DapaHF_n.25_censSC50<-Sim_EU_ASV_Est(sims=Sim_EU_500_DapaHF_n.25_censSC50,start=staASV)
EU_ASV_Est_EU_100_DapaHF_n.25_censSC50<-Sim_EU_ASV_Est(sims=Sim_EU_100_DapaHF_n.25_censSC50,start=staASV)
EU_ASV_Est_EU_50_DapaHF_n.25_censSC50<-Sim_EU_ASV_Est(sims=Sim_EU_50_DapaHF_n.25_censSC50,start=staASV)

#70
EU_ASV_Est_EU_500_DapaHF_n.25_censSC70<-Sim_EU_ASV_Est(sims=Sim_EU_500_DapaHF_n.25_censSC70,start=staASV)
EU_ASV_Est_EU_100_DapaHF_n.25_censSC70<-Sim_EU_ASV_Est(sims=Sim_EU_100_DapaHF_n.25_censSC70,start=staASV)
EU_ASV_Est_EU_50_DapaHF_n.25_censSC70<-Sim_EU_ASV_Est(sims=Sim_EU_50_DapaHF_n.25_censSC70,start=staASV)

####Estimation: PH data####
#NPPR
#b=0                                                                     #true value
#30                                                                      #censoring rate
Est_PH_500_DapaHF_0_censSC30<-Sim_est(sims=Sim_PH_500_DapaHF_0_censSC30)
Est_PH_100_DapaHF_0_censSC30<-Sim_est(sims=Sim_PH_100_DapaHF_0_censSC30)
Est_PH_50_DapaHF_0_censSC30<-Sim_est(sims=Sim_PH_50_DapaHF_0_censSC30)

#50
Est_PH_500_DapaHF_0_censSC50<-Sim_est(sims=Sim_PH_500_DapaHF_0_censSC50)
Est_PH_100_DapaHF_0_censSC50<-Sim_est(sims=Sim_PH_100_DapaHF_0_censSC50)
Est_PH_50_DapaHF_0_censSC50<-Sim_est(sims=Sim_PH_50_DapaHF_0_censSC50)

#70
Est_PH_500_DapaHF_0_censSC70<-Sim_est(sims=Sim_PH_500_DapaHF_0_censSC70)
Est_PH_100_DapaHF_0_censSC70<-Sim_est(sims=Sim_PH_100_DapaHF_0_censSC70)
Est_PH_50_DapaHF_0_censSC70<-Sim_est(sims=Sim_PH_50_DapaHF_0_censSC70)


#b=0.5
#30
Est_PH_500_DapaHF_.5_censSC30<-Sim_est(sims=Sim_PH_500_DapaHF_.5_censSC30)
Est_PH_100_DapaHF_.5_censSC30<-Sim_est(sims=Sim_PH_100_DapaHF_.5_censSC30)
Est_PH_50_DapaHF_.5_censSC30<-Sim_est(sims=Sim_PH_50_DapaHF_.5_censSC30)

#50
Est_PH_500_DapaHF_.5_censSC50<-Sim_est(sims=Sim_PH_500_DapaHF_.5_censSC50)
Est_PH_100_DapaHF_.5_censSC50<-Sim_est(sims=Sim_PH_100_DapaHF_.5_censSC50)
Est_PH_50_DapaHF_.5_censSC50<-Sim_est(sims=Sim_PH_50_DapaHF_.5_censSC50)

#70
Est_PH_500_DapaHF_.5_censSC70<-Sim_est(sims=Sim_PH_500_DapaHF_.5_censSC70)
Est_PH_100_DapaHF_.5_censSC70<-Sim_est(sims=Sim_PH_100_DapaHF_.5_censSC70)
Est_PH_50_DapaHF_.5_censSC70<-Sim_est(sims=Sim_PH_50_DapaHF_.5_censSC70)


#b=-0.5
#30
Est_PH_500_DapaHF_n.5_censSC30<-Sim_est(sims=Sim_PH_500_DapaHF_n.5_censSC30)
Est_PH_100_DapaHF_n.5_censSC30<-Sim_est(sims=Sim_PH_100_DapaHF_n.5_censSC30)
Est_PH_50_DapaHF_n.5_censSC30<-Sim_est(sims=Sim_PH_50_DapaHF_n.5_censSC30)

#50
Est_PH_500_DapaHF_n.5_censSC50<-Sim_est(sims=Sim_PH_500_DapaHF_n.5_censSC50)
Est_PH_100_DapaHF_n.5_censSC50<-Sim_est(sims=Sim_PH_100_DapaHF_n.5_censSC50)
Est_PH_50_DapaHF_n.5_censSC50<-Sim_est(sims=Sim_PH_50_DapaHF_n.5_censSC50)

#70
Est_PH_500_DapaHF_n.5_censSC70<-Sim_est(sims=Sim_PH_500_DapaHF_n.5_censSC70)
Est_PH_100_DapaHF_n.5_censSC70<-Sim_est(sims=Sim_PH_100_DapaHF_n.5_censSC70)
Est_PH_50_DapaHF_n.5_censSC70<-Sim_est(sims=Sim_PH_50_DapaHF_n.5_censSC70)


#b=0.25
#30
Est_PH_500_DapaHF_.25_censSC30<-Sim_est(sims=Sim_PH_500_DapaHF_.25_censSC30)
Est_PH_100_DapaHF_.25_censSC30<-Sim_est(sims=Sim_PH_100_DapaHF_.25_censSC30)
Est_PH_50_DapaHF_.25_censSC30<-Sim_est(sims=Sim_PH_50_DapaHF_.25_censSC30)

#50
Est_PH_500_DapaHF_.25_censSC50<-Sim_est(sims=Sim_PH_500_DapaHF_.25_censSC50)
Est_PH_100_DapaHF_.25_censSC50<-Sim_est(sims=Sim_PH_100_DapaHF_.25_censSC50)
Est_PH_50_DapaHF_.25_censSC50<-Sim_est(sims=Sim_PH_50_DapaHF_.25_censSC50)

#70
Est_PH_500_DapaHF_.25_censSC70<-Sim_est(sims=Sim_PH_500_DapaHF_.25_censSC70)
Est_PH_100_DapaHF_.25_censSC70<-Sim_est(sims=Sim_PH_100_DapaHF_.25_censSC70)
Est_PH_50_DapaHF_.25_censSC70<-Sim_est(sims=Sim_PH_50_DapaHF_.25_censSC70)


#b=-0.25
#30
Est_PH_500_DapaHF_n.25_censSC30<-Sim_est(sims=Sim_PH_500_DapaHF_n.25_censSC30)
Est_PH_100_DapaHF_n.25_censSC30<-Sim_est(sims=Sim_PH_100_DapaHF_n.25_censSC30)
Est_PH_50_DapaHF_n.25_censSC30<-Sim_est(sims=Sim_PH_50_DapaHF_n.25_censSC30)

#50
Est_PH_500_DapaHF_n.25_censSC50<-Sim_est(sims=Sim_PH_500_DapaHF_n.25_censSC50)
Est_PH_100_DapaHF_n.25_censSC50<-Sim_est(sims=Sim_PH_100_DapaHF_n.25_censSC50)
Est_PH_50_DapaHF_n.25_censSC50<-Sim_est(sims=Sim_PH_50_DapaHF_n.25_censSC50)

#70
Est_PH_500_DapaHF_n.25_censSC70<-Sim_est(sims=Sim_PH_500_DapaHF_n.25_censSC70)
Est_PH_100_DapaHF_n.25_censSC70<-Sim_est(sims=Sim_PH_100_DapaHF_n.25_censSC70)
Est_PH_50_DapaHF_n.25_censSC70<-Sim_est(sims=Sim_PH_50_DapaHF_n.25_censSC70)


####Pre-processing: NPPR estimator####
#PR data
#b=0                                            #true value
#30
sum(is.na(Est_EU_500_DapaHF_0_censSC30))        #check for NA, 1000-(results) are displayed in the tables on "numerical robustness"
range(Est_EU_500_DapaHF_0_censSC30)             #check for numerical problems
sum(is.na(Est_EU_100_DapaHF_0_censSC30))
range(Est_EU_100_DapaHF_0_censSC30)
sum(is.na(Est_EU_50_DapaHF_0_censSC30))
range(Est_EU_50_DapaHF_0_censSC30)

#50
sum(is.na(Est_EU_500_DapaHF_0_censSC50))
range(Est_EU_500_DapaHF_0_censSC50)
sum(is.na(Est_EU_100_DapaHF_0_censSC50))
range(Est_EU_100_DapaHF_0_censSC50)
sum(is.na(Est_EU_50_DapaHF_0_censSC50))
range(Est_EU_50_DapaHF_0_censSC50)

#70
sum(is.na(Est_EU_500_DapaHF_0_censSC70))
range(Est_EU_500_DapaHF_0_censSC70)
sum(is.na(Est_EU_100_DapaHF_0_censSC70))
range(Est_EU_100_DapaHF_0_censSC70)
sum(is.na(Est_EU_50_DapaHF_0_censSC70))#2                                                         
which(ifelse(is.na(Est_EU_50_DapaHF_0_censSC70)==T,1,0)==1)
Sim_EU_50_DapaHF_0_censSC70_nna<-Sim_EU_50_DapaHF_0_censSC70[!is.na(Est_EU_50_DapaHF_0_censSC70)] #exclude NA from Sims for CI
Est_EU_50_DapaHF_0_censSC70_nna<-Est_EU_50_DapaHF_0_censSC70[!is.na(Est_EU_50_DapaHF_0_censSC70)] #exclude NA from Ests for Output
range(Est_EU_50_DapaHF_0_censSC70_nna)


#b=0.5
#30
sum(is.na(Est_EU_500_DapaHF_.5_censSC30))
range(Est_EU_500_DapaHF_.5_censSC30)
sum(is.na(Est_EU_100_DapaHF_.5_censSC30))
range(Est_EU_100_DapaHF_.5_censSC30)
sum(is.na(Est_EU_50_DapaHF_.5_censSC30))
range(Est_EU_50_DapaHF_.5_censSC30)

#50
sum(is.na(Est_EU_500_DapaHF_.5_censSC50))
range(Est_EU_500_DapaHF_.5_censSC50)
sum(is.na(Est_EU_100_DapaHF_.5_censSC50))
range(Est_EU_100_DapaHF_.5_censSC50)
sum(is.na(Est_EU_50_DapaHF_.5_censSC50))
range(Est_EU_50_DapaHF_.5_censSC50)

#70
sum(is.na(Est_EU_500_DapaHF_.5_censSC70))
range(Est_EU_500_DapaHF_.5_censSC70)
sum(is.na(Est_EU_100_DapaHF_.5_censSC70))
range(Est_EU_100_DapaHF_.5_censSC70)
sum(is.na(Est_EU_50_DapaHF_.5_censSC70))#10
which(ifelse(is.na(Est_EU_50_DapaHF_.5_censSC70)==T,1,0)==1)
Sim_EU_50_DapaHF_.5_censSC70_nna<-Sim_EU_50_DapaHF_.5_censSC70[!is.na(Est_EU_50_DapaHF_.5_censSC70)]
Est_EU_50_DapaHF_.5_censSC70_nna<-Est_EU_50_DapaHF_.5_censSC70[!is.na(Est_EU_50_DapaHF_.5_censSC70)]
range(Est_EU_50_DapaHF_.5_censSC70_nna)


#b=-0.5
#30
sum(is.na(Est_EU_500_DapaHF_n.5_censSC30))
range(Est_EU_500_DapaHF_n.5_censSC30)
sum(is.na(Est_EU_100_DapaHF_n.5_censSC30))
range(Est_EU_100_DapaHF_n.5_censSC30)
sum(is.na(Est_EU_50_DapaHF_n.5_censSC30))
range(Est_EU_50_DapaHF_n.5_censSC30)

#50
sum(is.na(Est_EU_500_DapaHF_n.5_censSC50))
range(Est_EU_500_DapaHF_n.5_censSC50)
sum(is.na(Est_EU_100_DapaHF_n.5_censSC50))
range(Est_EU_100_DapaHF_n.5_censSC50)
sum(is.na(Est_EU_50_DapaHF_n.5_censSC50))
range(Est_EU_50_DapaHF_n.5_censSC50)

#70
sum(is.na(Est_EU_500_DapaHF_n.5_censSC70))
range(Est_EU_500_DapaHF_n.5_censSC70)
sum(is.na(Est_EU_100_DapaHF_n.5_censSC70))
range(Est_EU_100_DapaHF_n.5_censSC70)
sum(is.na(Est_EU_50_DapaHF_n.5_censSC70))#16
which(ifelse(is.na(Est_EU_50_DapaHF_n.5_censSC70)==T,1,0)==1)
Sim_EU_50_DapaHF_n.5_censSC70_nna<-Sim_EU_50_DapaHF_n.5_censSC70[!is.na(Est_EU_50_DapaHF_n.5_censSC70)]
Est_EU_50_DapaHF_n.5_censSC70_nna<-Est_EU_50_DapaHF_n.5_censSC70[!is.na(Est_EU_50_DapaHF_n.5_censSC70)]
range(Est_EU_50_DapaHF_n.5_censSC70_nna)


#b=0.25
#30
sum(is.na(Est_EU_500_DapaHF_.25_censSC30))
range(Est_EU_500_DapaHF_.25_censSC30)
sum(is.na(Est_EU_100_DapaHF_.25_censSC30))
range(Est_EU_100_DapaHF_.25_censSC30)
sum(is.na(Est_EU_50_DapaHF_.25_censSC30))
range(Est_EU_50_DapaHF_.25_censSC30)

#50
sum(is.na(Est_EU_500_DapaHF_.25_censSC50))
range(Est_EU_500_DapaHF_.25_censSC50)
sum(is.na(Est_EU_100_DapaHF_.25_censSC50))
range(Est_EU_100_DapaHF_.25_censSC50)
sum(is.na(Est_EU_50_DapaHF_.25_censSC50))
range(Est_EU_50_DapaHF_.25_censSC50)

#70
sum(is.na(Est_EU_500_DapaHF_.25_censSC70))
range(Est_EU_500_DapaHF_.25_censSC70)
sum(is.na(Est_EU_100_DapaHF_.25_censSC70))
range(Est_EU_100_DapaHF_.25_censSC70)
sum(is.na(Est_EU_50_DapaHF_.25_censSC70))#2
which(ifelse(is.na(Est_EU_50_DapaHF_.25_censSC70)==T,1,0)==1)
Sim_EU_50_DapaHF_.25_censSC70_nna<-Sim_EU_50_DapaHF_.25_censSC70[!is.na(Est_EU_50_DapaHF_.25_censSC70)]
Est_EU_50_DapaHF_.25_censSC70_nna<-Est_EU_50_DapaHF_.25_censSC70[!is.na(Est_EU_50_DapaHF_.25_censSC70)]
range(Est_EU_50_DapaHF_.25_censSC70_nna)


#b=-0.25
#30
sum(is.na(Est_EU_500_DapaHF_n.25_censSC30))
range(Est_EU_500_DapaHF_n.25_censSC30)
sum(is.na(Est_EU_100_DapaHF_n.25_censSC30))
range(Est_EU_100_DapaHF_n.25_censSC30)
sum(is.na(Est_EU_50_DapaHF_n.25_censSC30))
range(Est_EU_50_DapaHF_n.25_censSC30)

#50
sum(is.na(Est_EU_500_DapaHF_n.25_censSC50))
range(Est_EU_500_DapaHF_n.25_censSC50)
sum(is.na(Est_EU_100_DapaHF_n.25_censSC50))
range(Est_EU_100_DapaHF_n.25_censSC50)
sum(is.na(Est_EU_50_DapaHF_n.25_censSC50))
range(Est_EU_50_DapaHF_n.25_censSC50)

#70
sum(is.na(Est_EU_500_DapaHF_n.25_censSC70))
range(Est_EU_500_DapaHF_n.25_censSC70)
sum(is.na(Est_EU_100_DapaHF_n.25_censSC70))
range(Est_EU_100_DapaHF_n.25_censSC70)
sum(is.na(Est_EU_50_DapaHF_n.25_censSC70))#4
which(ifelse(is.na(Est_EU_50_DapaHF_n.25_censSC70)==T,1,0)==1)
Sim_EU_50_DapaHF_n.25_censSC70_nna<-Sim_EU_50_DapaHF_n.25_censSC70[!is.na(Est_EU_50_DapaHF_n.25_censSC70)]
Est_EU_50_DapaHF_n.25_censSC70_nna<-Est_EU_50_DapaHF_n.25_censSC70[!is.na(Est_EU_50_DapaHF_n.25_censSC70)]
range(Est_EU_50_DapaHF_n.25_censSC70_nna)



#PH data
#b=0                                                             
#30
sum(is.na(Est_PH_500_DapaHF_0_censSC30))
range(Est_PH_500_DapaHF_0_censSC30)
sum(is.na(Est_PH_100_DapaHF_0_censSC30))
range(Est_PH_100_DapaHF_0_censSC30)
sum(is.na(Est_PH_50_DapaHF_0_censSC30))
range(Est_PH_50_DapaHF_0_censSC30)

#50
sum(is.na(Est_PH_500_DapaHF_0_censSC50))
range(Est_PH_500_DapaHF_0_censSC50)
sum(is.na(Est_PH_100_DapaHF_0_censSC50))
range(Est_PH_100_DapaHF_0_censSC50)
sum(is.na(Est_PH_50_DapaHF_0_censSC50))
range(Est_PH_50_DapaHF_0_censSC50)

#70
sum(is.na(Est_PH_500_DapaHF_0_censSC70))
range(Est_PH_500_DapaHF_0_censSC70)
sum(is.na(Est_PH_100_DapaHF_0_censSC70))
range(Est_PH_100_DapaHF_0_censSC70)
sum(is.na(Est_PH_50_DapaHF_0_censSC70))#8
which(ifelse(is.na(Est_PH_50_DapaHF_0_censSC70)==T,1,0)==1)
Sim_PH_50_DapaHF_0_censSC70_nna<-Sim_PH_50_DapaHF_0_censSC70[!is.na(Est_PH_50_DapaHF_0_censSC70)]
Est_PH_50_DapaHF_0_censSC70_nna<-Est_PH_50_DapaHF_0_censSC70[!is.na(Est_PH_50_DapaHF_0_censSC70)]
range(Est_PH_50_DapaHF_0_censSC70_nna)


#b=0.5
#30
sum(is.na(Est_PH_500_DapaHF_.5_censSC30))
range(Est_PH_500_DapaHF_.5_censSC30)
sum(is.na(Est_PH_100_DapaHF_.5_censSC30))
range(Est_PH_100_DapaHF_.5_censSC30)
sum(is.na(Est_PH_50_DapaHF_.5_censSC30))
range(Est_PH_50_DapaHF_.5_censSC30)

#50
sum(is.na(Est_PH_500_DapaHF_.5_censSC50))
range(Est_PH_500_DapaHF_.5_censSC50)
sum(is.na(Est_PH_100_DapaHF_.5_censSC50))
range(Est_PH_100_DapaHF_.5_censSC50)
sum(is.na(Est_PH_50_DapaHF_.5_censSC50))
range(Est_PH_50_DapaHF_.5_censSC50)

#70
sum(is.na(Est_PH_500_DapaHF_.5_censSC70))
range(Est_PH_500_DapaHF_.5_censSC70)
sum(is.na(Est_PH_100_DapaHF_.5_censSC70))
range(Est_PH_100_DapaHF_.5_censSC70)
sum(is.na(Est_PH_50_DapaHF_.5_censSC70))#10
which(ifelse(is.na(Est_PH_50_DapaHF_.5_censSC70)==T,1,0)==1)
Sim_PH_50_DapaHF_.5_censSC70_nna<-Sim_PH_50_DapaHF_.5_censSC70[!is.na(Est_PH_50_DapaHF_.5_censSC70)]
Est_PH_50_DapaHF_.5_censSC70_nna<-Est_PH_50_DapaHF_.5_censSC70[!is.na(Est_PH_50_DapaHF_.5_censSC70)]
range(Est_PH_50_DapaHF_.5_censSC70_nna)


#b=-0.5
#30
sum(is.na(Est_PH_500_DapaHF_n.5_censSC30))
range(Est_PH_500_DapaHF_n.5_censSC30)
sum(is.na(Est_PH_100_DapaHF_n.5_censSC30))
range(Est_PH_100_DapaHF_n.5_censSC30)
sum(is.na(Est_PH_50_DapaHF_n.5_censSC30))
range(Est_PH_50_DapaHF_n.5_censSC30)

#50
sum(is.na(Est_PH_500_DapaHF_n.5_censSC50))
range(Est_PH_500_DapaHF_n.5_censSC50)
sum(is.na(Est_PH_100_DapaHF_n.5_censSC50))
range(Est_PH_100_DapaHF_n.5_censSC50)
sum(is.na(Est_PH_50_DapaHF_n.5_censSC50))
range(Est_PH_50_DapaHF_n.5_censSC50)

#70
sum(is.na(Est_PH_500_DapaHF_n.5_censSC70))
range(Est_PH_500_DapaHF_n.5_censSC70)
sum(is.na(Est_PH_100_DapaHF_n.5_censSC70))
range(Est_PH_100_DapaHF_n.5_censSC70)
sum(is.na(Est_PH_50_DapaHF_n.5_censSC70))#8
which(ifelse(is.na(Est_PH_50_DapaHF_n.5_censSC70)==T,1,0)==1)
Sim_PH_50_DapaHF_n.5_censSC70_nna<-Sim_PH_50_DapaHF_n.5_censSC70[!is.na(Est_PH_50_DapaHF_n.5_censSC70)]
Est_PH_50_DapaHF_n.5_censSC70_nna<-Est_PH_50_DapaHF_n.5_censSC70[!is.na(Est_PH_50_DapaHF_n.5_censSC70)]
range(Est_PH_50_DapaHF_n.5_censSC70_nna)


#b=0.25
#30
sum(is.na(Est_PH_500_DapaHF_.25_censSC30))
range(Est_PH_500_DapaHF_.25_censSC30)
sum(is.na(Est_PH_100_DapaHF_.25_censSC30))
range(Est_PH_100_DapaHF_.25_censSC30)
sum(is.na(Est_PH_50_DapaHF_.25_censSC30))
range(Est_PH_50_DapaHF_.25_censSC30)

#50
sum(is.na(Est_PH_500_DapaHF_.25_censSC50))
range(Est_PH_500_DapaHF_.25_censSC50)
sum(is.na(Est_PH_100_DapaHF_.25_censSC50))
range(Est_PH_100_DapaHF_.25_censSC50)
sum(is.na(Est_PH_50_DapaHF_.25_censSC50))
range(Est_PH_50_DapaHF_.25_censSC50)

#70
sum(is.na(Est_PH_500_DapaHF_.25_censSC70))
range(Est_PH_500_DapaHF_.25_censSC70)
sum(is.na(Est_PH_100_DapaHF_.25_censSC70))
range(Est_PH_100_DapaHF_.25_censSC70)
sum(is.na(Est_PH_50_DapaHF_.25_censSC70))#2
which(ifelse(is.na(Est_PH_50_DapaHF_.25_censSC70)==T,1,0)==1)
Sim_PH_50_DapaHF_.25_censSC70_nna<-Sim_PH_50_DapaHF_.25_censSC70[!is.na(Est_PH_50_DapaHF_.25_censSC70)]
Est_PH_50_DapaHF_.25_censSC70_nna<-Est_PH_50_DapaHF_.25_censSC70[!is.na(Est_PH_50_DapaHF_.25_censSC70)]
range(Est_PH_50_DapaHF_.25_censSC70_nna)


#b=-0.25
#30
sum(is.na(Est_PH_500_DapaHF_n.25_censSC30))
range(Est_PH_500_DapaHF_n.25_censSC30)
sum(is.na(Est_PH_100_DapaHF_n.25_censSC30))
range(Est_PH_100_DapaHF_n.25_censSC30)
sum(is.na(Est_PH_50_DapaHF_n.25_censSC30))
range(Est_PH_50_DapaHF_n.25_censSC30)

#50
sum(is.na(Est_PH_500_DapaHF_n.25_censSC50))
range(Est_PH_500_DapaHF_n.25_censSC50)
sum(is.na(Est_PH_100_DapaHF_n.25_censSC50))
range(Est_PH_100_DapaHF_n.25_censSC50)
sum(is.na(Est_PH_50_DapaHF_n.25_censSC50))
range(Est_PH_50_DapaHF_n.25_censSC50)

#70
sum(is.na(Est_PH_500_DapaHF_n.25_censSC70))
range(Est_PH_500_DapaHF_n.25_censSC70)
sum(is.na(Est_PH_100_DapaHF_n.25_censSC70))
range(Est_PH_100_DapaHF_n.25_censSC70)
sum(is.na(Est_PH_50_DapaHF_n.25_censSC70))#11
which(ifelse(is.na(Est_PH_50_DapaHF_n.25_censSC70)==T,1,0)==1)
Sim_PH_50_DapaHF_n.25_censSC70_nna<-Sim_PH_50_DapaHF_n.25_censSC70[!is.na(Est_PH_50_DapaHF_n.25_censSC70)]
Est_PH_50_DapaHF_n.25_censSC70_nna<-Est_PH_50_DapaHF_n.25_censSC70[!is.na(Est_PH_50_DapaHF_n.25_censSC70)]
range(Est_PH_50_DapaHF_n.25_censSC70_nna)



####Pre-processing: PPR model####
#b=0                                                                                           #true value
#30                                                                                            #censoring rate
Sim_EU_nna_test_log(EU_ASV_Est_EU_500_DapaHF_0_censSC30)                                       #check how often numerical problems occurred
EU_ASV_Est_EU_500_DapaHF_0_censSC30_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_500_DapaHF_0_censSC30)   #exclude and transform estimated values to RR
Sim_EU_nna_test_log(EU_ASV_Est_EU_100_DapaHF_0_censSC30)
EU_ASV_Est_EU_100_DapaHF_0_censSC30_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_100_DapaHF_0_censSC30)
Sim_EU_nna_test_log(EU_ASV_Est_EU_50_DapaHF_0_censSC30)
EU_ASV_Est_EU_50_DapaHF_0_censSC30_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_50_DapaHF_0_censSC30)

#50
Sim_EU_nna_test_log(EU_ASV_Est_EU_500_DapaHF_0_censSC50)
EU_ASV_Est_EU_500_DapaHF_0_censSC50_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_500_DapaHF_0_censSC50)
Sim_EU_nna_test_log(EU_ASV_Est_EU_100_DapaHF_0_censSC50)
EU_ASV_Est_EU_100_DapaHF_0_censSC50_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_100_DapaHF_0_censSC50)
Sim_EU_nna_test_log(EU_ASV_Est_EU_50_DapaHF_0_censSC50)
EU_ASV_Est_EU_50_DapaHF_0_censSC50_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_50_DapaHF_0_censSC50)

#70
Sim_EU_nna_test_log(EU_ASV_Est_EU_500_DapaHF_0_censSC70)
EU_ASV_Est_EU_500_DapaHF_0_censSC70_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_500_DapaHF_0_censSC70)
Sim_EU_nna_test_log(EU_ASV_Est_EU_100_DapaHF_0_censSC70)
EU_ASV_Est_EU_100_DapaHF_0_censSC70_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_100_DapaHF_0_censSC70)
Sim_EU_nna_test_log(EU_ASV_Est_EU_50_DapaHF_0_censSC70)
EU_ASV_Est_EU_50_DapaHF_0_censSC70_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_50_DapaHF_0_censSC70)


#b=0.5
#30
Sim_EU_nna_test_log(EU_ASV_Est_EU_500_DapaHF_.5_censSC30)
EU_ASV_Est_EU_500_DapaHF_.5_censSC30_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_500_DapaHF_.5_censSC30)
Sim_EU_nna_test_log(EU_ASV_Est_EU_100_DapaHF_.5_censSC30)
EU_ASV_Est_EU_100_DapaHF_.5_censSC30_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_100_DapaHF_.5_censSC30)
Sim_EU_nna_test_log(EU_ASV_Est_EU_50_DapaHF_.5_censSC30)
EU_ASV_Est_EU_50_DapaHF_.5_censSC30_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_50_DapaHF_.5_censSC30)

#50
Sim_EU_nna_test_log(EU_ASV_Est_EU_500_DapaHF_.5_censSC50)
EU_ASV_Est_EU_500_DapaHF_.5_censSC50_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_500_DapaHF_.5_censSC50)
Sim_EU_nna_test_log(EU_ASV_Est_EU_100_DapaHF_.5_censSC50)
EU_ASV_Est_EU_100_DapaHF_.5_censSC50_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_100_DapaHF_.5_censSC50)
Sim_EU_nna_test_log(EU_ASV_Est_EU_50_DapaHF_.5_censSC50)
EU_ASV_Est_EU_50_DapaHF_.5_censSC50_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_50_DapaHF_.5_censSC50)

#70
Sim_EU_nna_test_log(EU_ASV_Est_EU_500_DapaHF_.5_censSC70)
EU_ASV_Est_EU_500_DapaHF_.5_censSC70_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_500_DapaHF_.5_censSC70)
Sim_EU_nna_test_log(EU_ASV_Est_EU_100_DapaHF_.5_censSC70)
EU_ASV_Est_EU_100_DapaHF_.5_censSC70_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_100_DapaHF_.5_censSC70)
Sim_EU_nna_test_log(EU_ASV_Est_EU_50_DapaHF_.5_censSC70)
EU_ASV_Est_EU_50_DapaHF_.5_censSC70_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_50_DapaHF_.5_censSC70)


#b-=0.5
#30
Sim_EU_nna_test_log(EU_ASV_Est_EU_500_DapaHF_n.5_censSC30)
EU_ASV_Est_EU_500_DapaHF_n.5_censSC30_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_500_DapaHF_n.5_censSC30)
Sim_EU_nna_test_log(EU_ASV_Est_EU_100_DapaHF_n.5_censSC30)
EU_ASV_Est_EU_100_DapaHF_n.5_censSC30_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_100_DapaHF_n.5_censSC30)
Sim_EU_nna_test_log(EU_ASV_Est_EU_50_DapaHF_n.5_censSC30)
EU_ASV_Est_EU_50_DapaHF_n.5_censSC30_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_50_DapaHF_n.5_censSC30)

#50
Sim_EU_nna_test_log(EU_ASV_Est_EU_500_DapaHF_n.5_censSC50)
EU_ASV_Est_EU_500_DapaHF_n.5_censSC50_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_500_DapaHF_n.5_censSC50)
Sim_EU_nna_test_log(EU_ASV_Est_EU_100_DapaHF_n.5_censSC50)
EU_ASV_Est_EU_100_DapaHF_n.5_censSC50_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_100_DapaHF_n.5_censSC50)
Sim_EU_nna_test_log(EU_ASV_Est_EU_50_DapaHF_n.5_censSC50)
EU_ASV_Est_EU_50_DapaHF_n.5_censSC50_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_50_DapaHF_n.5_censSC50)

#70
Sim_EU_nna_test_log(EU_ASV_Est_EU_500_DapaHF_n.5_censSC70)
EU_ASV_Est_EU_500_DapaHF_n.5_censSC70_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_500_DapaHF_n.5_censSC70)
Sim_EU_nna_test_log(EU_ASV_Est_EU_100_DapaHF_n.5_censSC70)
EU_ASV_Est_EU_100_DapaHF_n.5_censSC70_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_100_DapaHF_n.5_censSC70)
Sim_EU_nna_test_log(EU_ASV_Est_EU_50_DapaHF_n.5_censSC70)
EU_ASV_Est_EU_50_DapaHF_n.5_censSC70_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_50_DapaHF_n.5_censSC70)


#b=0.25
#30
Sim_EU_nna_test_log(EU_ASV_Est_EU_500_DapaHF_.25_censSC30)
EU_ASV_Est_EU_500_DapaHF_.25_censSC30_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_500_DapaHF_.25_censSC30)
Sim_EU_nna_test_log(EU_ASV_Est_EU_100_DapaHF_.25_censSC30)
EU_ASV_Est_EU_100_DapaHF_.25_censSC30_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_100_DapaHF_.25_censSC30)
Sim_EU_nna_test_log(EU_ASV_Est_EU_50_DapaHF_.25_censSC30)
EU_ASV_Est_EU_50_DapaHF_.25_censSC30_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_50_DapaHF_.25_censSC30)

#50
Sim_EU_nna_test_log(EU_ASV_Est_EU_500_DapaHF_.25_censSC50)
EU_ASV_Est_EU_500_DapaHF_.25_censSC50_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_500_DapaHF_.25_censSC50)
Sim_EU_nna_test_log(EU_ASV_Est_EU_100_DapaHF_.25_censSC50)
EU_ASV_Est_EU_100_DapaHF_.25_censSC50_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_100_DapaHF_.25_censSC50)
Sim_EU_nna_test_log(EU_ASV_Est_EU_50_DapaHF_.25_censSC50)
EU_ASV_Est_EU_50_DapaHF_.25_censSC50_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_50_DapaHF_.25_censSC50)

#70
Sim_EU_nna_test_log(EU_ASV_Est_EU_500_DapaHF_.25_censSC70)
EU_ASV_Est_EU_500_DapaHF_.25_censSC70_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_500_DapaHF_.25_censSC70)
Sim_EU_nna_test_log(EU_ASV_Est_EU_100_DapaHF_.25_censSC70)
EU_ASV_Est_EU_100_DapaHF_.25_censSC70_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_100_DapaHF_.25_censSC70)
Sim_EU_nna_test_log(EU_ASV_Est_EU_50_DapaHF_.25_censSC70)
EU_ASV_Est_EU_50_DapaHF_.25_censSC70_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_50_DapaHF_.25_censSC70)


#b-=0.25
#30
Sim_EU_nna_test_log(EU_ASV_Est_EU_500_DapaHF_n.25_censSC30)
EU_ASV_Est_EU_500_DapaHF_n.25_censSC30_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_500_DapaHF_n.25_censSC30)
Sim_EU_nna_test_log(EU_ASV_Est_EU_100_DapaHF_n.25_censSC30)
EU_ASV_Est_EU_100_DapaHF_n.25_censSC30_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_100_DapaHF_n.25_censSC30)
Sim_EU_nna_test_log(EU_ASV_Est_EU_50_DapaHF_n.25_censSC30)
EU_ASV_Est_EU_50_DapaHF_n.25_censSC30_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_50_DapaHF_n.25_censSC30)

#50
Sim_EU_nna_test_log(EU_ASV_Est_EU_500_DapaHF_n.25_censSC50)
EU_ASV_Est_EU_500_DapaHF_n.25_censSC50_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_500_DapaHF_n.25_censSC50)
Sim_EU_nna_test_log(EU_ASV_Est_EU_100_DapaHF_n.25_censSC50)
EU_ASV_Est_EU_100_DapaHF_n.25_censSC50_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_100_DapaHF_n.25_censSC50)
Sim_EU_nna_test_log(EU_ASV_Est_EU_50_DapaHF_n.25_censSC50)
EU_ASV_Est_EU_50_DapaHF_n.25_censSC50_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_50_DapaHF_n.25_censSC50)

#70
Sim_EU_nna_test_log(EU_ASV_Est_EU_500_DapaHF_n.25_censSC70)
EU_ASV_Est_EU_500_DapaHF_n.25_censSC70_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_500_DapaHF_n.25_censSC70)
Sim_EU_nna_test_log(EU_ASV_Est_EU_100_DapaHF_n.25_censSC70)
EU_ASV_Est_EU_100_DapaHF_n.25_censSC70_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_100_DapaHF_n.25_censSC70)
Sim_EU_nna_test_log(EU_ASV_Est_EU_50_DapaHF_n.25_censSC70)
EU_ASV_Est_EU_50_DapaHF_n.25_censSC70_nna<-Sim_EU_nna_log(EU_ASV_Est_EU_50_DapaHF_n.25_censSC70)


####CI for NPPR estimator####
#PR data
#b=0                                                                                                   #true value
b1<-0
#30                                                                                                 #censoring rate
CIOut_EU_500_DapaHF_0_censSC30<-Sim_CI_fast(sim=Sim_EU_500_DapaHF_0_censSC30,b=b1,L=500,alpha=0.05)    #coverage for the given case
CIOut_EU_100_DapaHF_0_censSC30<-Sim_CI_fast(sim=Sim_EU_100_DapaHF_0_censSC30,b=b1,L=500,alpha=0.05)
CIOut_EU_50_DapaHF_0_censSC30<-Sim_CI_fast(sim=Sim_EU_50_DapaHF_0_censSC30,b=b1,L=500,alpha=0.05)

#50
CIOut_EU_500_DapaHF_0_censSC50<-Sim_CI_fast(sim=Sim_EU_500_DapaHF_0_censSC50,b=b1,L=500,alpha=0.05)
CIOut_EU_100_DapaHF_0_censSC50<-Sim_CI_fast(sim=Sim_EU_100_DapaHF_0_censSC50,b=b1,L=500,alpha=0.05)
CIOut_EU_50_DapaHF_0_censSC50<-Sim_CI_fast(sim=Sim_EU_50_DapaHF_0_censSC50,b=b1,L=500,alpha=0.05)

#70
CIOut_EU_500_DapaHF_0_censSC70<-Sim_CI_fast(sim=Sim_EU_500_DapaHF_0_censSC70,b=b1,L=500,alpha=0.05)
CIOut_EU_100_DapaHF_0_censSC70<-Sim_CI_fast(sim=Sim_EU_100_DapaHF_0_censSC70,b=b1,L=500,alpha=0.05)
CIOut_EU_50_DapaHF_0_censSC70<-Sim_CI_fast(sim=Sim_EU_50_DapaHF_0_censSC70_nna,b=b1,L=500,alpha=0.05)


#b=0.5
b1<-0.5
#30
CIOut_EU_500_DapaHF_.5_censSC30<-Sim_CI_fast(sim=Sim_EU_500_DapaHF_.5_censSC30,b=b1,L=500,alpha=0.05)
CIOut_EU_100_DapaHF_.5_censSC30<-Sim_CI_fast(sim=Sim_EU_100_DapaHF_.5_censSC30,b=b1,L=500,alpha=0.05)
CIOut_EU_50_DapaHF_.5_censSC30<-Sim_CI_fast(sim=Sim_EU_50_DapaHF_.5_censSC30,b=b1,L=500,alpha=0.05)

#50
CIOut_EU_500_DapaHF_.5_censSC50<-Sim_CI_fast(sim=Sim_EU_500_DapaHF_.5_censSC50,b=b1,L=500,alpha=0.05)
CIOut_EU_100_DapaHF_.5_censSC50<-Sim_CI_fast(sim=Sim_EU_100_DapaHF_.5_censSC50,b=b1,L=500,alpha=0.05)
CIOut_EU_50_DapaHF_.5_censSC50<-Sim_CI_fast(sim=Sim_EU_50_DapaHF_.5_censSC50,b=b1,L=500,alpha=0.05)

#70
CIOut_EU_500_DapaHF_.5_censSC70<-Sim_CI_fast(sim=Sim_EU_500_DapaHF_.5_censSC70,b=b1,L=500,alpha=0.05)
CIOut_EU_100_DapaHF_.5_censSC70<-Sim_CI_fast(sim=Sim_EU_100_DapaHF_.5_censSC70,b=b1,L=500,alpha=0.05)
CIOut_EU_50_DapaHF_.5_censSC70<-Sim_CI_fast(sim=Sim_EU_50_DapaHF_.5_censSC70_nna,b=b1,L=500,alpha=0.05)


#b=-0.5
b1<--0.5
#30
CIOut_EU_500_DapaHF_n.5_censSC30<-Sim_CI_fast(sim=Sim_EU_500_DapaHF_n.5_censSC30,b=b1,L=500,alpha=0.05)
CIOut_EU_100_DapaHF_n.5_censSC30<-Sim_CI_fast(sim=Sim_EU_100_DapaHF_n.5_censSC30,b=b1,L=500,alpha=0.05)
CIOut_EU_50_DapaHF_n.5_censSC30<-Sim_CI_fast(sim=Sim_EU_50_DapaHF_n.5_censSC30,b=b1,L=500,alpha=0.05)

#50
CIOut_EU_500_DapaHF_n.5_censSC50<-Sim_CI_fast(sim=Sim_EU_500_DapaHF_n.5_censSC50,b=b1,L=500,alpha=0.05)
CIOut_EU_100_DapaHF_n.5_censSC50<-Sim_CI_fast(sim=Sim_EU_100_DapaHF_n.5_censSC50,b=b1,L=500,alpha=0.05)
CIOut_EU_50_DapaHF_n.5_censSC50<-Sim_CI_fast(sim=Sim_EU_50_DapaHF_n.5_censSC50,b=b1,L=500,alpha=0.05)

#70
CIOut_EU_500_DapaHF_n.5_censSC70<-Sim_CI_fast(sim=Sim_EU_500_DapaHF_n.5_censSC70,b=b1,L=500,alpha=0.05)
CIOut_EU_100_DapaHF_n.5_censSC70<-Sim_CI_fast(sim=Sim_EU_100_DapaHF_n.5_censSC70,b=b1,L=500,alpha=0.05)
CIOut_EU_50_DapaHF_n.5_censSC70<-Sim_CI_fast(sim=Sim_EU_50_DapaHF_n.5_censSC70_nna,b=b1,L=500,alpha=0.05)


#b=0.25
b1<-0.25
#30
CIOut_EU_500_DapaHF_.25_censSC30<-Sim_CI_fast(sim=Sim_EU_500_DapaHF_.25_censSC30,b=b1,L=500,alpha=0.05)
CIOut_EU_100_DapaHF_.25_censSC30<-Sim_CI_fast(sim=Sim_EU_100_DapaHF_.25_censSC30,b=b1,L=500,alpha=0.05)
CIOut_EU_50_DapaHF_.25_censSC30<-Sim_CI_fast(sim=Sim_EU_50_DapaHF_.25_censSC30,b=b1,L=500,alpha=0.05)

#50
CIOut_EU_500_DapaHF_.25_censSC50<-Sim_CI_fast(sim=Sim_EU_500_DapaHF_.25_censSC50,b=b1,L=500,alpha=0.05)
CIOut_EU_100_DapaHF_.25_censSC50<-Sim_CI_fast(sim=Sim_EU_100_DapaHF_.25_censSC50,b=b1,L=500,alpha=0.05)
CIOut_EU_50_DapaHF_.25_censSC50<-Sim_CI_fast(sim=Sim_EU_50_DapaHF_.25_censSC50,b=b1,L=500,alpha=0.05)

#70
CIOut_EU_500_DapaHF_.25_censSC70<-Sim_CI_fast(sim=Sim_EU_500_DapaHF_.25_censSC70,b=b1,L=500,alpha=0.05)
CIOut_EU_100_DapaHF_.25_censSC70<-Sim_CI_fast(sim=Sim_EU_100_DapaHF_.25_censSC70,b=b1,L=500,alpha=0.05)
CIOut_EU_50_DapaHF_.25_censSC70<-Sim_CI_fast(sim=Sim_EU_50_DapaHF_.25_censSC70_nna,b=b1,L=500,alpha=0.05)


#b=-0.25
b1<--0.25
#30
CIOut_EU_500_DapaHF_n.25_censSC30<-Sim_CI_fast(sim=Sim_EU_500_DapaHF_n.25_censSC30,b=b1,L=500,alpha=0.05)
CIOut_EU_100_DapaHF_n.25_censSC30<-Sim_CI_fast(sim=Sim_EU_100_DapaHF_n.25_censSC30,b=b1,L=500,alpha=0.05)
CIOut_EU_50_DapaHF_n.25_censSC30<-Sim_CI_fast(sim=Sim_EU_50_DapaHF_n.25_censSC30,b=b1,L=500,alpha=0.05)

#50
CIOut_EU_500_DapaHF_n.25_censSC50<-Sim_CI_fast(sim=Sim_EU_500_DapaHF_n.25_censSC50,b=b1,L=500,alpha=0.05)
CIOut_EU_100_DapaHF_n.25_censSC50<-Sim_CI_fast(sim=Sim_EU_100_DapaHF_n.25_censSC50,b=b1,L=500,alpha=0.05)
CIOut_EU_50_DapaHF_n.25_censSC50<-Sim_CI_fast(sim=Sim_EU_50_DapaHF_n.25_censSC50,b=b1,L=500,alpha=0.05)

#70
CIOut_EU_500_DapaHF_n.25_censSC70<-Sim_CI_fast(sim=Sim_EU_500_DapaHF_n.25_censSC70,b=b1,L=500,alpha=0.05)
CIOut_EU_100_DapaHF_n.25_censSC70<-Sim_CI_fast(sim=Sim_EU_100_DapaHF_n.25_censSC70,b=b1,L=500,alpha=0.05)
CIOut_EU_50_DapaHF_n.25_censSC70<-Sim_CI_fast(sim=Sim_EU_50_DapaHF_n.25_censSC70_nna,b=b1,L=500,alpha=0.05)



#PH data
#b=0
b1<-0
#30
CIOut_PH_500_DapaHF_0_censSC30<-Sim_CI_fast(sim=Sim_PH_500_DapaHF_0_censSC30,b=b1,L=500,alpha=0.05)
CIOut_PH_100_DapaHF_0_censSC30<-Sim_CI_fast(sim=Sim_PH_100_DapaHF_0_censSC30,b=b1,L=500,alpha=0.05)
CIOut_PH_50_DapaHF_0_censSC30<-Sim_CI_fast(sim=Sim_PH_50_DapaHF_0_censSC30,b=b1,L=500,alpha=0.05)

#50
CIOut_PH_500_DapaHF_0_censSC50<-Sim_CI_fast(sim=Sim_PH_500_DapaHF_0_censSC50,b=b1,L=500,alpha=0.05)
CIOut_PH_100_DapaHF_0_censSC50<-Sim_CI_fast(sim=Sim_PH_100_DapaHF_0_censSC50,b=b1,L=500,alpha=0.05)
CIOut_PH_50_DapaHF_0_censSC50<-Sim_CI_fast(sim=Sim_PH_50_DapaHF_0_censSC50,b=b1,L=500,alpha=0.05)

#70
CIOut_PH_500_DapaHF_0_censSC70<-Sim_CI_fast(sim=Sim_PH_500_DapaHF_0_censSC70,b=b1,L=500,alpha=0.05)
CIOut_PH_100_DapaHF_0_censSC70<-Sim_CI_fast(sim=Sim_PH_100_DapaHF_0_censSC70,b=b1,L=500,alpha=0.05)
CIOut_PH_50_DapaHF_0_censSC70<-Sim_CI_fast(sim=Sim_PH_50_DapaHF_0_censSC70_nna,b=b1,L=500,alpha=0.05)


#b=0.5
b1<-0.5
#30
CIOut_PH_500_DapaHF_.5_censSC30<-Sim_CI_fast(sim=Sim_PH_500_DapaHF_.5_censSC30,b=b1,L=500,alpha=0.05)
CIOut_PH_100_DapaHF_.5_censSC30<-Sim_CI_fast(sim=Sim_PH_100_DapaHF_.5_censSC30,b=b1,L=500,alpha=0.05)
CIOut_PH_50_DapaHF_.5_censSC30<-Sim_CI_fast(sim=Sim_PH_50_DapaHF_.5_censSC30,b=b1,L=500,alpha=0.05)

#50
CIOut_PH_500_DapaHF_.5_censSC50<-Sim_CI_fast(sim=Sim_PH_500_DapaHF_.5_censSC50,b=b1,L=500,alpha=0.05)
CIOut_PH_100_DapaHF_.5_censSC50<-Sim_CI_fast(sim=Sim_PH_100_DapaHF_.5_censSC50,b=b1,L=500,alpha=0.05)
CIOut_PH_50_DapaHF_.5_censSC50<-Sim_CI_fast(sim=Sim_PH_50_DapaHF_.5_censSC50,b=b1,L=500,alpha=0.05)

#70
CIOut_PH_500_DapaHF_.5_censSC70<-Sim_CI_fast(sim=Sim_PH_500_DapaHF_.5_censSC70,b=b1,L=500,alpha=0.05)
CIOut_PH_100_DapaHF_.5_censSC70<-Sim_CI_fast(sim=Sim_PH_100_DapaHF_.5_censSC70,b=b1,L=500,alpha=0.05)
CIOut_PH_50_DapaHF_.5_censSC70<-Sim_CI_fast(sim=Sim_PH_50_DapaHF_.5_censSC70_nna,b=b1,L=500,alpha=0.05)


#b=-0.5
b1<--0.5
#30
CIOut_PH_500_DapaHF_n.5_censSC30<-Sim_CI_fast(sim=Sim_PH_500_DapaHF_n.5_censSC30,b=b1,L=500,alpha=0.05)
CIOut_PH_100_DapaHF_n.5_censSC30<-Sim_CI_fast(sim=Sim_PH_100_DapaHF_n.5_censSC30,b=b1,L=500,alpha=0.05)
CIOut_PH_50_DapaHF_n.5_censSC30<-Sim_CI_fast(sim=Sim_PH_50_DapaHF_n.5_censSC30,b=b1,L=500,alpha=0.05)

#50
CIOut_PH_500_DapaHF_n.5_censSC50<-Sim_CI_fast(sim=Sim_PH_500_DapaHF_n.5_censSC50,b=b1,L=500,alpha=0.05)
CIOut_PH_100_DapaHF_n.5_censSC50<-Sim_CI_fast(sim=Sim_PH_100_DapaHF_n.5_censSC50,b=b1,L=500,alpha=0.05)
CIOut_PH_50_DapaHF_n.5_censSC50<-Sim_CI_fast(sim=Sim_PH_50_DapaHF_n.5_censSC50,b=b1,L=500,alpha=0.05)

#70
CIOut_PH_500_DapaHF_n.5_censSC70<-Sim_CI_fast(sim=Sim_PH_500_DapaHF_n.5_censSC70,b=b1,L=500,alpha=0.05)
CIOut_PH_100_DapaHF_n.5_censSC70<-Sim_CI_fast(sim=Sim_PH_100_DapaHF_n.5_censSC70,b=b1,L=500,alpha=0.05)
CIOut_PH_50_DapaHF_n.5_censSC70<-Sim_CI_fast(sim=Sim_PH_50_DapaHF_n.5_censSC70_nna,b=b1,L=500,alpha=0.05)


#b=0.25
b1<-0.25
#30
CIOut_PH_500_DapaHF_.25_censSC30<-Sim_CI_fast(sim=Sim_PH_500_DapaHF_.25_censSC30,b=b1,L=500,alpha=0.05)
CIOut_PH_100_DapaHF_.25_censSC30<-Sim_CI_fast(sim=Sim_PH_100_DapaHF_.25_censSC30,b=b1,L=500,alpha=0.05)
CIOut_PH_50_DapaHF_.25_censSC30<-Sim_CI_fast(sim=Sim_PH_50_DapaHF_.25_censSC30,b=b1,L=500,alpha=0.05)

#50
CIOut_PH_500_DapaHF_.25_censSC50<-Sim_CI_fast(sim=Sim_PH_500_DapaHF_.25_censSC50,b=b1,L=500,alpha=0.05)
CIOut_PH_100_DapaHF_.25_censSC50<-Sim_CI_fast(sim=Sim_PH_100_DapaHF_.25_censSC50,b=b1,L=500,alpha=0.05)
CIOut_PH_50_DapaHF_.25_censSC50<-Sim_CI_fast(sim=Sim_PH_50_DapaHF_.25_censSC50,b=b1,L=500,alpha=0.05)

#70
CIOut_PH_500_DapaHF_.25_censSC70<-Sim_CI_fast(sim=Sim_PH_500_DapaHF_.25_censSC70,b=b1,L=500,alpha=0.05)
CIOut_PH_100_DapaHF_.25_censSC70<-Sim_CI_fast(sim=Sim_PH_100_DapaHF_.25_censSC70,b=b1,L=500,alpha=0.05)
CIOut_PH_50_DapaHF_.25_censSC70<-Sim_CI_fast(sim=Sim_PH_50_DapaHF_.25_censSC70_nna,b=b1,L=500,alpha=0.05)


#b=-0.25
b1<--0.25
#30
CIOut_PH_500_DapaHF_n.25_censSC30<-Sim_CI_fast(sim=Sim_PH_500_DapaHF_n.25_censSC30,b=b1,L=500,alpha=0.05)
CIOut_PH_100_DapaHF_n.25_censSC30<-Sim_CI_fast(sim=Sim_PH_100_DapaHF_n.25_censSC30,b=b1,L=500,alpha=0.05)
CIOut_PH_50_DapaHF_n.25_censSC30<-Sim_CI_fast(sim=Sim_PH_50_DapaHF_n.25_censSC30,b=b1,L=500,alpha=0.05)

#50
CIOut_PH_500_DapaHF_n.25_censSC50<-Sim_CI_fast(sim=Sim_PH_500_DapaHF_n.25_censSC50,b=b1,L=500,alpha=0.05)
CIOut_PH_100_DapaHF_n.25_censSC50<-Sim_CI_fast(sim=Sim_PH_100_DapaHF_n.25_censSC50,b=b1,L=500,alpha=0.05)
CIOut_PH_50_DapaHF_n.25_censSC50<-Sim_CI_fast(sim=Sim_PH_50_DapaHF_n.25_censSC50,b=b1,L=500,alpha=0.05)

#70
CIOut_PH_500_DapaHF_n.25_censSC70<-Sim_CI_fast(sim=Sim_PH_500_DapaHF_n.25_censSC70,b=b1,L=500,alpha=0.05)
CIOut_PH_100_DapaHF_n.25_censSC70<-Sim_CI_fast(sim=Sim_PH_100_DapaHF_n.25_censSC70,b=b1,L=500,alpha=0.05)
CIOut_PH_50_DapaHF_n.25_censSC70<-Sim_CI_fast(sim=Sim_PH_50_DapaHF_n.25_censSC70_nna,b=b1,L=500,alpha=0.05)


#Putting results together in one list
Results_PR.for.EU.CI<-list(
  CIOut_EU_500_DapaHF_0_censSC30,CIOut_EU_100_DapaHF_0_censSC30,CIOut_EU_50_DapaHF_0_censSC30,
  CIOut_EU_500_DapaHF_0_censSC50,CIOut_EU_100_DapaHF_0_censSC50,CIOut_EU_50_DapaHF_0_censSC50,
  CIOut_EU_500_DapaHF_0_censSC70,CIOut_EU_100_DapaHF_0_censSC70,CIOut_EU_50_DapaHF_0_censSC70,
  
  CIOut_EU_500_DapaHF_.5_censSC30,CIOut_EU_100_DapaHF_.5_censSC30,CIOut_EU_50_DapaHF_.5_censSC30,
  CIOut_EU_500_DapaHF_.5_censSC50,CIOut_EU_100_DapaHF_.5_censSC50,CIOut_EU_50_DapaHF_.5_censSC50,
  CIOut_EU_500_DapaHF_.5_censSC70,CIOut_EU_100_DapaHF_.5_censSC70,CIOut_EU_50_DapaHF_.5_censSC70,
  
  CIOut_EU_500_DapaHF_.25_censSC30,CIOut_EU_100_DapaHF_.25_censSC30,CIOut_EU_50_DapaHF_.25_censSC30,
  CIOut_EU_500_DapaHF_.25_censSC50,CIOut_EU_100_DapaHF_.25_censSC50,CIOut_EU_50_DapaHF_.25_censSC50,
  CIOut_EU_500_DapaHF_.25_censSC70,CIOut_EU_100_DapaHF_.25_censSC70,CIOut_EU_50_DapaHF_.25_censSC70,
  
  CIOut_EU_500_DapaHF_n.25_censSC30,CIOut_EU_100_DapaHF_n.25_censSC30,CIOut_EU_50_DapaHF_n.25_censSC30,
  CIOut_EU_500_DapaHF_n.25_censSC50,CIOut_EU_100_DapaHF_n.25_censSC50,CIOut_EU_50_DapaHF_n.25_censSC50,
  CIOut_EU_500_DapaHF_n.25_censSC70,CIOut_EU_100_DapaHF_n.25_censSC70,CIOut_EU_50_DapaHF_n.25_censSC70,
  
  CIOut_EU_500_DapaHF_n.5_censSC30,CIOut_EU_100_DapaHF_n.5_censSC30,CIOut_EU_50_DapaHF_n.5_censSC30,
  CIOut_EU_500_DapaHF_n.5_censSC50,CIOut_EU_100_DapaHF_n.5_censSC50,CIOut_EU_50_DapaHF_n.5_censSC50,
  CIOut_EU_500_DapaHF_n.5_censSC70,CIOut_EU_100_DapaHF_n.5_censSC70,CIOut_EU_50_DapaHF_n.5_censSC70
)

Results_PR.for.PH.CI<-list(
  CIOut_PH_500_DapaHF_0_censSC30,CIOut_PH_100_DapaHF_0_censSC30,CIOut_PH_50_DapaHF_0_censSC30,
  CIOut_PH_500_DapaHF_0_censSC50,CIOut_PH_100_DapaHF_0_censSC50,CIOut_PH_50_DapaHF_0_censSC50,
  CIOut_PH_500_DapaHF_0_censSC70,CIOut_PH_100_DapaHF_0_censSC70,CIOut_PH_50_DapaHF_0_censSC70,
  
  CIOut_PH_500_DapaHF_.5_censSC30,CIOut_PH_100_DapaHF_.5_censSC30,CIOut_PH_50_DapaHF_.5_censSC30,
  CIOut_PH_500_DapaHF_.5_censSC50,CIOut_PH_100_DapaHF_.5_censSC50,CIOut_PH_50_DapaHF_.5_censSC50,
  CIOut_PH_500_DapaHF_.5_censSC70,CIOut_PH_100_DapaHF_.5_censSC70,CIOut_PH_50_DapaHF_.5_censSC70,
  
  CIOut_PH_500_DapaHF_.25_censSC30,CIOut_PH_100_DapaHF_.25_censSC30,CIOut_PH_50_DapaHF_.25_censSC30,
  CIOut_PH_500_DapaHF_.25_censSC50,CIOut_PH_100_DapaHF_.25_censSC50,CIOut_PH_50_DapaHF_.25_censSC50,
  CIOut_PH_500_DapaHF_.25_censSC70,CIOut_PH_100_DapaHF_.25_censSC70,CIOut_PH_50_DapaHF_.25_censSC70,
  
  CIOut_PH_500_DapaHF_n.25_censSC30,CIOut_PH_100_DapaHF_n.25_censSC30,CIOut_PH_50_DapaHF_n.25_censSC30,
  CIOut_PH_500_DapaHF_n.25_censSC50,CIOut_PH_100_DapaHF_n.25_censSC50,CIOut_PH_50_DapaHF_n.25_censSC50,
  CIOut_PH_500_DapaHF_n.25_censSC70,CIOut_PH_100_DapaHF_n.25_censSC70,CIOut_PH_50_DapaHF_n.25_censSC70,
  
  CIOut_PH_500_DapaHF_n.5_censSC30,CIOut_PH_100_DapaHF_n.5_censSC30,CIOut_PH_50_DapaHF_n.5_censSC30,
  CIOut_PH_500_DapaHF_n.5_censSC50,CIOut_PH_100_DapaHF_n.5_censSC50,CIOut_PH_50_DapaHF_n.5_censSC50,
  CIOut_PH_500_DapaHF_n.5_censSC70,CIOut_PH_100_DapaHF_n.5_censSC70,CIOut_PH_50_DapaHF_n.5_censSC70
)

#combining results to data.frame
CI.PR.for.EU<-results_CI(Results_PR.for.EU.CI)
CI.PR.for.PH<-results_CI(Results_PR.for.PH.CI)



####CI for the PPR model####
#b=0                                                                                                               #true value
b1<-0
#30                                                                                                                #censoring rate                          
EU_ASV_CI_EU_500_DapaHF_0_censSC30<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_500_DapaHF_0_censSC30),b=b1,df=500) #coverage for the given case
EU_ASV_CI_EU_100_DapaHF_0_censSC30<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_100_DapaHF_0_censSC30),b=b1,df=100)
EU_ASV_CI_EU_50_DapaHF_0_censSC30<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_50_DapaHF_0_censSC30),b=b1,df=50)

#50
EU_ASV_CI_EU_500_DapaHF_0_censSC50<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_500_DapaHF_0_censSC50),b=b1,df=500)
EU_ASV_CI_EU_100_DapaHF_0_censSC50<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_100_DapaHF_0_censSC50),b=b1,df=100)
EU_ASV_CI_EU_50_DapaHF_0_censSC50<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_50_DapaHF_0_censSC50),b=b1,df=50)

#70
EU_ASV_CI_EU_500_DapaHF_0_censSC70<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_500_DapaHF_0_censSC70),b=b1,df=500)
EU_ASV_CI_EU_100_DapaHF_0_censSC70<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_100_DapaHF_0_censSC70),b=b1,df=100)
EU_ASV_CI_EU_50_DapaHF_0_censSC70<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_50_DapaHF_0_censSC70),b=b1,df=50)


#b=0.5
b1<-0.5
#30
EU_ASV_CI_EU_500_DapaHF_.5_censSC30<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_500_DapaHF_.5_censSC30),b=b1,df=500)
EU_ASV_CI_EU_100_DapaHF_.5_censSC30<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_100_DapaHF_.5_censSC30),b=b1,df=100)
EU_ASV_CI_EU_50_DapaHF_.5_censSC30<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_50_DapaHF_.5_censSC30),b=b1,df=50)

#50
EU_ASV_CI_EU_500_DapaHF_.5_censSC50<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_500_DapaHF_.5_censSC50),b=b1,df=500)
EU_ASV_CI_EU_100_DapaHF_.5_censSC50<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_100_DapaHF_.5_censSC50),b=b1,df=100)
EU_ASV_CI_EU_50_DapaHF_.5_censSC50<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_50_DapaHF_.5_censSC50),b=b1,df=50)

#70
EU_ASV_CI_EU_500_DapaHF_.5_censSC70<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_500_DapaHF_.5_censSC70),b=b1,df=500)
EU_ASV_CI_EU_100_DapaHF_.5_censSC70<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_100_DapaHF_.5_censSC70),b=b1,df=100)
EU_ASV_CI_EU_50_DapaHF_.5_censSC70<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_50_DapaHF_.5_censSC70),b=b1,df=50)


#b=-0.5
b1<--0.5
#30
EU_ASV_CI_EU_500_DapaHF_n.5_censSC30<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_500_DapaHF_n.5_censSC30),b=b1,df=500)
EU_ASV_CI_EU_100_DapaHF_n.5_censSC30<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_100_DapaHF_n.5_censSC30),b=b1,df=100)
EU_ASV_CI_EU_50_DapaHF_n.5_censSC30<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_50_DapaHF_n.5_censSC30),b=b1,df=50)

#50
EU_ASV_CI_EU_500_DapaHF_n.5_censSC50<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_500_DapaHF_n.5_censSC50),b=b1,df=500)
EU_ASV_CI_EU_100_DapaHF_n.5_censSC50<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_100_DapaHF_n.5_censSC50),b=b1,df=100)
EU_ASV_CI_EU_50_DapaHF_n.5_censSC50<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_50_DapaHF_n.5_censSC50),b=b1,df=50)

#70
EU_ASV_CI_EU_500_DapaHF_n.5_censSC70<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_500_DapaHF_n.5_censSC70),b=b1,df=500)
EU_ASV_CI_EU_100_DapaHF_n.5_censSC70<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_100_DapaHF_n.5_censSC70),b=b1,df=100)
EU_ASV_CI_EU_50_DapaHF_n.5_censSC70<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_50_DapaHF_n.5_censSC70),b=b1,df=50)


#b=0.25
b1<-0.25
#30
EU_ASV_CI_EU_500_DapaHF_.25_censSC30<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_500_DapaHF_.25_censSC30),b=b1,df=500)
EU_ASV_CI_EU_100_DapaHF_.25_censSC30<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_100_DapaHF_.25_censSC30),b=b1,df=100)
EU_ASV_CI_EU_50_DapaHF_.25_censSC30<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_50_DapaHF_.25_censSC30),b=b1,df=50)

#50
EU_ASV_CI_EU_500_DapaHF_.25_censSC50<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_500_DapaHF_.25_censSC50),b=b1,df=500)
EU_ASV_CI_EU_100_DapaHF_.25_censSC50<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_100_DapaHF_.25_censSC50),b=b1,df=100)
EU_ASV_CI_EU_50_DapaHF_.25_censSC50<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_50_DapaHF_.25_censSC50),b=b1,df=50)

#70
EU_ASV_CI_EU_500_DapaHF_.25_censSC70<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_500_DapaHF_.25_censSC70),b=b1,df=500)
EU_ASV_CI_EU_100_DapaHF_.25_censSC70<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_100_DapaHF_.25_censSC70),b=b1,df=100)
EU_ASV_CI_EU_50_DapaHF_.25_censSC70<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_50_DapaHF_.25_censSC70),b=b1,df=50)


#b=-0.25
b1<--0.25
#30
EU_ASV_CI_EU_500_DapaHF_n.25_censSC30<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_500_DapaHF_n.25_censSC30),b=b1,df=500)
EU_ASV_CI_EU_100_DapaHF_n.25_censSC30<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_100_DapaHF_n.25_censSC30),b=b1,df=100)
EU_ASV_CI_EU_50_DapaHF_n.25_censSC30<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_50_DapaHF_n.25_censSC30),b=b1,df=50)

#50
EU_ASV_CI_EU_500_DapaHF_n.25_censSC50<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_500_DapaHF_n.25_censSC50),b=b1,df=500)
EU_ASV_CI_EU_100_DapaHF_n.25_censSC50<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_100_DapaHF_n.25_censSC50),b=b1,df=100)
EU_ASV_CI_EU_50_DapaHF_n.25_censSC50<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_50_DapaHF_n.25_censSC50),b=b1,df=50)

#70
EU_ASV_CI_EU_500_DapaHF_n.25_censSC70<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_500_DapaHF_n.25_censSC70),b=b1,df=500)
EU_ASV_CI_EU_100_DapaHF_n.25_censSC70<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_100_DapaHF_n.25_censSC70),b=b1,df=100)
EU_ASV_CI_EU_50_DapaHF_n.25_censSC70<-Sim_EU_CI(ests=Sim_EU_nna_CI(EU_ASV_Est_EU_50_DapaHF_n.25_censSC70),b=b1,df=50)



#Putting results together in one list
Results_EU_CI<-list(
  EU_ASV_CI_EU_500_DapaHF_0_censSC30,EU_ASV_CI_EU_100_DapaHF_0_censSC30,EU_ASV_CI_EU_50_DapaHF_0_censSC30,
  EU_ASV_CI_EU_500_DapaHF_0_censSC50,EU_ASV_CI_EU_100_DapaHF_0_censSC50,EU_ASV_CI_EU_50_DapaHF_0_censSC50,
  EU_ASV_CI_EU_500_DapaHF_0_censSC70,EU_ASV_CI_EU_100_DapaHF_0_censSC70,EU_ASV_CI_EU_50_DapaHF_0_censSC70,
  
  EU_ASV_CI_EU_500_DapaHF_.5_censSC30,EU_ASV_CI_EU_100_DapaHF_.5_censSC30,EU_ASV_CI_EU_50_DapaHF_.5_censSC30,
  EU_ASV_CI_EU_500_DapaHF_.5_censSC50,EU_ASV_CI_EU_100_DapaHF_.5_censSC50,EU_ASV_CI_EU_50_DapaHF_.5_censSC50,
  EU_ASV_CI_EU_500_DapaHF_.5_censSC70,EU_ASV_CI_EU_100_DapaHF_.5_censSC70,EU_ASV_CI_EU_50_DapaHF_.5_censSC70,
  
  EU_ASV_CI_EU_500_DapaHF_.25_censSC30,EU_ASV_CI_EU_100_DapaHF_.25_censSC30,EU_ASV_CI_EU_50_DapaHF_.25_censSC30,
  EU_ASV_CI_EU_500_DapaHF_.25_censSC50,EU_ASV_CI_EU_100_DapaHF_.25_censSC50,EU_ASV_CI_EU_50_DapaHF_.25_censSC50,
  EU_ASV_CI_EU_500_DapaHF_.25_censSC70,EU_ASV_CI_EU_100_DapaHF_.25_censSC70,EU_ASV_CI_EU_50_DapaHF_.25_censSC70,
  
  EU_ASV_CI_EU_500_DapaHF_n.25_censSC30,EU_ASV_CI_EU_100_DapaHF_n.25_censSC30,EU_ASV_CI_EU_50_DapaHF_n.25_censSC30,
  EU_ASV_CI_EU_500_DapaHF_n.25_censSC50,EU_ASV_CI_EU_100_DapaHF_n.25_censSC50,EU_ASV_CI_EU_50_DapaHF_n.25_censSC50,
  EU_ASV_CI_EU_500_DapaHF_n.25_censSC70,EU_ASV_CI_EU_100_DapaHF_n.25_censSC70,EU_ASV_CI_EU_50_DapaHF_n.25_censSC70,
  
  EU_ASV_CI_EU_500_DapaHF_n.5_censSC30,EU_ASV_CI_EU_100_DapaHF_n.5_censSC30,EU_ASV_CI_EU_50_DapaHF_n.5_censSC30,
  EU_ASV_CI_EU_500_DapaHF_n.5_censSC50,EU_ASV_CI_EU_100_DapaHF_n.5_censSC50,EU_ASV_CI_EU_50_DapaHF_n.5_censSC50,
  EU_ASV_CI_EU_500_DapaHF_n.5_censSC70,EU_ASV_CI_EU_100_DapaHF_n.5_censSC70,EU_ASV_CI_EU_50_DapaHF_n.5_censSC70  
  )

#combining results to one data.frame
CI.EU<-results_CI(Results_EU_CI)


####Output NPPR estimator####
#PR data
#b=0                                                                                   #true value
b1<-0
#30                                                                                    #censoring rate
Out_log_EU_500_DapaHF_0_censSC30<-Sim_out_log(b_est=Est_EU_500_DapaHF_0_censSC30,b=b1) #data.frame including mean, Bias and MSE of the given case
Out_log_EU_100_DapaHF_0_censSC30<-Sim_out_log(b_est=Est_EU_100_DapaHF_0_censSC30,b=b1)
Out_log_EU_50_DapaHF_0_censSC30<-Sim_out_log(b_est=Est_EU_50_DapaHF_0_censSC30,b=b1)

#50
Out_log_EU_500_DapaHF_0_censSC50<-Sim_out_log(b_est=Est_EU_500_DapaHF_0_censSC50,b=b1)
Out_log_EU_100_DapaHF_0_censSC50<-Sim_out_log(b_est=Est_EU_100_DapaHF_0_censSC50,b=b1)
Out_log_EU_50_DapaHF_0_censSC50<-Sim_out_log(b_est=Est_EU_50_DapaHF_0_censSC50,b=b1)

#70
Out_log_EU_500_DapaHF_0_censSC70<-Sim_out_log(b_est=Est_EU_500_DapaHF_0_censSC70,b=b1)
Out_log_EU_100_DapaHF_0_censSC70<-Sim_out_log(b_est=Est_EU_100_DapaHF_0_censSC70,b=b1)
Out_log_EU_50_DapaHF_0_censSC70<-Sim_out_log(b_est=Est_EU_50_DapaHF_0_censSC70_nna,b=b1)


#b=0.5
b1<-0.5
#30
Out_log_EU_500_DapaHF_.5_censSC30<-Sim_out_log(b_est=Est_EU_500_DapaHF_.5_censSC30,b=b1)
Out_log_EU_100_DapaHF_.5_censSC30<-Sim_out_log(b_est=Est_EU_100_DapaHF_.5_censSC30,b=b1)
Out_log_EU_50_DapaHF_.5_censSC30<-Sim_out_log(b_est=Est_EU_50_DapaHF_.5_censSC30,b=b1)

#50
Out_log_EU_500_DapaHF_.5_censSC50<-Sim_out_log(b_est=Est_EU_500_DapaHF_.5_censSC50,b=b1)
Out_log_EU_100_DapaHF_.5_censSC50<-Sim_out_log(b_est=Est_EU_100_DapaHF_.5_censSC50,b=b1)
Out_log_EU_50_DapaHF_.5_censSC50<-Sim_out_log(b_est=Est_EU_50_DapaHF_.5_censSC50,b=b1)

#70
Out_log_EU_500_DapaHF_.5_censSC70<-Sim_out_log(b_est=Est_EU_500_DapaHF_.5_censSC70,b=b1)
Out_log_EU_100_DapaHF_.5_censSC70<-Sim_out_log(b_est=Est_EU_100_DapaHF_.5_censSC70,b=b1)
Out_log_EU_50_DapaHF_.5_censSC70<-Sim_out_log(b_est=Est_EU_50_DapaHF_.5_censSC70_nna,b=b1)


#b=-0.5
b1<--0.5
#30
Out_log_EU_500_DapaHF_n.5_censSC30<-Sim_out_log(b_est=Est_EU_500_DapaHF_n.5_censSC30,b=b1)
Out_log_EU_100_DapaHF_n.5_censSC30<-Sim_out_log(b_est=Est_EU_100_DapaHF_n.5_censSC30,b=b1)
Out_log_EU_50_DapaHF_n.5_censSC30<-Sim_out_log(b_est=Est_EU_50_DapaHF_n.5_censSC30,b=b1)

#50
Out_log_EU_500_DapaHF_n.5_censSC50<-Sim_out_log(b_est=Est_EU_500_DapaHF_n.5_censSC50,b=b1)
Out_log_EU_100_DapaHF_n.5_censSC50<-Sim_out_log(b_est=Est_EU_100_DapaHF_n.5_censSC50,b=b1)
Out_log_EU_50_DapaHF_n.5_censSC50<-Sim_out_log(b_est=Est_EU_50_DapaHF_n.5_censSC50,b=b1)

#70
Out_log_EU_500_DapaHF_n.5_censSC70<-Sim_out_log(b_est=Est_EU_500_DapaHF_n.5_censSC70,b=b1)
Out_log_EU_100_DapaHF_n.5_censSC70<-Sim_out_log(b_est=Est_EU_100_DapaHF_n.5_censSC70,b=b1)
Out_log_EU_50_DapaHF_n.5_censSC70<-Sim_out_log(b_est=Est_EU_50_DapaHF_n.5_censSC70_nna,b=b1)


#b=0.25
b1<-0.25
#30
Out_log_EU_500_DapaHF_.25_censSC30<-Sim_out_log(b_est=Est_EU_500_DapaHF_.25_censSC30,b=b1)
Out_log_EU_100_DapaHF_.25_censSC30<-Sim_out_log(b_est=Est_EU_100_DapaHF_.25_censSC30,b=b1)
Out_log_EU_50_DapaHF_.25_censSC30<-Sim_out_log(b_est=Est_EU_50_DapaHF_.25_censSC30,b=b1)

#50
Out_log_EU_500_DapaHF_.25_censSC50<-Sim_out_log(b_est=Est_EU_500_DapaHF_.25_censSC50,b=b1)
Out_log_EU_100_DapaHF_.25_censSC50<-Sim_out_log(b_est=Est_EU_100_DapaHF_.25_censSC50,b=b1)
Out_log_EU_50_DapaHF_.25_censSC50<-Sim_out_log(b_est=Est_EU_50_DapaHF_.25_censSC50,b=b1)

#70
Out_log_EU_500_DapaHF_.25_censSC70<-Sim_out_log(b_est=Est_EU_500_DapaHF_.25_censSC70,b=b1)
Out_log_EU_100_DapaHF_.25_censSC70<-Sim_out_log(b_est=Est_EU_100_DapaHF_.25_censSC70,b=b1)
Out_log_EU_50_DapaHF_.25_censSC70<-Sim_out_log(b_est=Est_EU_50_DapaHF_.25_censSC70_nna,b=b1)


#b=-0.25
b1<--0.25
#30
Out_log_EU_500_DapaHF_n.25_censSC30<-Sim_out_log(b_est=Est_EU_500_DapaHF_n.25_censSC30,b=b1)
Out_log_EU_100_DapaHF_n.25_censSC30<-Sim_out_log(b_est=Est_EU_100_DapaHF_n.25_censSC30,b=b1)
Out_log_EU_50_DapaHF_n.25_censSC30<-Sim_out_log(b_est=Est_EU_50_DapaHF_n.25_censSC30,b=b1)

#50
Out_log_EU_500_DapaHF_n.25_censSC50<-Sim_out_log(b_est=Est_EU_500_DapaHF_n.25_censSC50,b=b1)
Out_log_EU_100_DapaHF_n.25_censSC50<-Sim_out_log(b_est=Est_EU_100_DapaHF_n.25_censSC50,b=b1)
Out_log_EU_50_DapaHF_n.25_censSC50<-Sim_out_log(b_est=Est_EU_50_DapaHF_n.25_censSC50,b=b1)

#70
Out_log_EU_500_DapaHF_n.25_censSC70<-Sim_out_log(b_est=Est_EU_500_DapaHF_n.25_censSC70,b=b1)
Out_log_EU_100_DapaHF_n.25_censSC70<-Sim_out_log(b_est=Est_EU_100_DapaHF_n.25_censSC70,b=b1)
Out_log_EU_50_DapaHF_n.25_censSC70<-Sim_out_log(b_est=Est_EU_50_DapaHF_n.25_censSC70_nna,b=b1)



#PH data
#b=0
b1<-0
#30
Out_log_PH_500_DapaHF_0_censSC30<-Sim_out_log(b_est=Est_PH_500_DapaHF_0_censSC30,b=b1)
Out_log_PH_100_DapaHF_0_censSC30<-Sim_out_log(b_est=Est_PH_100_DapaHF_0_censSC30,b=b1)
Out_log_PH_50_DapaHF_0_censSC30<-Sim_out_log(b_est=Est_PH_50_DapaHF_0_censSC30,b=b1)

#50
Out_log_PH_500_DapaHF_0_censSC50<-Sim_out_log(b_est=Est_PH_500_DapaHF_0_censSC50,b=b1)
Out_log_PH_100_DapaHF_0_censSC50<-Sim_out_log(b_est=Est_PH_100_DapaHF_0_censSC50,b=b1)
Out_log_PH_50_DapaHF_0_censSC50<-Sim_out_log(b_est=Est_PH_50_DapaHF_0_censSC50,b=b1)

#70
Out_log_PH_500_DapaHF_0_censSC70<-Sim_out_log(b_est=Est_PH_500_DapaHF_0_censSC70,b=b1)
Out_log_PH_100_DapaHF_0_censSC70<-Sim_out_log(b_est=Est_PH_100_DapaHF_0_censSC70,b=b1)
Out_log_PH_50_DapaHF_0_censSC70<-Sim_out_log(b_est=Est_PH_50_DapaHF_0_censSC70_nna,b=b1)


#b=0.5
b1<-0.5
#30
Out_log_PH_500_DapaHF_.5_censSC30<-Sim_out_log(b_est=Est_PH_500_DapaHF_.5_censSC30,b=b1)
Out_log_PH_100_DapaHF_.5_censSC30<-Sim_out_log(b_est=Est_PH_100_DapaHF_.5_censSC30,b=b1)
Out_log_PH_50_DapaHF_.5_censSC30<-Sim_out_log(b_est=Est_PH_50_DapaHF_.5_censSC30,b=b1)

#50
Out_log_PH_500_DapaHF_.5_censSC50<-Sim_out_log(b_est=Est_PH_500_DapaHF_.5_censSC50,b=b1)
Out_log_PH_100_DapaHF_.5_censSC50<-Sim_out_log(b_est=Est_PH_100_DapaHF_.5_censSC50,b=b1)
Out_log_PH_50_DapaHF_.5_censSC50<-Sim_out_log(b_est=Est_PH_50_DapaHF_.5_censSC50,b=b1)

#70
Out_log_PH_500_DapaHF_.5_censSC70<-Sim_out_log(b_est=Est_PH_500_DapaHF_.5_censSC70,b=b1)
Out_log_PH_100_DapaHF_.5_censSC70<-Sim_out_log(b_est=Est_PH_100_DapaHF_.5_censSC70,b=b1)
Out_log_PH_50_DapaHF_.5_censSC70<-Sim_out_log(b_est=Est_PH_50_DapaHF_.5_censSC70_nna,b=b1)


#b=-0.5
b1<--0.5
#30
Out_log_PH_500_DapaHF_n.5_censSC30<-Sim_out_log(b_est=Est_PH_500_DapaHF_n.5_censSC30,b=b1)
Out_log_PH_100_DapaHF_n.5_censSC30<-Sim_out_log(b_est=Est_PH_100_DapaHF_n.5_censSC30,b=b1)
Out_log_PH_50_DapaHF_n.5_censSC30<-Sim_out_log(b_est=Est_PH_50_DapaHF_n.5_censSC30,b=b1)

#50
Out_log_PH_500_DapaHF_n.5_censSC50<-Sim_out_log(b_est=Est_PH_500_DapaHF_n.5_censSC50,b=b1)
Out_log_PH_100_DapaHF_n.5_censSC50<-Sim_out_log(b_est=Est_PH_100_DapaHF_n.5_censSC50,b=b1)
Out_log_PH_50_DapaHF_n.5_censSC50<-Sim_out_log(b_est=Est_PH_50_DapaHF_n.5_censSC50,b=b1)

#70
Out_log_PH_500_DapaHF_n.5_censSC70<-Sim_out_log(b_est=Est_PH_500_DapaHF_n.5_censSC70,b=b1)
Out_log_PH_100_DapaHF_n.5_censSC70<-Sim_out_log(b_est=Est_PH_100_DapaHF_n.5_censSC70,b=b1)
Out_log_PH_50_DapaHF_n.5_censSC70<-Sim_out_log(b_est=Est_PH_50_DapaHF_n.5_censSC70_nna,b=b1)


#b=0.25
b1<-0.25
#30
Out_log_PH_500_DapaHF_.25_censSC30<-Sim_out_log(b_est=Est_PH_500_DapaHF_.25_censSC30,b=b1)
Out_log_PH_100_DapaHF_.25_censSC30<-Sim_out_log(b_est=Est_PH_100_DapaHF_.25_censSC30,b=b1)
Out_log_PH_50_DapaHF_.25_censSC30<-Sim_out_log(b_est=Est_PH_50_DapaHF_.25_censSC30,b=b1)

#50
Out_log_PH_500_DapaHF_.25_censSC50<-Sim_out_log(b_est=Est_PH_500_DapaHF_.25_censSC50,b=b1)
Out_log_PH_100_DapaHF_.25_censSC50<-Sim_out_log(b_est=Est_PH_100_DapaHF_.25_censSC50,b=b1)
Out_log_PH_50_DapaHF_.25_censSC50<-Sim_out_log(b_est=Est_PH_50_DapaHF_.25_censSC50,b=b1)

#70
Out_log_PH_500_DapaHF_.25_censSC70<-Sim_out_log(b_est=Est_PH_500_DapaHF_.25_censSC70,b=b1)
Out_log_PH_100_DapaHF_.25_censSC70<-Sim_out_log(b_est=Est_PH_100_DapaHF_.25_censSC70,b=b1)
Out_log_PH_50_DapaHF_.25_censSC70<-Sim_out_log(b_est=Est_PH_50_DapaHF_.25_censSC70_nna,b=b1)


#b=-0.25
b1<--0.25
#30
Out_log_PH_500_DapaHF_n.25_censSC30<-Sim_out_log(b_est=Est_PH_500_DapaHF_n.25_censSC30,b=b1)
Out_log_PH_100_DapaHF_n.25_censSC30<-Sim_out_log(b_est=Est_PH_100_DapaHF_n.25_censSC30,b=b1)
Out_log_PH_50_DapaHF_n.25_censSC30<-Sim_out_log(b_est=Est_PH_50_DapaHF_n.25_censSC30,b=b1)

#50
Out_log_PH_500_DapaHF_n.25_censSC50<-Sim_out_log(b_est=Est_PH_500_DapaHF_n.25_censSC50,b=b1)
Out_log_PH_100_DapaHF_n.25_censSC50<-Sim_out_log(b_est=Est_PH_100_DapaHF_n.25_censSC50,b=b1)
Out_log_PH_50_DapaHF_n.25_censSC50<-Sim_out_log(b_est=Est_PH_50_DapaHF_n.25_censSC50,b=b1)

#70
Out_log_PH_500_DapaHF_n.25_censSC70<-Sim_out_log(b_est=Est_PH_500_DapaHF_n.25_censSC70,b=b1)
Out_log_PH_100_DapaHF_n.25_censSC70<-Sim_out_log(b_est=Est_PH_100_DapaHF_n.25_censSC70,b=b1)
Out_log_PH_50_DapaHF_n.25_censSC70<-Sim_out_log(b_est=Est_PH_50_DapaHF_n.25_censSC70_nna,b=b1)



#Putting results together in one list
Results_PR.for.EU.log<-list(
  Out_log_EU_500_DapaHF_0_censSC30,Out_log_EU_100_DapaHF_0_censSC30,Out_log_EU_50_DapaHF_0_censSC30,
  Out_log_EU_500_DapaHF_0_censSC50,Out_log_EU_100_DapaHF_0_censSC50,Out_log_EU_50_DapaHF_0_censSC50,
  Out_log_EU_500_DapaHF_0_censSC70,Out_log_EU_100_DapaHF_0_censSC70,Out_log_EU_50_DapaHF_0_censSC70,
  
  Out_log_EU_500_DapaHF_.5_censSC30,Out_log_EU_100_DapaHF_.5_censSC30,Out_log_EU_50_DapaHF_.5_censSC30,
  Out_log_EU_500_DapaHF_.5_censSC50,Out_log_EU_100_DapaHF_.5_censSC50,Out_log_EU_50_DapaHF_.5_censSC50,
  Out_log_EU_500_DapaHF_.5_censSC70,Out_log_EU_100_DapaHF_.5_censSC70,Out_log_EU_50_DapaHF_.5_censSC70,
  
  Out_log_EU_500_DapaHF_.25_censSC30,Out_log_EU_100_DapaHF_.25_censSC30,Out_log_EU_50_DapaHF_.25_censSC30,
  Out_log_EU_500_DapaHF_.25_censSC50,Out_log_EU_100_DapaHF_.25_censSC50,Out_log_EU_50_DapaHF_.25_censSC50,
  Out_log_EU_500_DapaHF_.25_censSC70,Out_log_EU_100_DapaHF_.25_censSC70,Out_log_EU_50_DapaHF_.25_censSC70,
  
  Out_log_EU_500_DapaHF_n.25_censSC30,Out_log_EU_100_DapaHF_n.25_censSC30,Out_log_EU_50_DapaHF_n.25_censSC30,
  Out_log_EU_500_DapaHF_n.25_censSC50,Out_log_EU_100_DapaHF_n.25_censSC50,Out_log_EU_50_DapaHF_n.25_censSC50,
  Out_log_EU_500_DapaHF_n.25_censSC70,Out_log_EU_100_DapaHF_n.25_censSC70,Out_log_EU_50_DapaHF_n.25_censSC70,
  
  Out_log_EU_500_DapaHF_n.5_censSC30,Out_log_EU_100_DapaHF_n.5_censSC30,Out_log_EU_50_DapaHF_n.5_censSC30,
  Out_log_EU_500_DapaHF_n.5_censSC50,Out_log_EU_100_DapaHF_n.5_censSC50,Out_log_EU_50_DapaHF_n.5_censSC50,
  Out_log_EU_500_DapaHF_n.5_censSC70,Out_log_EU_100_DapaHF_n.5_censSC70,Out_log_EU_50_DapaHF_n.5_censSC70
)

Results_PR.for.PH.log<-list(
  Out_log_PH_500_DapaHF_0_censSC30,Out_log_PH_100_DapaHF_0_censSC30,Out_log_PH_50_DapaHF_0_censSC30,
  Out_log_PH_500_DapaHF_0_censSC50,Out_log_PH_100_DapaHF_0_censSC50,Out_log_PH_50_DapaHF_0_censSC50,
  Out_log_PH_500_DapaHF_0_censSC70,Out_log_PH_100_DapaHF_0_censSC70,Out_log_PH_50_DapaHF_0_censSC70,
  
  Out_log_PH_500_DapaHF_.5_censSC30,Out_log_PH_100_DapaHF_.5_censSC30,Out_log_PH_50_DapaHF_.5_censSC30,
  Out_log_PH_500_DapaHF_.5_censSC50,Out_log_PH_100_DapaHF_.5_censSC50,Out_log_PH_50_DapaHF_.5_censSC50,
  Out_log_PH_500_DapaHF_.5_censSC70,Out_log_PH_100_DapaHF_.5_censSC70,Out_log_PH_50_DapaHF_.5_censSC70,
  
  Out_log_PH_500_DapaHF_.25_censSC30,Out_log_PH_100_DapaHF_.25_censSC30,Out_log_PH_50_DapaHF_.25_censSC30,
  Out_log_PH_500_DapaHF_.25_censSC50,Out_log_PH_100_DapaHF_.25_censSC50,Out_log_PH_50_DapaHF_.25_censSC50,
  Out_log_PH_500_DapaHF_.25_censSC70,Out_log_PH_100_DapaHF_.25_censSC70,Out_log_PH_50_DapaHF_.25_censSC70,
  
  Out_log_PH_500_DapaHF_n.25_censSC30,Out_log_PH_100_DapaHF_n.25_censSC30,Out_log_PH_50_DapaHF_n.25_censSC30,
  Out_log_PH_500_DapaHF_n.25_censSC50,Out_log_PH_100_DapaHF_n.25_censSC50,Out_log_PH_50_DapaHF_n.25_censSC50,
  Out_log_PH_500_DapaHF_n.25_censSC70,Out_log_PH_100_DapaHF_n.25_censSC70,Out_log_PH_50_DapaHF_n.25_censSC70,
  
  Out_log_PH_500_DapaHF_n.5_censSC30,Out_log_PH_100_DapaHF_n.5_censSC30,Out_log_PH_50_DapaHF_n.5_censSC30,
  Out_log_PH_500_DapaHF_n.5_censSC50,Out_log_PH_100_DapaHF_n.5_censSC50,Out_log_PH_50_DapaHF_n.5_censSC50,
  Out_log_PH_500_DapaHF_n.5_censSC70,Out_log_PH_100_DapaHF_n.5_censSC70,Out_log_PH_50_DapaHF_n.5_censSC70
)

#combining results (mean, bias, MSE) to one data.frame
R_PR.for.EU_log<-results_log(Results_PR.for.EU.log)
R_PR.for.PH_log<-results_log(Results_PR.for.PH.log)




####Output PPR model####
#b=0                                                                                                       #true value
b1<-0
#30                                                                                                        #censoring rate
EU_ASV_Out_log_EU_500_DapaHF_0_censSC30<-Sim_EU_out_log(ests=EU_ASV_Est_EU_500_DapaHF_0_censSC30_nna,b=b1) #data.frame including mean, bias and MSE of the given case
EU_ASV_Out_log_EU_100_DapaHF_0_censSC30<-Sim_EU_out_log(ests=EU_ASV_Est_EU_100_DapaHF_0_censSC30_nna,b=b1)
EU_ASV_Out_log_EU_50_DapaHF_0_censSC30<-Sim_EU_out_log(ests=EU_ASV_Est_EU_50_DapaHF_0_censSC30_nna,b=b1)

#50
EU_ASV_Out_log_EU_500_DapaHF_0_censSC50<-Sim_EU_out_log(ests=EU_ASV_Est_EU_500_DapaHF_0_censSC50_nna,b=b1)
EU_ASV_Out_log_EU_100_DapaHF_0_censSC50<-Sim_EU_out_log(ests=EU_ASV_Est_EU_100_DapaHF_0_censSC50_nna,b=b1)
EU_ASV_Out_log_EU_50_DapaHF_0_censSC50<-Sim_EU_out_log(ests=EU_ASV_Est_EU_50_DapaHF_0_censSC50_nna,b=b1)

#70
EU_ASV_Out_log_EU_500_DapaHF_0_censSC70<-Sim_EU_out_log(ests=EU_ASV_Est_EU_500_DapaHF_0_censSC70_nna,b=b1)
EU_ASV_Out_log_EU_100_DapaHF_0_censSC70<-Sim_EU_out_log(ests=EU_ASV_Est_EU_100_DapaHF_0_censSC70_nna,b=b1)
EU_ASV_Out_log_EU_50_DapaHF_0_censSC70<-Sim_EU_out_log(ests=EU_ASV_Est_EU_50_DapaHF_0_censSC70_nna,b=b1)


#b=0.5
b1<-0.5
#30
EU_ASV_Out_log_EU_500_DapaHF_.5_censSC30<-Sim_EU_out_log(ests=EU_ASV_Est_EU_500_DapaHF_.5_censSC30_nna,b=b1)
EU_ASV_Out_log_EU_100_DapaHF_.5_censSC30<-Sim_EU_out_log(ests=EU_ASV_Est_EU_100_DapaHF_.5_censSC30_nna,b=b1)
EU_ASV_Out_log_EU_50_DapaHF_.5_censSC30<-Sim_EU_out_log(ests=EU_ASV_Est_EU_50_DapaHF_.5_censSC30_nna,b=b1)

#50
EU_ASV_Out_log_EU_500_DapaHF_.5_censSC50<-Sim_EU_out_log(ests=EU_ASV_Est_EU_500_DapaHF_.5_censSC50_nna,b=b1)
EU_ASV_Out_log_EU_100_DapaHF_.5_censSC50<-Sim_EU_out_log(ests=EU_ASV_Est_EU_100_DapaHF_.5_censSC50_nna,b=b1)
EU_ASV_Out_log_EU_50_DapaHF_.5_censSC50<-Sim_EU_out_log(ests=EU_ASV_Est_EU_50_DapaHF_.5_censSC50_nna,b=b1)

#70
EU_ASV_Out_log_EU_500_DapaHF_.5_censSC70<-Sim_EU_out_log(ests=EU_ASV_Est_EU_500_DapaHF_.5_censSC70_nna,b=b1)
EU_ASV_Out_log_EU_100_DapaHF_.5_censSC70<-Sim_EU_out_log(ests=EU_ASV_Est_EU_100_DapaHF_.5_censSC70_nna,b=b1)
EU_ASV_Out_log_EU_50_DapaHF_.5_censSC70<-Sim_EU_out_log(ests=EU_ASV_Est_EU_50_DapaHF_.5_censSC70_nna,b=b1)


#b=-0.5
b1<--0.5
#30
EU_ASV_Out_log_EU_500_DapaHF_n.5_censSC30<-Sim_EU_out_log(ests=EU_ASV_Est_EU_500_DapaHF_n.5_censSC30_nna,b=b1)
EU_ASV_Out_log_EU_100_DapaHF_n.5_censSC30<-Sim_EU_out_log(ests=EU_ASV_Est_EU_100_DapaHF_n.5_censSC30_nna,b=b1)
EU_ASV_Out_log_EU_50_DapaHF_n.5_censSC30<-Sim_EU_out_log(ests=EU_ASV_Est_EU_50_DapaHF_n.5_censSC30_nna,b=b1)

#50
EU_ASV_Out_log_EU_500_DapaHF_n.5_censSC50<-Sim_EU_out_log(ests=EU_ASV_Est_EU_500_DapaHF_n.5_censSC50_nna,b=b1)
EU_ASV_Out_log_EU_100_DapaHF_n.5_censSC50<-Sim_EU_out_log(ests=EU_ASV_Est_EU_100_DapaHF_n.5_censSC50_nna,b=b1)
EU_ASV_Out_log_EU_50_DapaHF_n.5_censSC50<-Sim_EU_out_log(ests=EU_ASV_Est_EU_50_DapaHF_n.5_censSC50_nna,b=b1)

#70
EU_ASV_Out_log_EU_500_DapaHF_n.5_censSC70<-Sim_EU_out_log(ests=EU_ASV_Est_EU_500_DapaHF_n.5_censSC70_nna,b=b1)
EU_ASV_Out_log_EU_100_DapaHF_n.5_censSC70<-Sim_EU_out_log(ests=EU_ASV_Est_EU_100_DapaHF_n.5_censSC70_nna,b=b1)
EU_ASV_Out_log_EU_50_DapaHF_n.5_censSC70<-Sim_EU_out_log(ests=EU_ASV_Est_EU_50_DapaHF_n.5_censSC70_nna,b=b1)


#b=0.25
b1<-0.25
#30
EU_ASV_Out_log_EU_500_DapaHF_.25_censSC30<-Sim_EU_out_log(ests=EU_ASV_Est_EU_500_DapaHF_.25_censSC30_nna,b=b1)
EU_ASV_Out_log_EU_100_DapaHF_.25_censSC30<-Sim_EU_out_log(ests=EU_ASV_Est_EU_100_DapaHF_.25_censSC30_nna,b=b1)
EU_ASV_Out_log_EU_50_DapaHF_.25_censSC30<-Sim_EU_out_log(ests=EU_ASV_Est_EU_50_DapaHF_.25_censSC30_nna,b=b1)

#50
EU_ASV_Out_log_EU_500_DapaHF_.25_censSC50<-Sim_EU_out_log(ests=EU_ASV_Est_EU_500_DapaHF_.25_censSC50_nna,b=b1)
EU_ASV_Out_log_EU_100_DapaHF_.25_censSC50<-Sim_EU_out_log(ests=EU_ASV_Est_EU_100_DapaHF_.25_censSC50_nna,b=b1)
EU_ASV_Out_log_EU_50_DapaHF_.25_censSC50<-Sim_EU_out_log(ests=EU_ASV_Est_EU_50_DapaHF_.25_censSC50_nna,b=b1)

#70
EU_ASV_Out_log_EU_500_DapaHF_.25_censSC70<-Sim_EU_out_log(ests=EU_ASV_Est_EU_500_DapaHF_.25_censSC70_nna,b=b1)
EU_ASV_Out_log_EU_100_DapaHF_.25_censSC70<-Sim_EU_out_log(ests=EU_ASV_Est_EU_100_DapaHF_.25_censSC70_nna,b=b1)
EU_ASV_Out_log_EU_50_DapaHF_.25_censSC70<-Sim_EU_out_log(ests=EU_ASV_Est_EU_50_DapaHF_.25_censSC70_nna,b=b1)


#b=-0.25
b1<--0.25
#30
EU_ASV_Out_log_EU_500_DapaHF_n.25_censSC30<-Sim_EU_out_log(ests=EU_ASV_Est_EU_500_DapaHF_n.25_censSC30_nna,b=b1)
EU_ASV_Out_log_EU_100_DapaHF_n.25_censSC30<-Sim_EU_out_log(ests=EU_ASV_Est_EU_100_DapaHF_n.25_censSC30_nna,b=b1)
EU_ASV_Out_log_EU_50_DapaHF_n.25_censSC30<-Sim_EU_out_log(ests=EU_ASV_Est_EU_50_DapaHF_n.25_censSC30_nna,b=b1)

#50
EU_ASV_Out_log_EU_500_DapaHF_n.25_censSC50<-Sim_EU_out_log(ests=EU_ASV_Est_EU_500_DapaHF_n.25_censSC50_nna,b=b1)
EU_ASV_Out_log_EU_100_DapaHF_n.25_censSC50<-Sim_EU_out_log(ests=EU_ASV_Est_EU_100_DapaHF_n.25_censSC50_nna,b=b1)
EU_ASV_Out_log_EU_50_DapaHF_n.25_censSC50<-Sim_EU_out_log(ests=EU_ASV_Est_EU_50_DapaHF_n.25_censSC50_nna,b=b1)

#70
EU_ASV_Out_log_EU_500_DapaHF_n.25_censSC70<-Sim_EU_out_log(ests=EU_ASV_Est_EU_500_DapaHF_n.25_censSC70_nna,b=b1)
EU_ASV_Out_log_EU_100_DapaHF_n.25_censSC70<-Sim_EU_out_log(ests=EU_ASV_Est_EU_100_DapaHF_n.25_censSC70_nna,b=b1)
EU_ASV_Out_log_EU_50_DapaHF_n.25_censSC70<-Sim_EU_out_log(ests=EU_ASV_Est_EU_50_DapaHF_n.25_censSC70_nna,b=b1)



#Putting results together in one list
Results_ASV_EU.for.EU_log<-list(
  EU_ASV_Out_log_EU_500_DapaHF_0_censSC30,EU_ASV_Out_log_EU_100_DapaHF_0_censSC30,EU_ASV_Out_log_EU_50_DapaHF_0_censSC30,
  EU_ASV_Out_log_EU_500_DapaHF_0_censSC50,EU_ASV_Out_log_EU_100_DapaHF_0_censSC50,EU_ASV_Out_log_EU_50_DapaHF_0_censSC50,
  EU_ASV_Out_log_EU_500_DapaHF_0_censSC70,EU_ASV_Out_log_EU_100_DapaHF_0_censSC70,EU_ASV_Out_log_EU_50_DapaHF_0_censSC70,
  
  EU_ASV_Out_log_EU_500_DapaHF_.5_censSC30,EU_ASV_Out_log_EU_100_DapaHF_.5_censSC30,EU_ASV_Out_log_EU_50_DapaHF_.5_censSC30,
  EU_ASV_Out_log_EU_500_DapaHF_.5_censSC50,EU_ASV_Out_log_EU_100_DapaHF_.5_censSC50,EU_ASV_Out_log_EU_50_DapaHF_.5_censSC50,
  EU_ASV_Out_log_EU_500_DapaHF_.5_censSC70,EU_ASV_Out_log_EU_100_DapaHF_.5_censSC70,EU_ASV_Out_log_EU_50_DapaHF_.5_censSC70,
  
  EU_ASV_Out_log_EU_500_DapaHF_.25_censSC30,EU_ASV_Out_log_EU_100_DapaHF_.25_censSC30,EU_ASV_Out_log_EU_50_DapaHF_.25_censSC30,
  EU_ASV_Out_log_EU_500_DapaHF_.25_censSC50,EU_ASV_Out_log_EU_100_DapaHF_.25_censSC50,EU_ASV_Out_log_EU_50_DapaHF_.25_censSC50,
  EU_ASV_Out_log_EU_500_DapaHF_.25_censSC70,EU_ASV_Out_log_EU_100_DapaHF_.25_censSC70,EU_ASV_Out_log_EU_50_DapaHF_.25_censSC70,
  
  EU_ASV_Out_log_EU_500_DapaHF_n.25_censSC30,EU_ASV_Out_log_EU_100_DapaHF_n.25_censSC30,EU_ASV_Out_log_EU_50_DapaHF_n.25_censSC30,
  EU_ASV_Out_log_EU_500_DapaHF_n.25_censSC50,EU_ASV_Out_log_EU_100_DapaHF_n.25_censSC50,EU_ASV_Out_log_EU_50_DapaHF_n.25_censSC50,
  EU_ASV_Out_log_EU_500_DapaHF_n.25_censSC70,EU_ASV_Out_log_EU_100_DapaHF_n.25_censSC70,EU_ASV_Out_log_EU_50_DapaHF_n.25_censSC70,
  
  EU_ASV_Out_log_EU_500_DapaHF_n.5_censSC30,EU_ASV_Out_log_EU_100_DapaHF_n.5_censSC30,EU_ASV_Out_log_EU_50_DapaHF_n.5_censSC30,
  EU_ASV_Out_log_EU_500_DapaHF_n.5_censSC50,EU_ASV_Out_log_EU_100_DapaHF_n.5_censSC50,EU_ASV_Out_log_EU_50_DapaHF_n.5_censSC50,
  EU_ASV_Out_log_EU_500_DapaHF_n.5_censSC70,EU_ASV_Out_log_EU_100_DapaHF_n.5_censSC70,EU_ASV_Out_log_EU_50_DapaHF_n.5_censSC70  )

#combining results (mean, bias, MSE) to one data.frame
R_ASV_EU.for.EU_log<-results_log(Results_ASV_EU.for.EU_log)




####Small simulated study for the explanatory chart####
set.seed(23)

###Preparation 
#Censoring parameters
plotcens<-c(0,250)

#test for censoring rate
simulEU_SC_censrate(N=10000,param=ML_EU$par,rateC=plotcens,prob=0.5,b=0.5)

#simulate study
Plot_dataset<-Sim_EUSC(m=1,N=50,param=ML_EU$par,rateC =rateC_EU_DapaHF_0_censSC30,prob=0.5,b=0.75)

#extract simulated study from list (of 1)
Plot_dataset<-Plot_dataset[[1]]

#stratify into treatment groups
Plot_dataset_1<-subset(Plot_dataset,Plot_dataset$expressions==1)
Plot_dataset_0<-subset(Plot_dataset,Plot_dataset$expressions==0)

#extract the events
Plot_dataset_event_1<-subset(Plot_dataset_1$time,Plot_dataset_1$status==1)
Plot_dataset_event_0<-subset(Plot_dataset_0$time,Plot_dataset_0$status==1)
Plot_dataset_event<-c(Plot_dataset_event_1,Plot_dataset_event_0)         

#extract the censored observations
Plot_dataset_cens_1<-subset(Plot_dataset_1$time,Plot_dataset_1$status==0)
Plot_dataset_cens_0<-subset(Plot_dataset_0$time,Plot_dataset_0$status==0)
Plot_dataset_cens<-c(Plot_dataset_cens_1,Plot_dataset_cens_0)

#restrict to \tilde{T}
Plot_dataset_event<-subset(Plot_dataset_event,Plot_dataset_event>=max(min(Plot_dataset_event_1),min(Plot_dataset_event_0))&Plot_dataset_event<=min(max(Plot_dataset_event_1),max(Plot_dataset_event_0)))
min(max(Plot_dataset_1$time),max(Plot_dataset_0$time))


#explanatory chart on the RR scale
#KM-estimated CDF treatment
F_plot_1<-function(t){
  F_est(t,Plot_dataset_1,col = c(col_time=2,col_cens=3))
}
#control
F_plot_0<-function(t){
  F_est(t,Plot_dataset_0,col = c(col_time=2,col_cens=3))
}

#x-axis
X0_plot<-seq(0,max(Plot_dataset_0$time),by=0.001)
X1_plot<-seq(0,max(Plot_dataset_1$time),by=0.001)

#function RR
RR_plot<-function(t){
  F_plot_1(t)/F_plot_0(t)
}
#RR at every event time point in \tilde{T}
RR_ypoints<-as.numeric(lapply(Plot_dataset_event,RR_plot))

#SE functions
#treatment
W_plot_1<-function(t){
  KM_se(t,Plot_dataset_1,col=c(col_time=2,col_cens=3))
}
#control
W_plot_0<-function(t){
  KM_se(t,Plot_dataset_0,col=c(col_time=2,col_cens=3))
}

#weight function
Est_weights_plot<-function(t){
  1/(1/((1/F_plot_1(t)^2)*W_plot_1(t)^2+(1/F_plot_0(t)^2)*W_plot_0(t)^2))
}
#Weight sizes
weights_plot<-as.numeric(lapply(Plot_dataset_event,Est_weights_plot))
#ordered weight sizes
weights_plot_order<-order(weights_plot)

#plot on RR scale
dev.new(width=6, height=5, unit="cm")
par(mar=c(5,6,4,3)+.1)
plot(X0_plot,lapply(X0_plot,F_plot_0),col="red",type = "l",lty=1,xlab ="Time",ylab=c("Event probability","and RR"),xlim=c(0,100),ylim=c(0,1.2), cex.lab=1.5,ann=T,xaxt='n',yaxt="n")   #KM-estimated control
points(Plot_dataset_cens_0,lapply(Plot_dataset_cens_0,F_plot_0),col="red",pch=4)                                                                                                       #indicate censored observations
lines(X1_plot,lapply(X1_plot,F_plot_1),col="blue",type = "l",lty=1)                                                                                                                    #KM-estimated CDF treatment
points(Plot_dataset_cens_1,lapply(Plot_dataset_cens_1,F_plot_1),col="blue",pch=4)                                                                                                      #indicate censored observations
abline(v=max(min(Plot_dataset_event_1),min(Plot_dataset_event_0)),lty=3)                                                                                                               #indicating t_min
text(max(min(Plot_dataset_event_1),min(Plot_dataset_event_0))+3,-4,expression(tilde(t)[min]),cex=1)
abline(v=min(max(Plot_dataset_event_1),max(Plot_dataset_event_0)),lty=3)                                                                                                               #indicating t_max
text(min(max(Plot_dataset_event_1),max(Plot_dataset_event_0))
     +3,-4,expression(tilde(t)[max]),cex=1)
for(i in 1:length(Plot_dataset_event)){                                                                                                                                                #add \hat{beta}_t with size corresponding to weight
  x<-Plot_dataset_event[i]
  y<-RR_ypoints[i]
  p<-weights_plot_order[i]*0.1
  points(x,y,pch=3,cex=p,col="lightgreen")
}
text(x=90,y=0.85,expression(hat(F)[0]),col="red")
text(x=90,y=0.25,expression(hat(F)[1]),col="blue")
abline(h=exp(-0.7),lty=4,col="lightgreen")                                                                                                                                             #indicating true beta
text(x=95,y=exp(-0.7)+0.03,expression(paste(exp(-beta))),cex=1,col="lightgreen")
axis(1,line=0,at=Plot_dataset$time,col.ticks="black",col.axis="white",cex.axis=0.54)                                                                                                   #axis indicating all event/censored time points
axis(2,line=0,at=c(0,0.2,0.4,0.6,0.8,1),col.axis="black",col.ticks = "black")
axis(4,line=0,col.ticks = "lightgreen",col.axis="lightgreen")
legend(x=20,y=1.2, c(expression(paste(hat(F)[1](t)/hat(F)[0](t))),"Censored observation (control)","Censored observation (treatment)"),
       col=c("lightgreen","red","blue"),pch = c(3,4,4), cex=.9)



#chart on log scale
#log CDF treatment (negative)
nlogF_plot_1<-function(t){
  -log(F_est(t,Plot_dataset_1,col = c(col_time=2,col_cens=3)))
}
#control
logF_plot_0<-function(t){
  log(F_est(t,Plot_dataset_0,col = c(col_time=2,col_cens=3)))
}

#negative log RR
nlogRR_plot<-function(t){
  -log(F_plot_1(t)/F_plot_0(t))
}
nlogRR_ypoints<-as.numeric(lapply(Plot_dataset_event,nlogRR_plot))

#Plot the same as above on log scale
dev.new(width=6, height=5, unit="cm")
par(mar=c(5,6,4,3)+.1)
plot(X_plot,lapply(X_plot,logF_plot_0),col="red",type = "l",lty=1,xlab ="Time",ylab=expression(paste(beta , " and log. event probability")),xlim=c(0,100),ylim=c(-5,5), cex.lab=1.3,ann=T,xaxt='n')
points(Plot_dataset_cens_0,lapply(Plot_dataset_cens_0,logF_plot_0),col="red",pch=4)
lines(X_plot,lapply(X_plot,nlogF_plot_1),col="blue",type = "l",lty=1)
points(Plot_dataset_cens_1,lapply(Plot_dataset_cens_1,nlogF_plot_1),col="blue",pch=4)
abline(v=max(min(Plot_dataset_event_1),min(Plot_dataset_event_0)),lty=3)
text(max(min(Plot_dataset_event_1),min(Plot_dataset_event_0))+3,-4,expression(tilde(t)[min]),cex=1)
abline(v=min(max(Plot_dataset_event_1),max(Plot_dataset_event_0)),lty=3)
text(min(max(Plot_dataset_event_1),max(Plot_dataset_event_0))
     +3,-4,expression(tilde(t)[max]),cex=1)
for(i in 1:length(Plot_dataset_event)){
  x<-Plot_dataset_event[i]
  y<-nlogRR_ypoints[i]
  p<-weights_plot_order[i]*0.1
  points(x,y,pch=3,cex=p,col="lightgreen")
}
abline(h=0)
text(x=95,y=-1.5,expression(log(hat(F)[0])),col="red")
text(x=95,y=2,expression(-log(hat(F)[1])),col="blue")
abline(h=0.7,lty=4,col="lightgreen")
text(x=95,y=0.5,expression(paste(beta)),cex=1,col="lightgreen")
axis(1,line=0,at=Plot_dataset$time,col.ticks="black",col.axis="white",cex.axis=0.54)
axis(4,line=0,at=c(-4,-2,0,2,4),col.ticks = "lightgreen",col.axis="lightgreen")
legend(x=10,y=-3, c(expression(paste(-log(hat(F)[1](t)/hat(F)[0](t))," = ",-log(hat(F)[1](t))+log(hat(F)[0](t)))),"Censored observation (control)","Censored observation (treatment)"),
       col=c("lightgreen","red","blue"),pch = c(3,4,4), cex=.9)


