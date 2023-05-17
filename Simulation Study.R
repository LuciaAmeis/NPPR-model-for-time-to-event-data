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
ML_0<-optim(par=c(1.2,70),fn=negloglike_Weib,daten=subset(DapaHF,DapaHF$Arm==0),col=c(col_time=3,col_cens=1),lower = c(0, 0))



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
ML_EU<-optim(par=c(0.8,130),negloglikeEU,daten=subset(DapaHF,DapaHF$Arm==0),col=c(col_time=3,col_cens=1),method="Nelder-Mead")




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
    survival<-rweibull(n=Nn,shape=a1,scale = a2*(exp(b*pers_expr)^(1/a1)))
    
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
    
    survival<-rweibull(n=Nn,shape=a1,scale = a2*(exp(b*pers_expr)^(1/a1)))
    
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
sum(is.na(Est_PH_50_DapaHF_0_censSC70))#6
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
sum(is.na(Est_PH_50_DapaHF_.5_censSC70))#9
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
sum(is.na(Est_PH_50_DapaHF_n.5_censSC70))#11
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
sum(is.na(Est_PH_50_DapaHF_.25_censSC70))#7
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
sum(is.na(Est_PH_50_DapaHF_n.25_censSC70))#7
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

