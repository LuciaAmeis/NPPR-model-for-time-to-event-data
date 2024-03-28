#########################################################################
#                                                                       #
#    A non-parametric proportional risk model to assess a treatment     #
#                    effect in time-to-event data                       #
#                                                                       #
#         L. Ameis, O. Kuss, A. Hoyer and K. Moellenhoff                #
#                                                                       #
#                             PLOTS                                     #
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
#load workspace of simulations
load("./simulation/NPPR_simulation.RData")
load("./case_study/NPPR_case_study.RData")

#############
####Plots####
#############
#Figure 1 (and A)
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

#plot on RR scale (Figure 1)
dev.new(width=6, height=5, unit="cm")
postscript(file="./plot_generation/plot/Figure1.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
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
  points(x,y,pch=3,cex=p,col=rgb(0,0.6,0))
}
text(x=90,y=0.85,expression(hat(F)[0]),col="red")
text(x=90,y=0.25,expression(hat(F)[1]),col="blue")
abline(h=exp(-0.7),lty=4,col=rgb(0,0.6,0))                                                                                                                                             #indicating true beta
text(x=95,y=exp(-0.7)+0.03,expression(paste(exp(-beta))),cex=1,col=rgb(0,0.6,0))
axis(1,line=0,at=Plot_dataset$time,col.ticks="black",col.axis="white",cex.axis=0.54)                                                                                                   #axis indicating all event/censored time points
axis(2,line=0,at=c(0,0.2,0.4,0.6,0.8,1),col.axis="black",col.ticks = "black")
axis(4,line=0,col.ticks = rgb(0,0.6,0),col.axis=rgb(0,0.6,0))
legend(x=20,y=1.2, c(expression(paste(hat(F)[1](t)/hat(F)[0](t))),"Censored observation (control)","Censored observation (treatment)"),
       col=c(rgb(0,0.6,0),"red","blue"),pch = c(3,4,4), cex=.9)
dev.off()


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

#Plot the same as above on log scale (Figure A)
dev.new(width=6, height=5, unit="cm")
postscript(file="./plot_generation/plot/FigureA.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar=c(5,6,4,3)+.1)
plot(X0_plot,lapply(X0_plot,logF_plot_0),col="red",type = "l",lty=1,xlab ="Time",ylab=expression(paste(beta , " and log. event probability")),xlim=c(0,100),ylim=c(-5,5), cex.lab=1.3,ann=T,xaxt='n')
points(Plot_dataset_cens_0,lapply(Plot_dataset_cens_0,logF_plot_0),col="red",pch=4)
lines(X0_plot,lapply(X0_plot,nlogF_plot_1),col="blue",type = "l",lty=1)
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
  points(x,y,pch=3,cex=p,col=rgb(0,0.6,0))
}
abline(h=0)
text(x=95,y=-1.5,expression(log(hat(F)[0])),col="red")
text(x=95,y=2,expression(-log(hat(F)[1])),col="blue")
abline(h=0.7,lty=4,col=rgb(0,0.6,0))
text(x=95,y=0.5,expression(paste(beta)),cex=1,col=rgb(0,0.6,0))
axis(1,line=0,at=Plot_dataset$time,col.ticks="black",col.axis="white",cex.axis=0.54)
axis(4,line=0,at=c(-4,-2,0,2,4),col.ticks = rgb(0,0.6,0),col.axis=rgb(0,0.6,0))
legend(x=10,y=-3, c(expression(paste(-log(hat(F)[1](t)/hat(F)[0](t))," = ",-log(hat(F)[1](t))+log(hat(F)[0](t)))),"Censored observation (control)","Censored observation (treatment)"),
       col=c(rgb(0,0.6,0),"red","blue"),pch = c(3,4,4), cex=.9)
dev.off()


#Figure 2 and C
F_weib<-function(t,param){
  pweibull(t,shape=param[1],scale=param[2],lower.tail = TRUE, log.p = FALSE)
}


F_weib0<-function(t){
  F_weib(t,param=ML_0$par)
}

F_weib.5<-function(t){
  F_weib(t,param=c(ML_0$par[1],ML_0$par[2]*exp(0.5)^(1/ML_0$par[1])))
}
F_weibn.5<-function(t){
  F_weib(t,param=c(ML_0$par[1],ML_0$par[2]*exp(-0.5)^(1/ML_0$par[1])))
}
F_weib.25<-function(t){
  F_weib(t,param=c(ML_0$par[1],ML_0$par[2]*exp(0.25)^(1/ML_0$par[1])))
}
F_weibn.25<-function(t){
  F_weib(t,param=c(ML_0$par[1],ML_0$par[2]*exp(-0.25)^(1/ML_0$par[1])))
}


RR0<-function(t){
  F_weib0(t)/F_weib0(t)
}

RR.25<-function(t){
  F_weib.25(t)/F_weib0(t)
}

RRn.25<-function(t){
  F_weibn.25(t)/F_weib0(t)
}

RR.5<-function(t){
  F_weib.5(t)/F_weib0(t)
}

RRn.5<-function(t){
  F_weibn.5(t)/F_weib0(t)
}

F_EU0<-function(t){
  pEU(t,ML_EU$par[1],0,ML_EU$par[2])
}

F_EU.5<-function(t){
  pEU(t,ML_EU$par[1],0,ML_EU$par[2]*exp(-0.5)^(-1/ML_EU$par[1]))
}

F_EUn.5<-function(t){
  pEU(t,ML_EU$par[1],0,ML_EU$par[2]*exp(0.5)^(-1/ML_EU$par[1]))
}
F_EU.25<-function(t){
  pEU(t,ML_EU$par[1],0,ML_EU$par[2]*exp(-0.25)^(-1/ML_EU$par[1]))
}

F_EUn.25<-function(t){
  pEU(t,ML_EU$par[1],0,ML_EU$par[2]*exp(0.25)^(-1/ML_EU$par[1]))
}


HR_EU0<-function(t){
  h1<-  dEU(t,ML_EU$par[1],0,ML_EU$par[2])/(1-  pEU(t,ML_EU$par[1],0,ML_EU$par[2]))
  h0<-  dEU(t,ML_EU$par[1],0,ML_EU$par[2])/(1-  pEU(t,ML_EU$par[1],0,ML_EU$par[2]))
  h<-h1/h0
  h
}

HR_EU.5<-function(t){
  h1<-  dEU(t,ML_EU$par[1],0,ML_EU$par[2]*exp(-0.5)^(-1/ML_EU$par[1]))/(1-  pEU(t,ML_EU$par[1],0,ML_EU$par[2]*exp(-0.5)^(-1/ML_EU$par[1])))
  h0<-  dEU(t,ML_EU$par[1],0,ML_EU$par[2])/(1-  pEU(t,ML_EU$par[1],0,ML_EU$par[2]))
  h<-h1/h0
  h
}

HR_EUn.5<-function(t){
  h1<-  dEU(t,ML_EU$par[1],0,ML_EU$par[2]*exp(0.5)^(-1/ML_EU$par[1]))/(1-  pEU(t,ML_EU$par[1],0,ML_EU$par[2]*exp(0.5)^(-1/ML_EU$par[1])))
  h0<-  dEU(t,ML_EU$par[1],0,ML_EU$par[2])/(1-  pEU(t,ML_EU$par[1],0,ML_EU$par[2]))
  h<-h1/h0
  h
}

HR_EU.25<-function(t){
  h1<-  dEU(t,ML_EU$par[1],0,ML_EU$par[2]*exp(-0.25)^(-1/ML_EU$par[1]))/(1-  pEU(t,ML_EU$par[1],0,ML_EU$par[2]*exp(-0.25)^(-1/ML_EU$par[1])))
  h0<-  dEU(t,ML_EU$par[1],0,ML_EU$par[2])/(1-  pEU(t,ML_EU$par[1],0,ML_EU$par[2]))
  h<-h1/h0
  h
}

HR_EUn.25<-function(t){
  h1<-  dEU(t,ML_EU$par[1],0,ML_EU$par[2]*exp(0.25)^(-1/ML_EU$par[1]))/(1-  pEU(t,ML_EU$par[1],0,ML_EU$par[2]*exp(0.25)^(-1/ML_EU$par[1])))
  h0<-  dEU(t,ML_EU$par[1],0,ML_EU$par[2])/(1-  pEU(t,ML_EU$par[1],0,ML_EU$par[2]))
  h<-h1/h0
  h
}

OD_EU0<-function(t){
  fo1<-      pEU(t,ML_EU$par[1],0,ML_EU$par[2])/(1-  pEU(t,ML_EU$par[1],0,ML_EU$par[2]))
  fo0<-      pEU(t,ML_EU$par[1],0,ML_EU$par[2])/(1-  pEU(t,ML_EU$par[1],0,ML_EU$par[2]))
  fo<-fo1/fo0
}
OD_EU.5<-function(t){
  fo1<-      pEU(t,ML_EU$par[1],0,ML_EU$par[2]*exp(-0.5)^(-1/ML_EU$par[1]))/(1-  pEU(t,ML_EU$par[1],0,ML_EU$par[2]*exp(-0.5)^(-1/ML_EU$par[1])))
  fo0<-      pEU(t,ML_EU$par[1],0,ML_EU$par[2])/(1-  pEU(t,ML_EU$par[1],0,ML_EU$par[2]))
  fo<-fo1/fo0
}
OD_EUn.5<-function(t){
  fo1<-      pEU(t,ML_EU$par[1],0,ML_EU$par[2]*exp(0.5)^(-1/ML_EU$par[1]))/(1-  pEU(t,ML_EU$par[1],0,ML_EU$par[2]*exp(0.5)^(-1/ML_EU$par[1])))
  fo0<-      pEU(t,ML_EU$par[1],0,ML_EU$par[2])/(1-  pEU(t,ML_EU$par[1],0,ML_EU$par[2]))
  fo<-fo1/fo0
}
OD_EU.25<-function(t){
  fo1<-      pEU(t,ML_EU$par[1],0,ML_EU$par[2]*exp(-0.25)^(-1/ML_EU$par[1]))/(1-  pEU(t,ML_EU$par[1],0,ML_EU$par[2]*exp(-0.25)^(-1/ML_EU$par[1])))
  fo0<-      pEU(t,ML_EU$par[1],0,ML_EU$par[2])/(1-  pEU(t,ML_EU$par[1],0,ML_EU$par[2]))
  fo<-fo1/fo0
}
OD_EUn.25<-function(t){
  fo1<-      pEU(t,ML_EU$par[1],0,ML_EU$par[2]*exp(0.25)^(-1/ML_EU$par[1]))/(1-  pEU(t,ML_EU$par[1],0,ML_EU$par[2]*exp(0.25)^(-1/ML_EU$par[1])))
  fo0<-      pEU(t,ML_EU$par[1],0,ML_EU$par[2])/(1-  pEU(t,ML_EU$par[1],0,ML_EU$par[2]))
  fo<-fo1/fo0
}

F_LL<-function(t,param){
  pllogis(t,shape=param[1],scale=param[2],lower.tail = TRUE, log.p = FALSE)
}


F_LL0<-function(t){
  F_LL(t,param=ML_LL$par)
}

F_LL.5<-function(t){
  F_LL(t,param=c(ML_LL$par[1],ML_LL$par[2]*exp(0.5)^(1/ML_LL$par[1])))
}
F_LLn.5<-function(t){
  F_LL(t,param=c(ML_0$par[1],ML_LL$par[2]*exp(-0.5)^(1/ML_LL$par[1])))
}
F_LL.25<-function(t){
  F_LL(t,param=c(ML_0$par[1],ML_LL$par[2]*exp(0.25)^(1/ML_LL$par[1])))
}
F_LLn.25<-function(t){
  F_LL(t,param=c(ML_0$par[1],ML_LL$par[2]*exp(-0.25)^(1/ML_LL$par[1])))
}


RRLL0<-function(t){
  F_LL0(t)/F_LL0(t)
}

RRLL.25<-function(t){
  F_LL.25(t)/F_LL0(t)
}

RRLLn.25<-function(t){
  F_LLn.25(t)/F_LL0(t)
}

RRLL.5<-function(t){
  F_LL.5(t)/F_LL0(t)
}

RRLLn.5<-function(t){
  F_LLn.5(t)/F_LL0(t)
}

#Figure 2
dev.new(width=5, height=5, unit="cm")
postscript(file="./plot_generation/plot/Figure2.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar=c(5,6,4,1)+.1)
plot(seq(0,350,by=5),lapply(seq(0,350,by=5),RRn.5),col="orange",type = "p",pch=4,xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5,ann=T,xaxt='n')
lines(seq(0,350,by=5),lapply(seq(0,350,by=5),RRLLn.5),col="orange",type = "p",pch=20,xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),rep(exp(--0.5),times=length(seq(0,350,by=5))),col="orange",type = "l",xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),lapply(seq(0,350,by=5),RR.5),col="lightblue",type = "p",pch=4,xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),rep(exp(-0.5),times=length(seq(0,350,by=5))),col="lightblue",type = "l",xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),lapply(seq(0,350,by=5),RRLL.5),col="lightblue",type = "p",pch=20,xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),lapply(seq(0,350,by=5),RR.25),col="blue",type = "p",pch=4,xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),rep(exp(-0.25),times=length(seq(0,350,by=5))),col="blue",type = "l",xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),lapply(seq(0,350,by=5),RRLL.25),col="blue",type = "p",pch=20,xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),lapply(seq(0,350,by=5),RRn.25),col="red",type = "p",pch=4,xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),rep(exp(--0.25),times=length(seq(0,350,by=5))),col="red",type = "l",xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),lapply(seq(0,350,by=5),RRLLn.25),col="red",type = "p",pch=20,xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),lapply(seq(0,350,by=5),RR0),col="black",type = "p",pch=4,xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),rep(exp(0),times=length(seq(0,350,by=5))),col="black",type = "l",xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),lapply(seq(0,350,by=5),RRLL0),col="black",type = "p",pch=20,xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
axis(1,at=c(0,50,100,150,200),labels=c("0","50","100","150","200"))
legend(x="topright",c("Weib PH", "LL PO", "HR/OR=1","HR/OR=0.607","HR/OR=0.779","HR/OR=1.284","HR/OR=1.649"),
       col=c("black","black","black","lightblue","blue","red","orange"), lty=c(0,0,0,0,0,0,0),pch=c(4,20,15,15,15,15,15), cex=0.9,bg="white")
dev.off()

#Figure D
dev.new(width=5, height=5, unit="cm")
postscript(file="./plot_generation/plot/FigureD.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar=c(5,6,4,1)+.1)
plot(seq(0,350,by=5),lapply(seq(0,350,by=5),HR_EUn.5),col="orange",type = "p",pch=4,xlab ="Time (months)",ylab="HR and OR",xlim=c(0,100),ylim=c(0,17), cex.lab=1.5,ann=T,xaxt='n')
lines(seq(0,350,by=5),lapply(seq(0,350,by=5),OD_EUn.5),col="orange",type = "p",pch=20,xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),rep(exp(--0.5),times=length(seq(0,350,by=5))),col="orange",type = "l",xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),lapply(seq(0,350,by=5),HR_EU.5),col="lightblue",type = "p",pch=4,xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),rep(exp(-0.5),times=length(seq(0,350,by=5))),col="lightblue",type = "l",xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),lapply(seq(0,350,by=5),OD_EU.5),col="lightblue",type = "p",pch=20,xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),lapply(seq(0,350,by=5),HR_EU.25),col="blue",type = "p",pch=4,xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),rep(exp(-0.25),times=length(seq(0,350,by=5))),col="blue",type = "l",xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),lapply(seq(0,350,by=5),OD_EU.25),col="blue",type = "p",pch=20,xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),lapply(seq(0,350,by=5),HR_EUn.25),col="red",type = "p",pch=4,xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),rep(exp(--0.25),times=length(seq(0,350,by=5))),col="red",type = "l",xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),lapply(seq(0,350,by=5),OD_EUn.25),col="red",type = "p",pch=20,xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),lapply(seq(0,350,by=5),HR_EU0),col="black",type = "p",pch=4,xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),rep(exp(0),times=length(seq(0,350,by=5))),col="black",type = "l",xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
lines(seq(0,350,by=5),lapply(seq(0,350,by=5),OD_EU0),col="black",type = "p",pch=20,xlab ="Time (months)",ylab="RR",xlim=c(0,350),ylim=c(0.5,1.7), cex.lab=1.5)
axis(1,at=c(0,50,100),labels=c("0","50","100"))
legend(x="topright",c("HR", "OR", "RR=1","RR=0.607","RR=0.779","RR=1.284","RR=1.649"),
       col=c("black","black","black","lightblue","blue","red","orange"), lty=c(0,0,0,0,0,0,0),pch=c(4,20,15,15,15,15,15), cex=0.9,bg="white")
dev.off()


#Figure E
#b=0
dev.new(width=5, height=5, unit="cm")
postscript(file="./plot_generation/plot/FigureE_a.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
plot(seq(0,300,by=0.1),lapply(seq(0,300,by=0.1),F_EU0),col="black",type = "l",lty=1,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
lines(seq(0,300,by=5),lapply(seq(0,300,by=5),F_LL0),col="black",type = "p",pch=20,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
lines(seq(0,300,by=5),lapply(seq(0,300,by=5),F_weib0),col="black",type = "p",pch=4,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
legend(x="bottomright", legend=c("PPR", "Weib PH","LL PO"),
       col=c("black", "black","black"), lty=c(1,0,0),pch=c(0,4,20), cex=1.2)
dev.off()

#b=-0.5
dev.new(width=5, height=5, unit="cm")
postscript(file="./plot_generation/plot/FigureE_b.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
plot(seq(0,300,by=0.1),lapply(seq(0,300,by=0.1),F_EU0),col="red",type = "l",lty=1,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
lines(seq(0,300,by=0.1),lapply(seq(0,300,by=0.1),F_EUn.5),col="blue",type = "l",lty=1,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
lines(seq(0,300,by=5),lapply(seq(0,300,by=5),F_weib0),col="orange",type = "p",pch=4,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
lines(seq(0,300,by=5),lapply(seq(0,300,by=5),F_weibn.5),col="lightblue",type = "p",pch=4,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
lines(seq(0,300,by=5),lapply(seq(0,300,by=5),F_LL0),col=rgb(0.6,0,0.6),type = "p",pch=20,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
lines(seq(0,300,by=5),lapply(seq(0,300,by=5),F_LLn.5),col=rgb(0,0.6,0),type = "p",pch=20,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
legend(x="bottomright", legend=c("PPR control","PPR treatment", "Weib PH control","Weib PH treatment", "RR PPR","RR Weib PH"),
       col=c("red", "blue", "orange","lightblue",rgb(0.6,0,0.6),rgb(0,0.6,0)), lty=c(1,1,0,0,0,0),pch=c(0,0,4,4,20,20), cex=1.2)
dev.off()

#b=0.5
dev.new(width=5, height=5, unit="cm")
postscript(file="./plot_generation/plot/FigureE_c.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
plot(seq(0,350,by=0.1),lapply(seq(0,350,by=0.1),F_EU0),col="red",type = "l",lty=1,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
lines(seq(0,350,by=0.1),lapply(seq(0,350,by=0.1),F_EU.5),col="blue",type = "l",lty=1,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
lines(seq(0,350,by=5),lapply(seq(0,350,by=5),F_weib0),col="orange",type = "p",pch=4,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1),cex.lab=1.5)
lines(seq(0,350,by=5),lapply(seq(0,350,by=5),F_weib.5),col="lightblue",type = "p",pch=4,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
lines(seq(0,300,by=5),lapply(seq(0,300,by=5),F_LL0),col=rgb(0.6,0,0.6),type = "p",pch=20,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
lines(seq(0,300,by=5),lapply(seq(0,300,by=5),F_LL.5),col=rgb(0,0.6,0),type = "p",pch=20,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
legend(x="bottomright", legend=c("PPR control","PPR treatment", "Weib PH control","Weib PH treatment", "RR PPR","RR Weib PH"),
       col=c("red", "blue", "orange","lightblue",rgb(0.6,0,0.6),rgb(0,0.6,0)), lty=c(1,1,0,0,0,0),pch=c(0,0,4,4,20,20), cex=1.2)
dev.off()


#b=-0.25
dev.new(width=5, height=5, unit="cm")
postscript(file="./plot_generation/plot/FigureE_f.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
plot(seq(0,300,by=0.1),lapply(seq(0,300,by=0.1),F_EU0),col="red",type = "l",lty=1,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
lines(seq(0,300,by=0.1),lapply(seq(0,300,by=0.1),F_EUn.25),col="blue",type = "l",lty=1,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1),cex.lab=1.5)
lines(seq(0,300,by=5),lapply(seq(0,300,by=5),F_weib0),col="orange",type = "p",pch=4,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
lines(seq(0,300,by=5),lapply(seq(0,300,by=5),F_weibn.25),col="lightblue",type = "p",pch=4,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1),cex.lab=1.5)
lines(seq(0,300,by=5),lapply(seq(0,300,by=5),F_LL0),col=rgb(0.6,0,0.6),type = "p",pch=20,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
lines(seq(0,300,by=5),lapply(seq(0,300,by=5),F_LLn.25),col=rgb(0,0.6,0),type = "p",pch=20,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
legend(x="bottomright", legend=c("PPR control","PPR treatment", "Weib PH control","Weib PH treatment", "RR PPR","RR Weib PH"),
       col=c("red", "blue", "orange","lightblue",rgb(0.6,0,0.6),rgb(0,0.6,0)), lty=c(1,1,0,0,0,0),pch=c(0,0,4,4,20,20), cex=1.2)
dev.off()

#b=0.25
dev.new(width=5, height=5, unit="cm")
postscript(file="./plot_generation/plot/FigureE_g.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
plot(seq(0,300,by=0.1),lapply(seq(0,300,by=0.1),F_EU0),col="red",type = "l",lty=1,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
lines(seq(0,300,by=0.1),lapply(seq(0,300,by=0.1),F_EU.25),col="blue",type = "l",lty=1,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
lines(seq(0,300,by=5),lapply(seq(0,300,by=5),F_weib0),col="orange",type = "p",pch=4,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
lines(seq(0,300,by=5),lapply(seq(0,300,by=5),F_weib.25),col="lightblue",type = "p",pch=4,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1),cex.lab=1.5)
lines(seq(0,300,by=5),lapply(seq(0,300,by=5),F_LL0),col=rgb(0.6,0,0.6),type = "p",pch=20,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
lines(seq(0,300,by=5),lapply(seq(0,300,by=5),F_LL.25),col=rgb(0,0.6,0),type = "p",pch=20,xlab ="Time (months)",ylab="Event probability",xlim=c(0,250),ylim=c(0,1), cex.lab=1.5)
legend(x="bottomright", legend=c("PPR control","PPR treatment", "Weib PH control","Weib PH treatment", "RR PPR","RR Weib PH"),
       col=c("red", "blue", "orange","lightblue",rgb(0.6,0,0.6),rgb(0,0.6,0)), lty=c(1,1,0,0,0,0),pch=c(0,0,4,4,20,20), cex=1.2)
dev.off()



#Figure 3 and F-M
#function to filter the estimates out from the output of cox.ph
est_cox<-function(cox){
  b_est<-numeric(length(cox))
  for(i in 1:length(cox)){
    c<-cox[[i]]
    b_est[i]<--c$coefficients[1] #negative to compare with beta
  }
  return(b_est)
}

#PR case
#Figure F  (30% censoring)
dev.new(width=14, height=7, unit="cm",noRStudioGD = TRUE)
postscript(file="./plot_generation/plot/FigureF.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar=c(5,5,4,1)+.1,cex.lab=1.2)
boxplot(Est_PR_NPPR_nna[[1]],-log(Est_PR_PPR_nna[[1]]),est_cox(Est_PR_Cox_nna[[1]]),-log(Est_PR_LL_nna[[1]]),
        Est_PR_NPPR_nna[[10]],-log(Est_PR_PPR_nna[[10]]),est_cox(Est_PR_Cox_nna[[10]]),-log(Est_PR_LL_nna[[10]]),
        Est_PR_NPPR_nna[[28]],-log(Est_PR_PPR_nna[[28]]),est_cox(Est_PR_Cox_nna[[28]]),-log(Est_PR_LL_nna[[28]]),
        Est_PR_NPPR_nna[[37]],-log(Est_PR_PPR_nna[[37]]),est_cox(Est_PR_Cox_nna[[37]]),-log(Est_PR_LL_nna[[37]]),
        Est_PR_NPPR_nna[[19]],-log(Est_PR_PPR_nna[[19]]),est_cox(Est_PR_Cox_nna[[19]]),-log(Est_PR_LL_nna[[19]]),
        names=rep(c("NPPR","PPR","Cox","LL"),5),col=rep(c("lightblue","orange",rgb(0,0.6,0),rgb(0.6,0,0.6)),5),border =rep(c("black","black","black","black"),5),outcol=rep(c("lightblue","orange",rgb(0,0.6,0),rgb(0.6,0,0.6)),5),ylab=expression(paste(hat(beta))),xlab=expression(paste(beta)))
axis(side = 1, line = 0.7, at = c(1.5, 3.5,5.5,7.5,9.5), labels = c("0", "0.5", "0.25","-0.25","-0.5"), tick = F)
clip(0.5,4.5, -100, 100)
abline( h =0, col = "black",lwd = 2.5)
clip(4.5,8.5, -100, 100)
abline( h =0.5, col = "black",lwd = 2.5)
clip(8.5,12.5, -100, 100)
abline( h =0.25, col = "black",lwd = 2.5)
clip(12.5,16.5, -100, 100)
abline( h =-0.25, col = "black",lwd = 2.5)
clip(16.5,20.5, -100, 100)
abline( h =-0.5, col = "black",lwd = 2.5)
dev.off()


#Figure 3 (50% censoring)
dev.new(width=14, height=7, unit="cm",noRStudioGD = TRUE)
postscript(file="./plot_generation/plot/Figure3.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar=c(5,5,4,1)+.1,cex.lab=1.2)
boxplot(Est_PR_NPPR_nna[[4]],-log(Est_PR_PPR_nna[[4]]),est_cox(Est_PR_Cox_nna[[4]]),-log(Est_PR_LL_nna[[4]]),
        Est_PR_NPPR_nna[[13]],-log(Est_PR_PPR_nna[[13]]),est_cox(Est_PR_Cox_nna[[13]]),-log(Est_PR_LL_nna[[13]]),
        Est_PR_NPPR_nna[[31]],-log(Est_PR_PPR_nna[[31]]),est_cox(Est_PR_Cox_nna[[31]]),-log(Est_PR_LL_nna[[31]]),
        Est_PR_NPPR_nna[[40]],-log(Est_PR_PPR_nna[[40]]),est_cox(Est_PR_Cox_nna[[40]]),-log(Est_PR_LL_nna[[40]]),
        Est_PR_NPPR_nna[[22]],-log(Est_PR_PPR_nna[[22]]),est_cox(Est_PR_Cox_nna[[22]]),-log(Est_PR_LL_nna[[22]]),
        names=rep(c("NPPR","PPR","Cox","LL"),5),col=rep(c("lightblue","orange",rgb(0,0.6,0),rgb(0.6,0,0.6)),5),border =rep(c("black","black","black","black"),5),outcol=rep(c("lightblue","orange",rgb(0,0.6,0),rgb(0.6,0,0.6)),5),ylab=expression(paste(hat(beta))),xlab=expression(paste(beta)))
axis(side = 1, line = 0.7, at = c(1.5, 3.5,5.5,7.5,9.5), labels = c("0", "0.5", "0.25","-0.25","-0.5"), tick = F)
clip(0.5,4.5, -100, 100)
abline( h =0, col = "black",lwd = 2.5)
clip(4.5,8.5, -100, 100)
abline( h =0.5, col = "black",lwd = 2.5)
clip(8.5,12.5, -100, 100)
abline( h =0.25, col = "black",lwd = 2.5)
clip(12.5,16.5, -100, 100)
abline( h =-0.25, col = "black",lwd = 2.5)
clip(16.5,20.5, -100, 100)
abline( h =-0.5, col = "black",lwd = 2.5)
dev.off()

#Figure G (70% censoring)
dev.new(width=14, height=7, unit="cm",noRStudioGD = TRUE)
postscript(file="./plot_generation/plot/FigureG.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar=c(5,5,4,1)+.1,cex.lab=1.2)
boxplot(Est_PR_NPPR_nna[[7]],-log(Est_PR_PPR_nna[[7]]),est_cox(Est_PR_Cox_nna[[7]]),-log(Est_PR_LL_nna[[7]]),
        Est_PR_NPPR_nna[[16]],-log(Est_PR_PPR_nna[[16]]),est_cox(Est_PR_Cox_nna[[16]]),-log(Est_PR_LL_nna[[16]]),
        Est_PR_NPPR_nna[[34]],-log(Est_PR_PPR_nna[[34]]),est_cox(Est_PR_Cox_nna[[34]]),-log(Est_PR_LL_nna[[34]]),
        Est_PR_NPPR_nna[[43]],-log(Est_PR_PPR_nna[[43]]),est_cox(Est_PR_Cox_nna[[43]]),-log(Est_PR_LL_nna[[43]]),
        Est_PR_NPPR_nna[[25]],-log(Est_PR_PPR_nna[[25]]),est_cox(Est_PR_Cox_nna[[25]]),-log(Est_PR_LL_nna[[25]]),
        names=rep(c("NPPR","PPR","Cox","LL"),5),col=rep(c("lightblue","orange",rgb(0,0.6,0),rgb(0.6,0,0.6)),5),border =rep(c("black","black","black","black"),5),outcol=rep(c("lightblue","orange",rgb(0,0.6,0),rgb(0.6,0,0.6)),5),ylab=expression(paste(hat(beta))),xlab=expression(paste(beta)))
axis(side = 1, line = 0.7, at = c(1.5, 3.5,5.5,7.5,9.5), labels = c("0", "0.5", "0.25","-0.25","-0.5"), tick = F)
clip(0.5,4.5, -100, 100)
abline( h =0, col = "black",lwd = 2.5)
clip(4.5,8.5, -100, 100)
abline( h =0.5, col = "black",lwd = 2.5)
clip(8.5,12.5, -100, 100)
abline( h =0.25, col = "black",lwd = 2.5)
clip(12.5,16.5, -100, 100)
abline( h =-0.25, col = "black",lwd = 2.5)
clip(16.5,20.5, -100, 100)
abline( h =-0.5, col = "black",lwd = 2.5)
dev.off()


#PH case
#Figure H (30% censoring)
dev.new(width=10, height=10, unit="cm",noRStudioGD = TRUE)
postscript(file="./plot_generation/plot/FigureH.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar=c(5,5,4,1)+.1,cex.lab=1.2)
boxplot(Est_PH_NPPR_nna[[1]],est_cox(Est_PH_Cox_nna[[1]]),
        Est_PH_NPPR_nna[[10]],est_cox(Est_PH_Cox_nna[[10]]),
        Est_PH_NPPR_nna[[28]],est_cox(Est_PH_Cox_nna[[28]]),
        Est_PH_NPPR_nna[[37]],est_cox(Est_PH_Cox_nna[[37]]),
        Est_PH_NPPR_nna[[19]],est_cox(Est_PH_Cox_nna[[19]]),
        names=rep(c("NPPR","Cox"),5),col=rep(c("lightblue",rgb(0,0.6,0)),5),border =rep(c("black","black"),5),outcol=rep(c("lightblue",rgb(0,0.6,0)),5),ylab=expression(paste(hat(beta))),xlab=expression(paste(beta)))
axis(side = 1, line = 0.7, at = c(1.5, 3.5,5.5,7.5,9.5), labels = c("0", "0.5", "0.25","-0.25","-0.5"), tick = F)
clip(0.25,2.75, -100, 100)
abline( h =0, col = "black",lwd = 2.5)
clip(2.25,4.75, -100, 100)
abline( h =0.5, col = "black",lwd = 2.5)
clip(4.25,6.75, -100, 100)
abline( h =0.25, col = "black",lwd = 2.5)
clip(6.25,8.75, -100, 100)
abline( h =-0.25, col = "black",lwd = 2.5)
clip(8.25,10.75, -100, 100)
abline( h =-0.5, col = "black",lwd = 2.5)
dev.off()

#Figure I (50% censoring)
dev.new(width=10, height=10, unit="cm",noRStudioGD = TRUE)
postscript(file="./plot_generation/plot/FigureI.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar=c(5,5,4,1)+.1,cex.lab=1.2)
boxplot(Est_PH_NPPR_nna[[4]],est_cox(Est_PH_Cox_nna[[4]]),
        Est_PH_NPPR_nna[[13]],est_cox(Est_PH_Cox_nna[[13]]),
        Est_PH_NPPR_nna[[31]],est_cox(Est_PH_Cox_nna[[31]]),
        Est_PH_NPPR_nna[[40]],est_cox(Est_PH_Cox_nna[[40]]),
        Est_PH_NPPR_nna[[22]],est_cox(Est_PH_Cox_nna[[22]]),
        names=rep(c("NPPR","Cox"),5),col=rep(c("lightblue",rgb(0,0.6,0)),5),border =rep(c("black","black"),5),outcol=rep(c("lightblue",rgb(0,0.6,0)),5),ylab=expression(paste(hat(beta))),xlab=expression(paste(beta)))
axis(side = 1, line = 0.7, at = c(1.5, 3.5,5.5,7.5,9.5), labels = c("0", "0.5", "0.25","-0.25","-0.5"), tick = F)
clip(0.25,2.75, -100, 100)
abline( h =0, col = "black",lwd = 2.5)
clip(2.25,4.75, -100, 100)
abline( h =0.5, col = "black",lwd = 2.5)
clip(4.25,6.75, -100, 100)
abline( h =0.25, col = "black",lwd = 2.5)
clip(6.25,8.75, -100, 100)
abline( h =-0.25, col = "black",lwd = 2.5)
clip(8.25,10.75, -100, 100)
abline( h =-0.5, col = "black",lwd = 2.5)
dev.off()

#Figure J (70% censoring)
dev.new(width=10, height=10, unit="cm",noRStudioGD = TRUE)
postscript(file="./plot_generation/plot/FigureJ.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar=c(5,5,4,1)+.1,cex.lab=1.2)
boxplot(Est_PH_NPPR_nna[[7]],est_cox(Est_PH_Cox_nna[[7]]),
        Est_PH_NPPR_nna[[16]],est_cox(Est_PH_Cox_nna[[16]]),
        Est_PH_NPPR_nna[[34]],est_cox(Est_PH_Cox_nna[[34]]),
        Est_PH_NPPR_nna[[43]],est_cox(Est_PH_Cox_nna[[43]]),
        Est_PH_NPPR_nna[[25]],est_cox(Est_PH_Cox_nna[[25]]),
        names=rep(c("NPPR","Cox"),5),col=rep(c("lightblue",rgb(0,0.6,0)),5),border =rep(c("black","black"),5),outcol=rep(c("lightblue",rgb(0,0.6,0)),5),ylab=expression(paste(hat(beta))),xlab=expression(paste(beta)))
axis(side = 1, line = 0.7, at = c(1.5, 3.5,5.5,7.5,9.5), labels = c("0", "0.5", "0.25","-0.25","-0.5"), tick = F)
clip(0.25,2.75, -100, 100)
abline( h =0, col = "black",lwd = 2.5)
clip(2.25,4.75, -100, 100)
abline( h =0.5, col = "black",lwd = 2.5)
clip(4.25,6.75, -100, 100)
abline( h =0.25, col = "black",lwd = 2.5)
clip(6.25,8.75, -100, 100)
abline( h =-0.25, col = "black",lwd = 2.5)
clip(8.25,10.75, -100, 100)
abline( h =-0.5, col = "black",lwd = 2.5)
dev.off()

#PO case
#Figure K (30% censoring)
dev.new(width=10, height=10, unit="cm",noRStudioGD = TRUE)
postscript(file="./plot_generation/plot/FigureK.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar=c(5,5,4,1)+.1,cex.lab=1.2)
boxplot(Est_PO_NPPR_nna[[1]],-log(Est_PO_LL_nna[[1]]),
        Est_PO_NPPR_nna[[10]],-log(Est_PO_LL_nna[[10]]),
        Est_PO_NPPR_nna[[28]],-log(Est_PO_LL_nna[[28]]),
        Est_PO_NPPR_nna[[37]],-log(Est_PO_LL_nna[[37]]),
        Est_PO_NPPR_nna[[19]],-log(Est_PO_LL_nna[[19]]),
        names=rep(c("NPPR","LL"),5),col=rep(c("lightblue",rgb(0.6,0,0.6)),5),border =rep(c("black","black"),5),outcol=rep(c("lightblue",rgb(0.6,0,0.6)),5),ylab=expression(paste(hat(beta))),xlab=expression(paste(beta)))
axis(side = 1, line = 0.7, at = c(1.5, 3.5,5.5,7.5,9.5), labels = c("0", "0.5", "0.25","-0.25","-0.5"), tick = F)
clip(0.25,2.75, -100, 100)
abline( h =0, col = "black",lwd = 2.5)
clip(2.25,4.75, -100, 100)
abline( h =0.5, col = "black",lwd = 2.5)
clip(4.25,6.75, -100, 100)
abline( h =0.25, col = "black",lwd = 2.5)
clip(6.25,8.75, -100, 100)
abline( h =-0.25, col = "black",lwd = 2.5)
clip(8.25,10.75, -100, 100)
abline( h =-0.5, col = "black",lwd = 2.5)
dev.off()

#Figure L (50% censoring)
dev.new(width=10, height=10, unit="cm",noRStudioGD = TRUE)
postscript(file="./plot_generation/plot/FigureL.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar=c(5,5,4,1)+.1,cex.lab=1.2)
boxplot(Est_PO_NPPR_nna[[4]],-log(Est_PO_LL_nna[[4]]),
        Est_PO_NPPR_nna[[13]],-log(Est_PO_LL_nna[[13]]),
        Est_PO_NPPR_nna[[31]],-log(Est_PO_LL_nna[[31]]),
        Est_PO_NPPR_nna[[40]],-log(Est_PO_LL_nna[[40]]),
        Est_PO_NPPR_nna[[22]],-log(Est_PO_LL_nna[[22]]),
        names=rep(c("NPPR","LL"),5),col=rep(c("lightblue",rgb(0.6,0,0.6)),5),border =rep(c("black","black"),5),outcol=rep(c("lightblue",rgb(0.6,0,0.6)),5),ylab=expression(paste(hat(beta))),xlab=expression(paste(beta)))
axis(side = 1, line = 0.7, at = c(1.5, 3.5,5.5,7.5,9.5), labels = c("0", "0.5", "0.25","-0.25","-0.5"), tick = F)
clip(0.25,2.75, -100, 100)
abline( h =0, col = "black",lwd = 2.5)
clip(2.25,4.75, -100, 100)
abline( h =0.5, col = "black",lwd = 2.5)
clip(4.25,6.75, -100, 100)
abline( h =0.25, col = "black",lwd = 2.5)
clip(6.25,8.75, -100, 100)
abline( h =-0.25, col = "black",lwd = 2.5)
clip(8.25,10.75, -100, 100)
abline( h =-0.5, col = "black",lwd = 2.5)
dev.off()

#Figure M (70% censoring)
dev.new(width=10, height=10, unit="cm",noRStudioGD = TRUE)
postscript(file="./plot_generation/plot/FigureM.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar=c(5,5,4,1)+.1,cex.lab=1.2)
boxplot(Est_PO_NPPR_nna[[7]],-log(Est_PO_LL_nna[[7]]),
        Est_PO_NPPR_nna[[16]],-log(Est_PO_LL_nna[[16]]),
        Est_PO_NPPR_nna[[34]],-log(Est_PO_LL_nna[[34]]),
        Est_PO_NPPR_nna[[43]],-log(Est_PO_LL_nna[[43]]),
        Est_PO_NPPR_nna[[25]],-log(Est_PO_LL_nna[[25]]),
        names=rep(c("NPPR","LL"),5),col=rep(c("lightblue",rgb(0.6,0,0.6)),5),border =rep(c("black","black"),5),outcol=rep(c("lightblue",rgb(0.6,0,0.6)),5),ylab=expression(paste(hat(beta))),xlab=expression(paste(beta)))
axis(side = 1, line = 0.7, at = c(1.5, 3.5,5.5,7.5,9.5), labels = c("0", "0.5", "0.25","-0.25","-0.5"), tick = F)
clip(0.25,2.75, -100, 100)
abline( h =0, col = "black",lwd = 2.5)
clip(2.25,4.75, -100, 100)
abline( h =0.5, col = "black",lwd = 2.5)
clip(4.25,6.75, -100, 100)
abline( h =0.25, col = "black",lwd = 2.5)
clip(6.25,8.75, -100, 100)
abline( h =-0.25, col = "black",lwd = 2.5)
clip(8.25,10.75, -100, 100)
abline( h =-0.5, col = "black",lwd = 2.5)


#Figure 4-5 and O (event: death of all causes)
DapaHF_ac_1<-subset(DapaHF_ac,DapaHF_ac$Treat_bin==1)
DapaHF_ac_0<-subset(DapaHF_ac,DapaHF_ac$Treat_bin==0)

DapaHF_ac_1_event<-subset(DapaHF_ac_1,DapaHF_ac_1$Event==1)
DapaHF_ac_0_event<-subset(DapaHF_ac_0,DapaHF_ac_0$Event==1)


DapaHF_ac_1_cens<-subset(DapaHF_ac_1,DapaHF_ac_1$Event==0)
DapaHF_ac_0_cens<-subset(DapaHF_ac_0,DapaHF_ac_0$Event==0)


DapaHF_ac_event<-T_jump(DapaHF_ac_1_event$SurvivalTimeMonths,DapaHF_ac_0_event$SurvivalTimeMonths)


F_Dapa_ac_1<-function(t){
  F_est(t,DapaHF_ac_1,col=c(col_time=3,col_cens=2))
}

F_Dapa_ac_0<-function(t){
  F_est(t,DapaHF_ac_0,col=c(col_time=3,col_cens=2))
}

W_Dapa_ac_1<-function(t){
  KM_se(t,DapaHF_ac_1,col=c(col_time=3,col_cens=2))
}
W_Dapa_ac_0<-function(t){
  KM_se(t,DapaHF_ac_0,col=c(col_time=3,col_cens=2))
}

Est_timepoint_DapaHF_ac<-function(t){
  -(log(F_Dapa_ac_1(t)/F_Dapa_ac_0(t)))
}
Est_weights_DapaHF_ac<-function(t){
  1/( 1/((1/F_Dapa_ac_1(t)^2)*W_Dapa_ac_1(t)^2+(1/F_Dapa_ac_0(t)^2)*W_Dapa_ac_0(t)^2))
}

t_Dapa_ac_max<-max(DapaHF_ac$SurvivalTimeMonths)
t_Dapa_ac_max0<-max(DapaHF_ac_0$SurvivalTimeMonths)
t_Dapa_ac_max1<-max(DapaHF_ac_1$SurvivalTimeMonths)

#Figure 4
dev.new(width=5, height=5, unit="cm")
postscript(file="./plot_generation/plot/Figure4.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
plot(seq(0,t_Dapa_ac_max1,by=0.1),lapply(seq(0,t_Dapa_ac_max1,by=0.1),F_Dapa_ac_1),col="blue",type = "l",xlab ="Time (months)",ylab="Event probability",xlim=c(0,t_Dapa_ac_max),ylim=c(0,0.3),cex.lab=1.5)
points(DapaHF_ac_1_cens$SurvivalTimeMonths,lapply(DapaHF_ac_1_cens$SurvivalTimeMonths,F_Dapa_ac_1),col="blue",pch=4)
lines(seq(0,t_Dapa_ac_max0,by=0.1),lapply(seq(0,t_Dapa_ac_max0,by=0.1),F_Dapa_ac_0),col="red",type = "l",xlab ="Time (months)",ylab="Event probability",xlim=c(0,t_Dapa_ac_max),ylim=c(0,0.3),cex.lab=1.5)
points(DapaHF_ac_0_cens$SurvivalTimeMonths,lapply(DapaHF_ac_0_cens$SurvivalTimeMonths,F_Dapa_ac_0),col="red",pch=4)
legend(x="topleft", legend=c("Placebo","Placebo censored", "Dapagliflozin","Dapagliflozin censored"),
       col=c("red","red", "blue", "blue"), lty=c(1,0,1,0),pch=c(1,4,1,4), cex=1.2)
dev.off()


#Figure 5a
dev.new(width=6, height=5, unit="cm")
postscript(file="./plot_generation/plot/Figure5_a.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar=c(5,6,4,1)+.1)
plot(DapaHF_ac_event,lapply(DapaHF_ac_event,Est_timepoint_DapaHF_ac),col="blue",type = "p",xlab ="Time (months)",ylab=expression(beta),xlim=c(0,t_Dapa_ac_max),cex.lab=1.5)
abline( h =estimate_ac, col = "black")
dev.off()

#Figure 5b
dev.new(width=5, height=5, unit="cm")
postscript(file="./plot_generation/plot/Figure5_b.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar=c(5,6,4,1)+.1)
plot(DapaHF_ac_event,lapply(DapaHF_ac_event,Est_weights_DapaHF_ac),col="red",xlab ="Time (months)",ylab="Weights",xlim=c(0,t_Dapa_ac_max),cex.lab=1.5)
dev.off()


rd_Dapahf_ac<-function(t){
  rd<-(1-exp(-estimate_ac))*F_Dapa_ac_0(t)
  return(rd)
}
needtotreat_Dapahf_ac<-function(t){
  nnt<-1/rd_Dapahf_ac(t)
  return(nnt)
}
max(as.numeric(lapply(seq(min(DapaHF_ac_0_event$SurvivalTimeMonths),t_Dapa_ac_max1,by=0.1),needtotreat_Dapahf_ac)))

#Figure L
dev.new(width=6, height=5, unit="cm")
postscript(file="./plot_generation/plot/FigureO.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar=c(5,6,4,1)+.1)
plot(seq(min(DapaHF_ac_0_event$SurvivalTimeMonths),t_Dapa_ac_max1,by=0.1),lapply(seq(min(DapaHF_ac_0_event$SurvivalTimeMonths),t_Dapa_ac_max1,by=0.1),needtotreat_Dapahf_ac),col="lightgreen",type = "l",xlab ="Time (months)",ylab="Number needed to treat",xlim=c(0,t_Dapa_ac_max),ylim=c(0,1500),cex.lab=1.5)
dev.off()



#Figure C (event: heart failure)
DapaHF_1<-subset(DapaHF,DapaHF$Arm==1)
DapaHF_0<-subset(DapaHF,DapaHF$Arm==0)

DapaHF_1_event<-subset(DapaHF_1,DapaHF_1$Event==1)
DapaHF_0_event<-subset(DapaHF_0,DapaHF_0$Event==1)


DapaHF_1_cens<-subset(DapaHF_1,DapaHF_1$Event==0)
DapaHF_0_cens<-subset(DapaHF_0,DapaHF_0$Event==0)


DapaHF_event<-T_jump(DapaHF_1_event$SurvivalTimeMonths,DapaHF_0_event$SurvivalTimeMonths)


F_Dapa_1<-function(t){
  F_est(t,DapaHF_1,col=c(col_time=3,col_cens=1))
}

F_Dapa_0<-function(t){
  F_est(t,DapaHF_0,col=c(col_time=3,col_cens=1))
}

W_Dapa_1<-function(t){
  KM_se(t,DapaHF_1,col=c(col_time=3,col_cens=1))
}
W_Dapa_0<-function(t){
  KM_se(t,DapaHF_0,col=c(col_time=3,col_cens=1))
}

Est_timepoint_DapaHF<-function(t){
  -(log(F_Dapa_1(t)/F_Dapa_0(t)))
}
Est_weights_DapaHF<-function(t){
  1/( 1/((1/F_Dapa_1(t)^2)*W_Dapa_1(t)^2+(1/F_Dapa_0(t)^2)*W_Dapa_0(t)^2))
}

t_Dapa_max<-max(DapaHF$SurvivalTimeMonths)
t_Dapa_max0<-max(DapaHF_0$SurvivalTimeMonths)
t_Dapa_max1<-max(DapaHF_1$SurvivalTimeMonths)

#Figure C a)
dev.new(width=5, height=5, unit="cm")
postscript(file="./plot_generation/plot/FigureC_a.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
plot(seq(0,t_Dapa_max1,by=0.1),lapply(seq(0,t_Dapa_max1,by=0.1),F_Dapa_1),col="blue",type = "l",xlab ="Time (months)",ylab="Event probability",xlim=c(0,t_Dapa_max),ylim=c(0,0.3),cex.lab=1.5)
points(DapaHF_1_cens$SurvivalTimeMonths,lapply(DapaHF_1_cens$SurvivalTimeMonths,F_Dapa_1),col="blue",pch=4)
lines(seq(0,t_Dapa_max0,by=0.1),lapply(seq(0,t_Dapa_max0,by=0.1),F_Dapa_0),col="red",type = "l",xlab ="Time (months)",ylab="Event probability",xlim=c(0,t_Dapa_max),ylim=c(0,0.3),cex.lab=1.5)
points(DapaHF_0_cens$SurvivalTimeMonths,lapply(DapaHF_0_cens$SurvivalTimeMonths,F_Dapa_0),col="red",pch=4)
legend(x="topleft", legend=c("Placebo","Placebo censored", "Dapagliflozin","Dapagliflozin censored"),
       col=c("red","red", "blue", "blue"), lty=c(1,0,1,0),pch=c(1,4,1,4), cex=1.2)
dev.off()

#Figure C b)
dev.new(width=6, height=5, unit="cm")
postscript(file="./plot_generation/plot/FigureC_b.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar=c(5,6,4,1)+.1)
plot(DapaHF_event,lapply(DapaHF_event,Est_timepoint_DapaHF),col="blue",type = "p",xlab ="Time (months)",ylab=expression(beta),xlim=c(0,t_Dapa_max),cex.lab=1.5)
abline( h =estimate, col = "black")
dev.off()

#Figure C c)
dev.new(width=5, height=5, unit="cm")
postscript(file="./plot_generation/plot/FigureC_c.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar=c(5,6,4,1)+.1)
plot(DapaHF_event,lapply(DapaHF_event,Est_weights_DapaHF),col="red",xlab ="Time (months)",ylab="Weights",xlim=c(0,t_Dapa_max),cex.lab=1.5)
dev.off()

rd_Dapahf<-function(t){
  rd<-(1-exp(-estimate))*F_Dapa_0(t)
  return(rd)
}
needtotreat_Dapahf<-function(t){
  nnt<-1/rd_Dapahf(t)
  return(nnt)
}
max(as.numeric(lapply(seq(min(DapaHF_0_event$SurvivalTimeMonths),t_Dapa_max1,by=0.1),needtotreat_Dapahf)))

#Figure C d)
dev.new(width=6, height=5, unit="cm")
postscript(file="./plot_generation/plot/FigureC_d.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar=c(5,6,4,1)+.1)
plot(seq(min(DapaHF_0_event$SurvivalTimeMonths),t_Dapa_max1,by=0.1),lapply(seq(min(DapaHF_0_event$SurvivalTimeMonths),t_Dapa_max1,by=0.1),needtotreat_Dapahf),col="lightgreen",type = "l",xlab ="Time (months)",ylab="Number needed to treat",xlim=c(0,t_Dapa_max),ylim=c(0,1500),cex.lab=1.5)
dev.off()


#Figure N Boxplot nocens
dev.new(width=14, height=7, unit="cm",noRStudioGD = TRUE)
postscript(file="./plot_generation/plot/FigureN.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
par(mar=c(5,5,4,1)+.1,cex.lab=1.2)
boxplot(Est_PR_NPPR_nocens[[10]],Est_PR_NPPR_nocens[[11]],Est_PR_NPPR_nocens[[12]],
        est_cox(Est_PH_Cox_nocens[[10]]),est_cox(Est_PH_Cox_nocens[[11]]),est_cox(Est_PH_Cox_nocens[[12]]),
        -log(Est_PO_LL_nocens_nna[[10]]),-log(Est_PO_LL_nocens_nna[[11]]),-log(Est_PO_LL_nocens_nna[[12]]),
        names=c("NPPR 500","NPPR 100","NPPR 50","Cox 500","Cox 100","Cox 50","LL 500","LL 100","LL 50"),
        col=c("lightblue","lightblue","lightblue",rgb(0,0.6,0),rgb(0,0.6,0),rgb(0,0.6,0),rgb(0.6,0,0.6),rgb(0.6,0,0.6),rgb(0.6,0,0.6)),
        border =rep(c("black"),9),
        outcol=c("lightblue","lightblue","lightblue",rgb(0,0.6,0),rgb(0,0.6,0),rgb(0,0.6,0),rgb(0.6,0,0.6),rgb(0.6,0,0.6),rgb(0.6,0,0.6)),
        ylab=expression(paste(hat(beta))))
abline( h =0.25, col = "black",lwd = 2.5)
dev.off()


#Figure B Cox RR
Cox_RR.example<-function(t){
  EU<-(0.009*t)^0.859
  num<-1-(1-EU)^exp(0.5)
  rr<-num/EU
  return(rr)
}

dev.new(width=6, height=5, unit="cm")
postscript(file="./plot_generation/plot/FigureB.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
plot(1:100,lapply(1:100,Cox_RR.example),type="l",xlab="Time",ylab="RR")
dev.off()