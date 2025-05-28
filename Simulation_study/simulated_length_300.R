#This code was used to simulate 20 curves of length 300 following the MAR(1,1) process described
#in the Simulation Study.

# Set directory to the folder with the current file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("../fmd_functions.R")#It is necessary to upload the 
#file with all the necessary functions

#Mixed causal-noncausal autoregressive process
library("stabledist")
library("fda")
#Function to generate from mixed causal-noncausal autoregressive process of order 1
MAR11=function(time,alpha,beta,scale,location,
               psi_noncausal,phi_causal,
               u_terminal,yt_init){
  epsilon=u=yt=mar=c()
  yt[1]=yt_init
  u[time]=u_terminal
  for(t in 2:time){
    epsilon[time-t+1]=rstable(1,alpha=alpha,beta=beta,gamma=scale,delta=location)
    u[time-t+1]=psi_noncausal*u[time-(t-2)]+epsilon[time-t+1]
    yt[t]=phi_causal*yt[t-1]
  }
  mar=yt+u
  newlist=list(noncausal=u,causal=yt,mixed=mar)
  return(newlist)
  
}



#separate seed for each simulated curve
setseeds_new=c(5080,5190,5255,5455,5625,5635,5775,5860,5880,5900,5915,
               5925,5960,6210,6265,6295,6355,6370,6640,6500)
mar2b=matrix(NA,300,20)
for(i in 1:20){
  set.seed(setseeds_new[i])
  mar2b[,i]=MAR11(time=300,alpha=1.6,beta=1,scale=1,location=0,
                  psi_noncausal=0.9,phi_causal=0.7,u_terminal=5,
                  yt_init=5)$mixed
}
matplot(mar2b,type="l",xlab="Time", ylab="",ylim=c(-10,90))
simulated_data=mar2b



#Smoothing with local polynomial regression
data_smoothed=list()
data_smoothed[[1]]=matrix(NA,nrow=nrow(simulated_data),ncol=ncol(simulated_data))
data_smoothed[[2]]=matrix(NA,nrow=nrow(simulated_data),ncol=ncol(simulated_data))
for(i in 1:ncol(simulated_data)){
  data_smoothed[[1]][,i]=smooth_locpoly_single(simulated_data[,i],degree=3,bandwidth=2,v0=TRUE,v1=TRUE)$v0
  data_smoothed[[2]][,i]=smooth_locpoly_single(simulated_data[,i],degree=3,bandwidth=2,v0=TRUE,v1=TRUE)$v1
}
names(data_smoothed)=c("data_smoothed","data_smoothed_derivative")
simulated300_data=data_smoothed


save(simulated300_data, file = './sim300_fmd/data_smoothed.RData')
