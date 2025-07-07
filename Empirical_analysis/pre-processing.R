#This code shows how we detrend the parts corresponding to the same time period
#of time series by substracting the linear trend, 
#and how we implement the hierarchical clustering using the Spearman coefficient
#in order to select the less correlated input curves (the less correlated stocks).

#At the end we detrend and smooth those curves throughout their entire length and save them as input curves for
#functional motif discovery with ProbKMA.

#The same procedure has been done with all the other curves that have been excluded from functional motif disvery
#due to the correlation matters, but that were employed for functional motif search.

#For the correlation matters we consider the training interval in common
#for all the stocks from April 3, 2014 to December 31, 2020.



# Set directory to the folder with the current file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("../fmd_functions.R")#It is necessary to upload the 
#file with all the necessary functions

library("dendextend")
#Function for hierarchical clustering using the spearman coefficient
clustering_spearman <- function(data,file_name,h){
  correlation <- cor(data,method="spearman");
  abs.corr=as.dist(1-abs(correlation))
  fit <- hclust(abs.corr)
  pdf(paste0("clustering_spearman",file_name,".pdf"),width=12,height=12)
  par(mar=c(4,1,1,14)+0.1)
  dend=as.dendrogram(fit,hang=0.05)
  dend_cut=cutree(dend,h=h)[order.dendrogram(dend)]
  dend_col=rep(1,length(dend_cut))
  cluster=which(table(dend_cut)>1)
  for(k in seq_along(cluster)){
    dend_col[dend_cut==cluster[k]]=k+2
  }
  labels_colors(dend)=dend_col
  plot(dend,xlab= "1-|Spearman's correlation|", ylab="",main=NULL,cex=0.7,cex.lab=1.5,horiz=TRUE,xlim=c(1,0))
  lines(c(h,h),c(-2,nrow(correlation)+1),col='red',lty='dashed',lwd=2)
  dev.off()
  return(list(fit = fit, cluster = cluster))
}

library("readxl")
#Uploading the file excel with the stock prices
titles=read_excel("SP500_wrds.xlsx")
names_stocks=c("AAPL","MSFT","AMZN","GOOGL","GOOG","TSLA","BRKb","NVDA",
               "JNJ","UNH","FB","XOM","JPM","VISA","PG","CVX",
               "HD","MA","PFE","ABBV","BAC","KO","LLY","AVGO","PEP","MRK",
               "TMO","VZ","COST","ABT","ADBE","DIS","ACN","CMCSA","CSCO",
               "MCD","CRM","WMT","INTC","WFC","AMD","LIN","DHR","PM","BMY",
               "QCOM","NEE","TXN","NKE","COP","SP500")
list_entire=list()
for(i in 1:length(names_stocks)){
  list_entire[[i]]=data.frame(read_excel("SP500_wrds.xlsx",sheet=names_stocks[i])$Date,read_excel("SP500_wrds.xlsx",sheet=names_stocks[i])$Price)
}
names(list_entire)=names_stocks

#Giving the names to the columns
names_col=function(x){
  colnames(x)=c("Date","Price")
  return(x)
}
list_entire_fin=lapply(list_entire,names_col)



#Comparing the dates between SP500 and the other stocks
SP_dates=list()
only_sp_dates=list()#here are the dates that are in SP500 but not in the other stocks
library(dplyr)
library(anytime)
for(i in 1:length(list_entire_fin)){
  if(length(anydate(list_entire_fin[[i]]$Date,tz="UTC"))!=length(anydate(list_entire_fin$SP500$Date,tz="UTC"))){
    SP_dates[[i]]=data.frame(Date=((list_entire_fin$SP500[anydate(list_entire_fin$SP500$Date,tz="UTC")>=anydate(list_entire_fin[[i]]$Date,tz="UTC")[1],])$Date),
                             Price=((list_entire_fin$SP500[anydate(list_entire_fin$SP500$Date,tz="UTC")>=anydate(list_entire_fin[[i]]$Date,tz="UTC")[1],])$Price))
    
    #removing the prices for the date not present in SP500 but present in the other stocks
    list_entire_fin[[i]]=list_entire_fin[[i]][which(list_entire_fin[[i]]$Date%in%SP_dates[[i]]$Date),]
    #imputig the prices present in SP500 but absent in the other stocks
    
    list_entire_fin[[i]]= SP_dates[[i]] %>%
      dplyr::select(Date) %>%
      left_join(list_entire_fin[[i]], by = 'Date') %>%
      mutate(Price = zoo::na.approx(Price))
    
  } }


#For the correlation purposes we consider
#the time interval from April 3, 2014 to December 31, 2020.
#Because GOOG is the stock that has less available data and it starts on April 3, 2014
common=function(x){
  x=x[x$Date>="2014-04-03",]#Year, month, day
}


list_common=lapply(list_entire_fin,common)

training=function(x){
  x=x[x$Date<="2020-12-31",]
}

list_entire_training=lapply(list_entire_fin,training)

list_common_training=lapply(list_common,training)


#Detrending the common part substracting each title's trend separately
library("lubridate")#for decimal_date              #year-month-day
dates_stocks=c(rep(decimal_date(as.Date("2014-04-03")),length(list_entire_fin)))
list_common_training_ts=list()
for(i in 1:length(names_stocks)){
  list_common_training_ts[[i]]=ts(list_common_training[[i]]$Price,start = dates_stocks[[i]], frequency = 250)
}
names(list_common_training_ts)=names_stocks


#Function to substract the linear trend
detrend_calculate=function(series){
  reg=lm(series ~ seq(1:length(series)))
  detrended = series - reg$fitted.values
  return(detrended)
}



titles_common_detrended=lapply(list_common_training_ts,detrend_calculate)
titles_common_detrended_2=data.frame(titles_common_detrended)

# Clustering based on Spearman correlation
res_spearman = clustering_spearman(titles_common_detrended_2,"_h=0.3_regression_wrds",h=0.3)
pdf('./plot_spearman.pdf', height = 12, width = 18)
plot(res_spearman$fit)
dev.off()
plot(res_spearman$fit)

#With h=0.3 we consider 30 stocks

#The final dataset 
final_titles_common=titles_common_detrended_2



BAC_norm=(final_titles_common$BAC-min(final_titles_common$BAC))/(max(final_titles_common$BAC)-min(final_titles_common$BAC))

JPM_norm=(final_titles_common$JPM-min(final_titles_common$JPM))/(max(final_titles_common$JPM)-min(final_titles_common$JPM))

WFC_norm=(final_titles_common$WFC-min(final_titles_common$WFC))/(max(final_titles_common$WFC)-min(final_titles_common$WFC))


plot(BAC_norm,type="l",ylab="Normalized price",xlab="Time")
lines(JPM_norm,col="red")
lines(WFC_norm,col="green")
legend("topright", legend=c("BAC", "JPM","WFC"),
       col=c("black", "red","green"),lty=1, cex=0.8)
#remove the most correlated stocks
final_titles_common$ABBV=NULL
final_titles_common$TXN=NULL
final_titles_common$CSCO=NULL
final_titles_common$BAC=NULL
final_titles_common$WFC=NULL
final_titles_common$XOM=NULL
final_titles_common$GOOG=NULL
final_titles_common$AMD=NULL
final_titles_common$TSLA=NULL
final_titles_common$CRM=NULL
final_titles_common$ADBE=NULL
final_titles_common$LIN=NULL
final_titles_common$ABT=NULL
final_titles_common$MA=NULL
final_titles_common$VISA=NULL
final_titles_common$COST=NULL
final_titles_common$DHR=NULL
final_titles_common$TMO=NULL
final_titles_common$QCOM=NULL
final_titles_common$NEE=NULL
final_titles_common$PG=NULL

######################################
######################################
######################################
#Save the uncorrelated and smoothed stocks
list_entire_uncorrelated=list(list_entire_fin$AAPL,list_entire_fin$MSFT,
                              list_entire_fin$AMZN,list_entire_fin$GOOGL,
                              list_entire_fin$BRKb,list_entire_fin$NVDA,
                              list_entire_fin$JNJ,list_entire_fin$UNH,
                              list_entire_fin$FB,list_entire_fin$JPM,
                              list_entire_fin$CVX,list_entire_fin$HD,
                              list_entire_fin$PFE,list_entire_fin$KO,
                              list_entire_fin$LLY,list_entire_fin$AVGO,
                              list_entire_fin$PEP,list_entire_fin$MRK,
                              list_entire_fin$VZ,list_entire_fin$DIS,
                              list_entire_fin$ACN,list_entire_fin$CMCSA,
                              list_entire_fin$MCD,list_entire_fin$WMT,
                              list_entire_fin$INTC,list_entire_fin$PM,
                              list_entire_fin$BMY,list_entire_fin$NKE,
                              list_entire_fin$COP,list_entire_fin$SP500)
names(list_entire_uncorrelated)=names(final_titles_common)
#Giving the names to the columns
names_col=function(x){
  colnames(x)=c("Date","Price")
  return(x)
}
list_entire_fin_uncorrelated=lapply(list_entire_uncorrelated,names_col)

#Comparing the dates between SP500 and the other stocks
SP_dates=list()
only_sp_dates=list()#here are the dates that are in SP500 but not in the other stocks
library(dplyr)
library(anytime)
for(i in 1:length(list_entire_fin_uncorrelated)){
  if(length(anydate(list_entire_fin_uncorrelated[[i]]$Date,tz="UTC"))!=length(anydate(list_entire_fin_uncorrelated$SP500$Date,tz="UTC"))){
    SP_dates[[i]]=data.frame(Date=((list_entire_fin_uncorrelated$SP500[anydate(list_entire_fin_uncorrelated$SP500$Date,tz="UTC")>=anydate(list_entire_fin_uncorrelated[[i]]$Date,tz="UTC")[1],])$Date),
                             Price=((list_entire_fin_uncorrelated$SP500[anydate(list_entire_fin_uncorrelated$SP500$Date,tz="UTC")>=anydate(list_entire_fin_uncorrelated[[i]]$Date,tz="UTC")[1],])$Price))
    
    #removing the prices for the date not present in SP500 but present in the other stocks
    list_entire_fin_uncorrelated[[i]]=list_entire_fin_uncorrelated[[i]][which(list_entire_fin_uncorrelated[[i]]$Date%in%SP_dates[[i]]$Date),]
    #imputig the prices present in SP500 but absent in the other stocks
    
    list_entire_fin_uncorrelated[[i]]= SP_dates[[i]] %>%
      dplyr::select(Date) %>%
      left_join(list_entire_fin_uncorrelated[[i]], by = 'Date') %>%
      mutate(Price = zoo::na.approx(Price))
    
  } }


#Subsetting the dates until Dec 31, 2020
training=function(x){
  x=x[x$Date<="2020-12-31",]
}

test=function(x){
  x=x[x$Date>"2020-12-31",]
}

list_training_uncorrelated=lapply(list_entire_fin_uncorrelated,training)


library("lubridate")
Dates=c(rep("2009-03-10",4),
        rep("2009-03-10",4),
        "2012-05-18",rep("2009-03-10",6),"2009-08-06",
        rep("2009-03-10",14)
)
list_training_ts_uncorr=list()
for(i in 1:length(list_training_uncorrelated)){
  list_training_ts_uncorr[[i]]=ts(list_training_uncorrelated[[i]]$Price,start=decimal_date(as.Date(Dates[i])),frequency=250)
}

names(list_training_ts_uncorr)=names(list_training_uncorrelated)


list_test=lapply(list_entire_fin_uncorrelated,test)

titles_of_interest_detrend=lapply(list_training_ts_uncorr,detrend_calculate)


names(titles_of_interest_detrend)=names(final_titles_common)



# smoothing detrended data

bandwidth_smooth = 10

data_smoothed=list()
for(i in 1:length(titles_of_interest_detrend)){
  data_smoothed[[i]]=smooth_locpoly_single(titles_of_interest_detrend[[i]],degree=3,bandwidth=bandwidth_smooth,v0=TRUE,v1=TRUE)
}


names(data_smoothed)=names(titles_of_interest_detrend)


save(data_smoothed,file="data_wrds.Rdata")

par(mfrow=c(2,1),mar=c(4,4,2,2))
plot(titles_of_interest_detrend$AAPL,type="l",main="AAPL",xlab="Time",ylab="Price")
plot(data_smoothed$AAPL$v0,type="l",main="AAPL",xlab="Time",ylab="Price")






list_training=lapply(list_entire_fin,training)

library("lubridate")
Dates=c(rep("2009-03-10",4),
        "2014-04-03","2010-06-29","2009-03-10","2009-03-10",
        "2009-03-10","2009-03-10","2012-05-18","2009-03-10",
        rep("2009-03-10",4),
        "2009-03-10","2009-03-10","2009-03-10","2013-01-02",
        "2009-03-10","2009-03-10","2009-03-10","2009-08-06",
        rep("2009-03-10",27)
)
list_training_ts=list()
for(i in 1:length(list_training)){
  list_training_ts[[i]]=ts(list_training[[i]]$Price,start=decimal_date(as.Date(Dates[i])),frequency=250)
}

names(list_training_ts)=names_stocks








list_test=lapply(list_entire_fin,test)

#detrend substracting each series trend

titles_of_interest_detrend=lapply(list_training_ts,detrend_calculate)


names(titles_of_interest_detrend)=c("AAPL","MSFT","AMZN","GOOGL","GOOG","TSLA","BRKb","NVDA",
                                    "JNJ","UNH","FB","XOM","JPM","VISA","PG","CVX",
                                    "HD","MA","PFE","ABBV","BAC","KO","LLY","AVGO","PEP","MRK",
                                    "TMO","VZ","COST","ABT","ADBE","DIS","ACN","CMCSA","CSCO",
                                    "MCD","CRM","WMT","INTC","WFC","AMD","LIN","DHR","PM","BMY",
                                    "QCOM","NEE","TXN","NKE","COP","SP500")


data_smoothed_all_stocks=list()
for(i in 1:length(titles_of_interest_detrend)){
  data_smoothed_all_stocks[[i]]=smooth_locpoly_single(titles_of_interest_detrend[[i]],degree=3,bandwidth=bandwidth_smooth,v0=TRUE,v1=TRUE)
}


names(data_smoothed_all_stocks)=names(titles_of_interest_detrend)



save(data_smoothed_all_stocks,file="wrds_all_data.Rdata")

