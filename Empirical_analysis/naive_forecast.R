#This code was used to compute the naive forecast

# Set directory to the folder with the current file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library("readxl")
#Uploading the file excel with the stock prices
titles=read_excel("SP500_wrds.xlsx")
names_titles=c("AAPL","MSFT","AMZN","GOOGL","GOOG","TSLA","BRKb","NVDA",
                     "JNJ","UNH","FB","XOM","JPM","VISA","PG","CVX",
                     "HD","MA","PFE","ABBV","BAC","KO","LLY","AVGO","PEP","MRK",
                     "TMO","VZ","COST","ABT","ADBE","DIS","ACN","CMCSA","CSCO",
                     "MCD","CRM","WMT","INTC","WFC","AMD","LIN","DHR","PM","BMY",
                     "QCOM","NEE","TXN","NKE","COP","SP500")
list_entire=list()
for(i in 1:length(names_titles)){
  list_entire[[i]]=data.frame(read_excel("SP500_wrds.xlsx",sheet=names_titles[i])$Date,read_excel("SP500_wrds.xlsx",sheet=names_titles[i])$Price)
}
names(list_entire)=names_titles

#Giving the names to the columns
names_col=function(x){
  colnames(x)=c("Date","Price")
  return(x)
}
list_entire_fin=lapply(list_entire,names_col)

#Comparing the dates between SP500 and the other stocks
SP_dates=list()
only_sp_dates=list()#here are the dates that are in SP500 but not in the other stocks
#'%!in%' <- function(x,y)!('%in%'(x,y)) 
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





#Transforming the date into an appropriate format and subsetting until Dec 31, 2020
training=function(x){
  x=x[x$Date<="2020-12-31",]
}

test=function(x){
  x=x[x$Date>"2020-12-31",]
}

list_training=lapply(list_entire_fin,training)


list_test=lapply(list_entire_fin,test)



rmse=function(series,predictions)
{
  sqrt(mean((series-predictions)^2))#(series-predictions) or error/residuals
}

mae=function(series, predictions)
{
  mean(abs(series-predictions))
}


#Naive forecast where the price at time t is equal to the price at time (t-1)
forecast_all=final_forecast_all=list()
for(i in 1:length(list_training)){
  forecast_all[[i]]=list()
  forecast_all[[i]][[1]]=tail(list_training[[i]]$Price,1)
  for(t in 1:length(list_test[[i]]$Price)){
    forecast_all[[i]][[t+1]]=list_test[[i]]$Price[t]
  }
  final_forecast_all[[i]]=head(unlist(forecast_all[[i]]),-1)
}

#Calculate the accuracy measures of the prices
rmse_all=mae_all=list()
for(i in (1:length(list_test))){
  rmse_all[[i]]=rmse(list_test[[i]]$Price,final_forecast_all[[i]])
  mae_all[[i]]=mae(list_test[[i]]$Price,final_forecast_all[[i]])
}
names(rmse_all)=names(mae_all)=names(list_test)

accuracy=data.frame(stocks=names(list_test),RMSE=unlist(rmse_all),MAE=unlist(mae_all))
rownames(accuracy)=names(list_test)


for(i in 1:length(list_training)){
  plot(list_test[[i]]$Date,list_test[[i]]$Price,type="l",ylim=range(c(list_test[[i]]$Price),c(final_forecast_all[[i]])),main=names(list_test)[[i]])
  lines(list_test[[i]]$Date,final_forecast_all[[i]],col="red")
}
dev.off()


#Calculating the observed returns
real_returns=list_test
for(i in 1:length(list_test)){
  for(t in 2:length(list_test[[1]]$Price)){
    real_returns[[i]][t,2]=(list_test[[i]][t,2]-list_test[[i]][t-1,2])/(list_test[[i]][t-1,2])
    
  }
}
#remove the first row because it contains prices
real_returns=lapply(real_returns, function(x) x[-1, ])


#Calculating predicted returns using the predicted price at time t and the real price at time (t-1)
predictions_returns=final_forecast_all
for(i in 1:length(list_test)){
  for(t in 2:length(list_test[[1]]$Price)){
    predictions_returns[[i]][t]=(final_forecast_all[[i]][t]-list_test[[i]]$Price[t-1])/(list_test[[i]]$Price[t-1])
  }
  
}

predictions_returns=lapply(predictions_returns, function(x) x[-1])

#Calculate the accuracy measures of the returns
rmse_all=mae_all=list()
for(i in 1:length(list_test)){
  rmse_all[[i]]=rmse(real_returns[[i]]$Price,predictions_returns[[i]])
  mae_all[[i]]=mae(real_returns[[i]]$Price,predictions_returns[[i]])
}
names(rmse_all)=names(mae_all)=names(list_test)

accuracy_returns=data.frame(stocks=names(list_test),RMSE=unlist(rmse_all),MAE=unlist(mae_all))
rownames(accuracy_returns)=names(list_test)


for(i in 1:length(list_test)){
  plot(real_returns[[i]]$Date,real_returns[[i]]$Price,type="l",ylim=range(c(real_returns[[i]]$Price),c(predictions_returns[[i]])),main=names(list_test)[[i]])
  lines(real_returns[[i]]$Date,predictions_returns[[i]],col="red")
}
dev.off()
