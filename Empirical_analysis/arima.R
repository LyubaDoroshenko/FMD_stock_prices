#This code models and predicts stock prices with the help of ARIMA model

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
 
library(anytime)
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





#Transforming the date into an appropriate format and subsetting until Dec 31, 2020

training=function(x){
  x=x[x$Date<="2020-12-31",]
}

test=function(x){
  x=x[x$Date>"2020-12-31",]
}
#Selecting the training period for stock prices
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

names(list_training_ts)=names_titles

#Selecting the test period for stock prices
list_test=lapply(list_entire_fin,test)




rmse=function(series,predictions)
{
  sqrt(mean((series-predictions)^2))#(series-predictions) or error/residuals
}

mae=function(series, predictions)
{
  mean(abs(series-predictions))
}

library("forecast")
#ARIMA
#auto.arima returns best ARIMA model according to either AIC, AICc or BIC value.
#The function conducts search over possible model within the provided order constraints.
#We choose the max order 10 of Autoregressive and Moving Average processes 
#with differences d=1, restricting search to non-seasonal models.
arima_fit=function(x){
  arima_fit=auto.arima(x,max.p=10,max.q=10,d=1,approximation=TRUE,seasonal=FALSE,stepwise=FALSE)#d=1
}

fitted_all=lapply(list_training_ts,arima_fit)

for(i in 1:length(list_training_ts)){
  print(names(list_training_ts)[[i]])
  print(fitted_all[[i]])
}

#Forecasting one unit at a time and refitting the ARIMA model at each time
#considering the real observations for the refit.
new_fitted_all=new_forecast_all=final_forecast_all=vector(mode="list",length=length(list_test))
first_forecast_all=c()
for(i in 1:length(list_test)){
  new_fitted_all[[i]]=list()
  new_forecast_all[[i]]=list()
  for(j in 1:(length(list_test$AAPL$Price))){
    new_fitted_all[[i]][[j]]=Arima(ts(c(fitted_all[[i]]$x,list_test[[i]]$Price[1:j]),start=decimal_date(as.Date(list_test[[i]]$Date[1])),frequency=250),model=fitted_all[[i]])
    new_forecast_all[[i]][[j]]=fitted(new_fitted_all[[i]][[j]])[(length(list_entire_fin[[i]]$Price)-length(list_test[[i]]$Price)+j)]
  }
  final_forecast_all[[i]]=ts(c(unlist(new_forecast_all[[i]])),start = decimal_date(as.Date("2021-01-04")), frequency = 250)
}

names(final_forecast_all)=names(list_test)


rmse_all=mae_all=arima_forecast=list()
for(i in (1:length(list_test))){
  rmse_all[[i]]=rmse(list_test[[i]]$Price,final_forecast_all[[i]])
  mae_all[[i]]=mae(list_test[[i]]$Price,final_forecast_all[[i]])
}
names(rmse_all)=names(mae_all)=names(list_test)

accuracy_arima=data.frame(stocks=names(list_test),RMSE=unlist(rmse_all),MAE=unlist(mae_all))
rownames(accuracy_arima)=names(list_test)


for(i in 1:length(list_training)){
  plot(list_test[[i]]$Date,list_test[[i]]$Price,type="l",ylim=range(c(list_test[[i]]$Price),c(final_forecast_all[[i]])),main=names(list_test)[[i]])
  lines(list_test[[i]]$Date,final_forecast_all[[i]],col="red")
}
dev.off()


#Calculating the returns of the real stock prices
real_returns=list_test
for(i in 1:length(list_test)){
  for(t in 2:length(list_test[[1]]$Price)){
    real_returns[[i]][t,2]=(list_test[[i]][t,2]-list_test[[i]][t-1,2])/(list_test[[i]][t-1,2])
    
  }
}
#Remove the first row because it contains prices since we can obtain the returns
#starting from time 2
real_returns=lapply(real_returns, function(x) x[-1, ])


#Calculating the returns from the predicted stock prices, where we consider the real prices at time (t-1) and predicted prices at time t.
predictions_returns=final_forecast_all
for(i in 1:length(list_test)){
  for(t in 2:length(list_test[[1]]$Price)){
    predictions_returns[[i]][t]=(final_forecast_all[[i]][t]-list_test[[i]][t-1,2])/(list_test[[i]][t-1,2])
  }
  
}

predictions_returns=lapply(predictions_returns, function(x) x[-1])


#Calculating the accuracy measures on the predicted returns
rmse_all=mae_all=list()
for(i in 1:length(list_test)){
  rmse_all[[i]]=rmse(real_returns[[i]]$Price,predictions_returns[[i]])
  mae_all[[i]]=mae(real_returns[[i]]$Price,predictions_returns[[i]])
}
names(rmse_all)=names(mae_all)=names(list_test)

accuracy=data.frame(stocks=names(list_test),RMSE=unlist(rmse_all),MAE=unlist(mae_all))
rownames(accuracy)=names(list_test)


for(i in 1:length(list_test)){
  plot(real_returns[[i]]$Date,real_returns[[i]]$Price,type="l",ylim=range(c(real_returns[[i]]$Price),c(predictions_returns[[i]])),main=names(list_test)[[i]])
  lines(real_returns[[i]]$Date,predictions_returns[[i]],col="red")
}
dev.off()
