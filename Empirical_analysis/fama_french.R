#We used this code to make fama&french model based predictions
# Set directory to the folder with the current file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#SP500_wrds.xlsx is the file that contains the daily prices of the analyzed stocks.
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

#Comparing the dates between SP500 and the other stocks, i.e. imputing the days present in SP500 but missing in the other stocks
#and removing the dates missing in SP500 but present in the other stocks.
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
    #imputing the prices present in SP500 but absent in the other stocks
    list_entire_fin[[i]]= SP_dates[[i]] %>%
      dplyr::select(Date) %>%
      left_join(list_entire_fin[[i]], by = 'Date') %>%
      mutate(Price = zoo::na.approx(Price))
    
    
  } }


list_entire_returns=list_entire_fin
for(i in 1:length(list_entire_fin)){
  for(t in 2:nrow(list_entire_fin[[i]])){
    list_entire_returns[[i]][t,2]=(list_entire_fin[[i]][t,2]-list_entire_fin[[i]][t-1,2])/(list_entire_fin[[i]][t-1,2])
  }
}
#Remove the first row because it contains prices at time=1 since we
#can calculate returns starting from t=2
list_entire_fin_returns=lapply(list_entire_returns, function(x) x[-1, ])



#Transforming the date into an appropriate format and subsetting until Dec 31, 2020
training=function(x){
  x=x[x$Date<="2020-12-31",]
}

test=function(x){
  x=x[x$Date>"2020-12-31",]
}
list_training=lapply(list_entire_fin,training)

list_test=lapply(list_entire_fin,test)

list_test_returns=lapply(list_entire_fin_returns,test)

list_training_returns_fin=lapply(list_entire_fin_returns,training)

names(list_test)=names(list_training)=names(list_test_returns)=names(list_training_returns_fin)=names(list_entire)

#The function to calculate the Root Mean Squared Errors 
rmse=function(series,predictions)
{
  sqrt(mean((series-predictions)^2)
}

library("forecast")

#FAMA and FRENCH
#Uploading Fama and French 5 factors
factors=read.csv("fama_french_5_factors.csv",header=T)
factors$Date=as.Date(strptime(factors$Date, format = "%Y%m%d"))
library("lubridate")
#5 Fama and French factors
factors$RF=ts(factors$RF,start = decimal_date(as.Date("2009-03-10")), frequency = 250)
factors$Mkt.RF=ts(factors$Mkt.RF,start = decimal_date(as.Date("2009-03-10")), frequency = 250)
factors$SMB=ts(factors$SMB,start = decimal_date(as.Date("2009-03-10")), frequency = 250)
factors$HML=ts(factors$HML,start = decimal_date(as.Date("2009-03-10")), frequency = 250)
factors$RMW=ts(factors$RMW,start = decimal_date(as.Date("2009-03-10")), frequency = 250)
factors$CMA=ts(factors$CMA,start = decimal_date(as.Date("2009-03-10")), frequency = 250)

factors_list=list(RF=data.frame(Date=factors$Date,Price=factors$RF),
                  Mkt.RF=data.frame(Date=factors$Date,Price=factors$Mkt.RF),
                  SMB=data.frame(Date=factors$Date,Price=factors$SMB),
                  HML=data.frame(Date=factors$Date,Price=factors$HML),
                  RMW=data.frame(Date=factors$Date,Price=factors$RMW),
                  CMA=data.frame(Date=factors$Date,Price=factors$CMA))
names(factors_list)=c("RF","Mkt.RF","SMB","HML","RMW","CMA")



library("dplyr")

new_list=list()
#I will impute the missing values for factors
for(i in 1:length(factors_list)){
  new_list[[i]]=data.frame(Date=anydate(list_entire_fin_returns$SP500$Date,tz="UTC"), Price=list_entire_fin_returns$SP500$Price)
  factors_list[[i]]= new_list[[i]] %>%
    dplyr::select(Date) %>%
    left_join(factors_list[[i]], by = 'Date') %>%
    mutate(Price = zoo::na.approx(Price))
}

factors_new=data.frame(Date=factors_list$RF$Date,RF=factors_list$RF$Price,
                       Mkt.RF=factors_list$Mkt.RF$Price,SMB=factors_list$SMB$Price,
                       HML=factors_list$HML$Price,RMW=factors_list$RMW$Price,
                       CMA=factors_list$CMA$Price)



factors_training=lapply(factors_list,training)
factors_test=lapply(factors_list,test)

#common dates
common_dates_training=subsets_common_training=list()
dates_not_in_common=list()
for(i in 1:length(list_training_returns_fin)){
  common_dates_training[[i]]=intersect(as.factor(factors_new$Date), as.factor(list_training_returns_fin[[i]]$Date))
  subsets_common_training[[i]]=subset(list_training_returns_fin[[i]],as.factor(Date) %in% common_dates_training[[i]])
  
}

setdiff(as.factor(list_training_returns_fin[[1]]$Date),as.factor(common_dates_training[[1]]))

#Subset the 5 factors for training period
#Convert returns to excess returns to use the Fama&French model
subset_RF_training=subset_Mkt.RF_training=subset_SMB_training=
  subset_HML_training=subset_RMW_training=subset_CMA_training=list()
for(i in 1:length(list_training_returns_fin)){
  subset_RF_training[[i]]=subset(factors_new$RF,as.factor(factors_new$Date) %in% common_dates_training[[i]])
  subset_Mkt.RF_training[[i]]=subset(factors_new$Mkt.RF,as.factor(factors_new$Date) %in% common_dates_training[[i]])
  subset_SMB_training[[i]]=subset(factors_new$SMB,as.factor(factors_new$Date) %in% common_dates_training[[i]])
  subset_HML_training[[i]]=subset(factors_new$HML,as.factor(factors_new$Date) %in% common_dates_training[[i]])
  subset_RMW_training[[i]]=subset(factors_new$RMW,as.factor(factors_new$Date) %in% common_dates_training[[i]])
  subset_CMA_training[[i]]=subset(factors_new$CMA,as.factor(factors_new$Date) %in% common_dates_training[[i]])
  
}
names(subset_RF_training)=names(subset_Mkt.RF_training)=names(subset_SMB_training)=
  names(subset_HML_training)=names(subset_RMW_training)=names( subset_CMA_training)=
  names(list_training_returns_fin)



###############################################
###############################################
#Considering the lagged factors with the maximum lag=5 for subsetting matters
training_subsets=excess_returns=subset_RF_training_new=list()
for(i in 1:length(list_training_returns_fin)){
  subset_RF_training_new[[i]]=subset_RF_training[[i]][6:length(subset_RF_training[[i]])]
  
  training_subsets[[i]]=data.frame(Mkt.RF_lag1=subset_Mkt.RF_training[[i]][5:(length(subset_Mkt.RF_training[[i]])-1)],
                                   SMB_lag1=subset_SMB_training[[i]][5:(length(subset_SMB_training[[i]])-1)],
                                   HML_lag1=subset_HML_training[[i]][5:(length(subset_HML_training[[i]])-1)],
                                   RMW_lag1=subset_RMW_training[[i]][5:(length(subset_RMW_training[[i]])-1)],
                                   CMA_lag1=subset_CMA_training[[i]][5:(length(subset_CMA_training[[i]])-1)],
                                   Mkt.RF_lag2=subset_Mkt.RF_training[[i]][4:(length(subset_Mkt.RF_training[[i]])-2)],
                                   SMB_lag2=subset_SMB_training[[i]][4:(length(subset_SMB_training[[i]])-2)],
                                   HML_lag2=subset_HML_training[[i]][4:(length(subset_HML_training[[i]])-2)],
                                   RMW_lag2=subset_RMW_training[[i]][4:(length(subset_RMW_training[[i]])-2)],
                                   CMA_lag2=subset_CMA_training[[i]][4:(length(subset_CMA_training[[i]])-2)],
                                   Mkt.RF_lag3=subset_Mkt.RF_training[[i]][3:(length(subset_Mkt.RF_training[[i]])-3)],
                                   SMB_lag3=subset_SMB_training[[i]][3:(length(subset_SMB_training[[i]])-3)],
                                   HML_lag3=subset_HML_training[[i]][3:(length(subset_HML_training[[i]])-3)],
                                   RMW_lag3=subset_RMW_training[[i]][3:(length(subset_RMW_training[[i]])-3)],
                                   CMA_lag3=subset_CMA_training[[i]][3:(length(subset_CMA_training[[i]])-3)],
                                   Mkt.RF_lag4=subset_Mkt.RF_training[[i]][2:(length(subset_Mkt.RF_training[[i]])-4)],
                                   SMB_lag4=subset_SMB_training[[i]][2:(length(subset_SMB_training[[i]])-4)],
                                   HML_lag4=subset_HML_training[[i]][2:(length(subset_HML_training[[i]])-4)],
                                   RMW_lag4=subset_RMW_training[[i]][2:(length(subset_RMW_training[[i]])-4)],
                                   CMA_lag4=subset_CMA_training[[i]][2:(length(subset_CMA_training[[i]])-4)],
                                   Mkt.RF_lag5=subset_Mkt.RF_training[[i]][1:(length(subset_Mkt.RF_training[[i]])-5)],
                                   SMB_lag5=subset_SMB_training[[i]][1:(length(subset_SMB_training[[i]])-5)],
                                   HML_lag5=subset_HML_training[[i]][1:(length(subset_HML_training[[i]])-5)],
                                   RMW_lag5=subset_RMW_training[[i]][1:(length(subset_RMW_training[[i]])-5)],
                                   CMA_lag5=subset_CMA_training[[i]][1:(length(subset_CMA_training[[i]])-5)]
  )                                
  
  excess_returns[[i]]=list_training_returns_fin[[i]]$Price[6:length(list_training_returns_fin[[i]]$Price)]-subset_RF_training_new[[i]]
  
}

library("leaps")
#Choosing the factors and the lagged factors to consider with the help of 
#model selection by exhaustive search
regsubsets_full=regsubsets_summary=list()
for(i in 1:length(list_training_returns_fin)){
  regsubsets_full[[i]] = regsubsets(excess_returns[[i]] ~ ., data = training_subsets[[i]], nvmax = 25,method="exhaustive")
  regsubsets_summary[[i]] = summary(regsubsets_full[[i]])
}
regsubsets_summary
plot(regsubsets_full[[i]],scale="bic")
#Considering Adjusted R2, Cp and BIC and AIC
aic=params=list()
for(i in 1:length(list_training_returns_fin)){
  aic[[i]]=matrix(NA,nrow=25,ncol=1)
  params[[i]]=seq(1:25)
  for(j in 1:25)
    aic[[i]][j]=regsubsets_summary[[i]]$bic[j]-params[[i]][j]*log(regsubsets_full[[i]]$nn)+
    2*params[[i]][j]
}
#Here I am adding the AIC
criteria=list()
for(i in 1:length(list_training_returns_fin)){
  criteria[[i]]=data.frame(AdjR2=which.max(regsubsets_summary[[i]]$adjr2),
                           Cp=which.min(regsubsets_summary[[i]]$cp),
                           BIC=which.min(regsubsets_summary[[i]]$bic),
                           AIC=which.min(as.vector(aic[[i]])))
  
}




#Factors for predictions
factors_pred=list()
RF_pred=list()
for(i in 1:length(list_training_returns_fin)){
  factors_pred[[i]]=list()
  RF_pred[[i]]=list()
  for(t in 1:length(list_test_returns[[1]]$Price)){
    RF_pred[[i]][[t]]=c(subset_RF_training[[i]],factors_test$RF$Price)[(length(subset_RF_training[[i]])+t-1)]
    factors_pred[[i]][[t]]=data.frame(Mkt.RF_lag1=c(subset_Mkt.RF_training[[i]],factors_test$Mkt.RF$Price)[(length(subset_Mkt.RF_training[[i]])+t-1)],
                                      SMB_lag1=c(subset_SMB_training[[i]],factors_test$SMB$Price)[(length(subset_SMB_training[[i]])+t-1)],
                                      HML_lag1=c(subset_HML_training[[i]],factors_test$HML$Price)[(length(subset_HML_training[[i]])+t-1)],
                                      RMW_lag1=c(subset_RMW_training[[i]],factors_test$RMW$Price)[(length(subset_RMW_training[[i]])+t-1)],
                                      CMA_lag1=c(subset_CMA_training[[i]],factors_test$CMA$Price)[(length(subset_CMA_training[[i]])+t-1)],
                                      Mkt.RF_lag2=c(subset_Mkt.RF_training[[i]],factors_test$Mkt.RF$Price)[(length(subset_Mkt.RF_training[[i]])+t-2)],
                                      SMB_lag2=c(subset_SMB_training[[i]],factors_test$SMB$Price)[(length(subset_SMB_training[[i]])+t-2)],
                                      HML_lag2=c(subset_HML_training[[i]],factors_test$HML$Price)[(length(subset_HML_training[[i]])+t-2)],
                                      RMW_lag2=c(subset_RMW_training[[i]],factors_test$RMW$Price)[(length(subset_RMW_training[[i]])+t-2)],
                                      CMA_lag2=c(subset_CMA_training[[i]],factors_test$CMA$Price)[(length(subset_CMA_training[[i]])+t-2)],
                                      Mkt.RF_lag3=c(subset_Mkt.RF_training[[i]],factors_test$Mkt.RF$Price)[(length(subset_Mkt.RF_training[[i]])+t-3)],
                                      SMB_lag3=c(subset_SMB_training[[i]],factors_test$SMB$Price)[(length(subset_SMB_training[[i]])+t-3)],
                                      HML_lag3=c(subset_HML_training[[i]],factors_test$HML$Price)[(length(subset_HML_training[[i]])+t-3)],
                                      RMW_lag3=c(subset_RMW_training[[i]],factors_test$RMW$Price)[(length(subset_RMW_training[[i]])+t-3)],
                                      CMA_lag3=c(subset_CMA_training[[i]],factors_test$CMA$Price)[(length(subset_CMA_training[[i]])+t-3)],
                                      Mkt.RF_lag4=c(subset_Mkt.RF_training[[i]],factors_test$Mkt.RF$Price)[(length(subset_Mkt.RF_training[[i]])+t-1)],
                                      SMB_lag4=c(subset_SMB_training[[i]],factors_test$SMB$Price)[(length(subset_SMB_training[[i]])+t-4)],
                                      HML_lag4=c(subset_HML_training[[i]],factors_test$HML$Price)[(length(subset_HML_training[[i]])+t-4)],
                                      RMW_lag4=c(subset_RMW_training[[i]],factors_test$RMW$Price)[(length(subset_RMW_training[[i]])+t-4)],
                                      CMA_lag4=c(subset_CMA_training[[i]],factors_test$CMA$Price)[(length(subset_CMA_training[[i]])+t-4)],
                                      Mkt.RF_lag5=c(subset_Mkt.RF_training[[i]],factors_test$Mkt.RF$Price)[(length(subset_Mkt.RF_training[[i]])+t-1)],
                                      SMB_lag5=c(subset_SMB_training[[i]],factors_test$SMB$Price)[(length(subset_SMB_training[[i]])+t-5)],
                                      HML_lag5=c(subset_HML_training[[i]],factors_test$HML$Price)[(length(subset_HML_training[[i]])+t-5)],
                                      RMW_lag5=c(subset_RMW_training[[i]],factors_test$RMW$Price)[(length(subset_RMW_training[[i]])+t-5)],
                                      CMA_lag5=c(subset_CMA_training[[i]],factors_test$CMA$Price)[(length(subset_CMA_training[[i]])+t-5)])
  }
  #RF will be taken at time t-1
}



library("ISLR")

# Iterate over each stock
coef_adj_R=coef_Cp=coef_bic=coef_aic=list()
test_mat=pred_adj_R=pred_Cp=pred_bic=pred_aic=list()

for(i in 1:length(list_training_returns_fin)){
  # Extract the vector of predictors in the best fit model on i predictors
  coef_adj_R[[i]] = coef(regsubsets_full[[i]], id = criteria[[i]]$AdjR2)
  coef_Cp[[i]]=coef(regsubsets_full[[i]], id = criteria[[i]]$Cp)
  coef_bic[[i]]=coef(regsubsets_full[[i]], id = criteria[[i]]$BIC)
  coef_aic[[i]]=coef(regsubsets_full[[i]], id = criteria[[i]]$AIC)
  
  #here for each t
  test_mat[[i]]=pred_adj_R[[i]]=pred_Cp[[i]]=pred_bic[[i]]=pred_aic[[i]]=list()
  for(t in 1:length(list_test_returns[[1]]$Price)){
    test_mat[[i]][[t]]=model.matrix(excess_returns[[i]][1]~ .,data=factors_pred[[i]][[t]])
    # Make predictions using matrix multiplication of the test matrix and the coefficients vector
    pred_adj_R[[i]][[t]]=test_mat[[i]][[t]][,names(coef_adj_R[[i]])]%*%coef_adj_R[[i]]
    pred_Cp[[i]][[t]]=test_mat[[i]][[t]][,names(coef_Cp[[i]])]%*%coef_Cp[[i]]
    pred_bic[[i]][[t]]=test_mat[[i]][[t]][,names(coef_bic[[i]])]%*%coef_bic[[i]]
    pred_aic[[i]][[t]]=test_mat[[i]][[t]][,names(coef_aic[[i]])]%*%coef_aic[[i]]
  }
  
}




#Predictions
predictions_lag=list()
predictions_lag_returns_adj_R=predictions_lag_returns_Cp=predictions_lag_returns_bic=predictions_lag_returns_aic=list()
ff_fin_predictions_adj_R=ff_fin_predictions_Cp=ff_fin_predictions_bic=ff_fin_predictions_aic=list()
rmse_returns_adj_R=rmse_returns_Cp=rmse_returns_bic=rmse_returns_aic=list()
for(i in 1:length(list_training_returns_fin)){
  predictions_lag_returns_adj_R[[i]]=predictions_lag_returns_Cp[[i]]=
    predictions_lag_returns_bic[[i]]=predictions_lag_returns_aic[[i]]=list()
  for(t in 1:length(list_test_returns[[1]]$Price)){
    #now calculating the returns adding back the risk free rates
    predictions_lag_returns_adj_R[[i]][[t]]=pred_adj_R[[i]][[t]]+RF_pred[[i]][[t]]
    predictions_lag_returns_Cp[[i]][[t]]=pred_Cp[[i]][[t]]+RF_pred[[i]][[t]]
    predictions_lag_returns_bic[[i]][[t]]=pred_bic[[i]][[t]]+RF_pred[[i]][[t]]
    predictions_lag_returns_aic[[i]][[t]]=pred_aic[[i]][[t]]+RF_pred[[i]][[t]]
  }
  ff_fin_predictions_adj_R[[i]]=unlist(predictions_lag_returns_adj_R[[i]])
  ff_fin_predictions_Cp[[i]]=unlist(predictions_lag_returns_Cp[[i]])
  ff_fin_predictions_bic[[i]]=unlist(predictions_lag_returns_bic[[i]])
  ff_fin_predictions_aic[[i]]=unlist(predictions_lag_returns_aic[[i]])
  
  
  rmse_returns_adj_R[[i]] = rmse(ff_fin_predictions_adj_R[[i]],list_test_returns[[i]]$Price)
  rmse_returns_Cp[[i]]=rmse(ff_fin_predictions_Cp[[i]],list_test_returns[[i]]$Price)
  rmse_returns_bic[[i]]=rmse(ff_fin_predictions_bic[[i]],list_test_returns[[i]]$Price)
  rmse_returns_aic[[i]]=rmse(ff_fin_predictions_aic[[i]],list_test_returns[[i]]$Price)
}

accuracy_returns=data.frame(stocks=names(list_test),RMSE_adj_R=unlist(rmse_returns_adj_R),
                            RMSE_Cp=unlist(rmse_returns_Cp),
                            RMSE_bic=unlist(rmse_returns_bic),
                            RMSE_aic=unlist(rmse_returns_aic))
rownames(accuracy_returns)=names(list_test)

for(i in 1:length(list_test)){
  plot(list_test_returns[[i]]$Date,list_test_returns[[i]]$Price,type="l")
  lines(list_test_returns[[i]]$Date,ff_fin_predictions_adj_R[[i]],type="l",col="red")
}

for(i in 1:length(list_test)){
  plot(list_test_returns[[i]]$Date,list_test_returns[[i]]$Price,type="l")
  lines(list_test_returns[[i]]$Date,ff_fin_predictions_Cp[[i]],type="l",col="red")
}

       for(i in 1:length(list_test)){
  plot(list_test_returns[[i]]$Date,list_test_returns[[i]]$Price,type="l")
  lines(list_test_returns[[i]]$Date,ff_fin_predictions_bic[[i]],type="l",col="red")
}




#Going back to prices using the real prices at t-1
ff_fin_pred_prices_adj_R=ff_fin_pred_prices_Cp=ff_fin_pred_prices_bic=list()
true_prices=vector(mode="list",length=length(list_training_returns_fin))
rmse_adj_R=rmse_Cp=rmse_bic=list()
for(i in 1:length(list_training_returns_fin)){
  ff_fin_pred_prices_adj_R[[i]]=ff_fin_pred_prices_Cp[[i]]=ff_fin_pred_prices_bic[[i]]=matrix(NA,nrow=length(list_test[[i]]$Price),ncol=1)
  true_prices[[i]]=matrix(NA,nrow=length(list_test[[i]]$Price),ncol=1)
  true_prices[[i]][1,1]=tail(list_training[[i]]$Price,1)
  ff_fin_pred_prices_adj_R[[i]][1,1]=ff_fin_predictions_adj_R[[i]][1]*true_prices[[i]][1,1]+true_prices[[i]][1,1]
  ff_fin_pred_prices_Cp[[i]][1,1]=ff_fin_predictions_Cp[[i]][1]*true_prices[[i]][1,1]+true_prices[[i]][1,1]
  ff_fin_pred_prices_bic[[i]][1,1]=ff_fin_predictions_bic[[i]][1]*true_prices[[i]][1,1]+true_prices[[i]][1,1]
  
  for(t in 2:(length(list_test[[i]]$Price))){
    true_prices[[i]][t,1]=list_test[[i]]$Price[t-1]
    ff_fin_pred_prices_adj_R[[i]][t,1]=ff_fin_predictions_adj_R[[i]][t]*true_prices[[i]][t,1]+
      true_prices[[i]][t,1]
    ff_fin_pred_prices_Cp[[i]][t,1]=ff_fin_predictions_Cp[[i]][t]*true_prices[[i]][t,1]+
      true_prices[[i]][t,1]
    ff_fin_pred_prices_bic[[i]][t,1]=ff_fin_predictions_bic[[i]][t]*true_prices[[i]][t,1]+
      true_prices[[i]][t,1]
  }
  
  # Calculate the RMSE
  rmse_adj_R[[i]] = rmse(ff_fin_pred_prices_adj_R[[i]],list_test[[i]]$Price)
  rmse_Cp[[i]]=rmse(ff_fin_pred_prices_Cp[[i]],list_test[[i]]$Price)
  rmse_bic[[i]]=rmse(ff_fin_pred_prices_bic[[i]],list_test[[i]]$Price)
}
names(ff_fin_pred_prices_adj_R)=names(ff_fin_pred_prices_Cp)=names(ff_fin_pred_prices_bic)=
  names(list_training_returns_fin)


names(rmse_adj_R)=names(rmse_Cp)=names(rmse_bic)=names(list_test)


accuracy_ff_adj_R=data.frame(stocks=names(list_test),RMSE=unlist(rmse_adj_R))
rownames(accuracy_ff_adj_R)=names(list_test)

accuracy_ff_Cp=data.frame(stocks=names(list_test),RMSE=unlist(rmse_Cp))
rownames(accuracy_ff_Cp)=names(list_test)

accuracy_ff_bic=data.frame(stocks=names(list_test),RMSE=unlist(rmse_bic))
rownames(accuracy_ff_bic)=names(list_test)

