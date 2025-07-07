library("renv")
renv::use_python(python="/home/lyuba/myenv/bin/python",
                 name="/home/lyuba/myenv",type="virtualenv")


library("reticulate")

library(keras)
library(tensorflow)

#------------------------------------------------------------------------------------------

library(reticulate)

# Set directory to the folder with the current file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


set.seed(12345)
library("readxl")
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


library(anytime)
#Comparing the dates between SP500 and the other stocks
SP_dates=list()
only_sp_dates=list()#here are the dates that are in SP500 but not in the other stocks
library(zoo)
library(dplyr)
library(anytime)
for(i in 1:length(list_entire_fin)){
  if(length(anydate(list_entire_fin[[i]]$Date,tz="UTC"))!=length(anydate(list_entire_fin$SP500$Date,tz="UTC"))){
    SP_dates[[i]]=data.frame(Date=((list_entire_fin$SP500[anydate(list_entire_fin$SP500$Date,tz="UTC")>=anydate(list_entire_fin[[i]]$Date,tz="UTC")[1],])$Date),
                             Price=((list_entire_fin$SP500[anydate(list_entire_fin$SP500$Date,tz="UTC")>=anydate(list_entire_fin[[i]]$Date,tz="UTC")[1],])$Price))
    
    #removing the prices for the date not present in SP500 but present in the other stocks
    list_entire_fin[[i]]=list_entire_fin[[i]][which(list_entire_fin[[i]]$Date%in%SP_dates[[i]]$Date),]
    #imputig the prices present in Sp500 but absent in the other stocks
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

##Here we consider the first 5 two years as a training period
##for the nested cross validation. Each year has approximately 250 observations.
training_cv_first=list()
for(i in 1:length(list_entire_fin)){
  training_cv_first[[i]]=list_entire_fin[[i]][1:(250*5),]
}
##How many intervals we will consider for nested cross validation
chunks=list()
rest=list()

##we predict 12 months
for(i in 1:length(list_entire_fin)){
  chunks[[i]]=(1:floor((nrow(list_training[[i]])-nrow(training_cv_first[[i]]))/250))
  rest[[i]]=nrow(list_training[[i]])-(tail(chunks[[i]],1)*250+nrow(training_cv_first[[i]]))
  #the last chunk to predict will be shorter
}

to_predict=vector(mode="list",length=length(list_entire_fin))
for(i in 1:length(list_entire_fin)){
  to_predict[[i]]=list()
  for(j in 1:length(chunks[[i]])){
    to_predict[[i]][[j]]=250
  }
  to_predict[[i]][[(length(chunks[[i]])+1)]]=rest[[i]]
  to_predict[[i]]=unlist(to_predict[[i]])
}

training_cv=test_cv=list()
for(i in 1:length(list_entire_fin)){
  training_cv[[i]]=list()
  test_cv[[i]]=list()
  training_cv[[i]][[1]]=training_cv_first[[i]]$Price
  for(j in 1:length(chunks[[i]])){
    training_cv[[i]][[j+1]]=c(training_cv[[i]][[1]],
                              list_training[[i]]$Price[(length(training_cv[[i]][[1]])+1):
                                                         ((length(training_cv[[i]][[1]])+1)+chunks[[i]][j]*250-1)])
    test_cv[[i]][[j]]=list_training[[i]]$Price[(length(training_cv[[i]][[j]])+1):((length(training_cv[[i]][[j]])+1)+1*250-1)]
  }
  test_cv[[i]][[(length(chunks[[i]])+1)]]=list_training[[i]]$Price[(length(list_training[[i]]$Price)-rest[[i]]+1):
                                                                     length(list_training[[i]]$Price)]
  
}


rmse=function(series,predictions)
{
  sqrt(mean((series-predictions)^2))#(series-predictions) or error/residuals
}

mae=function(series, predictions)
{
  mean(abs(series-predictions))
}


#we train the lstm on the lists of training_cv and then we predict the number of
#observations in the list to_predict
library(anytime)
#Taking the first differences
Series=diffed=list()
for(i in 1:length(list_entire_fin)){
  Series[[i]]=diffed[[i]]=list()
  for(j in 1:(length(training_cv[[i]])-1)){
    Series[[i]][[j]]=training_cv[[i]][[j+1]]#training plus test
  }
  Series[[i]][[length(training_cv[[i]])]]=list_training[[i]]$Price 
  #the last test period is smaller than half a year
  diffed[[i]]=lapply(Series[[i]],diff)
}



N_to_rescale=c()
n_to_rescale=list()
for(i in 1:length(list_entire_fin)){
  N_to_rescale[i]=length(diffed[[i]][[length(diffed[[i]])]])
  n_to_rescale[[i]]=list()
  for(j in 1:length(to_predict[[i]])){
    n_to_rescale[[i]][[j]]=length(training_cv[[i]][[j]])
  }
  n_to_rescale[[i]]=unlist(n_to_rescale[[i]])
}


lag_transform <- function(x, k= 1){
  
  lagged =  c(rep(NA, k), x[1:(length(x)-k)])
  DF = as.data.frame(cbind(lagged, x))
  colnames(DF) <- c( paste0('x-', k), 'x')
  DF[is.na(DF)] <- 0
  return(DF)
}


split_sequence=function(sequence,n_steps){
  matrix_x_y=matrix(NA,nrow=(length(sequence)-n_steps),ncol=(n_steps+1))
  matrix_x=matrix(NA,nrow=(length(sequence)-n_steps),ncol=n_steps)
  matrix_y=matrix(NA,nrow=(length(sequence)-n_steps),ncol=1)
  
  end_ix=c()
  seq_X=seq_Y=list()
  #split a univariate sequence into samples
  for(i in 1:length(sequence)){
    # find the end of this pattern
    end_ix[i]=i+n_steps
    # check if we are beyond the sequence
    if(end_ix[i] > length(sequence)){
      break
    }
    # gather input and output parts of the pattern
    seq_X[[i]]=sequence[i:(end_ix[i]-1)]
    seq_Y[[i]]=sequence[end_ix[i]]
  }
  matrix_x=matrix(unlist(seq_X),ncol=n_steps,byrow=TRUE)
  
  
  matrix_y=matrix(unlist(seq_Y),ncol=1,byrow=TRUE)
  
  
  matrix_x_y=cbind(matrix_x,matrix_y)
  
  
  my_list=list(x=matrix_x,y=matrix_y,x_y=matrix_x_y)
  return(my_list)
}


n_steps=c(40,50,60,20)

supervised=list()
for(i in 1:length(diffed)){
  supervised[[i]]=list()
  for(s in 1:length(n_steps)){
    supervised[[i]][[s]]=list()
    for(j in 1:length(to_predict[[i]])){
      supervised[[i]][[s]][[j]]=split_sequence(diffed[[i]][[j]],n_steps=n_steps[s])$x_y
    }
  }
}
#Split dataset into training and testing sets
N=n=c()
for(i in 1:length(supervised)){
  N[[i]]=n[[i]]=list()
  for(s in 1:length(n_steps)){
    N[[i]][[s]]=n[[i]][[s]]=list()
    for(j in 1:length(to_predict[[i]])){
      N[[i]][[s]][[j]] = nrow(supervised[[i]][[s]][[j]])
      n[[i]][[s]][[j]]=N[[i]][[s]][[j]]-to_predict[[i]][[j]]
    }
  }
}


train=test=list()
for(i in 1:length(supervised)){
  train[[i]]=test[[i]]=list()
  for(s in 1:length(n_steps)){
    train[[i]][[s]]=test[[i]][[s]]=list()
    for(j in 1:length(to_predict[[i]])){
      train[[i]][[s]][[j]] = supervised[[i]][[s]][[j]][1:n[[i]][[s]][[j]], ]
      test[[i]][[s]][[j]]=supervised[[i]][[s]][[j]][(n[[i]][[s]][[j]]+1):N[[i]][[s]][[j]],]
    }
  }
}


#Normalize the data
# scale data
scale_data = function(train, test, feature_range = c(0, 1)) {
  x = train
  fr_min = feature_range[1]
  fr_max = feature_range[2]
  std_train = ((x - min(x) ) / (max(x) - min(x)  ))
  std_test  = ((test - min(x) ) / (max(x) - min(x)  ))
  
  scaled_train = std_train *(fr_max -fr_min) + fr_min
  scaled_test = std_test *(fr_max -fr_min) + fr_min
  
  return( list(scaled_train = scaled_train, scaled_test = scaled_test ,
               scaler= c(min =min(x), max = max(x))) )
  
}


Scaled=y_train=x_train=y_test=x_test=list()
for(i in 1:length(list_entire_fin)){
  Scaled[[i]]=y_train[[i]]=x_train[[i]]=y_test[[i]]=x_test[[i]]=list()
  for(s in 1:length(n_steps)){
    Scaled[[i]][[s]]=y_train[[i]][[s]]=x_train[[i]][[s]]=y_test[[i]][[s]]=x_test[[i]][[s]]=list()
    for(j in 1:length(to_predict[[i]])){
      Scaled[[i]][[s]][[j]] = scale_data(train[[i]][[s]][[j]],test[[i]][[s]][[j]], c(-1, 1))
      y_train[[i]][[s]][[j]] = Scaled[[i]][[s]][[j]]$scaled_train[, (n_steps[s]+1)]
      x_train[[i]][[s]][[j]] = Scaled[[i]][[s]][[j]]$scaled_train[, (1:n_steps[s])]
      
      y_test[[i]][[s]][[j]]=Scaled[[i]][[s]][[j]]$scaled_test[, (n_steps[s]+1)]
      x_test[[i]][[s]][[j]]=Scaled[[i]][[s]][[j]]$scaled_test[,(1:n_steps[s])]
      
    }
  }
}

#The following code will be required to revert the predicted values to the original scale.
## inverse-transform
invert_scaling = function(scaled, scaler, feature_range = c(0, 1)){
  min = scaler[1]
  max = scaler[2]
  t = length(scaled)
  mins = feature_range[1]
  maxs = feature_range[2]
  inverted_dfs = numeric(t)
  
  for( i in 1:t){
    X = (scaled[i]- mins)/(maxs - mins)
    rawValues = X *(max - min) + min
    inverted_dfs[i] <- rawValues
  }
  return(inverted_dfs)
}

#############################
#############################
#############################
#############################
#############################
#Define the model
# Reshape the input to 3-dim
for(i in 1:length(list_entire_fin)){
  for(s in 1:length(n_steps)){
    for(j in 1:length(to_predict[[i]])){
      dim(x_train[[i]][[s]][[j]])=c(nrow(x_train[[i]][[s]][[j]]),ncol(x_train[[i]][[s]][[j]]),1)
    }
  }
}


X_shape2=X_shape3=list()
for(i in 1:length(list_entire_fin)){
  X_shape2[[i]]=X_shape3[[i]]=list()
  for(s in 1:length(n_steps)){
    X_shape2[[i]][[s]]=X_shape3[[i]][[s]]=list()
    for(j in 1:length(to_predict[[i]])){
      X_shape2[[i]][[s]][[j]]=dim(x_train[[i]][[s]][[j]])[2]
      X_shape3[[i]][[s]][[j]]=dim(x_train[[i]][[s]][[j]])[3]
    }
  }
}



batch_size = 32                
x_train_new=x_train
y_train_new=y_train
for(i in 1:length(list_entire_fin)){
  for(s in 1:length(n_steps)){
    for(j in 1:length(to_predict[[i]])){
      if(nrow(x_train[[i]][[s]][[j]])%%batch_size!=0){
        x_train_new[[i]][[s]][[j]]=x_train[[i]][[s]][[j]][-(1:(nrow(x_train[[i]][[s]][[j]])%%batch_size)),,1]
        y_train_new[[i]][[s]][[j]]=y_train[[i]][[s]][[j]][-(1:(length(y_train[[i]][[s]][[j]])%%batch_size))]
      }
      else{
        x_train_new[[i]][[s]][[j]]=x_train[[i]][[s]][[j]]
        y_train_new[[i]][[s]][[j]]=y_train[[i]][[s]][[j]]
      }
    }
  }
}





set.seed(12345)
###Here I am adding a dropout grid
dropouts = seq(0.1,0.9,0.1)
x_test2=x_test
for(i in 1:length(list_test)){
  for(s in 1:length(n_steps)){
    for(j in 1:length(to_predict[[i]])){
      x_test2[[i]][[s]][[j]]=rbind(tail(x_train[[i]][[s]][[j]][,,1],(batch_size-1)),x_test[[i]][[s]][[j]])
    }
  }
}

L=list()
for(i in 1:length(list_entire_fin)){
  L[[i]]=list()
  for(j in 1:length(to_predict[[i]])){
    L[[i]][[j]]=c(seq(1:to_predict[[i]][j]))
  }
}


###TO NOT RUN
###In this part of the code we perform a nested cross-validation.
###This part of the code is particularly computationally expensive. Uncomment this part of the code in order to perform a nested cross-validation.
#scaler=predictions=X=yhat_original=yhat_no_diff=yhat=list()
#diff_wise=sums=rmse_fin=sums_fin=list()
#for (i in 1:length(list_entire_fin)){
#  scaler[[i]]=predictions[[i]]=X[[i]]=yhat_original[[i]]=yhat_no_diff[[i]]=yhat[[i]]=list()
#  diff_wise[[i]]=sums[[i]]=sums_fin[[i]]=rmse_fin[[i]]=list()
#  for(s in 1:length(n_steps)){
#    scaler[[i]][[s]]=predictions[[i]][[s]]=X[[i]][[s]]=yhat_original[[i]][[s]]=yhat_no_diff[[i]][[s]]=yhat[[i]][[s]]=list()
#    diff_wise[[i]][[s]]=sums[[i]][[s]]=sums_fin[[i]][[s]]=rmse_fin[[i]][[s]]=list()
#    for(d in 1:length(dropouts)){
#      predictions[[i]][[s]][[d]]=list()
#      yhat_original[[i]][[s]][[d]]=yhat_no_diff[[i]][[s]][[d]]=yhat[[i]][[s]][[d]]=list()
#      diff_wise[[i]][[s]][[d]]=sums[[i]][[s]][[d]]=rmse_fin[[i]][[s]][[d]]=list()
#      for(j in 1:length(to_predict[[i]])){#
#        scaler[[i]][[s]][[j]] = Scaled[[i]][[s]][[j]]$scaler
#        predictions[[i]][[s]][[d]][[j]]=X[[i]][[s]][[j]]=list()
#        yhat_original[[i]][[s]][[d]][[j]]=yhat_no_diff[[i]][[s]][[d]][[j]]=yhat[[i]][[s]][[d]][[j]]=list()
#        
#        model <- keras_model_sequential()
#        model %>%
#          layer_lstm(
#            units=10,
#            batch_input_shape = c(batch_size,
#                                  X_shape2[[i]][[s]][[j]], X_shape3[[i]][[s]][[j]]),
#            stateful = TRUE,
#            dropout = dropouts[[d]]
#          ) %>%
#          layer_dense(units=10)%>%#new hidden layer with 10 neurons
#          layer_dense(units = 1, activation = "linear")   
#        #Compile the model
#        model %>% compile(
#          loss = 'mean_squared_error',
#          optimizer = optimizer_adam(learning_rate = 0.02),
#          metrics = c("MSE")
#        )
#        
#        model %>% fit(
#          x_train_new[[i]][[s]][[j]],
#          y_train_new[[i]][[s]][[j]],
#          epochs = 500,
#          batch_size = batch_size,
#          verbose = 1,
#          shuffle = FALSE,
#          callbacks=list(callback_early_stopping(monitor="MSE",
#                                                 min_delta=0,
#                                                 patience=20,
#                                                 verbose=0,
#                                                 mode="min",
#                                                 restore_best_weights = TRUE)
#
#          ))
#        
#        for(l in 1:length(L[[i]][[j]])){
#          X[[i]][[s]][[j]][[l]]=x_test2[[i]][[s]][[j]][(l:(l+(batch_size-1))),]
#          dim(X[[i]][[s]][[j]][[l]]) = c(batch_size,n_steps[s],1)
#          yhat_original[[i]][[s]][[d]][[j]][[l]] = tail((model %>% predict(X[[i]][[s]][[j]][[l]],batch_size=batch_size)),1)#, batch_size=batch_size
#          #invert scaling
#          yhat_no_diff[[i]][[s]][[d]][[j]][[l]] = invert_scaling(yhat_original[[i]][[s]][[d]][[j]][[l]], scaler[[i]][[s]][[j]],  c(-1, 1))
#          # invert differencing
#          yhat[[i]][[s]][[d]][[j]][[l]]  = yhat_no_diff[[i]][[s]][[d]][[j]][[l]] + Series[[i]][[j]][(n_to_rescale[[i]][j]+l)]#we consider Xt-1 from 2976:3478;the total
#          # store                                                           #number of training+test is 3479. n=2975.  
#          predictions[[i]][[s]][[d]][[j]][[l]] <- yhat[[i]][[s]][[d]][[j]][[l]]
#        }#closing l loop
#        predictions[[i]][[s]][[d]][[j]]=unlist(predictions[[i]][[s]][[d]][[j]])
#        
#        diff_wise[[i]][[s]][[d]][[j]]=(predictions[[i]][[s]][[d]][[j]]-test_cv[[i]][[j]])^2
#        sums[[i]][[s]][[d]][[j]]=Reduce("+",diff_wise[[i]][[s]][[d]][[j]])
#      }#closing j loop
#      sums_fin[[i]][[s]][[d]]=Reduce("+",sums[[i]][[s]][[d]])
#      rmse_fin[[i]][[s]][[d]]=sqrt(sums_fin[[i]][[s]][[d]]/sum(to_predict[[i]]))
#      
#    }#closing d loop
#  }#closing s loop
#}#closing i loop

#rmse_matrix=rmse_min=list()
#rmse_min_matrix=matrix(NA,nrow=length(list_entire_fin),ncol=4)
#for(i in 1:length(list_entire_fin)){
#  rmse_matrix[[i]]=matrix(NA,nrow=length(n_steps),ncol=length(dropouts))
#  rmse_min[[i]]=matrix(NA,nrow=1,ncol=2)
#  for(s in 1:length(n_steps)){
#    for(d in 1:length(dropouts)){
#      rmse_matrix[[i]][s,d]=rmse_fin[[i]][[s]][[d]]
#    }
#    
#  }
#  rmse_min[[i]][1,]=which(rmse_matrix[[i]]==min(rmse_matrix[[i]]),arr.ind=TRUE)[1,]    
#  #the first element is s, i.e. timesteps, corresponds to row of the min value of
#  #matrix rmse_matrix
#  #the second element is d, i.e. dropouts, corresponds to column of the min value of
#  #matrix rmse_matrix
#  rmse_min_matrix[i,]=c(names(list_entire_fin)[i],rmse_matrix[[i]][rmse_min[[i]]],n_steps[rmse_min[[i]][1,1]],dropouts[rmse_min[[i]][1,2]])
#}



#rmse_min_matrix=data.frame(rmse_min_matrix)
#colnames(rmse_min_matrix)=c("Stocks","RMSE","timesteps","dropout")


#rmse_dropout_01_40=rmse_dropout_02_40=rmse_dropout_03_40=rmse_dropout_04_40=
#  rmse_dropout_05_40=rmse_dropout_06_40=rmse_dropout_07_40=rmse_dropout_08_40=
#  rmse_dropout_09_40=list()

#rmse_dropout_01_50=rmse_dropout_02_50=rmse_dropout_03_50=rmse_dropout_04_50=
#  rmse_dropout_05_50=rmse_dropout_06_50=rmse_dropout_07_50=rmse_dropout_08_50=
#  rmse_dropout_09_50=list()

#rmse_dropout_01_60=rmse_dropout_02_60=rmse_dropout_03_60=rmse_dropout_04_60=
#  rmse_dropout_05_60=rmse_dropout_06_60=rmse_dropout_07_60=rmse_dropout_08_60=
#  rmse_dropout_09_60=list()

#rmse_dropout_01_20=rmse_dropout_02_20=rmse_dropout_03_20=rmse_dropout_04_20=
#  rmse_dropout_05_20=rmse_dropout_06_20=rmse_dropout_07_20=rmse_dropout_08_20=
#  rmse_dropout_09_20=list()

#for(i in 1:length(list_entire_fin)){
#  rmse_dropout_01_40[[i]]=rmse_dropout_02_40[[i]]=rmse_dropout_03_40[[i]]=rmse_dropout_04_40[[i]]=
#    rmse_dropout_05_40[[i]]=rmse_dropout_06_40[[i]]=rmse_dropout_07_40[[i]]=rmse_dropout_08_40[[i]]=
#    rmse_dropout_09_40[[i]]=list()
  
#  rmse_dropout_01_50[[i]]=rmse_dropout_02_50[[i]]=rmse_dropout_03_50[[i]]=rmse_dropout_04_50[[i]]=
#    rmse_dropout_05_50[[i]]=rmse_dropout_06_50[[i]]=rmse_dropout_07_50[[i]]=rmse_dropout_08_50[[i]]=
#    rmse_dropout_09_50[[i]]=list()
  
#  rmse_dropout_01_60[[i]]=rmse_dropout_02_60[[i]]=rmse_dropout_03_60[[i]]=rmse_dropout_04_60[[i]]=
#    rmse_dropout_05_60[[i]]=rmse_dropout_06_60[[i]]=rmse_dropout_07_60[[i]]=rmse_dropout_08_60[[i]]=
#    rmse_dropout_09_60[[i]]=list()
  
#  rmse_dropout_01_20[[i]]=rmse_dropout_02_20[[i]]=rmse_dropout_03_20[[i]]=rmse_dropout_04_20[[i]]=
#    rmse_dropout_05_20[[i]]=rmse_dropout_06_20[[i]]=rmse_dropout_07_20[[i]]=rmse_dropout_08_20[[i]]=
#    rmse_dropout_09_20[[i]]=list()
  
  
#  rmse_dropout_01_40[[i]]=rmse_fin[[i]][[1]][[1]]
#  rmse_dropout_02_40[[i]]=rmse_fin[[i]][[1]][[2]]
#  rmse_dropout_03_40[[i]]=rmse_fin[[i]][[1]][[3]]
#  rmse_dropout_04_40[[i]]=rmse_fin[[i]][[1]][[4]]
#  rmse_dropout_05_40[[i]]=rmse_fin[[i]][[1]][[5]]
#  rmse_dropout_06_40[[i]]=rmse_fin[[i]][[1]][[6]]
#  rmse_dropout_07_40[[i]]=rmse_fin[[i]][[1]][[7]]
#  rmse_dropout_08_40[[i]]=rmse_fin[[i]][[1]][[8]]
#  rmse_dropout_09_40[[i]]=rmse_fin[[i]][[1]][[9]]
  
#  rmse_dropout_01_50[[i]]=rmse_fin[[i]][[2]][[1]]
#  rmse_dropout_02_50[[i]]=rmse_fin[[i]][[2]][[2]]
#  rmse_dropout_03_50[[i]]=rmse_fin[[i]][[2]][[3]]
#  rmse_dropout_04_50[[i]]=rmse_fin[[i]][[2]][[4]]
#  rmse_dropout_05_50[[i]]=rmse_fin[[i]][[2]][[5]]
#  rmse_dropout_06_50[[i]]=rmse_fin[[i]][[2]][[6]]
#  rmse_dropout_07_50[[i]]=rmse_fin[[i]][[2]][[7]]
#  rmse_dropout_08_50[[i]]=rmse_fin[[i]][[2]][[8]]
#  rmse_dropout_09_50[[i]]=rmse_fin[[i]][[2]][[9]]
  
#  rmse_dropout_01_60[[i]]=rmse_fin[[i]][[3]][[1]]
#  rmse_dropout_02_60[[i]]=rmse_fin[[i]][[3]][[2]]
#  rmse_dropout_03_60[[i]]=rmse_fin[[i]][[3]][[3]]
#  rmse_dropout_04_60[[i]]=rmse_fin[[i]][[3]][[4]]
#  rmse_dropout_05_60[[i]]=rmse_fin[[i]][[3]][[5]]
#  rmse_dropout_06_60[[i]]=rmse_fin[[i]][[3]][[6]]
#  rmse_dropout_07_60[[i]]=rmse_fin[[i]][[3]][[7]]
#  rmse_dropout_08_60[[i]]=rmse_fin[[i]][[3]][[8]]
#  rmse_dropout_09_60[[i]]=rmse_fin[[i]][[3]][[9]]
  
#  rmse_dropout_01_20[[i]]=rmse_fin[[i]][[4]][[1]]
#  rmse_dropout_02_20[[i]]=rmse_fin[[i]][[4]][[2]]
#  rmse_dropout_03_20[[i]]=rmse_fin[[i]][[4]][[3]]
#  rmse_dropout_04_20[[i]]=rmse_fin[[i]][[4]][[4]]
#  rmse_dropout_05_20[[i]]=rmse_fin[[i]][[4]][[5]]
#  rmse_dropout_06_20[[i]]=rmse_fin[[i]][[4]][[6]]
#  rmse_dropout_07_20[[i]]=rmse_fin[[i]][[4]][[7]]
#  rmse_dropout_08_20[[i]]=rmse_fin[[i]][[4]][[8]]
#  rmse_dropout_09_20[[i]]=rmse_fin[[i]][[4]][[9]]
#}
#accuracy=data.frame(stocks=names(list_entire_fin),RMSE_01_40=unlist(rmse_dropout_01_40),
#                    RMSE_02_40=unlist(rmse_dropout_02_40),RMSE_03_40=unlist(rmse_dropout_03_40),
#                    RMSE_04_40=unlist(rmse_dropout_04_40),RMSE_05_40=unlist(rmse_dropout_05_40),
#                    RMSE_06_40=unlist(rmse_dropout_06_40),RMSE_07_40=unlist(rmse_dropout_07_40),
#                    RMSE_08_40=unlist(rmse_dropout_08_40),RMSE_09_40=unlist(rmse_dropout_09_40),
                    
#                    RMSE_01_50=unlist(rmse_dropout_01_50),
#                    RMSE_02_50=unlist(rmse_dropout_02_50),RMSE_03_50=unlist(rmse_dropout_03_50),
#                    RMSE_04_50=unlist(rmse_dropout_04_50),RMSE_05_50=unlist(rmse_dropout_05_50),
#                    RMSE_06_50=unlist(rmse_dropout_06_50),RMSE_07_50=unlist(rmse_dropout_07_50),
#                    RMSE_08_50=unlist(rmse_dropout_08_50),RMSE_09_50=unlist(rmse_dropout_09_50),
                    
#                    RMSE_01_60=unlist(rmse_dropout_01_60),
#                    RMSE_02_60=unlist(rmse_dropout_02_60),RMSE_03_60=unlist(rmse_dropout_03_60),
#                    RMSE_04_60=unlist(rmse_dropout_04_60),RMSE_05_60=unlist(rmse_dropout_05_60),
#                    RMSE_06_60=unlist(rmse_dropout_06_60),RMSE_07_60=unlist(rmse_dropout_07_60),
#                    RMSE_08_60=unlist(rmse_dropout_08_60),RMSE_09_60=unlist(rmse_dropout_09_60),
                    
#                    RMSE_01_20=unlist(rmse_dropout_01_20),
#                    RMSE_02_20=unlist(rmse_dropout_02_20),RMSE_03_20=unlist(rmse_dropout_03_20),
#                    RMSE_04_20=unlist(rmse_dropout_04_20),RMSE_05_20=unlist(rmse_dropout_05_20),
#                    RMSE_06_20=unlist(rmse_dropout_06_20),RMSE_07_20=unlist(rmse_dropout_07_20),
#                    RMSE_08_20=unlist(rmse_dropout_08_20),RMSE_09_20=unlist(rmse_dropout_09_20))
#rownames(accuracy)=names(list_entire_fin)




##########################
#########Now doing the final predictions using the combinations
#########of dropout values and time steps that
#########generated the minimum values of RMSE in the nested cross-validation
library(anytime)

#Transform data to stationary
Series=list()
for(i in 1:length(list_entire_fin)){
  Series[[i]]=list_entire_fin[[i]]$Price
}
diffed=lapply(Series,diff)


n_to_rescale=N_to_rescale=c()
for(i in 1:length(list_entire_fin)){
  N_to_rescale[i]=length(diffed[[i]])
  n_to_rescale[i]=length(diffed[[i]])-length(list_test[[1]]$Price)
}


lag_transform <- function(x, k= 1){
  
  lagged =  c(rep(NA, k), x[1:(length(x)-k)])
  DF = as.data.frame(cbind(lagged, x))
  colnames(DF) <- c( paste0('x-', k), 'x')
  DF[is.na(DF)] <- 0
  return(DF)
}


split_sequence=function(sequence,n_steps){

  matrix_x_y=matrix(NA,nrow=(length(sequence)-n_steps),ncol=(n_steps+1))
  matrix_x=matrix(NA,nrow=(length(sequence)-n_steps),ncol=n_steps)
  matrix_y=matrix(NA,nrow=(length(sequence)-n_steps),ncol=1)
  
  end_ix=c()
  seq_X=seq_Y=list()
  #split a univariate sequence into samples
  for(i in 1:length(sequence)){
    # find the end of this pattern
    end_ix[i]=i+n_steps
    # check if we are beyond the sequence
    if(end_ix[i] > length(sequence)){
      break
    }
    # gather input and output parts of the pattern
    seq_X[[i]]=sequence[i:(end_ix[i]-1)]
    seq_Y[[i]]=sequence[end_ix[i]]
  }
  matrix_x=matrix(unlist(seq_X),ncol=n_steps,byrow=TRUE)
  
  
  matrix_y=matrix(unlist(seq_Y),ncol=1,byrow=TRUE)
  
  
  matrix_x_y=cbind(matrix_x,matrix_y)
  
  
  my_list=list(x=matrix_x,y=matrix_y,x_y=matrix_x_y)
  return(my_list)
}


supervised=list()
for(i in 1:length(diffed)){
  supervised[[i]]=split_sequence(diffed[[i]],n_steps=n_steps[rmse_min[[i]][1,1]])$x_y
}

#Split dataset into training and testing sets
N=n=c()
for(i in 1:length(supervised)){
  N[i] = nrow(supervised[[i]])
  n[i]=N[i]-length(list_test[[1]]$Price)
}

train=test=list()
for(i in 1:length(supervised)){
  train[[i]] = supervised[[i]][1:n[i], ]
  test[[i]]=supervised[[i]][(n[i]+1):N[i],]
}



#Normalize the data
# scale data
scale_data = function(train, test, feature_range = c(0, 1)) {
  x = train
  fr_min = feature_range[1]
  fr_max = feature_range[2]
  std_train = ((x - min(x) ) / (max(x) - min(x)  ))
  std_test  = ((test - min(x) ) / (max(x) - min(x)  ))
  
  scaled_train = std_train *(fr_max -fr_min) + fr_min
  scaled_test = std_test *(fr_max -fr_min) + fr_min
  
  return( list(scaled_train = scaled_train, scaled_test = scaled_test ,
               scaler= c(min =min(x), max = max(x))) )
  
}



Scaled=y_train=x_train=y_test=x_test=list()
for(i in 1:length(list_entire_fin)){
  Scaled[[i]] = scale_data(train[[i]],test[[i]], feature_range=c(-1, 1))
  y_train[[i]] = Scaled[[i]]$scaled_train[, (n_steps[rmse_min[[i]][1,1]]+1)]
  x_train[[i]] = Scaled[[i]]$scaled_train[, (1:n_steps[rmse_min[[i]][1,1]])]
  
  y_test[[i]]=Scaled[[i]]$scaled_test[, (n_steps[rmse_min[[i]][1,1]]+1)]
  x_test[[i]]=Scaled[[i]]$scaled_test[,(1:n_steps[rmse_min[[i]][1,1]])]
  
}

#The following code will be required to revert the predicted values to the original scale.
## inverse-transform
invert_scaling = function(scaled, scaler, feature_range = c(0, 1)){
  min = scaler[1]
  max = scaler[2]
  t = length(scaled)
  mins = feature_range[1]
  maxs = feature_range[2]
  inverted_dfs = numeric(t)
  
  for( i in 1:t){
    X = (scaled[i]- mins)/(maxs - mins)
    rawValues = X *(max - min) + min
    inverted_dfs[i] <- rawValues
  }
  return(inverted_dfs)
}

#############################
#############################
#############################
#############################
#############################
#Define the model
# Reshape the input to 3-dim
for(i in 1:length(list_entire_fin)){
  dim(x_train[[i]])=c(nrow(x_train[[i]]),ncol(x_train[[i]]),1)
}


X_shape2=X_shape3=list()
for(i in 1:length(list_entire_fin)){
  X_shape2[[i]]=dim(x_train[[i]])[2]
  X_shape3[[i]]=dim(x_train[[i]])[3]
}




batch_size = 32                
x_train_new=x_train
y_train_new=y_train
for(i in 1:length(list_entire_fin)){
  if(nrow(x_train[[i]])%%batch_size!=0){
    x_train_new[[i]]=x_train[[i]][-(1:(nrow(x_train[[i]])%%batch_size)),,1]
    y_train_new[[i]]=y_train[[i]][-(1:(length(y_train[[i]])%%batch_size))]
  }
  else{
    x_train_new=x_train
    y_train_new=y_train
  }
}


set.seed(12345)

x_test2=x_test
for(i in 1:length(list_test)){
  x_test2[[i]]=rbind(tail(x_train[[i]][,,1],(batch_size-1)),x_test[[i]])
}


###TO NOT RUN
###In this part of the code we perform the out-of-sample predictions using the best model chosen by a nested cross-validation.
###This part of the code is particularly computationally expensive. Uncomment this part of the code in order to perform the out-of-sample predictions.

#L=c(seq(1:length(list_test[[1]]$Price)))
#scaler=predictions_final=X=yhat_original=yhat_no_diff=yhat=list()
#for (i in 1:length(list_entire_fin)){
#  scaler[[i]] = Scaled[[i]]$scaler
#  predictions_final[[i]] = X[[i]]=list()
#  yhat_original[[i]]=yhat_no_diff[[i]]=yhat[[i]]=list()
  
#  model <- keras_model_sequential()
#  model %>%
#    layer_lstm(
#      units=10,
#      batch_input_shape = c(batch_size,
#                            X_shape2[[i]], X_shape3[[i]]),
#      stateful = TRUE,
#      dropout = dropouts[rmse_min[[i]][1,2]]
#    ) %>%
#    layer_dense(units=10)%>%#new hidden layer with 10 neurons
#    layer_dense(units = 1, activation = "linear")
#  #Compile the model
#  model %>% compile(
#    loss = 'mean_squared_error',
#    optimizer = optimizer_adam(learning_rate = 0.02),
#    metrics = c("MSE")
#  )
  
#  model %>% fit(
#    x_train_new[[i]],
#    y_train_new[[i]],
#    epochs = 500,
#    batch_size = batch_size,
#    verbose = 1,
#    shuffle = FALSE,
#    callbacks=list(callback_early_stopping(monitor="MSE",
#                                           min_delta=0,
#                                           patience=20,
#                                           verbose=0,
#                                           mode="min",
#                                           restore_best_weights = TRUE)
#
#    ))
#    for(j in 1:length(L)){
#    X[[i]][[j]]=x_test2[[i]][(j:(j+(batch_size-1))),]
#    dim(X[[i]][[j]]) = c(batch_size,n_steps[rmse_min[[i]][1,1]],1)
#    yhat_original[[i]][[j]] = tail((model %>% predict(X[[i]][[j]],batch_size=batch_size)),1)
#    # invert scaling
#    
#    yhat_no_diff[[i]][[j]] = invert_scaling(yhat_original[[i]][[j]], scaler[[i]],  c(-1, 1))
#    # invert differencing
#    yhat[[i]][[j]]  = yhat_no_diff[[i]][[j]] + Series[[i]][(n_to_rescale[i]+j)]#we consider Xt-1 from 2976:3478;the total
#    # store                                                           #number of training+test is 3479. n=2975.  
#    predictions_final[[i]][[j]] <- yhat[[i]][[j]]
#  }
#  predictions_final[[i]]=unlist(predictions_final[[i]])
#}

load("./lstm_data/model2_lin_predictions.Rdata",model2_lin_predictions <- new.env())
predictions_final=model2_lin_predictions$model2_lin_predictions


rmse=function(series,predictions)
{
  sqrt(mean((series-predictions)^2))#(series-predictions) or error/residuals
}

mae=function(series, predictions)
{
  mean(abs(series-predictions))
}



rmse_all=list()
for(i in (1:length(list_entire_fin))){
  rmse_all[[i]]=rmse(list_test[[i]]$Price,predictions_final[[i]])
  
}


accuracy_final=data.frame(stocks=names(list_entire_fin),RMSE=unlist(rmse_all))
rownames(accuracy_final)=names(list_entire_fin)


for(i in (1:length(list_entire_fin))){
  plot(list_test[[i]]$Date,list_test[[i]]$Price,type="l",ylim=range(c(list_test[[i]]$Price),c(predictions_final[[i]])),main=names(list_test)[[i]])
  lines(list_test[[i]]$Date,predictions_final[[i]],col="red")  
}
dev.off()



real_returns=list_test
for(i in (1:length(list_entire_fin))){
  for(t in 2:length(list_test[[1]]$Price)){
    real_returns[[i]][t,2]=(list_test[[i]][t,2]-list_test[[i]][t-1,2])/(list_test[[i]][t-1,2])
    
  }
}
#remove the first row because it contains prices
real_returns=lapply(real_returns, function(x) x[-1, ])




#Do the same with the returns
#Predicted returns
predictions_returns=list()
for(i in (1:length(list_entire_fin))){
  predictions_returns[[i]]=predictions_final[[i]]
  for(t in 2:length(list_test[[1]]$Price)){
    predictions_returns[[i]][t]=(predictions_final[[i]][t]-list_test[[i]][t-1,2])/(list_test[[i]][t-1,2])
  }
  
}

predictions_returns=lapply(predictions_returns, function(x) x[-1])

rmse_returns=list()
for(i in (1:length(list_entire_fin))){
  rmse_returns[[i]]=rmse(real_returns[[i]]$Price,predictions_returns[[i]])
}


accuracy_returns=data.frame(stocks=names(list_test),RMSE=unlist(rmse_returns))

rownames(accuracy_returns)=names(list_test)

