# This code was used to obtain the candidate motifs
# via functional motif discovery using the partially random initialization
# as well as random initialization 
# After the motifs are discovered, a motif search on all curves is performed


# Set directory to the folder with the current file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


source("../fmd_functions.R")#It is necessary to upload the 
#file with all the necessary functions


load("data_wrds.Rdata",smoothed_part2_env <- new.env())


# create directory for saving results
dir.create("stock_prices_motifs",showWarnings=FALSE)


data_smoothed=smoothed_part2_env$data_smoothed



ncurves=length(data_smoothed)

load("wrds_all_data.Rdata",smoothed_all_stocks<- new.env())
data_smoothed_all_stocks=smoothed_all_stocks$data_smoothed_all_stocks

#######################
###     RUN FMD     ###
#######################
set.seed(13333)
# use Sobolev-like distance d_0.5
diss = 'd0_d1_L2'
alpha = 0.9

max_gap = 0 # no gaps allowed
iter4elong =1 #5-10 # perform elongation
trials_elong =150 # try all possible elongations 60-70
c_max = 150 # maximum motif length 60-70


### run probKMA multiple times (2x3x10=60 times)

initializations=10
n_init = initializations # number of random initializations to try

V_initt=list()
V_initt[[1]]=list()         #2 clusters
V_initt[[1]][[1]]=list() #2 clusters of length 40

V_initt[[1]][[2]]=list() #2 clusters of length 50

V_initt[[1]][[3]]=list() #2 clusters of length 60


V_initt[[2]]=list()           #3 clusters
V_initt[[2]][[1]]=list()   #3 clusters of length 40

V_initt[[2]][[2]]=list()   #3 clusters of length 50

V_initt[[2]][[3]]=list()   #3 clusters of length 60



for(i in 1:4){
  a40=motifs_init(40,10,40,80,100)
  b40=motifs_init(40,10,40,80,100)
  V_initt[[1]][[1]][[i]]=list(list(v0=(matrix(a40$motif,ncol=1)),v1=(matrix(a40$motif_derivative,ncol=1))),
                              list(v0=(matrix(b40$motif,ncol=1)),v1=(matrix(b40$motif_derivative,ncol=1))))
  a50=motifs_init(50,10,40,80,100)
  b50=motifs_init(50,10,40,80,100)
  V_initt[[1]][[2]][[i]]=list(list(v0=(matrix(a50$motif,ncol=1)),v1=(matrix(a50$motif_derivative,ncol=1))),
                              list(v0=(matrix(b50$motif,ncol=1)),v1=(matrix(b50$motif_derivative,ncol=1))))
  a60=motifs_init(60,10,40,80,100)
  b60=motifs_init(60,10,40,80,100)
  V_initt[[1]][[3]][[i]]=list(list(v0=(matrix(a60$motif,ncol=1)),v1=(matrix(a60$motif_derivative,ncol=1))),
                              list(v0=(matrix(b60$motif,ncol=1)),v1=(matrix(b60$motif_derivative,ncol=1))))
  
  
  aa40=motifs_init(40,10,40,80,100)
  bb40=motifs_init(40,10,40,80,100)
  cc40=motifs_init(40,10,40,80,100)
  V_initt[[2]][[1]][[i]]=list(list(v0=matrix(aa40$motif,ncol=1),v1=matrix(aa40$motif_derivative,ncol=1)),
                              list(v0=matrix(bb40$motif,ncol=1),v1=matrix(bb40$motif_derivative,ncol=1)),
                              list(v0=matrix(cc40$motif,ncol=1),v1=matrix(cc40$motif_derivative,ncol=1)))
  aa50=motifs_init(50,10,40,80,100)
  bb50=motifs_init(50,10,40,80,100)
  cc50=motifs_init(50,10,40,80,100)
  V_initt[[2]][[2]][[i]]=list(list(v0=matrix(aa50$motif,ncol=1),v1=matrix(aa50$motif_derivative,ncol=1)),
                              list(v0=matrix(bb50$motif,ncol=1),v1=matrix(bb50$motif_derivative,ncol=1)),
                              list(v0=matrix(cc50$motif,ncol=1),v1=matrix(cc50$motif_derivative,ncol=1)))
  aa60=motifs_init(60,10,40,80,100)
  bb60=motifs_init(60,10,40,80,100)
  cc60=motifs_init(60,10,40,80,100)
  V_initt[[2]][[3]][[i]]=list(list(v0=matrix(aa60$motif,ncol=1),v1=matrix(aa60$motif_derivative,ncol=1)),
                              list(v0=matrix(bb60$motif,ncol=1),v1=matrix(bb60$motif_derivative,ncol=1)),
                              list(v0=matrix(cc60$motif,ncol=1),v1=matrix(cc60$motif_derivative,ncol=1)))
  
  
}

for(i in 5:8){
  a_2_40=motifs_init_rev(40,10,40,80,100)
  b_2_40=motifs_init_rev(40,10,40,80,100)
  V_initt[[1]][[1]][[i]]=list(list(v0=(matrix(a_2_40$motif,ncol=1)),v1=(matrix(a_2_40$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_2_40$motif,ncol=1)),v1=(matrix(b_2_40$motif_derivative,ncol=1))))
  a_2_50=motifs_init_rev(50,10,40,80,100)
  b_2_50=motifs_init_rev(50,10,40,80,100)
  V_initt[[1]][[2]][[i]]=list(list(v0=(matrix(a_2_50$motif,ncol=1)),v1=(matrix(a_2_50$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_2_50$motif,ncol=1)),v1=(matrix(b_2_50$motif_derivative,ncol=1))))
  a_2_60=motifs_init_rev(60,10,40,80,100)
  b_2_60=motifs_init_rev(60,10,40,80,100)
  V_initt[[1]][[3]][[i]]=list(list(v0=(matrix(a_2_60$motif,ncol=1)),v1=(matrix(a_2_60$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_2_60$motif,ncol=1)),v1=(matrix(b_2_60$motif_derivative,ncol=1))))
  
  
  aa_2_40=motifs_init_rev(40,10,40,80,100)
  bb_2_40=motifs_init_rev(40,10,40,80,100)
  cc_2_40=motifs_init_rev(40,10,40,80,100)
  V_initt[[2]][[1]][[i]]=list(list(v0=matrix(aa_2_40$motif,ncol=1),v1=matrix(aa_2_40$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_2_40$motif,ncol=1),v1=matrix(bb_2_40$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_2_40$motif,ncol=1),v1=matrix(cc_2_40$motif_derivative,ncol=1)))
  aa_2_50=motifs_init_rev(50,10,40,80,100)
  bb_2_50=motifs_init_rev(50,10,40,80,100)
  cc_2_50=motifs_init_rev(50,10,40,80,100)
  V_initt[[2]][[2]][[i]]=list(list(v0=matrix(aa_2_50$motif,ncol=1),v1=matrix(aa_2_50$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_2_50$motif,ncol=1),v1=matrix(bb_2_50$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_2_50$motif,ncol=1),v1=matrix(cc_2_50$motif_derivative,ncol=1)))
  aa_2_60=motifs_init_rev(60,10,40,80,100)
  bb_2_60=motifs_init_rev(60,10,40,80,100)
  cc_2_60=motifs_init_rev(60,10,40,80,100)
  V_initt[[2]][[3]][[i]]=list(list(v0=matrix(aa_2_60$motif,ncol=1),v1=matrix(aa_2_60$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_2_60$motif,ncol=1),v1=matrix(bb_2_60$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_2_60$motif,ncol=1),v1=matrix(cc_2_60$motif_derivative,ncol=1)))
  
  
}

plot(motifs_init(40,10,40,80,100)$motif,type="l",ylab="",xlab="Motif length")
lines(motifs_init(40,10,40,80,100)$motif,col="red")
lines(motifs_init(40,10,40,80,100)$motif,col="green")
lines(motifs_init(40,10,40,80,100)$motif,col="blue")




#Increasing lines
for(i in 9:9){
  a_3_40=motifs_line(40,10,40,80,100)
  b_3_40=motifs_line(40,10,40,80,100)
  V_initt[[1]][[1]][[i]]=list(list(v0=(matrix(a_3_40$motif,ncol=1)),v1=(matrix(a_3_40$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_3_40$motif,ncol=1)),v1=(matrix(b_3_40$motif_derivative,ncol=1))))
  a_3_50=motifs_line(50,10,40,80,100)
  b_3_50=motifs_line(50,10,40,80,100)
  V_initt[[1]][[2]][[i]]=list(list(v0=(matrix(a_3_50$motif,ncol=1)),v1=(matrix(a_3_50$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_3_50$motif,ncol=1)),v1=(matrix(b_3_50$motif_derivative,ncol=1))))
  a_3_60=motifs_line(60,10,40,80,100)
  b_3_60=motifs_line(60,10,40,80,100)
  V_initt[[1]][[3]][[i]]=list(list(v0=(matrix(a_3_60$motif,ncol=1)),v1=(matrix(a_3_60$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_3_60$motif,ncol=1)),v1=(matrix(b_3_60$motif_derivative,ncol=1))))
  
  
  aa_3_40=motifs_line(40,10,40,80,100)
  bb_3_40=motifs_line(40,10,40,80,100)
  cc_3_40=motifs_line(40,10,40,80,100)
  V_initt[[2]][[1]][[i]]=list(list(v0=matrix(aa_3_40$motif,ncol=1),v1=matrix(aa_3_40$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_3_40$motif,ncol=1),v1=matrix(bb_3_40$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_3_40$motif,ncol=1),v1=matrix(cc_3_40$motif_derivative,ncol=1)))
  aa_3_50=motifs_line(50,10,40,80,100)
  bb_3_50=motifs_line(50,10,40,80,100)
  cc_3_50=motifs_line(50,10,40,80,100)
  V_initt[[2]][[2]][[i]]=list(list(v0=matrix(aa_3_50$motif,ncol=1),v1=matrix(aa_3_50$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_3_50$motif,ncol=1),v1=matrix(bb_3_50$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_3_50$motif,ncol=1),v1=matrix(cc_3_50$motif_derivative,ncol=1)))
  aa_3_60=motifs_line(60,10,40,80,100)
  bb_3_60=motifs_line(60,10,40,80,100)
  cc_3_60=motifs_line(60,10,40,80,100)
  V_initt[[2]][[3]][[i]]=list(list(v0=matrix(aa_3_60$motif,ncol=1),v1=matrix(aa_3_60$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_3_60$motif,ncol=1),v1=matrix(bb_3_60$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_3_60$motif,ncol=1),v1=matrix(cc_3_60$motif_derivative,ncol=1)))
  
  
}

#Decreasing lines
for(i in 10:10){
  a_4_40=motifs_line_rev(40,10,40,80,100)
  b_4_40=motifs_line_rev(40,10,40,80,100)
  V_initt[[1]][[1]][[i]]=list(list(v0=(matrix(a_4_40$motif,ncol=1)),v1=(matrix(a_4_40$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_4_40$motif,ncol=1)),v1=(matrix(b_4_40$motif_derivative,ncol=1))))
  a_4_50=motifs_line_rev(50,10,40,80,100)
  b_4_50=motifs_line_rev(50,10,40,80,100)
  V_initt[[1]][[2]][[i]]=list(list(v0=(matrix(a_4_50$motif,ncol=1)),v1=(matrix(a_4_50$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_4_50$motif,ncol=1)),v1=(matrix(b_4_50$motif_derivative,ncol=1))))
  a_4_60=motifs_line_rev(60,10,40,80,100)
  b_4_60=motifs_line_rev(60,10,40,80,100)
  V_initt[[1]][[3]][[i]]=list(list(v0=(matrix(a_4_60$motif,ncol=1)),v1=(matrix(a_4_60$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_4_60$motif,ncol=1)),v1=(matrix(b_4_60$motif_derivative,ncol=1))))
  
  
  aa_4_40=motifs_line_rev(40,10,40,80,100)
  bb_4_40=motifs_line_rev(40,10,40,80,100)
  cc_4_40=motifs_line_rev(40,10,40,80,100)
  V_initt[[2]][[1]][[i]]=list(list(v0=matrix(aa_4_40$motif,ncol=1),v1=matrix(aa_4_40$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_4_40$motif,ncol=1),v1=matrix(bb_4_40$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_4_40$motif,ncol=1),v1=matrix(cc_4_40$motif_derivative,ncol=1)))
  aa_4_50=motifs_line_rev(50,10,40,80,100)
  bb_4_50=motifs_line_rev(50,10,40,80,100)
  cc_4_50=motifs_line_rev(50,10,40,80,100)
  V_initt[[2]][[2]][[i]]=list(list(v0=matrix(aa_4_50$motif,ncol=1),v1=matrix(aa_4_50$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_4_50$motif,ncol=1),v1=matrix(bb_4_50$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_4_50$motif,ncol=1),v1=matrix(cc_4_50$motif_derivative,ncol=1)))
  aa_4_60=motifs_line_rev(60,10,40,80,100)
  bb_4_60=motifs_line_rev(60,10,40,80,100)
  cc_4_60=motifs_line_rev(60,10,40,80,100)
  V_initt[[2]][[3]][[i]]=list(list(v0=matrix(aa_4_60$motif,ncol=1),v1=matrix(aa_4_60$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_4_60$motif,ncol=1),v1=matrix(bb_4_60$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_4_60$motif,ncol=1),v1=matrix(cc_4_60$motif_derivative,ncol=1)))
  
  
}


V_init=V_initt
K = c(2, 3) # number of clusters to try
c=c(40,50,60)
n_init = 10 # number of partially random initializations to try
# NOTE: rename "results" folder to re-run everything
#(TIME CONSUMING)

set.seed(13333)
files = list.files('./stock_prices_motifs')

Y0=Y1=list()
for(i in 1:length(data_smoothed)){
  Y0[[i]]=matrix(data_smoothed[[i]]$v0,ncol=1)
  Y1[[i]]=matrix(data_smoothed[[i]]$v1,ncol=1)#Just the uncorrelated stocks
  
}

Y0_all=Y1_all=list()
for(i in 1:length(data_smoothed_all_stocks)){
  Y0_all[[i]]=matrix(data_smoothed_all_stocks[[i]]$v0,ncol=1)
  Y1_all[[i]]=matrix(data_smoothed_all_stocks[[i]]$v1,ncol=1)#All the stocks
  
}

###############################
########## NOT RUN ############
###############################
# It is necessary to load the saved motifs candidates in the file motifs_candidate.RData in order to avoid a time consuming sequential candidate motifs discovery.
# Otherwise a sequential functional motif discovery can be performed by running the commented code.


# The following part of code is very time consuming, especially if run in sequential mode
# (worker_number = 1)
#if('motifs_candidate.RData' %in% files){
  #if candidate motifs already present, load them

load('./stock_prices_motifs/motifs_candidate.RData')

#}else{
# find candidate motifs with partially random initialization

#find_candidate_motifs_results = find_candidate_motifs(Y0, Y1, K, c, n_init,V_init=V_init,
#                                                        name = './stock_prices_motifs/motifs', names_var = 'x(t)',
#                                                        probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = 1000,
#                                                                               iter4elong = iter4elong, trials_elong = trials_elong, max_gap = max_gap,
#                                                                               return_options = TRUE, return_init = TRUE,
#                                                                               diss = diss, alpha = alpha,transformed=TRUE),
#                                                        plot = TRUE, worker_number = 1)



# find candidate motifs with random initialization
n_init = 20 # number of partially random plus random initialization to try

#find_candidate_motifs_results = find_candidate_motifs(Y0, Y1, K, c, n_init,V_init=NULL,
#                                                      name = './stock_prices_motifs/motifs', names_var = 'x(t)',
#                                                      probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = 1000,
#                                                                             iter4elong = iter4elong, trials_elong = trials_elong, max_gap = max_gap,
#                                                                             return_options = TRUE, return_init = TRUE,
#                                                                             diss = diss, alpha = alpha,transformed=TRUE),
#                                                      plot = TRUE, worker_number = 1)
#save(find_candidate_motifs_results, file = './stock_prices_motifs/motifs_candidate.RData')
#}


### filter candidate motifs based on silhouette average and size
silhouette_average = Reduce(rbind, Reduce(rbind,
                                          find_candidate_motifs_results$silhouette_average_sd))[ , 1] # retrieve silhouette average for all candidate motifs

filter_candidate_motifs_results = filter_candidate_motifs(find_candidate_motifs_results,
                                                          sil_threshold = quantile(silhouette_average, 0.7),
                                                          size_threshold = 5)

### cluster candidate motifs based on their distance and select radii
cluster_candidate_motifs_results = cluster_candidate_motifs(filter_candidate_motifs_results, k_knn = 6,votes_knn_Rm=0.5,
                                                            votes_knn_Rm_finding = 0.15,
                                                            motif_overlap = 0.8)

### plot cluster candidate motifs results
pdf('./stock_prices_motifs/clustering_candidate_motifs.pdf', height = 12, width = 9)
cluster_candidate_motifs_plot(cluster_candidate_motifs_results, ask = FALSE,R_all=0.0002)#R_all=0.0008
dev.off()


#Use all the curves (also the ones that were excluded due to the correlation matters)
#for motifs search
cluster_candidate_motifs_results$Y0=Y0_all
cluster_candidate_motifs_results$Y1=Y1_all

### search selected motifs
motifs_search_results = motifs_search(cluster_candidate_motifs_results,R_all=0.0002,
                                      use_real_occurrences = FALSE, length_diff = 0.3,
                                      different_R_m_finding = TRUE)


save(motifs_search_results, file = './stock_prices_motifs/motifs_search_results.RData')



library("anytime")
library("lubridate")
library("dplyr")
#We create the ticks for the dates of each stock
#The file SP500_wrds.xlsx contains the daily prices of the stocks
library("readxl")
CMCSA_date=anydate(read_excel("SP500_wrds.xlsx",sheet="CMCSA")$Date,tz="UTC")
CMCSA_price=read_excel("SP500_wrds.xlsx",sheet="CMCSA")$Price
CMCSA_dates=data.frame(CMCSA_price,year(CMCSA_date),month(CMCSA_date),day(CMCSA_date))
names(CMCSA_dates)=c("price","year","month","day")   

#we manually create the ticks on the first available day of January for each year
years_names=c("2010","2011","2012",
              "2013","2014","2015","2016",
              "2017","2018","2019","2020")

TSLA_date=anydate(read_excel("SP500_wrds.xlsx",sheet="TSLA")$Date,tz="UTC")
TSLA_price=read_excel("SP500_wrds.xlsx",sheet="TSLA")$Price
TSLA_dates=data.frame(TSLA_price,year(TSLA_date),month(TSLA_date),day(TSLA_date))
names(TSLA_dates)=c("price","year","month","day")   

years_TSLA_names=c("2011","2012",
                   "2013","2014","2015","2016",
                   "2017","2018","2019","2020")


GOOG_date=anydate(read_excel("SP500_wrds.xlsx",sheet="GOOG")$Date,tz="UTC")
GOOG_price=read_excel("SP500_wrds.xlsx",sheet="GOOG")$Price
GOOG_dates=data.frame(GOOG_price,year(GOOG_date),month(GOOG_date),day(GOOG_date))
names(GOOG_dates)=c("price","year","month","day")   

years_GOOG_names=c("2015","2016","2017","2018","2019",
                   "2020")

FB_date=anydate(read_excel("SP500_wrds.xlsx",sheet="FB")$Date,tz="UTC")
FB_price=read_excel("SP500_wrds.xlsx",sheet="FB")$Price
FB_dates=data.frame(FB_price,year(FB_date),month(FB_date),day(FB_date))
names(FB_dates)=c("price","year","month","day")   

years_FB_names=c("2013","2014","2015","2016",
                 "2017","2018","2019","2020")


ABBV_date=anydate(read_excel("SP500_wrds.xlsx",sheet="ABBV")$Date,tz="UTC")
ABBV_price=read_excel("SP500_wrds.xlsx",sheet="ABBV")$Price
ABBV_dates=data.frame(ABBV_price,year(ABBV_date),month(ABBV_date),day(ABBV_date))
names(ABBV_dates)=c("price","year","month","day")   

years_ABBV_names=c("2013","2014","2015","2016",
                   "2017","2018","2019","2020")

AVGO_date=anydate(read_excel("SP500_wrds.xlsx",sheet="AVGO")$Date,tz="UTC")
AVGO_price=read_excel("SP500_wrds.xlsx",sheet="AVGO")$Price
AVGO_dates=data.frame(AVGO_price,year(AVGO_date),month(AVGO_date),day(AVGO_date))
names(AVGO_dates)=c("price","year","month","day")   

years_AVGO_names=c("2010","2011","2012",
                   "2013","2014","2015","2016",
                   "2017","2018","2019","2020")




#I will take a tick on January 2-4 of each year
years=c(2010:2020)
CMCSA_ticks=c()
for(i in 1:length(years)){
  if(years[i] == 2010){
    CMCSA_ticks[i]=which(CMCSA_dates$year == years[i] & CMCSA_dates$month == 1 & CMCSA_dates$day == 4)
  }
  if(years[i]==2016){
    CMCSA_ticks[i]=which(CMCSA_dates$year == years[i] & CMCSA_dates$month == 1 & CMCSA_dates$day == 4)
  }
  if(years[i] <= 2012 &years[i]>=2011){
    CMCSA_ticks[i]=which(CMCSA_dates$year == years[i] & CMCSA_dates$month == 1 & CMCSA_dates$day == 3)
  }
  if(years[i] <= 2015 &years[i]>=2013){
    CMCSA_ticks[i]=which(CMCSA_dates$year == years[i] & CMCSA_dates$month == 1 & CMCSA_dates$day == 2)
  }
  if(years[i] == 2017){
    CMCSA_ticks[i]=which(CMCSA_dates$year == years[i] & CMCSA_dates$month == 1 & CMCSA_dates$day == 3)
  }
  
  if(years[i] <= 2020 &years[i]>=2018){
    CMCSA_ticks[i]=which(CMCSA_dates$year == years[i] & CMCSA_dates$month == 1 & CMCSA_dates$day == 2)
  }
}

ticks_all=list()
years_all=list()
for(i in 1:length(data_smoothed_all_stocks)){
  ticks_all[[i]]=CMCSA_ticks
  years_all[[i]]=years_names
}


years_TSLA=c(2011:2020)
TSLA_ticks=c()
for(i in 1:length(years_TSLA)){
  if(years_TSLA[[i]]==2016){
    TSLA_ticks[i]=which(TSLA_dates$year == years_TSLA[i] & TSLA_dates$month == 1 & TSLA_dates$day == 4)
  }
  if(years_TSLA[[i]]==2017){
    TSLA_ticks[i]=which(TSLA_dates$year == years_TSLA[i] & TSLA_dates$month == 1 & TSLA_dates$day == 4)
  }
  if(years_TSLA[i] <= 2012 &years_TSLA[i]>=2011){
    TSLA_ticks[i]=which(TSLA_dates$year == years_TSLA[i] & TSLA_dates$month == 1 & TSLA_dates$day == 3)
  }
  if(years_TSLA[i] <= 2015 &years_TSLA[i]>=2013){
    TSLA_ticks[i]=which(TSLA_dates$year == years_TSLA[i] & TSLA_dates$month == 1 & TSLA_dates$day == 2)
  }
  if(years_TSLA[i]<=2020&years_TSLA[i]>=2018){
    TSLA_ticks[i]=which(TSLA_dates$year == years_TSLA[i] & TSLA_dates$month == 1 & TSLA_dates$day == 2)
  }
  
}

years_GOOG=c(2015:2020)
GOOG_ticks=c()
for(i in 1:length(years_GOOG)){
  if(years_GOOG[i]==2015){
    GOOG_ticks[i]=which(GOOG_dates$year == years_GOOG[i] & GOOG_dates$month == 1 & GOOG_dates$day == 2)
  }
  
  if(years_GOOG[i]==2016){
    GOOG_ticks[i]=which(GOOG_dates$year == years_GOOG[i] & GOOG_dates$month == 1 & GOOG_dates$day == 4)
  }
  
  if(years_GOOG[i]==2017){
    GOOG_ticks[i]=which(GOOG_dates$year == years_GOOG[i] & GOOG_dates$month == 1 & GOOG_dates$day == 3)
  }
  
  if(years_GOOG[i] <= 2019 &years_GOOG[i]>=2018){
    GOOG_ticks[i]=which(GOOG_dates$year == years_GOOG[i] & GOOG_dates$month == 1 & GOOG_dates$day == 2)
  }
  
  if(years_GOOG[i]==2020){
    GOOG_ticks[i]=which(GOOG_dates$year == years_GOOG[i] & GOOG_dates$month == 1 & GOOG_dates$day == 3)
  }
  
}

years_FB=c(2013:2020)
FB_ticks=c()
for(i in 1:length(years_FB)){
  if(years_FB[i] <= 2015 &years_FB[i]>=2013){
    FB_ticks[i]=which(FB_dates$year == years_FB[i] & FB_dates$month == 1 & FB_dates$day == 2)
  }
  if(years_FB[i] == 2016){
    FB_ticks[i]=which(FB_dates$year == years_FB[i] & FB_dates$month == 1 & FB_dates$day == 4)
  }
  if(years_FB[i] == 2017){
    FB_ticks[i]=which(FB_dates$year == years_FB[i] & FB_dates$month == 1 & FB_dates$day == 3)
  }
  if(years_FB[i] <= 2020 &years_FB[i]>=2018){
    FB_ticks[i]=which(FB_dates$year == years_FB[i] & FB_dates$month == 1 & FB_dates$day == 2)
  }
  
}


years_ABBV=c(2013:2020)
ABBV_ticks=c()
for(i in 1:length(years_ABBV)){
  
  if(years_ABBV[i] <= 2015 &years_ABBV[i]>=2013){
    ABBV_ticks[i]=which(ABBV_dates$year == years_ABBV[i] & ABBV_dates$month == 1 & ABBV_dates$day == 2)
  }
  
  if(years_ABBV[i] == 2016){
    ABBV_ticks[i]=which(ABBV_dates$year == years_ABBV[i] & ABBV_dates$month == 1 & ABBV_dates$day == 4)
  }
  
  if(years_ABBV[i] == 2017){
    ABBV_ticks[i]=which(ABBV_dates$year == years_ABBV[i] & ABBV_dates$month == 1 & ABBV_dates$day == 3)
  }
  
  if(years_ABBV[i]<=2020&years_ABBV[i]>=2018){
    ABBV_ticks[i]=which(ABBV_dates$year == years_ABBV[i] & ABBV_dates$month == 1 & ABBV_dates$day == 2)
  }
  
}




years_AVGO=c(2010:2020)
AVGO_ticks=c()
for(i in 1:length(years_AVGO)){
  if(years_AVGO[i]==2010){
    AVGO_ticks[i]=which(AVGO_dates$year == years_AVGO[i] & AVGO_dates$month == 1 & AVGO_dates$day == 4)
  }
  
  if(years_AVGO[i] <= 2012 &years_AVGO[i]>=2011){
    AVGO_ticks[i]=which(AVGO_dates$year == years_AVGO[i] & AVGO_dates$month == 1 & AVGO_dates$day == 3)
  }
  
  
  if(years_AVGO[i] <= 2015 &years_AVGO[i]>=2013){
    AVGO_ticks[i]=which(AVGO_dates$year == years_AVGO[i] & AVGO_dates$month == 1 & AVGO_dates$day == 2)
  }
  
  if(years_AVGO[i]==2016){
    AVGO_ticks[i]=which(AVGO_dates$year == years_AVGO[i] & AVGO_dates$month == 1 & AVGO_dates$day == 4)
  }
  
  if(years_AVGO[i]==2017){
    AVGO_ticks[i]=which(AVGO_dates$year == years_AVGO[i] & AVGO_dates$month == 1 & AVGO_dates$day == 3)
  }
  
  if(years_AVGO[i]<=2020&years_AVGO[i]>=2018){
    AVGO_ticks[i]=which(AVGO_dates$year == years_AVGO[i] & AVGO_dates$month == 1 & AVGO_dates$day == 2)
  }
  
}





ticks_all[[6]]=TSLA_ticks
years_all[[6]]=years_TSLA_names

ticks_all[[11]]=FB_ticks
years_all[[11]]=years_FB_names


ticks_all[[5]]=GOOG_ticks
years_all[[5]]=years_GOOG_names

ticks_all[[20]]=ABBV_ticks
years_all[[20]]=years_ABBV_names

ticks_all[[24]]=AVGO_ticks
years_all[[24]]=years_AVGO_names


#normalized with time
pdf('./stock_prices_motifs/motifs_search_results.pdf', height = 13, width = 38)
motifs_search_plot_norm_time(motifs_search_results, ylab = 'x(t)', freq_threshold = 5,transformed=TRUE,data_smoothed=data_smoothed_all_stocks)
dev.off()

#From here I take my observations for predictions
library("readxl")
#Uploading the file excel with the stock prices
stocks=read_excel("SP500_wrds.xlsx")
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
library(zoo)
for(i in 1:length(list_entire_fin)){
  if(length(anydate(list_entire_fin[[i]]$Date,tz="UTC"))!=length(anydate(list_entire_fin$SP500$Date,tz="UTC"))){
    SP_dates[[i]]=data.frame(Date=((list_entire_fin$SP500[anydate(list_entire_fin$SP500$Date,tz="UTC")>=anydate(list_entire_fin[[i]]$Date,tz="UTC")[1],])$Date),
                             Price=((list_entire_fin$SP500[anydate(list_entire_fin$SP500$Date,tz="UTC")>=anydate(list_entire_fin[[i]]$Date,tz="UTC")[1],])$Price))
    
    #removing the prices for the date not present in SP500 but present in the other stocks
    list_entire_fin[[i]]=list_entire_fin[[i]][which(list_entire_fin[[i]]$Date%in%SP_dates[[i]]$Date),]
    #imputing the prices present in Sp500 but absent in the other stocks
    list_entire_fin[[i]]= SP_dates[[i]] %>%
      dplyr::select(Date) %>%
      left_join(list_entire_fin[[i]], by = 'Date') %>%
      mutate(Price = zoo::na.approx(Price))
    
  } }


#Transforming the date into an appropriate format and subsetting until Dec 31, 2020
library(anytime)
training=function(x){
  x=x[x$Date<="2020-12-31",]
}

test=function(x){
  x=x[x$Date>"2020-12-31",]
}

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




#The function for detrending the times series by substracting the linear trend
detrend_calculate=function(series){
  reg=lm(series ~ seq(1:length(series)))
  detrended = series - reg$fitted.values
  return(detrended)
}


#This is my list of stocks in the test period
list_test=lapply(list_entire_fin,test)

names(list_test)=names(list_training)


#detrend substracting each series trend
stocks_of_interest_detrend=lapply(list_training_ts,detrend_calculate)
names(stocks_of_interest_detrend)=c("AAPL","MSFT","AMZN","GOOGL","GOOG","TSLA","BRKb","NVDA",
                                    "JNJ","UNH","FB","XOM","JPM","VISA","PG","CVX",
                                    "HD","MA","PFE","ABBV","BAC","KO","LLY","AVGO","PEP","MRK",
                                    "TMO","VZ","COST","ABT","ADBE","DIS","ACN","CMCSA","CSCO",
                                    "MCD","CRM","WMT","INTC","WFC","AMD","LIN","DHR","PM","BMY",
                                    "QCOM","NEE","TXN","NKE","COP","SP500")

bandwidth_smooth = 10

###################
###################
###################
#Making the predictions. Start here
require(KernSmooth)
forecast_motifs=motifs_forecast(motifs_search_results=motifs_search_results,list_test=list_test,
                                list_training_ts=list_training_ts,closest=3,overlap=0.4,bandwidth_smooth=bandwidth_smooth)$forecast

save(forecast_motifs, file = './stock_prices_motifs/forecast_motifs.RData')

rmse=function(series,predictions)
{
  sqrt(mean((series-predictions)^2))
}


rmse_all=list()
for(i in 1:length(list_test)){
  rmse_all[[i]]=rmse(list_test[[i]]$Price,forecast_motifs$forecast[[i]])
}
names(rmse_all)=names(list_test)

accuracy=data.frame(stocks=names(list_test),RMSE=unlist(rmse_all))
rownames(accuracy)=names(list_test)


#Real returns
real_returns=list_test
for(i in 1:length(list_test)){
  real_returns[[i]]$Price[2:nrow(list_test[[i]])] = 
    (list_test[[i]]$Price[2:nrow(list_test[[i]])] - list_test[[i]]$Price[2:nrow(list_test[[i]])-1]) / list_test[[i]]$Price[2:nrow(list_test[[i]])-1]
  #compute return for the first test time
  real_returns[[i]]$Price[1] = (list_test[[i]]$Price[1] - list_training[[i]]$Price[nrow(list_training[[i]])]) / list_training[[i]]$Price[nrow(list_training[[i]])]
  names(real_returns[[i]])[2] = "Return"
}

#Predicted returns
predictions_returns=forecast_motifs$forecast
for(i in 1:length(list_test)){
  predictions_returns[[i]][2:nrow(list_test[[i]])] = 
    (forecast_motifs$forecast[[i]][2:nrow(list_test[[i]])] - list_test[[i]]$Price[2:nrow(list_test[[i]])-1]) / list_test[[i]]$Price[2:nrow(list_test[[i]])-1]
  #compute return for the first test time
  predictions_returns[[i]][1] = (forecast_motifs$forecast[[i]][1] - list_training[[i]]$Price[nrow(list_training[[i]])]) / list_training[[i]]$Price[nrow(list_training[[i]])]
}


rmse_all_returns=list()
for(i in 1:length(list_test)){
  rmse_all_returns[[i]]=rmse(real_returns[[i]]$Return,predictions_returns[[i]])
}
names(rmse_all_returns)=names(list_test)

accuracy_returns=data.frame(stocks=names(list_test),RMSE=unlist(rmse_all_returns))
rownames(accuracy_returns)=names(list_test)






