#This code was used to simulate the 20 curves of length 300 and to do the 
#functional motif discovery and search using a partially random initialization
#as well as random initialization



# Set directory to the folder with the current file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# create directory for saving results
dir.create("sim_300_fmd",showWarnings=FALSE)


library(combinat)
library(parallel)
library(class)
library(dendextend)


#Load smoothed curves simulated with Mixed causal-noncausal autoregressive process
load("./sim300_fmd/data_smoothed.RData",smoothed_env <- new.env())

data_smoothed=smoothed_env$simulated300_data$data_smoothed
data_smoothed_derivative=smoothed_env$simulated300_data$data_smoothed_derivative


source("../fmd_functions.R")#It is necessary to upload the 
#file with all the necessary functions



ncurves=ncol(data_smoothed)
#######################
###     RUN FMD     ###
#######################
set.seed(13333)
# use Sobolev-like distance d_0.9
diss = 'd0_d1_L2'
alpha = 0.9

n_init=10
max_gap = 0 # no gaps allowed
iter4elong =1 # perform elongation
trials_elong = 50 # try all possible elongations 
c_max = 50 # maximum motif length 


### run probKMA multiple times (2x3x10=60 times)
K = c(2, 3) # number of clusters to try
c = c(15,20,30) # minimum motif lengths to try 


n_init = 10 # number of random initializations to try

V_initt=list()
V_initt[[1]]=list()         #2 clusters
V_initt[[1]][[1]]=list() #2 clusters of length 15

V_initt[[1]][[2]]=list() #2 clusters of length 20

V_initt[[1]][[3]]=list() #2 clusters of length 30


V_initt[[2]]=list()           #3 clusters
V_initt[[2]][[1]]=list()   #3 clusters of length 15

V_initt[[2]][[2]]=list()   #3 clusters of length 20

V_initt[[2]][[3]]=list()   #3 clusters of length 30



for(i in 1:4){
  a15=motifs_init(15,10,40,80,100)
  b15=motifs_init(15,10,40,80,100)
  V_initt[[1]][[1]][[i]]=list(list(v0=(matrix(a15$motif,ncol=1)),v1=(matrix(a15$motif_derivative,ncol=1))),
                              list(v0=(matrix(b15$motif,ncol=1)),v1=(matrix(b15$motif_derivative,ncol=1))))
  a20=motifs_init(20,10,40,80,100)
  b20=motifs_init(20,10,40,80,100)
  V_initt[[1]][[2]][[i]]=list(list(v0=(matrix(a20$motif,ncol=1)),v1=(matrix(a20$motif_derivative,ncol=1))),
                              list(v0=(matrix(b20$motif,ncol=1)),v1=(matrix(b20$motif_derivative,ncol=1))))
  a30=motifs_init(30,10,40,80,100)
  b30=motifs_init(30,10,40,80,100)
  V_initt[[1]][[3]][[i]]=list(list(v0=(matrix(a30$motif,ncol=1)),v1=(matrix(a30$motif_derivative,ncol=1))),
                              list(v0=(matrix(b30$motif,ncol=1)),v1=(matrix(b30$motif_derivative,ncol=1))))
  
  
  aa15=motifs_init(15,10,40,80,100)
  bb15=motifs_init(15,10,40,80,100)
  cc15=motifs_init(15,10,40,80,100)
  V_initt[[2]][[1]][[i]]=list(list(v0=matrix(aa15$motif,ncol=1),v1=matrix(aa15$motif_derivative,ncol=1)),
                              list(v0=matrix(bb15$motif,ncol=1),v1=matrix(bb15$motif_derivative,ncol=1)),
                              list(v0=matrix(cc15$motif,ncol=1),v1=matrix(cc15$motif_derivative,ncol=1)))
  aa20=motifs_init(20,10,40,80,100)
  bb20=motifs_init(20,10,40,80,100)
  cc20=motifs_init(20,10,40,80,100)
  V_initt[[2]][[2]][[i]]=list(list(v0=matrix(aa20$motif,ncol=1),v1=matrix(aa20$motif_derivative,ncol=1)),
                              list(v0=matrix(bb20$motif,ncol=1),v1=matrix(bb20$motif_derivative,ncol=1)),
                              list(v0=matrix(cc20$motif,ncol=1),v1=matrix(cc20$motif_derivative,ncol=1)))
  aa30=motifs_init(30,10,40,80,100)
  bb30=motifs_init(30,10,40,80,100)
  cc30=motifs_init(30,10,40,80,100)
  V_initt[[2]][[3]][[i]]=list(list(v0=matrix(aa30$motif,ncol=1),v1=matrix(aa30$motif_derivative,ncol=1)),
                              list(v0=matrix(bb30$motif,ncol=1),v1=matrix(bb30$motif_derivative,ncol=1)),
                              list(v0=matrix(cc30$motif,ncol=1),v1=matrix(cc30$motif_derivative,ncol=1)))
  
  
}

for(i in 5:8){
  a_2_15=motifs_init_rev(15,10,40,80,100)
  b_2_15=motifs_init_rev(15,10,40,80,100)
  V_initt[[1]][[1]][[i]]=list(list(v0=(matrix(a_2_15$motif,ncol=1)),v1=(matrix(a_2_15$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_2_15$motif,ncol=1)),v1=(matrix(b_2_15$motif_derivative,ncol=1))))
  a_2_20=motifs_init_rev(20,10,40,80,100)
  b_2_20=motifs_init_rev(20,10,40,80,100)
  V_initt[[1]][[2]][[i]]=list(list(v0=(matrix(a_2_20$motif,ncol=1)),v1=(matrix(a_2_20$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_2_20$motif,ncol=1)),v1=(matrix(b_2_20$motif_derivative,ncol=1))))
  a_2_30=motifs_init_rev(30,10,40,80,100)
  b_2_30=motifs_init_rev(30,10,40,80,100)
  V_initt[[1]][[3]][[i]]=list(list(v0=(matrix(a_2_30$motif,ncol=1)),v1=(matrix(a_2_30$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_2_30$motif,ncol=1)),v1=(matrix(b_2_30$motif_derivative,ncol=1))))
  
  
  aa_2_15=motifs_init_rev(15,10,40,80,100)
  bb_2_15=motifs_init_rev(15,10,40,80,100)
  cc_2_15=motifs_init_rev(15,10,40,80,100)
  V_initt[[2]][[1]][[i]]=list(list(v0=matrix(aa_2_15$motif,ncol=1),v1=matrix(aa_2_15$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_2_15$motif,ncol=1),v1=matrix(bb_2_15$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_2_15$motif,ncol=1),v1=matrix(cc_2_15$motif_derivative,ncol=1)))
  aa_2_20=motifs_init_rev(20,10,40,80,100)
  bb_2_20=motifs_init_rev(20,10,40,80,100)
  cc_2_20=motifs_init_rev(20,10,40,80,100)
  V_initt[[2]][[2]][[i]]=list(list(v0=matrix(aa_2_20$motif,ncol=1),v1=matrix(aa_2_20$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_2_20$motif,ncol=1),v1=matrix(bb_2_20$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_2_20$motif,ncol=1),v1=matrix(cc_2_20$motif_derivative,ncol=1)))
  aa_2_30=motifs_init_rev(30,10,40,80,100)
  bb_2_30=motifs_init_rev(30,10,40,80,100)
  cc_2_30=motifs_init_rev(30,10,40,80,100)
  V_initt[[2]][[3]][[i]]=list(list(v0=matrix(aa_2_30$motif,ncol=1),v1=matrix(aa_2_30$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_2_30$motif,ncol=1),v1=matrix(bb_2_30$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_2_30$motif,ncol=1),v1=matrix(cc_2_30$motif_derivative,ncol=1)))
  
  
}



#Increasing lines
for(i in 9:9){
  a_3_15=motifs_line(15,10,40,80,100)
  b_3_15=motifs_line(15,10,40,80,100)
  V_initt[[1]][[1]][[i]]=list(list(v0=(matrix(a_3_15$motif,ncol=1)),v1=(matrix(a_3_15$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_3_15$motif,ncol=1)),v1=(matrix(b_3_15$motif_derivative,ncol=1))))
  a_3_20=motifs_line(20,10,40,80,100)
  b_3_20=motifs_line(20,10,40,80,100)
  V_initt[[1]][[2]][[i]]=list(list(v0=(matrix(a_3_20$motif,ncol=1)),v1=(matrix(a_3_20$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_3_20$motif,ncol=1)),v1=(matrix(b_3_20$motif_derivative,ncol=1))))
  a_3_30=motifs_line(30,10,40,80,100)
  b_3_30=motifs_line(30,10,40,80,100)
  V_initt[[1]][[3]][[i]]=list(list(v0=(matrix(a_3_30$motif,ncol=1)),v1=(matrix(a_3_30$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_3_30$motif,ncol=1)),v1=(matrix(b_3_30$motif_derivative,ncol=1))))
  
  
  aa_3_15=motifs_line(15,10,40,80,100)
  bb_3_15=motifs_line(15,10,40,80,100)
  cc_3_15=motifs_line(15,10,40,80,100)
  V_initt[[2]][[1]][[i]]=list(list(v0=matrix(aa_3_15$motif,ncol=1),v1=matrix(aa_3_15$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_3_15$motif,ncol=1),v1=matrix(bb_3_15$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_3_15$motif,ncol=1),v1=matrix(cc_3_15$motif_derivative,ncol=1)))
  aa_3_20=motifs_line(20,10,40,80,100)
  bb_3_20=motifs_line(20,10,40,80,100)
  cc_3_20=motifs_line(20,10,40,80,100)
  V_initt[[2]][[2]][[i]]=list(list(v0=matrix(aa_3_20$motif,ncol=1),v1=matrix(aa_3_20$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_3_20$motif,ncol=1),v1=matrix(bb_3_20$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_3_20$motif,ncol=1),v1=matrix(cc_3_20$motif_derivative,ncol=1)))
  aa_3_30=motifs_line(30,10,40,80,100)
  bb_3_30=motifs_line(30,10,40,80,100)
  cc_3_30=motifs_line(30,10,40,80,100)
  V_initt[[2]][[3]][[i]]=list(list(v0=matrix(aa_3_30$motif,ncol=1),v1=matrix(aa_3_30$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_3_30$motif,ncol=1),v1=matrix(bb_3_30$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_3_30$motif,ncol=1),v1=matrix(cc_3_30$motif_derivative,ncol=1)))
  
  
}

#Decreasing lines
for(i in 10:10){
  a_4_15=motifs_line_rev(15,10,40,80,100)
  b_4_15=motifs_line_rev(15,10,40,80,100)
  V_initt[[1]][[1]][[i]]=list(list(v0=(matrix(a_4_15$motif,ncol=1)),v1=(matrix(a_4_15$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_4_15$motif,ncol=1)),v1=(matrix(b_4_15$motif_derivative,ncol=1))))
  a_4_20=motifs_line_rev(20,10,40,80,100)
  b_4_20=motifs_line_rev(20,10,40,80,100)
  V_initt[[1]][[2]][[i]]=list(list(v0=(matrix(a_4_20$motif,ncol=1)),v1=(matrix(a_4_20$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_4_20$motif,ncol=1)),v1=(matrix(b_4_20$motif_derivative,ncol=1))))
  a_4_30=motifs_line_rev(30,10,40,80,100)
  b_4_30=motifs_line_rev(30,10,40,80,100)
  V_initt[[1]][[3]][[i]]=list(list(v0=(matrix(a_4_30$motif,ncol=1)),v1=(matrix(a_4_30$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_4_30$motif,ncol=1)),v1=(matrix(b_4_30$motif_derivative,ncol=1))))
  
  
  aa_4_15=motifs_line_rev(15,10,40,80,100)
  bb_4_15=motifs_line_rev(15,10,40,80,100)
  cc_4_15=motifs_line_rev(15,10,40,80,100)
  V_initt[[2]][[1]][[i]]=list(list(v0=matrix(aa_4_15$motif,ncol=1),v1=matrix(aa_4_15$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_4_15$motif,ncol=1),v1=matrix(bb_4_15$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_4_15$motif,ncol=1),v1=matrix(cc_4_15$motif_derivative,ncol=1)))
  aa_4_20=motifs_line_rev(20,10,40,80,100)
  bb_4_20=motifs_line_rev(20,10,40,80,100)
  cc_4_20=motifs_line_rev(20,10,40,80,100)
  V_initt[[2]][[2]][[i]]=list(list(v0=matrix(aa_4_20$motif,ncol=1),v1=matrix(aa_4_20$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_4_20$motif,ncol=1),v1=matrix(bb_4_20$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_4_20$motif,ncol=1),v1=matrix(cc_4_20$motif_derivative,ncol=1)))
  aa_4_30=motifs_line_rev(30,10,40,80,100)
  bb_4_30=motifs_line_rev(30,10,40,80,100)
  cc_4_30=motifs_line_rev(30,10,40,80,100)
  V_initt[[2]][[3]][[i]]=list(list(v0=matrix(aa_4_30$motif,ncol=1),v1=matrix(aa_4_30$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_4_30$motif,ncol=1),v1=matrix(bb_4_30$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_4_30$motif,ncol=1),v1=matrix(cc_4_30$motif_derivative,ncol=1)))
  
  
}


V_init=V_initt

files = list.files('./sim300_fmd')#sim300_not_random is the name of the folder where
#the candidate motifs need to be present or where the candidate
#motifs will be saved if they are not there
Y0=Y1=list()
for(i in 1:ncol(data_smoothed)){
  Y0[[i]]=matrix(data_smoothed[,i],ncol=1)
  Y1[[i]]=matrix(data_smoothed_derivative[,i],ncol=1)
}







###############################
########## NOT RUN ############
###############################
if('motifs_candidate.RData' %in% files){
  # candidate motifs already present, load them
  load('./sim300_fmd/motifs_candidate.RData')
}else{
  # The following part of code is very time consuming, especially if run in sequential mode
  # (worker_number = 1)
  # find candidate motifs with partially random initialization# find candidate motifs with partially random initialization
  find_candidate_motifs_results = find_candidate_motifs(Y0, Y1, K, c, n_init,V_init=V_init,
                                                        name = './sim300_fmd/len300', names_var = 'x(t)',
                                                        probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = 1000,
                                                                               iter4elong = iter4elong, trials_elong = trials_elong, max_gap = max_gap,
                                                                               return_options = TRUE, return_init = TRUE,
                                                                               diss = diss, alpha = alpha,transformed=TRUE),
                                                        plot = TRUE, worker_number = 1)
  # find candidate motifs with random initialization
  n_init = 20 # number of partially random plus random initialization to try
  find_candidate_motifs_results = find_candidate_motifs(Y0, Y1, K, c, n_init,
                                                        name = './sim300_fmd/len300', names_var = 'x(t)',V_init=NULL,
                                                        probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = 1000,
                                                                               iter4elong = iter4elong, trials_elong = trials_elong, max_gap = max_gap,
                                                                               return_options = TRUE, return_init = TRUE,
                                                                               diss = diss, alpha = alpha,transformed=TRUE),
                                                        plot = TRUE, worker_number = 1)
  save(find_candidate_motifs_results, file = './sim300_fmd/motifs_candidate.RData')
}

#load('./sim300_fmd/motifs_candidate.RData')


### filter candidate motifs based on silhouette average and size
silhouette_average = Reduce(rbind, Reduce(rbind,
                                          find_candidate_motifs_results$silhouette_average_sd))[ , 1] # retrieve silhouette average for all candidate motifs
filter_candidate_motifs_results = filter_candidate_motifs(find_candidate_motifs_results,
                                                          sil_threshold = quantile(silhouette_average, 0.7),
                                                          size_threshold = 5)

### cluster candidate motifs based on their distance and select radii
cluster_candidate_motifs_results = cluster_candidate_motifs(filter_candidate_motifs_results,
                                                            motif_overlap = 0.8)


### plot cluster candidate motifs results
pdf('./sim300_fmd/clustering_candidate_motifs.pdf', height = 12, width = 9)
cluster_candidate_motifs_plot(cluster_candidate_motifs_results, ask = FALSE)# ,R_all=0.0002
dev.off()

### search selected motifs
motifs_search_results = motifs_search(cluster_candidate_motifs_results, #R_all = 0.005, 
                                      use_real_occurrences = FALSE, length_diff = 0.3)#length_diff=+Inf


### plot FMD results (NB: no threshold of frequencies of motif found!)
pdf('./sim300_fmd/motifs_search_results.pdf', height = 10, width = 20)
motifs_search_plot_norm(motifs_search_results, ylab = 'x(t)', freq_threshold = 5,
                        transformed=TRUE)
dev.off()

save(find_candidate_motifs_results, silhouette_average,
     filter_candidate_motifs_results,
     cluster_candidate_motifs_results, motifs_search_results,
     file='./sim300_fmd/results_all.RData')


#Motifs analysis for fmd 
column_names=c("Motif","Length","Frequency","Radius")
motifs_analysis=data.frame(c(1:19),motifs_search_results$V_length,
                           motifs_search_results$V_frequencies,
                           motifs_search_results$R_motifs)
colnames(motifs_analysis)=column_names


#Considering only the motifs with more than 5 occurrences:
filtered_motifs_analysis=motifs_analysis[motifs_analysis$Frequency >= 5, ]
# Renumber the Motif column
filtered_motifs_analysis$Motif=seq_len(nrow(filtered_motifs_analysis))

motifs_analysis_ordered=filtered_motifs_analysis[order(filtered_motifs_analysis$Length),]


library("writexl")
motifs_analysis_table=data.frame(Motif=motifs_analysis_ordered$Motif,
                                 Length=motifs_analysis_ordered$Length,
                                 Frequency=motifs_analysis_ordered$Frequency,
                                 Radius=motifs_analysis_ordered$Radius)

write_xlsx(motifs_analysis_ordered,"./sim300_fmd/simulated_motifs_analysis.xlsx")



#Let's multiply the radius for 1000
motifs_analysis_table_latex=data.frame(Motif=motifs_analysis_ordered$Motif,
                                 Length=motifs_analysis_ordered$Length,
                                 Frequency=motifs_analysis_ordered$Frequency,
                                 Radius=motifs_analysis_ordered$Radius*1000)
library("xtable")
print(xtable(motifs_analysis_table_latex),digits=c(5,5,5,5),include.rownames = FALSE)
