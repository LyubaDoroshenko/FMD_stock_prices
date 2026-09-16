#This code was used to simulate the 20 curves of length 300 and to do the 
#functional motif discovery and search using a partially random initialization
#as well as random initialization



# Set directory to the folder with the current file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# create directory for saving results
dir.create("sim300_only_random_fmd",showWarnings=FALSE)


library(combinat)
library(parallel)
library(class)
library(dendextend)


#Load smoothed curves simulated with Mixed causal-noncausal autoregressive process
load("./data_smoothed.RData",smoothed_env <- new.env())

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


files = list.files('./sim300_only_random_fmd/')#sim300_not_random is the name of the folder where
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
  load('./sim300_only_random_fmd/motifs_candidate.RData')
}else{
  # The following part of code is very time consuming, especially if run in sequential mode
  # (worker_number = 1)
  # find candidate motifs with random initialization
  n_init = 10 # number of random initialization to try
  find_candidate_motifs_results = find_candidate_motifs(Y0, Y1, K, c, n_init,
                                                        name = './sim300_only_random_fmd/len300', names_var = 'x(t)',V_init=NULL,
                                                        probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = 1000,
                                                                               iter4elong = iter4elong, trials_elong = trials_elong, max_gap = max_gap,
                                                                               return_options = TRUE, return_init = TRUE,
                                                                               diss = diss, alpha = alpha,transformed=TRUE),
                                                        plot = TRUE, worker_number = 1)
  save(find_candidate_motifs_results, file = './sim300_only_random_fmd/motifs_candidate.RData')
}

#load('./sim300_only_random_fmd/motifs_candidate.RData')


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
pdf('./sim300_only_random_fmd/clustering_candidate_motifs.pdf', height = 12, width = 9)
cluster_candidate_motifs_plot(cluster_candidate_motifs_results, ask = FALSE)#, R_all=0.0008)
dev.off()

### search selected motifs
motifs_search_results = motifs_search(cluster_candidate_motifs_results, #R_all = 0.0008, 
                                      use_real_occurrences = FALSE, length_diff = 0.3)#length_diff=+Inf


### plot FMD results (NB: no threshold of frequencies of motif found!)
pdf('./sim300_only_random_fmd/motifs_search_results.pdf', height = 10, width = 20)
motifs_search_plot_norm(motifs_search_results, ylab = 'x(t)', freq_threshold = 5,
                        transformed=TRUE)
dev.off()

save(find_candidate_motifs_results, silhouette_average,
     filter_candidate_motifs_results,
     cluster_candidate_motifs_results, motifs_search_results,
     file='./sim300_only_random_fmd/results_all.RData')







#Motifs analysis for fmd 
column_names=c("Motif","Length","Frequency","Radius")
motifs_analysis=data.frame(c(1:length(motifs_search_results$V_length)),motifs_search_results$V_length,
                           motifs_search_results$V_frequencies,
                           motifs_search_results$R_motifs*1000)
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

write_xlsx(motifs_analysis_ordered,"./sim300_only_random_fmd/simulated_motifs_analysis.xlsx")



#Let's multiply the radius for 1000
motifs_analysis_table_latex=data.frame(Motif=motifs_analysis_ordered$Motif,
                                       Length=motifs_analysis_ordered$Length,
                                       Frequency=motifs_analysis_ordered$Frequency,
                                       Radius=motifs_analysis_ordered$Radius)
library("xtable")
print(xtable(motifs_analysis_table_latex),digits=c(5,5,5,5),include.rownames = FALSE)



pdf('./sim300_only_random_fmd/motifs_search_results_for_paper_all_motifs.pdf', height = 4, width = 6)
motifs_search_plot_for_paper(motifs_search_results, freq_threshold = 5,
                             transformed=TRUE)
dev.off()

pdf('./sim300_only_random_fmd/motifs_search_results_for_paper.pdf', height = 4, width = 6)
motifs_search_plot_for_paper(motifs_search_results, index_plot_in_curves = c(8,2), freq_threshold = 5,
                             transformed=TRUE)
dev.off()

