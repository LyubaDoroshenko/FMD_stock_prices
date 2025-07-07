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
dir.create("stock_prices_motifs_shorter",showWarnings=FALSE)


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
V_initt[[1]][[1]]=list() #2 clusters of length 20

V_initt[[1]][[2]]=list() #2 clusters of length 30

V_initt[[1]][[3]]=list() #2 clusters of length 40

V_initt[[1]][[4]]=list() #2 clusters of length 50

V_initt[[1]][[5]]=list() #2 clusters of length 60

V_initt[[2]]=list()           #3 clusters
V_initt[[2]][[1]]=list()   #3 clusters of length 20

V_initt[[2]][[2]]=list()   #3 clusters of length 30

V_initt[[2]][[3]]=list()   #3 clusters of length 40

V_initt[[2]][[4]]=list()   #3 clusters of length 50

V_initt[[2]][[5]]=list()   #3 clusters of length 60

for(i in 1:4){
  a20=motifs_init(20,10,40,80,100)
  b20=motifs_init(20,10,40,80,100)
  V_initt[[1]][[1]][[i]]=list(list(v0=(matrix(a20$motif,ncol=1)),v1=(matrix(a20$motif_derivative,ncol=1))),
                              list(v0=(matrix(b20$motif,ncol=1)),v1=(matrix(b20$motif_derivative,ncol=1))))
  
  a30=motifs_init(30,10,40,80,100)
  b30=motifs_init(30,10,40,80,100)
  V_initt[[1]][[2]][[i]]=list(list(v0=(matrix(a30$motif,ncol=1)),v1=(matrix(a30$motif_derivative,ncol=1))),
                              list(v0=(matrix(b30$motif,ncol=1)),v1=(matrix(b30$motif_derivative,ncol=1))))
  
  a40=motifs_init(40,10,40,80,100)
  b40=motifs_init(40,10,40,80,100)
  V_initt[[1]][[3]][[i]]=list(list(v0=(matrix(a40$motif,ncol=1)),v1=(matrix(a40$motif_derivative,ncol=1))),
                              list(v0=(matrix(b40$motif,ncol=1)),v1=(matrix(b40$motif_derivative,ncol=1))))
  a50=motifs_init(50,10,40,80,100)
  b50=motifs_init(50,10,40,80,100)
  V_initt[[1]][[4]][[i]]=list(list(v0=(matrix(a50$motif,ncol=1)),v1=(matrix(a50$motif_derivative,ncol=1))),
                              list(v0=(matrix(b50$motif,ncol=1)),v1=(matrix(b50$motif_derivative,ncol=1))))
  a60=motifs_init(60,10,40,80,100)
  b60=motifs_init(60,10,40,80,100)
  V_initt[[1]][[5]][[i]]=list(list(v0=(matrix(a60$motif,ncol=1)),v1=(matrix(a60$motif_derivative,ncol=1))),
                              list(v0=(matrix(b60$motif,ncol=1)),v1=(matrix(b60$motif_derivative,ncol=1))))
  
  
  
  aa20=motifs_init(20,10,40,80,100)
  bb20=motifs_init(20,10,40,80,100)
  cc20=motifs_init(20,10,40,80,100)
  V_initt[[2]][[1]][[i]]=list(list(v0=matrix(aa20$motif,ncol=1),v1=matrix(aa20$motif_derivative,ncol=1)),
                              list(v0=matrix(bb20$motif,ncol=1),v1=matrix(bb20$motif_derivative,ncol=1)),
                              list(v0=matrix(cc20$motif,ncol=1),v1=matrix(cc20$motif_derivative,ncol=1)))
  
  aa30=motifs_init(30,10,40,80,100)
  bb30=motifs_init(30,10,40,80,100)
  cc30=motifs_init(30,10,40,80,100)
  V_initt[[2]][[2]][[i]]=list(list(v0=matrix(aa30$motif,ncol=1),v1=matrix(aa30$motif_derivative,ncol=1)),
                              list(v0=matrix(bb30$motif,ncol=1),v1=matrix(bb30$motif_derivative,ncol=1)),
                              list(v0=matrix(cc30$motif,ncol=1),v1=matrix(cc30$motif_derivative,ncol=1)))
  
  aa40=motifs_init(40,10,40,80,100)
  bb40=motifs_init(40,10,40,80,100)
  cc40=motifs_init(40,10,40,80,100)
  V_initt[[2]][[3]][[i]]=list(list(v0=matrix(aa40$motif,ncol=1),v1=matrix(aa40$motif_derivative,ncol=1)),
                              list(v0=matrix(bb40$motif,ncol=1),v1=matrix(bb40$motif_derivative,ncol=1)),
                              list(v0=matrix(cc40$motif,ncol=1),v1=matrix(cc40$motif_derivative,ncol=1)))
  aa50=motifs_init(50,10,40,80,100)
  bb50=motifs_init(50,10,40,80,100)
  cc50=motifs_init(50,10,40,80,100)
  V_initt[[2]][[4]][[i]]=list(list(v0=matrix(aa50$motif,ncol=1),v1=matrix(aa50$motif_derivative,ncol=1)),
                              list(v0=matrix(bb50$motif,ncol=1),v1=matrix(bb50$motif_derivative,ncol=1)),
                              list(v0=matrix(cc50$motif,ncol=1),v1=matrix(cc50$motif_derivative,ncol=1)))
  aa60=motifs_init(60,10,40,80,100)
  bb60=motifs_init(60,10,40,80,100)
  cc60=motifs_init(60,10,40,80,100)
  V_initt[[2]][[5]][[i]]=list(list(v0=matrix(aa60$motif,ncol=1),v1=matrix(aa60$motif_derivative,ncol=1)),
                              list(v0=matrix(bb60$motif,ncol=1),v1=matrix(bb60$motif_derivative,ncol=1)),
                              list(v0=matrix(cc60$motif,ncol=1),v1=matrix(cc60$motif_derivative,ncol=1)))
  
}


for(i in 5:8){
  a_2_20=motifs_init_rev(20,10,40,80,100)
  b_2_20=motifs_init_rev(20,10,40,80,100)
  V_initt[[1]][[1]][[i]]=list(list(v0=(matrix(a_2_20$motif,ncol=1)),v1=(matrix(a_2_20$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_2_20$motif,ncol=1)),v1=(matrix(b_2_20$motif_derivative,ncol=1))))
  a_2_30=motifs_init_rev(30,10,40,80,100)
  b_2_30=motifs_init_rev(30,10,40,80,100)
  V_initt[[1]][[2]][[i]]=list(list(v0=(matrix(a_2_30$motif,ncol=1)),v1=(matrix(a_2_30$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_2_30$motif,ncol=1)),v1=(matrix(b_2_30$motif_derivative,ncol=1))))
  a_2_40=motifs_init_rev(40,10,40,80,100)
  b_2_40=motifs_init_rev(40,10,40,80,100)
  V_initt[[1]][[3]][[i]]=list(list(v0=(matrix(a_2_40$motif,ncol=1)),v1=(matrix(a_2_40$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_2_40$motif,ncol=1)),v1=(matrix(b_2_40$motif_derivative,ncol=1))))
  a_2_50=motifs_init_rev(50,10,40,80,100)
  b_2_50=motifs_init_rev(50,10,40,80,100)
  V_initt[[1]][[4]][[i]]=list(list(v0=(matrix(a_2_50$motif,ncol=1)),v1=(matrix(a_2_50$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_2_50$motif,ncol=1)),v1=(matrix(b_2_50$motif_derivative,ncol=1))))
  a_2_60=motifs_init_rev(60,10,40,80,100)
  b_2_60=motifs_init_rev(60,10,40,80,100)
  V_initt[[1]][[5]][[i]]=list(list(v0=(matrix(a_2_60$motif,ncol=1)),v1=(matrix(a_2_60$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_2_60$motif,ncol=1)),v1=(matrix(b_2_60$motif_derivative,ncol=1))))
  
  
  
  aa_2_20=motifs_init_rev(20,10,40,80,100)
  bb_2_20=motifs_init_rev(20,10,40,80,100)
  cc_2_20=motifs_init_rev(20,10,40,80,100)
  V_initt[[2]][[1]][[i]]=list(list(v0=matrix(aa_2_20$motif,ncol=1),v1=matrix(aa_2_20$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_2_20$motif,ncol=1),v1=matrix(bb_2_20$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_2_20$motif,ncol=1),v1=matrix(cc_2_20$motif_derivative,ncol=1)))
  aa_2_30=motifs_init_rev(30,10,40,80,100)
  bb_2_30=motifs_init_rev(30,10,40,80,100)
  cc_2_30=motifs_init_rev(30,10,40,80,100)
  V_initt[[2]][[2]][[i]]=list(list(v0=matrix(aa_2_30$motif,ncol=1),v1=matrix(aa_2_30$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_2_30$motif,ncol=1),v1=matrix(bb_2_30$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_2_30$motif,ncol=1),v1=matrix(cc_2_30$motif_derivative,ncol=1)))
  aa_2_40=motifs_init_rev(40,10,40,80,100)
  bb_2_40=motifs_init_rev(40,10,40,80,100)
  cc_2_40=motifs_init_rev(40,10,40,80,100)
  V_initt[[2]][[3]][[i]]=list(list(v0=matrix(aa_2_40$motif,ncol=1),v1=matrix(aa_2_40$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_2_40$motif,ncol=1),v1=matrix(bb_2_40$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_2_40$motif,ncol=1),v1=matrix(cc_2_40$motif_derivative,ncol=1)))
  aa_2_50=motifs_init_rev(50,10,40,80,100)
  bb_2_50=motifs_init_rev(50,10,40,80,100)
  cc_2_50=motifs_init_rev(50,10,40,80,100)
  V_initt[[2]][[4]][[i]]=list(list(v0=matrix(aa_2_50$motif,ncol=1),v1=matrix(aa_2_50$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_2_50$motif,ncol=1),v1=matrix(bb_2_50$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_2_50$motif,ncol=1),v1=matrix(cc_2_50$motif_derivative,ncol=1)))
  aa_2_60=motifs_init_rev(60,10,40,80,100)
  bb_2_60=motifs_init_rev(60,10,40,80,100)
  cc_2_60=motifs_init_rev(60,10,40,80,100)
  V_initt[[2]][[5]][[i]]=list(list(v0=matrix(aa_2_60$motif,ncol=1),v1=matrix(aa_2_60$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_2_60$motif,ncol=1),v1=matrix(bb_2_60$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_2_60$motif,ncol=1),v1=matrix(cc_2_60$motif_derivative,ncol=1)))
  
  
}

#Increasing lines
for(i in 9:9){
  a_3_20=motifs_line(20,10,40,80,100)
  b_3_20=motifs_line(20,10,40,80,100)
  V_initt[[1]][[1]][[i]]=list(list(v0=(matrix(a_3_20$motif,ncol=1)),v1=(matrix(a_3_20$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_3_20$motif,ncol=1)),v1=(matrix(b_3_20$motif_derivative,ncol=1))))
  
  a_3_30=motifs_line(30,10,40,80,100)
  b_3_30=motifs_line(30,10,40,80,100)
  V_initt[[1]][[2]][[i]]=list(list(v0=(matrix(a_3_30$motif,ncol=1)),v1=(matrix(a_3_30$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_3_30$motif,ncol=1)),v1=(matrix(b_3_30$motif_derivative,ncol=1))))
  
  a_3_40=motifs_line(40,10,40,80,100)
  b_3_40=motifs_line(40,10,40,80,100)
  V_initt[[1]][[3]][[i]]=list(list(v0=(matrix(a_3_40$motif,ncol=1)),v1=(matrix(a_3_40$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_3_40$motif,ncol=1)),v1=(matrix(b_3_40$motif_derivative,ncol=1))))
  a_3_50=motifs_line(50,10,40,80,100)
  b_3_50=motifs_line(50,10,40,80,100)
  V_initt[[1]][[4]][[i]]=list(list(v0=(matrix(a_3_50$motif,ncol=1)),v1=(matrix(a_3_50$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_3_50$motif,ncol=1)),v1=(matrix(b_3_50$motif_derivative,ncol=1))))
  a_3_60=motifs_line(60,10,40,80,100)
  b_3_60=motifs_line(60,10,40,80,100)
  V_initt[[1]][[5]][[i]]=list(list(v0=(matrix(a_3_60$motif,ncol=1)),v1=(matrix(a_3_60$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_3_60$motif,ncol=1)),v1=(matrix(b_3_60$motif_derivative,ncol=1))))
  
  
  aa_3_20=motifs_line(20,10,40,80,100)
  bb_3_20=motifs_line(20,10,40,80,100)
  cc_3_20=motifs_line(20,10,40,80,100)
  V_initt[[2]][[1]][[i]]=list(list(v0=matrix(aa_3_20$motif,ncol=1),v1=matrix(aa_3_20$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_3_20$motif,ncol=1),v1=matrix(bb_3_20$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_3_20$motif,ncol=1),v1=matrix(cc_3_20$motif_derivative,ncol=1)))
  
  aa_3_30=motifs_line(30,10,40,80,100)
  bb_3_30=motifs_line(30,10,40,80,100)
  cc_3_30=motifs_line(30,10,40,80,100)
  V_initt[[2]][[2]][[i]]=list(list(v0=matrix(aa_3_30$motif,ncol=1),v1=matrix(aa_3_30$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_3_30$motif,ncol=1),v1=matrix(bb_3_30$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_3_30$motif,ncol=1),v1=matrix(cc_3_30$motif_derivative,ncol=1)))
  
  aa_3_40=motifs_line(40,10,40,80,100)
  bb_3_40=motifs_line(40,10,40,80,100)
  cc_3_40=motifs_line(40,10,40,80,100)
  V_initt[[2]][[3]][[i]]=list(list(v0=matrix(aa_3_40$motif,ncol=1),v1=matrix(aa_3_40$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_3_40$motif,ncol=1),v1=matrix(bb_3_40$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_3_40$motif,ncol=1),v1=matrix(cc_3_40$motif_derivative,ncol=1)))
  aa_3_50=motifs_line(50,10,40,80,100)
  bb_3_50=motifs_line(50,10,40,80,100)
  cc_3_50=motifs_line(50,10,40,80,100)
  V_initt[[2]][[4]][[i]]=list(list(v0=matrix(aa_3_50$motif,ncol=1),v1=matrix(aa_3_50$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_3_50$motif,ncol=1),v1=matrix(bb_3_50$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_3_50$motif,ncol=1),v1=matrix(cc_3_50$motif_derivative,ncol=1)))
  aa_3_60=motifs_line(60,10,40,80,100)
  bb_3_60=motifs_line(60,10,40,80,100)
  cc_3_60=motifs_line(60,10,40,80,100)
  V_initt[[2]][[5]][[i]]=list(list(v0=matrix(aa_3_60$motif,ncol=1),v1=matrix(aa_3_60$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_3_60$motif,ncol=1),v1=matrix(bb_3_60$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_3_60$motif,ncol=1),v1=matrix(cc_3_60$motif_derivative,ncol=1)))
  
  
}

#Decreasing lines
for(i in 10:10){
  a_4_20=motifs_line_rev(20,10,40,80,100)
  b_4_20=motifs_line_rev(20,10,40,80,100)
  V_initt[[1]][[1]][[i]]=list(list(v0=(matrix(a_4_20$motif,ncol=1)),v1=(matrix(a_4_20$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_4_20$motif,ncol=1)),v1=(matrix(b_4_20$motif_derivative,ncol=1))))
  
  a_4_30=motifs_line_rev(30,10,40,80,100)
  b_4_30=motifs_line_rev(30,10,40,80,100)
  V_initt[[1]][[2]][[i]]=list(list(v0=(matrix(a_4_30$motif,ncol=1)),v1=(matrix(a_4_30$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_4_30$motif,ncol=1)),v1=(matrix(b_4_30$motif_derivative,ncol=1))))
  
  a_4_40=motifs_line_rev(40,10,40,80,100)
  b_4_40=motifs_line_rev(40,10,40,80,100)
  V_initt[[1]][[3]][[i]]=list(list(v0=(matrix(a_4_40$motif,ncol=1)),v1=(matrix(a_4_40$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_4_40$motif,ncol=1)),v1=(matrix(b_4_40$motif_derivative,ncol=1))))
  a_4_50=motifs_line_rev(50,10,40,80,100)
  b_4_50=motifs_line_rev(50,10,40,80,100)
  V_initt[[1]][[4]][[i]]=list(list(v0=(matrix(a_4_50$motif,ncol=1)),v1=(matrix(a_4_50$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_4_50$motif,ncol=1)),v1=(matrix(b_4_50$motif_derivative,ncol=1))))
  a_4_60=motifs_line_rev(60,10,40,80,100)
  b_4_60=motifs_line_rev(60,10,40,80,100)
  V_initt[[1]][[5]][[i]]=list(list(v0=(matrix(a_4_60$motif,ncol=1)),v1=(matrix(a_4_60$motif_derivative,ncol=1))),
                              list(v0=(matrix(b_4_60$motif,ncol=1)),v1=(matrix(b_4_60$motif_derivative,ncol=1))))
  
  
  
  aa_4_20=motifs_line_rev(20,10,40,80,100)
  bb_4_20=motifs_line_rev(20,10,40,80,100)
  cc_4_20=motifs_line_rev(20,10,40,80,100)
  V_initt[[2]][[1]][[i]]=list(list(v0=matrix(aa_4_20$motif,ncol=1),v1=matrix(aa_4_20$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_4_20$motif,ncol=1),v1=matrix(bb_4_20$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_4_20$motif,ncol=1),v1=matrix(cc_4_20$motif_derivative,ncol=1)))
  
  aa_4_30=motifs_line_rev(30,10,40,80,100)
  bb_4_30=motifs_line_rev(30,10,40,80,100)
  cc_4_30=motifs_line_rev(30,10,40,80,100)
  V_initt[[2]][[2]][[i]]=list(list(v0=matrix(aa_4_30$motif,ncol=1),v1=matrix(aa_4_30$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_4_30$motif,ncol=1),v1=matrix(bb_4_30$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_4_30$motif,ncol=1),v1=matrix(cc_4_30$motif_derivative,ncol=1)))
  
  aa_4_40=motifs_line_rev(40,10,40,80,100)
  bb_4_40=motifs_line_rev(40,10,40,80,100)
  cc_4_40=motifs_line_rev(40,10,40,80,100)
  V_initt[[2]][[3]][[i]]=list(list(v0=matrix(aa_4_40$motif,ncol=1),v1=matrix(aa_4_40$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_4_40$motif,ncol=1),v1=matrix(bb_4_40$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_4_40$motif,ncol=1),v1=matrix(cc_4_40$motif_derivative,ncol=1)))
  aa_4_50=motifs_line_rev(50,10,40,80,100)
  bb_4_50=motifs_line_rev(50,10,40,80,100)
  cc_4_50=motifs_line_rev(50,10,40,80,100)
  V_initt[[2]][[4]][[i]]=list(list(v0=matrix(aa_4_50$motif,ncol=1),v1=matrix(aa_4_50$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_4_50$motif,ncol=1),v1=matrix(bb_4_50$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_4_50$motif,ncol=1),v1=matrix(cc_4_50$motif_derivative,ncol=1)))
  aa_4_60=motifs_line_rev(60,10,40,80,100)
  bb_4_60=motifs_line_rev(60,10,40,80,100)
  cc_4_60=motifs_line_rev(60,10,40,80,100)
  V_initt[[2]][[5]][[i]]=list(list(v0=matrix(aa_4_60$motif,ncol=1),v1=matrix(aa_4_60$motif_derivative,ncol=1)),
                              list(v0=matrix(bb_4_60$motif,ncol=1),v1=matrix(bb_4_60$motif_derivative,ncol=1)),
                              list(v0=matrix(cc_4_60$motif,ncol=1),v1=matrix(cc_4_60$motif_derivative,ncol=1)))
  
  
}


V_init=V_initt
K = c(2, 3) # number of clusters to try
c=c(20,30,40,50,60)
n_init = 10 # number of partially random initializations to try
# NOTE: rename "results" folder to re-run everything
#(TIME CONSUMING)

set.seed(13333)
files = list.files('./stock_prices_motifs_shorter')

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

# The following part of code is very time consuming, especially if run in sequential mode
# (worker_number = 1)
if('motifs_candidate.RData' %in% files){
  #if candidate motifs already present, load them
  load('./stock_prices_motifs_shorter/motifs_candidate.RData')
}else{
  # find candidate motifs with partially random initialization
find_candidate_motifs_results = find_candidate_motifs(Y0, Y1, K, c, n_init,V_init=V_init,
                                                      name = './stock_prices_motifs_shorter/motifs', names_var = 'x(t)',
                                                      probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = 1000,
                                                                             iter4elong = iter4elong, trials_elong = trials_elong, max_gap = max_gap,
                                                                             return_options = TRUE, return_init = TRUE,
                                                                             diss = diss, alpha = alpha,transformed=TRUE),
                                                      plot = TRUE, worker_number = 1)



# find candidate motifs with random initialization
n_init = 20 # number of partially random plus random initialization to try
find_candidate_motifs_results = find_candidate_motifs(Y0, Y1, K, c, n_init,V_init=NULL,
                                                      name = './stock_prices_motifs_shorter/motifs', names_var = 'x(t)',
                                                      probKMA_options = list(c_max = c_max, standardize = FALSE, iter_max = 1000,
                                                                             iter4elong = iter4elong, trials_elong = trials_elong, max_gap = max_gap,
                                                                             return_options = TRUE, return_init = TRUE,
                                                                             diss = diss, alpha = alpha,transformed=TRUE),
                                                      plot = TRUE, worker_number = 1)
save(find_candidate_motifs_results, file = './stock_prices_motifs_shorter/motifs_candidate.RData')
}

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

#save(cluster_candidate_motifs_results, file = './stock_prices_motifs_shorter/cluster_candidate_motifs_results.RData')


### plot cluster candidate motifs results
pdf('./stock_prices_motifs_shorter/clustering_candidate_motifs.pdf', height = 12, width = 9)
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
pdf('./stock_prices_motifs_shorter/motifs_search_results.pdf', height = 18, width = 30)
motifs_search_plot_norm_time(motifs_search_results, ylab = 'x(t)', freq_threshold = 5,transformed=TRUE,data_smoothed=data_smoothed_all_stocks)
dev.off()



#Motifs analysis for fmd including shorter motifs
column_names=c("Motif","Length","Frequency","Radius")
motifs_analysis=data.frame(c(1:68),motifs_search_results$V_length,
                           motifs_search_results$V_frequencies,
                           motifs_search_results$R_motifs)
colnames(motifs_analysis)=column_names


#Considering only the motifs with more than 5 occurrences:
#filtered_motifs_analysis=motifs_analysis[motifs_analysis$Frequency >= 5, ]
# Renumber the Motif column
#filtered_motifs_analysis$Motif=seq_len(nrow(filtered_motifs_analysis))

motifs_analysis_ordered=motifs_analysis[order(motifs_analysis$Length),]


library("writexl")
motifs_analysis_table=data.frame(Motif=motifs_analysis_ordered$Motif,
                                 Length=motifs_analysis_ordered$Length,
                                 Frequency=motifs_analysis_ordered$Frequency,
                                 Radius=motifs_analysis_ordered$Radius)

write_xlsx(motifs_analysis_table,"./stock_prices_motifs_shorter/simulated_motifs_analysis.xlsx")



#Let's multiply the radius for 1000
motifs_analysis_table_latex=data.frame(Motif=motifs_analysis_ordered$Motif,
                                       Length=motifs_analysis_ordered$Length,
                                       Frequency=motifs_analysis_ordered$Frequency,
                                       Radius=motifs_analysis_ordered$Radius*100000)
library("xtable")
print(xtable(motifs_analysis_table_latex),digits=c(5,5,5,5),include.rownames = FALSE)


