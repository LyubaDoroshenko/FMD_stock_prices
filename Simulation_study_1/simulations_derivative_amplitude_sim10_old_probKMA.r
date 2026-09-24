require(fda)
require(combinat)
require(Hmisc)

# Set directory to the folder with the current file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("../fmd_functions.R")


folder='N20_K2-10-10_derivative_amplitude_old_probKMA/'
sim_simulation=paste0('simulation',1:10)
len_simulation=c(200,300,400,500)
sd_simulation=c(0.1,0.5,1,2)

N=20
nmotifs=2
norder=3
dist_knots=10
len_motifs=60
freq_motifs=list(c1=c(1,0),c2=c(1,0),c3=c(1,0),c4=c(1,0),c5=c(1,0),c6=c(1,0),
                 c7=c(0,1),c8=c(0,1),c9=c(0,1),c10=c(0,1),c11=c(0,1),c12=c(0,1),
                 c13=c(1,1),c14=c(1,1),
                 c15=c(2,0),c16=c(2,0),
                 c17=c(0,2),c18=c(0,2),
                 c19=c(0,0),c20=c(0,0))



#######################
###     RUN FMD     ###
#######################

standardize=FALSE
transformed=FALSE
iter_max=1000
max_gap=0 # no gaps allowed
return_options=TRUE
return_init=TRUE
diss='d1_L2'
alpha=1

sim=sim_simulation[10]
# # run functional motif discovery 10 times
# for(run in 1:5){
#   
#   
#   # different curve lengths
#   j=1
#   len=len_simulation[j]
#   
#   # different levels of noise
#   for(jj in seq_along(sd_simulation)){
#     set.seed(jj)
#     sd_noise=sd_simulation[jj]
#     load(paste0(folder,sim,'/len',len,'_sd',sd_noise,'.RData'))
#     
#     # evaluate curves
#     Y0=evaluate_curves(curves$fd_curves,by=1,Lfdobj=0)
#     Y1=evaluate_curves(curves$fd_curves,by=1,Lfdobj=1)
#     
#     iter4elong=1 # perform elongation
#     trials_elong=length(Y0[[1]]) # try all possible elongations
#     c_max=len_motifs+1+10 
#     
#     K=c(2,3)
#     c=len_motifs+1-c(0,10,20)
#     n_init=10
#     worker_number=32
#     
#     files=list.files(paste0(paste0(folder,sim)))
#     if(paste0('len',len,'_sd',sd_noise,'_run',run,'_candidate.RData') %in% files){
#       # candidate motifs already present, load them
#       load(paste0(folder,sim,'/len',len,'_sd',sd_noise,'_run',run,'_candidate.RData'))
#     }else{
#       message(paste0('len',len,'_sd',sd_noise))
#       start=proc.time()
#       # find candidate motifs
#       find_candidate_motifs_results=find_candidate_motifs(Y0,Y1,K,c,n_init, 
#                                                           name=paste0(folder,sim,'/len',len,'_sd',sd_noise,'_run',run),
#                                                           probKMA_options=list(c_max=c_max,standardize=standardize,transformed=transformed,iter_max=iter_max,
#                                                                                iter4elong=iter4elong,trials_elong=trials_elong,max_gap=max_gap,
#                                                                                return_options=return_options,return_init=return_init,
#                                                                                diss=diss,alpha=alpha,worker_number=1),
#                                                           plot=TRUE,worker_number=worker_number)
#       end=proc.time()
#       time_find=end-start
#       save(find_candidate_motifs_results,time_find,file=paste0(folder,sim,'/len',len,'_sd',sd_noise,'_run',run,'_candidate.RData'))
#     }
#     
#     start=proc.time()
#     # retrieve silhouette average for all candidate motifs, in order to select a threshold based on their distribution (0.9 quantile)
#     silhouette_average=Reduce(rbind,Reduce(rbind,find_candidate_motifs_results$silhouette_average_sd))[,1]
#     # filter candidate motifs based on silhouette average and size
#     filter_candidate_motifs_results=filter_candidate_motifs(find_candidate_motifs_results,
#                                                             sil_threshold=quantile(silhouette_average,0.9),
#                                                             size_threshold=5)
#     # cluster candidate motifs based on their distance and select radii
#     cluster_candidate_motifs_results=cluster_candidate_motifs(filter_candidate_motifs_results,
#                                                               motif_overlap=0.6,worker_number=worker_number)
#     # plot cluster candidate motifs results
#     pdf(paste0(folder,sim,'/len',len,'_sd',sd_noise,'_run',run,'_clustering_candidate_motifs.pdf'),height=12,width=9)
#     cluster_candidate_motifs_plot(cluster_candidate_motifs_results,ask=FALSE)
#     dev.off()
#     # search selected motifs
#     motifs_search_results=motifs_search(cluster_candidate_motifs_results,
#                                         use_real_occurrences=FALSE,length_diff=+Inf,worker_number=worker_number)
#     
#     # plot FMD results (NB: no threshold of frequencies of motif found!)
#     pdf(paste0(folder,sim,'/len',len,'_sd',sd_noise,'_run',run,'_results.pdf'),height=7,width=17)
#     motifs_search_plot_norm(motifs_search_results,ylab='x(t)',freq_threshold=1,transformed=TRUE)
#     dev.off()
#     end=proc.time()
#     time=time_find+(end-start)
#     
#     save(find_candidate_motifs_results,silhouette_average,filter_candidate_motifs_results,
#          cluster_candidate_motifs_results,motifs_search_results,time_find,time,
#          file=paste0(folder,sim,'/len',len,'_sd',sd_noise,'_run',run,'.RData'))
#   }
# }


#######################
### EVALUATE RESULT ###
#######################

matching='frequency' # matching between motifs found and real motifs, based on 'frequency' (false negative) or 'distance'
overlap_D=0.5 # minimum overlap (in percentage) between motifs found and real motifs, in order to compute their distance
overlap_freq=0.5 # minimum overlap (in percentage) between occurrences found and real occurrences, in order to compute true/false positives

##### arrays to store results #####
# number of motifs found
nmotifs_found=array(data=NA,dim=c(length(sim_simulation),length(len_simulation),length(sd_simulation),10),
                    dimnames=list(sim_simulation,len_simulation,sd_simulation,1:10))
# running time
time_run=array(data=NA,dim=c(length(sim_simulation),length(len_simulation),length(sd_simulation),10),
               dimnames=list(sim_simulation,len_simulation,sd_simulation,1:10))
# frequency (true positives) of each motif found (only the ones matched to real motifs)
freq_motifs_found=array(data=NA,dim=c(length(sim_simulation),length(len_simulation),length(sd_simulation),10,nmotifs),
                        dimnames=list(sim_simulation,len_simulation,sd_simulation,1:10,1:nmotifs))
# wrong frequency (false positives) of each motif found (only the ones matched to real motifs)
wrong_freq_motifs_found=array(data=NA,dim=c(length(sim_simulation),length(len_simulation),length(sd_simulation),10,nmotifs),
                              dimnames=list(sim_simulation,len_simulation,sd_simulation,1:10,1:nmotifs))
# length of each motif found (only the ones matched to real motifs)
len_motifs_found=array(data=NA,dim=c(length(sim_simulation),length(len_simulation),length(sd_simulation),10,nmotifs),
                       dimnames=list(sim_simulation,len_simulation,sd_simulation,1:10,1:nmotifs))
# distance between each real motif and its matched motif found
MV_dist=array(data=NA,dim=c(length(sim_simulation),length(len_simulation),length(sd_simulation),10,nmotifs),
              dimnames=list(sim_simulation,len_simulation,sd_simulation,1:10,1:nmotifs))
# alignment (shift) between each real motif and its matched motif found
MV_shift=array(data=NA,dim=c(length(sim_simulation),length(len_simulation),length(sd_simulation),10,nmotifs),
               dimnames=list(sim_simulation,len_simulation,sd_simulation,1:10,1:nmotifs))
# permutation of motifs found in order to match real motifs
motifs_found_permut=array(data=NA,dim=c(length(sim_simulation),length(len_simulation),length(sd_simulation),10,nmotifs),
                          dimnames=list(sim_simulation,len_simulation,sd_simulation,1:10,1:nmotifs))
# distance between each real motif and the average of all its occurrences 
# (this is usually different from 0 since the noise is Gaussian in the spline coefficients and not in the curves)
MM_average_dist=array(data=NA,dim=c(length(sim_simulation),length(len_simulation),length(sd_simulation),nmotifs),
                      dimnames=list(sim_simulation,len_simulation,sd_simulation,1:nmotifs))
# distance between each real motif and the average of all its occurrences, considering only their starting part (from position 1 to len_occ)
# this is to see how the distance varies with the length of the considered region
M_occurrences_mean_dist=array(data=NA,dim=c(length(sim_simulation),length(len_simulation),length(sd_simulation),nmotifs,61),
                              dimnames=list(sim_simulation,len_simulation,sd_simulation,1:nmotifs,1:61))

##### store results and plot #####
if(alpha==1){
  use0=FALSE
  use1=TRUE
}else{
  use0=TRUE
  use1=TRUE
}


sim=sim_simulation[10]
# different curve lengths
j=1
len=len_simulation[j]
# different levels of noise
for(jj in seq_along(sd_simulation)){
  sd_noise=sd_simulation[jj]
  load(paste0(folder,sim,'/len',len,'_sd',sd_noise,'.RData'))
  M=vector('list',nmotifs)
  M_average=vector('list',nmotifs)
  M_occurrences=vector('list',nmotifs)
  for(motif in 1:nmotifs){
    breaks=seq(0,(length(curves$coeff_motifs[[motif]])-norder+1+2*(norder-1))*dist_knots,by=dist_knots)-(norder-1)*dist_knots
    basis=create.bspline.basis(norder=norder,breaks=breaks)
    curve=fd(coef=c(rep(0,norder-1),curves$coeff_motifs[[motif]],rep(0,norder-1)),basisobj=basis)
    # real motif (0=curve, 1=derivative)
    M[[motif]]$m0=eval.fd(seq(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=1),curve,Lfdobj=0)
    M[[motif]]$m1=eval.fd(seq(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=1),curve,Lfdobj=1)
    if(transformed){
      M_min=min(M[[motif]]$m0)
      M_max=max(M[[motif]]$m0)
      M[[motif]]$m0=(M[[motif]]$m0-M_min)/(M_max-M_min)
      M[[motif]]$m1=M[[motif]]$m1/(M_max-M_min)
    }
    # real motif occurrences in each curve (0=curve, 1=derivative)
    xy0=mapply(function(fd_curve,motifs_in_curves_i){
      x=lapply(motifs_in_curves_i$pos_motifs[motifs_in_curves_i$id_motifs==motif]-1,function(start) seq(start*dist_knots,(start+length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=1))
      if(is.null(x))
        return(NULL)
      xy=lapply(x,function(x) cbind(x,eval.fd(x,fd_curve,Lfdobj=0)))
      return(xy)
    },curves$fd_curves,curves$motifs_in_curves,SIMPLIFY=FALSE)
    xy1=mapply(function(fd_curve,motifs_in_curves_i){
      x=lapply(motifs_in_curves_i$pos_motifs[motifs_in_curves_i$id_motifs==motif]-1,function(start) seq(start*dist_knots,(start+length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=1))
      if(is.null(x))
        return(NULL)
      xy=lapply(x,function(x) cbind(x,eval.fd(x,fd_curve,Lfdobj=1)))
      return(xy)
    },curves$fd_curves,curves$motifs_in_curves,SIMPLIFY=FALSE)
    M_occurrences[[motif]]=mapply(function(m0,m1) list(m0=as.matrix(m0[,2]),m1=as.matrix(m1[,2])),
                                  unlist(xy0,recursive=FALSE),unlist(xy1,recursive=FALSE),SIMPLIFY=FALSE)
    # average of real motif occurrences (0=curve, 1=derivative)
    if(transformed){
      M_average[[motif]]$m0=as.matrix(colMeans(Reduce(rbind,lapply(xy0,
                                                                   function(xy){
                                                                     if(!is.null(xy))
                                                                       Reduce(rbind,lapply(xy,function(xy) (xy[,2]-min(xy[,2]))/(max(xy[,2])-min(xy[,2]))))
                                                                   }))))
      M_average[[motif]]$m1=as.matrix(colMeans(Reduce(rbind,mapply(function(xy,xy_0){
        if(!is.null(xy))
          Reduce(rbind,mapply(function(xy,xy_0) xy[,2]/(max(xy_0[,2])-min(xy_0[,2])),xy,xy_0,SIMPLIFY = FALSE))
      },xy1,xy0, SIMPLIFY = FALSE))))
    }else{
      M_average[[motif]]$m0=as.matrix(colMeans(Reduce(rbind,lapply(xy0,
                                                                   function(xy){
                                                                     if(!is.null(xy))
                                                                       Reduce(rbind,lapply(xy,function(xy) xy[,2]))
                                                                   }))))
      M_average[[motif]]$m1=as.matrix(colMeans(Reduce(rbind,lapply(xy1,
                                                                   function(xy){
                                                                     if(!is.null(xy))
                                                                       Reduce(rbind,lapply(xy,function(xy) xy[,2]))
                                                                   }))))
    }
    # distance between real motif and the average of all its occurrences 
    MM_average_dist[sim,j,jj,motif]=.find_min_diss(M[[motif]],M_average[[motif]],transform_y = transformed, transform_v = transformed,
                                                   alpha=alpha,w=1,c_k=length(M[[motif]]$m0),d=1,use0=use0,use1=use1)[2]
    # distance between real motif and the average of all its occurrences, considering only their starting part (from position 1 to len_occ)
    for(len_occ in 1:61){
      M_occurrences_mean_dist[sim,j,jj,motif,len_occ]=mean(Reduce(rbind,lapply(lapply(M_occurrences[[motif]],function(m) list(m0=as.matrix(m$m0[1:len_occ,]),m1=as.matrix(m$m1[1:len_occ,]))),
                                                                               .find_min_diss,v=list(m0=as.matrix(M_average[[motif]]$m0[1:len_occ,]),m1=as.matrix(M_average[[motif]]$m1[1:len_occ,])),transform_y = transformed, transform_v = transformed,
                                                                               alpha=alpha,w=1,c_k=len_occ,d=1,use0=use0,use1=use1))[,2])
    }
  }
  for(run in 1:5){
    if(paste0('len',len,'_sd',sd_noise,'_run',run,'.RData') %in% list.files(paste0(folder,sim))){
      load(paste0(folder,sim,'/len',len,'_sd',sd_noise,'_run',run,'.RData'))
      # number of motifs found (NB: no threshold of frequencies!)
      nmotifs_found[sim,j,jj,run]=length(motifs_search_results$V_frequencies)
      # running time
      time_run[sim,j,jj,run]=time[3]
      # motifs found (v0=curve, v1=derivatie)
      V=mapply(function(v0,v1) list(v0=v0,v1=v1),
               motifs_search_results$V0[1:nmotifs_found[sim,j,jj,run]],
               motifs_search_results$V1[1:nmotifs_found[sim,j,jj,run]],SIMPLIFY=FALSE)
      # compute distances (and shift) between each of the motifs found and each of the real motifs
      # require that their intersection is at least overlap_D % of the shortest in each pair
      VV=rep(V,each=nmotifs)
      VV_lengths=unlist(lapply(VV,function(vv) nrow(vv$v0)))
      MM=rep(M,nmotifs_found[sim,j,jj,run])
      MM_lengths=unlist(lapply(MM,function(mm) nrow(mm$m0)))
      swap=MM_lengths<VV_lengths # swap order between MM and VV to find their minimum distance (align the shortest to the longest)
      MM_swap=MM
      MM_swap[swap]=VV[swap]
      VV[swap]=MM[swap]
      SD=mapply(.find_min_diss,MM_swap,VV,round(pmin(MM_lengths,VV_lengths)*overlap_D),
                MoreArgs=list(alpha=filter_candidate_motifs_results$alpha,transform_y = transformed, transform_v = transformed,w=filter_candidate_motifs_results$w,d=1,use0=use0,use1=use1),SIMPLIFY=TRUE)
      MV_D=matrix(SD[2,],nrow=nmotifs,ncol=nmotifs_found[sim,j,jj,run]) # distance
      MV_S=matrix(SD[1,],nrow=nmotifs,ncol=nmotifs_found[sim,j,jj,run]) # shift
      if(nmotifs_found[sim,j,jj,run]<nmotifs){ # add NA if the motifs found are less than the real motifs
        MV_D=cbind(MV_D,matrix(NA,nrow=nmotifs,ncol=nmotifs-nmotifs_found[sim,j,jj,run]))
        MV_S=cbind(MV_S,matrix(NA,nrow=nmotifs,ncol=nmotifs-nmotifs_found[sim,j,jj,run]))
      }
      # count the number of occurrences of each of the motifs found that match the actual occurrences of the real motifs
      # require that their intersection is at least overlap_freq of the shortest in each pair
      real_motifs_in_curves=Reduce(rbind,lapply(seq_along(curves$motifs_in_curves),
                                                function(i_curve){
                                                  if(is.null(curves$motifs_in_curves[[i_curve]])){
                                                    c()
                                                  }else{
                                                    cbind(i_curve,(curves$motifs_in_curves[[i_curve]]$pos_motifs-1)*dist_knots,curves$motifs_in_curves[[i_curve]]$id_motifs)
                                                  }
                                                }))
      colnames(real_motifs_in_curves)=c('curve','shift','motif')
      real_motifs_in_curves=as.data.frame(real_motifs_in_curves)
      real_freq=rep(NA,nmotifs)
      MV_freq=matrix(NA,nrow=nmotifs,ncol=max(nmotifs_found[sim,j,jj,run],nmotifs)) # real frequencies (true positives)
      for(motif in 1:nmotifs){
        real_motif_in_curves=real_motifs_in_curves[real_motifs_in_curves$motif==motif,1:2]
        real_freq[motif]=nrow(real_motif_in_curves)
        real_length=length(M[[motif]]$m0)
        for(motif_found in 1:nmotifs_found[sim,j,jj,run]){
          found_motif_in_curves=matrix(motifs_search_results$V_occurrences[[motif_found]][,1:2],ncol=2)
          found_length=motifs_search_results$V_length[motif_found]
          min_inter_length=round(min(real_length,found_length)*overlap_freq)
          real_occurrences_found=apply(real_motif_in_curves,1,
                                       function(real){
                                         TRUE %in% unlist(lapply(found_motif_in_curves[found_motif_in_curves[,1]==real[1],2],
                                                                 function(x) length(intersect(seq(real[2],real[2]+real_length),seq(x,x+found_length)))>=min_inter_length))
                                       })
          MV_freq[motif,motif_found]=sum(real_occurrences_found)
        }
      }
      MV_freq_wrong=matrix(motifs_search_results$V_frequencies,nrow=nmotifs,ncol=max(nmotifs_found[sim,j,jj,run],nmotifs),byrow=TRUE)-MV_freq # wrong frequencies (false positives)
      # try all possible matching between motifs found and real motifs
      tot_MV_freq=unlist(lapply(permn(1:min(max(nmotifs_found[sim,j,jj,run],nmotifs),8)),function(index) sum(real_freq-diag(MV_freq[,index]),na.rm=TRUE)))
      # try all possible matching between motifs found and real motifs
      tot_MV_D=unlist(lapply(permn(1:min(max(nmotifs_found[sim,j,jj,run],nmotifs),8)),function(index) sum(diag(MV_D[,index]),na.rm=TRUE)))
      # select the permutation that produce less false negative
      min_permut_freq=which(min(tot_MV_freq)==tot_MV_freq)
      if(length(min_permut_freq)>1){
        min_permut_freq=min_permut_freq[which.min(tot_MV_D[min_permut_freq])]
      }
      # select the permutation that produce a smaller distance
      min_permut_D=which.min(tot_MV_D) 
      if(min_permut_freq!=min_permut_D){
        message(sim,' len',len,' sd',sd_noise,' run',run,' - frequency and distance select different matching')
      }
      if(matching=='frequency'){
        min_permut=min_permut_freq
      }else if(matching=='distance'){
        min_permut=min_permut_D
      }else{
        stop('matching should be frequency or distance')
      }
      # permutation
      motifs_found_permut[sim,j,jj,run,]=permn(1:min(max(nmotifs_found[sim,j,jj,run],nmotifs),8))[[min_permut]][1:nmotifs] 
      # frequency of motifs found (true positives)
      freq_motifs_found[sim,j,jj,run,]=diag(MV_freq[,motifs_found_permut[sim,j,jj,run,]]) 
      # wrong frequency of motifs found (false positives)
      wrong_freq_motifs_found[sim,j,jj,run,]=diag(MV_freq_wrong[,motifs_found_permut[sim,j,jj,run,]]) 
      # distance between found and real motifs
      MV_dist[sim,j,jj,run,]=diag(MV_D[,motifs_found_permut[sim,j,jj,run,]]) 
      if(nmotifs_found[sim,j,jj,run]<nmotifs){ # add NA if the motifs found are less than the real motifs
        motifs_found_permut[sim,j,jj,run,is.na(MV_dist[sim,j,jj,run,])]=NA
      }
      # shift that produces the smallest distance
      MV_shift[sim,j,jj,run,]=diag(((MV_S-1)*(-swap)+(MV_S-1)*(!swap))[,motifs_found_permut[sim,j,jj,run,]])
      # length of motifs found
      len_motifs_found[sim,j,jj,run,]=motifs_search_results$V_length[motifs_found_permut[sim,j,jj,run,]]
      
      if(transformed){
        ylim=c(0,1)
      }else{
        ylim=c(-15,20)
      }
      pdf(paste0(folder,sim,'/len',len,'_sd',sd_noise,'_run',run,'_motifs_found.pdf'),10,10)
      for(motif in 1:nmotifs){
        plot(1,type='n',xlab='t',ylab='x(t)',xlim=c(-10,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots+10),ylim=ylim,
             main=paste('Motif',motif),cex.main=2,cex.axis=1.5,cex.lab=1.5)
        breaks=seq(0,(length(curves$coeff_motifs[[motif]])-norder+1+2*(norder-1))*dist_knots,by=dist_knots)-(norder-1)*dist_knots
        basis=create.bspline.basis(norder=norder,breaks=breaks)
        curve=fd(coef=c(rep(0,norder-1),curves$coeff_motifs[[motif]],rep(0,norder-1)),basisobj=basis)
        if(!is.na(motifs_found_permut[sim,j,jj,run,motif]))
          lines(1:nrow(motifs_search_results$V0[[motifs_found_permut[sim,j,jj,run,motif]]])-1+MV_shift[sim,j,jj,run,motif],motifs_search_results$V0[[motifs_found_permut[sim,j,jj,run,motif]]],col=1+motif,lty=2,lwd=4)
        curve_motif=eval.fd(seq(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=0.5),curve,Lfdobj=0)
        if(transformed){
          curve_motif=(curve_motif-min(curve_motif))/(max(curve_motif)-min(curve_motif))
        }
        lines(seq(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=0.5),curve_motif,col=1,lwd=2)
        lines(seq(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=1),M_average[[motif]]$m0,col=1,lwd=2,lty=3)
        legend('topleft',legend=c('Real','Sample mean','Estimated'),col=c(1,1,motif+1),lty=c(1,3,2),lwd=c(2,2,4),cex=1.5)
      }
      dev.off()
      
      pdf(paste0(folder,sim,'/len',len,'_sd',sd_noise,'_run',run,'_curves_with_motifs.pdf'),10,5*N)
      par(mfrow=c(N,1))
      freq=freq_motifs_found[sim,j,jj,run,]+wrong_freq_motifs_found[sim,j,jj,run,]
      freq[is.na(freq)]=0
      motifs_in_curves=cbind(rep(1:nmotifs,freq),Reduce(rbind,motifs_search_results$V_occurrences[motifs_found_permut[sim,j,jj,run,]]))
      for(i in 1:N){
        plot(seq(0,curves$fd_curves[[i]]$basis$rangeval[2],by=0.5),eval.fd(seq(0,curves$fd_curves[[i]]$basis$rangeval[2],by=0.5),curves$fd_curves[[i]],Lfdobj=0),type='l',ylim=c(-35,50),xlab='t',ylab='x(t)',main=paste('Random curve',i),cex.main=2,cex.axis=1.5,cex.lab=1.5)
        index=which(motifs_in_curves[,2]==i)
        for(ind in index){
          id_motif=motifs_in_curves[ind,1]
          pos_motif=motifs_in_curves[ind,3]
          len_motif=motifs_search_results$V_length[motifs_found_permut[sim,j,jj,run,id_motif]]
          lines(seq(pos_motif,min(pos_motif-1+len_motif,len),by=0.5),
                eval.fd(seq(pos_motif,min(pos_motif-1+len_motif,len),by=0.5),curves$fd_curves[[i]],Lfdobj=0),
                col=1+id_motif,lwd=3)
        }
      }
      dev.off()
    }
  }
  
  pdf(paste0(folder,sim,'/len',len,'_sd',sd_noise,'_motifs_found.pdf'),10,10)
  boxplot(nmotifs_found[sim,j,jj,],ylim=range(c(nmotifs_found[sim,j,jj,],nmotifs),na.rm=TRUE),main='Number of motifs',cex.main=2,cex.axis=1.5)
  abline(h=nmotifs)
  stripchart(nmotifs_found[sim,j,jj,],vertical=TRUE,method="jitter",pch=19,cex=1.5,add=TRUE,col=1)
  layout(matrix(c(1,1,1,1,1,1,1,1,2,3,4,5),nrow=3,byrow=TRUE))
  for(motif in 1:nmotifs){
    plot(1,type='n',xlab='t',ylab='x(t)',xlim=c(-10,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots+10),ylim=ylim,
         main=paste('Motif',motif),cex.main=3,cex.axis=1.5,cex.lab=1.5)
    breaks=seq(0,(length(curves$coeff_motifs[[motif]])-norder+1+2*(norder-1))*dist_knots,by=dist_knots)-(norder-1)*dist_knots
    basis=create.bspline.basis(norder=norder,breaks=breaks)
    curve=fd(coef=c(rep(0,norder-1),curves$coeff_motifs[[motif]],rep(0,norder-1)),basisobj=basis)
    for(run in 1:10){
      if(paste0('len',len,'_sd',sd_noise,'_run',run,'.RData') %in% list.files(paste0(folder,sim))){
        load(paste0(folder,sim,'/len',len,'_sd',sd_noise,'_run',run,'.RData'))
        if(!is.na(motifs_found_permut[sim,j,jj,run,motif]))
          lines(1:nrow(motifs_search_results$V0[[motifs_found_permut[sim,j,jj,run,motif]]])-1+MV_shift[sim,j,jj,run,motif],motifs_search_results$V0[[motifs_found_permut[sim,j,jj,run,motif]]],col=1+motif,lty=2,lwd=2)
      }
    }
    curve_motif=eval.fd(seq(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=0.5),curve,Lfdobj=0)
    if(transformed){
      curve_motif=(curve_motif-min(curve_motif))/(max(curve_motif)-min(curve_motif))
    }
    lines(seq(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=0.5),curve_motif,col=1,lwd=2)
    lines(seq(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=1),M_average[[motif]]$m0,col=1,lwd=2,lty=3)
    legend('topright',legend=c('Real','Sample mean','Estimated'),col=c(1,1,motif+1),lty=c(1,3,2),lwd=c(2,2,2),cex=1.5)
    boxplot(MV_dist[sim,j,jj,,motif],main='Distance',ylim=c(0,max(c(MV_dist[sim,j,jj,,],MM_average_dist[sim,j,jj,motif]),na.rm=TRUE)),cex.main=2,cex.axis=1.5)
    abline(h=MM_average_dist[sim,j,jj,motif])
    stripchart(MV_dist[sim,j,jj,,motif],vertical=TRUE,method="jitter",pch=19,cex=1.5,add=TRUE,col=1+motif)
    boxplot(len_motifs_found[sim,j,jj,,motif]-1,main='Length',ylim=c(0,max(c(len_motifs_found[sim,j,jj,,]-1,len_motifs),na.rm=TRUE)),cex.main=2,cex.axis=1.5)
    abline(h=len_motifs)
    stripchart(len_motifs_found[sim,j,jj,,motif]-1,vertical=TRUE,method="jitter",pch=19,cex=1.5,add=TRUE,col=1+motif)
    boxplot(freq_motifs_found[sim,j,jj,,motif],ylim=c(0,max(c(freq_motifs_found[sim,j,jj,,],colSums(Reduce(rbind,freq_motifs))),na.rm=TRUE)),main='True positives',cex.main=2,cex.axis=1.5)
    abline(h=colSums(Reduce(rbind,freq_motifs))[motif])
    stripchart(freq_motifs_found[sim,j,jj,,motif],vertical=TRUE,method="jitter",pch=19,cex=1.5,add=TRUE,col=1+motif)
    boxplot(wrong_freq_motifs_found[sim,j,jj,,motif],ylim=c(0,max(c(wrong_freq_motifs_found[sim,j,jj,,],colSums(Reduce(rbind,freq_motifs))),na.rm=TRUE)),main='False positives',cex.main=2,cex.axis=1.5)
    abline(h=0)
    stripchart(wrong_freq_motifs_found[sim,j,jj,,motif],vertical=TRUE,method="jitter",pch=19,cex=1.5,add=TRUE,col=1+motif)
  }
  dev.off()
}

pdf(paste0(folder,sim,'/len',len,'_motifs_found.pdf'),10,10)
boxplot(t(nmotifs_found[sim,j,,]),ylim=range(c(nmotifs_found[sim,j,,],nmotifs),na.rm=TRUE),xlab='Noise sd',cex.lab=1.5,main='Number of motifs',cex.main=2,cex.axis=1.5)
abline(h=nmotifs)
stripchart(as.vector(t(nmotifs_found[sim,j,,]))~rep(sd_simulation,each=10),vertical=TRUE,method="jitter",pch=19,cex=1.5,add=TRUE,col=1)
layout(matrix(c(1,1,1,1,1,1,1,1,2,3,4,5),nrow=3,byrow=TRUE))
for(motif in 1:nmotifs){
  plot(1,type='n',xlab='t',ylab='x(t)',xlim=c(-10,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots+10),ylim=ylim,
       main=paste('Motif',motif),cex.main=3,cex.axis=1.5,cex.lab=1.5)
  breaks=seq(0,(length(curves$coeff_motifs[[motif]])-norder+1+2*(norder-1))*dist_knots,by=dist_knots)-(norder-1)*dist_knots
  basis=create.bspline.basis(norder=norder,breaks=breaks)
  curve=fd(coef=c(rep(0,norder-1),curves$coeff_motifs[[motif]],rep(0,norder-1)),basisobj=basis)
  col=c(rgb(colorRamp(c('white',motif+1))((seq_along(sd_simulation)/length(sd_simulation)))[-1,],alpha=150,max=255),motif+1)
  for(jj in seq_along(sd_simulation)){
    sd_noise=sd_simulation[jj]
    for(run in 1:10){
      if(paste0('len',len,'_sd',sd_noise,'_run',run,'.RData') %in% list.files(paste0(folder,sim))){
        load(paste0(folder,sim,'/len',len,'_sd',sd_noise,'_run',run,'.RData'))
        if(!is.na(motifs_found_permut[sim,j,jj,run,motif]))
          lines(1:nrow(motifs_search_results$V0[[motifs_found_permut[sim,j,jj,run,motif]]])-1+MV_shift[sim,j,jj,run,motif],motifs_search_results$V0[[motifs_found_permut[sim,j,jj,run,motif]]],col=col[jj],type='b',pch=jj,cex=1.5,lty=2,lwd=2)
      }
    }
  }
  curve_plot=eval.fd(seq(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=0.5),curve,Lfdobj=0)
  if(transformed){
    curve_plot=(curve_plot-min(curve_plot))/(max(curve_plot)-min(curve_plot))
  }
  lines(seq(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=0.5),curve_plot,col=1,lwd=2)
  legend('topright',legend=c('Real',paste0('Estimated - sd=',sd_simulation)),
         col=c(1,col),lty=c(1,rep(2,length(sd_simulation))),lwd=c(2,rep(2,length(sd_simulation))),
         pch=c(NA,seq_along(sd_simulation)),cex=1.5)
  boxplot(t(MV_dist[sim,j,,,motif]),main='Distance',ylim=c(0,max(c(MV_dist[sim,j,,,],MM_average_dist[sim,j,,motif]),na.rm=TRUE)),xlab='Noise sd',cex.lab=1.5,cex.main=2,cex.axis=1.5,las=3)
  lines(rep(seq_along(sd_simulation),each=2)+c(-0.5,0.5),rep(MM_average_dist[sim,j,,motif],each=2),type='s')
  stripchart(as.vector(t(MV_dist[sim,j,,,motif]))~rep(sd_simulation,each=10),vertical=TRUE,method="jitter",pch=seq_along(sd_simulation),cex=1.5,add=TRUE,col=col)
  boxplot(t(len_motifs_found[sim,j,,,motif]-1),main='Length',ylim=c(0,max(c(len_motifs_found[sim,j,,,]-1,len_motifs),na.rm=TRUE)),xlab='Noise sd',cex.lab=1.5,cex.main=2,cex.axis=1.5,las=3)
  abline(h=len_motifs)
  stripchart(as.vector(t(len_motifs_found[sim,j,,,motif]-1))~rep(sd_simulation,each=10),vertical=TRUE,method="jitter",pch=seq_along(sd_simulation),cex=1.5,add=TRUE,col=col)
  boxplot(t(freq_motifs_found[sim,j,,,motif]),ylim=c(0,max(c(freq_motifs_found[sim,j,,,],colSums(Reduce(rbind,freq_motifs))),na.rm=TRUE)),main='True positives',xlab='Noise sd',cex.lab=1.5,cex.main=2,cex.axis=1.5,las=3)
  abline(h=colSums(Reduce(rbind,freq_motifs))[motif])
  stripchart(as.vector(t(freq_motifs_found[sim,j,,,motif]))~rep(sd_simulation,each=10),vertical=TRUE,method="jitter",pch=seq_along(sd_simulation),cex=1.5,add=TRUE,col=col)
  boxplot(t(wrong_freq_motifs_found[sim,j,,,motif]),ylim=c(0,max(c(wrong_freq_motifs_found[sim,j,,,],colSums(Reduce(rbind,freq_motifs))),na.rm=TRUE)),main='False positives',xlab='Noise sd',cex.lab=1.5,cex.main=2,cex.axis=1.5,las=3)
  abline(h=0)
  stripchart(as.vector(t(wrong_freq_motifs_found[sim,j,,,motif]))~rep(sd_simulation,each=10),vertical=TRUE,method="jitter",pch=seq_along(sd_simulation),cex=1.5,add=TRUE,col=col)
}
dev.off()


##################
### PLOT PAPER ###
##################

########## Figures about simulation motifs and curves ###########
for(jj in 1:4){
  sd_noise=sd_simulation[jj]
  load(paste0(folder,sim,'/len',len,'_sd',sd_noise,'.RData'))
  colors=c('red','blue')
  
  pdf(paste0(folder,sim,'/ex_motifs_ampl_with_noise',gsub("\\.","",as.character(sd_noise)),'.pdf'),height=4,width=7)
  par(mar=c(2,5,3,1)+0.1)
  for(motif in 1:nmotifs){
    xy=mapply(function(fd_curve,motifs_in_curves_i){
      x=lapply(motifs_in_curves_i$pos_motifs[motifs_in_curves_i$id_motifs==motif]-1,function(start) seq(start*dist_knots,(start+length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=0.5))
      if(is.null(x))
        return(NULL)
      xy=lapply(x,function(x) cbind(x,eval.fd(x,fd_curve,Lfdobj=0)))
      return(xy)
    },curves$fd_curves,curves$motifs_in_curves,SIMPLIFY=FALSE)
    plot(1,type='n',xlab='',ylab='v(t)',xlim=c(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots),ylim=c(-55,68),yaxt='n',
         main=substitute(bold("Motif"~m*","~sigma*"="*s),list(m=as.character(motif),s=as.character(sd_noise))),cex.main=2,cex.axis=1.5,cex.lab=1.7)
    lapply(xy,function(xy){
      if(!is.null(xy))
        lapply(xy,function(xy) lines(xy[,1]-xy[1,1],xy[,2],col=colors[motif],lty=2,lwd=2))
    })
    breaks=seq(0,(length(curves$coeff_motifs[[motif]])-norder+1+2*(norder-1))*dist_knots,by=dist_knots)-(norder-1)*dist_knots
    basis=create.bspline.basis(norder=norder,breaks=breaks)
    curve=fd(coef=c(rep(0,norder-1),curves$coeff_motifs[[motif]],rep(0,norder-1)),basisobj=basis)
    #lines(seq(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=0.5),eval.fd(seq(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=0.5),curve,Lfdobj=0),lwd=3)
    pos=ifelse(motif==1,'topright','topleft')
    axis(2,at=seq(-55,68,5),labels=c(NA,-50,NA,-40,NA,-30,NA,-20,NA,-10,NA,0,NA,10,NA,20,NA,30,NA,40,NA,50,NA,60,NA),cex.axis=1.5,las=2)
    #legend(pos,legend=c('Motif','Occurrences'),col=c('black',colors[motif]),bty='n',lwd=c(3,2),lty=c(1,2),cex=1.5)
  }
  dev.off()
  
  pdf(paste0(folder,sim,'/ex_curves_and_motifs_ampl_with_noise',gsub("\\.","",as.character(sd_noise)),'.pdf'),height=4,width=14)
  par(mar=c(2,5,3,1)+0.1)
  i=1
  plot(seq(0,curves$fd_curves[[i]]$basis$rangeval[2],by=0.5),eval.fd(seq(0,curves$fd_curves[[i]]$basis$rangeval[2],by=0.5),curves$fd_curves[[i]],Lfdobj=0),
       type='l',ylim=c(-55,68),yaxt='n',xlab='',ylab='x(t)',main=substitute(bold("Curves with motifs, "~sigma*"="*s),list(s=as.character(sd_noise))),lty=i,cex.main=2,cex.axis=1.5,cex.lab=1.7)
  mapply(function(id_motif,pos_motif){
    lines(seq((pos_motif-1)*dist_knots,(pos_motif-1+length(curves$coeff_motifs[[id_motif]])-norder+1)*dist_knots,by=0.5),eval.fd(seq((pos_motif-1)*dist_knots,(pos_motif-1+length(curves$coeff_motifs[[id_motif]])-norder+1)*dist_knots,by=0.5),curves$fd_curves[[i]],Lfdobj=0),
          col=colors[id_motif],lty=i,lwd=2)
  },curves$motifs_in_curves[[i]]$id_motifs,curves$motifs_in_curves[[i]]$pos_motifs)
  axis(2,at=seq(-55,68,5),labels=c(NA,-50,NA,-40,NA,-30,NA,-20,NA,-10,NA,0,NA,10,NA,20,NA,30,NA,40,NA,50,NA,60,NA),cex.axis=1.5,las=2)
  for(i in 2:N){
    lines(seq(0,curves$fd_curves[[i]]$basis$rangeval[2],by=0.5),eval.fd(seq(0,curves$fd_curves[[i]]$basis$rangeval[2],by=0.5),curves$fd_curves[[i]],Lfdobj=0),lty=i)
    mapply(function(id_motif,pos_motif){
      lines(seq((pos_motif-1)*dist_knots,(pos_motif-1+length(curves$coeff_motifs[[id_motif]])-norder+1)*dist_knots,by=0.5),eval.fd(seq((pos_motif-1)*dist_knots,(pos_motif-1+length(curves$coeff_motifs[[id_motif]])-norder+1)*dist_knots,by=0.5),curves$fd_curves[[i]],Lfdobj=0),
            col=colors[id_motif],lty=i,lwd=2)
    },curves$motifs_in_curves[[i]]$id_motifs,curves$motifs_in_curves[[i]]$pos_motifs)
  }
  dev.off()
}


########## Figures about simulation results ##########
pdf(paste0(folder,sim,"/ex_results_ampl_",len,".pdf"),height=7,width=11)
par(mar=c(4,5,3,1)+0.1)
layout(matrix(c(1,1,1,1,1,1,2,2,3,3,4,4,5,5,6,6),nrow=2,byrow=TRUE))
colors=c('red','blue')
for(motif in 1:nmotifs){
  plot(1,type='n',xlab='',ylab='v(t)',xlim=c(-10,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots+10),ylim=ylim,
       main=paste('Motif',motif),cex.main=2,cex.axis=1.8,cex.lab=2,las=2)
  breaks=seq(0,(length(curves$coeff_motifs[[motif]])-norder+1+2*(norder-1))*dist_knots,by=dist_knots)-(norder-1)*dist_knots
  basis=create.bspline.basis(norder=norder,breaks=breaks)
  curve=fd(coef=c(rep(0,norder-1),curves$coeff_motifs[[motif]],rep(0,norder-1)),basisobj=basis)
  col=c(rgb(colorRamp(c('white',colors[motif]))((seq_along(sd_simulation)/length(sd_simulation)))[-1,],alpha=150,max=255),colors[motif])
  for(jj in seq_along(sd_simulation)){
    sd_noise=sd_simulation[jj]
    for(run in 1:10){
      if(paste0('len',len,'_sd',sd_noise,'_run',run,'.RData') %in% list.files(paste0(folder,sim))){
        load(paste0(folder,sim,'/len',len,'_sd',sd_noise,'_run',run,'.RData'))
        if(!is.na(motifs_found_permut[sim,j,jj,run,motif]))
          lines(1:nrow(motifs_search_results$V0[[motifs_found_permut[sim,j,jj,run,motif]]])-1+MV_shift[sim,j,jj,run,motif],motifs_search_results$V0[[motifs_found_permut[sim,j,jj,run,motif]]],
                col=col[jj],type='b',pch=jj,cex=1.5,lty=2,lwd=2)
      }
    }
  }
  curve_motif=eval.fd(seq(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=0.5),curve,Lfdobj=0)
  if(transformed){
    curve_motif=(curve_motif-min(curve_motif))/(max(curve_motif)-min(curve_motif))
  }
  lines(seq(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=0.5),curve_motif,col=1,lwd=2)
  plot.new()
  legend('left',legend=c('True motif',
                         as.expression(substitute(paste('Estimated,'~sigma*"="*s),list(s=as.character(sd_simulation[1])))),
                         as.expression(substitute(paste('Estimated,'~sigma*"="*s),list(s=as.character(sd_simulation[2])))),
                         as.expression(substitute(paste('Estimated,'~sigma*"="*s),list(s=as.character(sd_simulation[3])))),
                         as.expression(substitute(paste('Estimated,'~sigma*"="*s),list(s=as.character(sd_simulation[4]))))),
         inset=c(-0.1,0),xpd=TRUE,bty="n",col=c(1,col),lty=c(1,rep(2,length(sd_simulation))),lwd=c(2,rep(2,length(sd_simulation))),
         pch=c(NA,seq_along(sd_simulation)),cex=1.8)
  options(scipen = 0)
  boxplot(t(MV_dist[sim,j,,,motif]),main='Distance',ylim=c(0,max(c(MV_dist[sim,j,,,],MM_average_dist[sim,j,,motif]),na.rm=TRUE)),
          xlab=expression(sigma),cex.lab=1.8,cex.main=2,cex.axis=1.8)
  lines(rep(seq_along(sd_simulation),each=2)+c(-0.5,0.5),rep(MM_average_dist[sim,j,,motif],each=2),type='s')
  stripchart(as.vector(t(MV_dist[sim,j,,,motif]))~rep(sd_simulation,each=10),vertical=TRUE,method="jitter",pch=seq_along(sd_simulation),
             cex=1.5,add=TRUE,col=col)
  boxplot(t(len_motifs_found[sim,j,,,motif]-1),main='Length',ylim=c(0,max(c(len_motifs_found[sim,j,,,]-1,len_motifs),na.rm=TRUE)),
          xlab=expression(sigma),cex.lab=1.8,cex.main=2,cex.axis=1.8,las=2)
  abline(h=len_motifs)
  stripchart(as.vector(t(len_motifs_found[sim,j,,,motif]-1))~rep(sd_simulation,each=10),vertical=TRUE,method="jitter",pch=seq_along(sd_simulation),cex=1.5,add=TRUE,col=col)
  boxplot(t(freq_motifs_found[sim,j,,,motif]),ylim=c(0,max(c(freq_motifs_found[sim,j,,,],colSums(Reduce(rbind,freq_motifs))),na.rm=TRUE)),
          main='True positives',xlab=expression(sigma),cex.lab=1.8,cex.main=2,cex.axis=1.8,las=2)
  abline(h=colSums(Reduce(rbind,freq_motifs))[motif])
  stripchart(as.vector(t(freq_motifs_found[sim,j,,,motif]))~rep(sd_simulation,each=10),vertical=TRUE,method="jitter",pch=seq_along(sd_simulation),cex=1.5,add=TRUE,col=col)
  boxplot(t(wrong_freq_motifs_found[sim,j,,,motif]),ylim=c(0,max(c(wrong_freq_motifs_found[sim,j,,,],colSums(Reduce(rbind,freq_motifs))),na.rm=TRUE)),
          main='False positives',xlab=expression(sigma),cex.lab=1.8,cex.main=2,cex.axis=1.8,las=2)
  abline(h=0)
  stripchart(as.vector(t(wrong_freq_motifs_found[sim,j,,,motif]))~rep(sd_simulation,each=10),vertical=TRUE,method="jitter",pch=seq_along(sd_simulation),cex=1.5,add=TRUE,col=col)
}
dev.off()






