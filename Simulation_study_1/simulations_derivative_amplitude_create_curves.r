require(fda)
require(combinat)
require(Hmisc)

# Set directory to the folder with the current file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("random_curves.r")
source("../fmd_functions.R")
#source("find_candidate_motifs_same_init.r")


folder='N20_K2-10-10_derivative_amplitude_new_probKMA/'
folder_sobolev='N20_K2-10-10_same_init_and_background/'
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
### GENERATE CURVES ###
#######################

set.seed(2026)
sim=sim_simulation[4]
  dir.create(paste0(folder,sim))
  load(paste0(folder_sobolev,sim,'/base_curves_no_noise.RData'))
  curves_no_noise=generate_curves(N=N,coeff_motifs=curves_no_noise$coeff_motifs,motifs_in_curves=curves_no_noise$motifs_in_curves,
                         dist_knots=dist_knots,len=200,norder=norder,
                         sd_noise=0,distrib='beta',
                         only_der=TRUE,coeff_min_only_der=-10,coeff_max_only_der=10,
                         amplitude=TRUE,coeff_min_amplitude=0.1,coeff_max_amplitude=5)
  save(N,norder,dist_knots,nmotifs,len_motifs,freq_motifs,curves_no_noise,file=paste0(folder,sim,'/base_curves_no_noise.RData'))
  
  # different levels of noise
  for(jj in seq_along(sd_simulation)){
    sd_noise=sd_simulation[jj]
    load(paste0(folder_sobolev,sim,'/base_curves_sd',sd_noise,'.RData'))
    curves_noise_old=curves_noise
    curves_noise=curves_no_noise
    curves_noise$fd_curves=mapply(function(fd_curves_i,motifs_in_curves_i,fd_curves_old_i){
                                    len_motifs=unlist(lapply(curves_noise$coeff_motifs,length))[motifs_in_curves_i$id_motifs]
                                    pos_coeff_motifs=unlist(mapply(function(a,b) seq(a)+b,len_motifs,motifs_in_curves_i$pos_motifs-1,SIMPLIFY=FALSE))
                                    fd_curves_i$coefs[pos_coeff_motifs,]=fd_curves_i$coefs[pos_coeff_motifs,]+fd_curves_old_i$coefs[pos_coeff_motifs,]-unlist(curves_noise$coeff_motifs[motifs_in_curves_i$id_motifs])
                                    return(fd_curves_i)
                                  },curves_noise$fd_curves,curves_noise$motifs_in_curves,curves_noise_old$fd_curves,SIMPLIFY=FALSE)
    save(N,norder,dist_knots,nmotifs,len_motifs,freq_motifs,sd_noise,curves_noise,file=paste0(folder,sim,'/base_curves_sd',sd_noise,'.RData'))
  }
  
  # different curve lengths and levels of noise
  for(j in seq_along(len_simulation)){
    len=len_simulation[j]
    curves_len=generate_curves(N=N,dist_knots=dist_knots,len=len,norder=norder,motifs_in_curves=curves_no_noise$motifs_in_curves,
                               coeff_motifs=curves_no_noise$coeff_motifs,nmotifs=nmotifs,len_motifs=len_motifs,freq_motifs=freq_motifs,
                               sd_noise=0,distrib='beta')
    # different levels of noise
    for(jj in seq_along(sd_simulation)){
      sd_noise=sd_simulation[jj]
      curves=curves_len
      load(paste0(folder,sim,'/base_curves_sd',sd_noise,'.RData'))
      curves$fd_curves=mapply(function(fd_curves_i,fd_curves_noise_i){
                                fd_curves_i$coefs[1:length(fd_curves_noise_i$coefs)]=fd_curves_noise_i$coefs
                                return(fd_curves_i)
                              },curves$fd_curves,curves_noise$fd_curves,SIMPLIFY=FALSE)
      save(N,norder,dist_knots,len,nmotifs,len_motifs,freq_motifs,sd_noise,curves,file=paste0(folder,sim,'/len',len,'_sd',sd_noise,'.RData'))
      
      # plot curves with motifs, motifs and curves
      pdf(file=paste0(folder,sim,'/len',len,'_sd',sd_noise,'.pdf'),10,5*N)
      par(mfrow=c(N,1))
      for(i in 1:N){
        plot(seq(0,curves$fd_curves[[i]]$basis$rangeval[2],by=0.5),eval.fd(seq(0,curves$fd_curves[[i]]$basis$rangeval[2],by=0.5),curves$fd_curves[[i]],Lfdobj=0),type='l',ylim=c(-25,35),xlab='t',ylab='x(t)',main=paste('Random curve',i),cex.main=2,cex.axis=1.5,cex.lab=1.5)
        mapply(function(id_motif,pos_motif){
                 lines(seq((pos_motif-1)*dist_knots,(pos_motif-1+length(curves$coeff_motifs[[id_motif]])-norder+1)*dist_knots,by=0.5),eval.fd(seq((pos_motif-1)*dist_knots,(pos_motif-1+length(curves$coeff_motifs[[id_motif]])-norder+1)*dist_knots,by=0.5),curves$fd_curves[[i]],Lfdobj=0),col=1+id_motif,lwd=3)
               },curves$motifs_in_curves[[i]]$id_motifs,curves$motifs_in_curves[[i]]$pos_motifs)
      }
      dev.off()
      pdf(file=paste0(folder,sim,'/len',len,'_sd',sd_noise,'_motifs.pdf'),10,10)
      for(motif in 1:nmotifs){
        xy=mapply(function(fd_curve,motifs_in_curves_i){
                    x=lapply(motifs_in_curves_i$pos_motifs[motifs_in_curves_i$id_motifs==motif]-1,function(start) seq(start*dist_knots,(start+length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=0.5))
                    if(is.null(x))
                      return(NULL)
                    xy=lapply(x,function(x) cbind(x,eval.fd(x,fd_curve,Lfdobj=0)))
                    return(xy)
                  },curves$fd_curves,curves$motifs_in_curves,SIMPLIFY=FALSE)
        plot(1,type='n',xlab='t',ylab='x(t)',xlim=c(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots),ylim=c(-25,35),
             main=paste('Motif',motif),cex.main=2,cex.axis=1.5,cex.lab=1.5)
        lapply(xy,function(xy){
                    if(!is.null(xy))
                      lapply(xy,function(xy) lines(xy[,1]-xy[1,1],xy[,2]))
                  })
        breaks=seq(0,(length(curves$coeff_motifs[[motif]])-norder+1+2*(norder-1))*dist_knots,by=dist_knots)-(norder-1)*dist_knots
        basis=create.bspline.basis(norder=norder,breaks=breaks)
        curve=fd(coef=c(rep(0,norder-1),curves$coeff_motifs[[motif]],rep(0,norder-1)),basisobj=basis)
        lines(seq(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=0.5),eval.fd(seq(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=0.5),curve,Lfdobj=0),col=1+motif,lwd=3)
        curve_average=colMeans(Reduce(rbind,lapply(xy,
                                                   function(xy){
                                                     if(!is.null(xy))
                                                       Reduce(rbind,lapply(xy,function(xy) xy[,2]))
                                                   })))
        lines(seq(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=0.5),curve_average,col=1+motif,lwd=3,lty=2)
        legend('topleft',legend=c('Simulated','Sample mean','Real'),col=c(1,1+motif,1+motif),lwd=c(1,3,3),lty=c(1,2,1))
      }
      dev.off()
      pdf(file=paste0(folder,sim,'/len',len,'_sd',sd_noise,'_curves.pdf'),10,5*N)
      par(mfrow=c(N,1))
      for(i in 1:N){
        plot(seq(0,curves$fd_curves[[i]]$basis$rangeval[2],by=0.5),eval.fd(seq(0,curves$fd_curves[[i]]$basis$rangeval[2],by=0.5),curves$fd_curves[[i]],Lfdobj=0),type='l',ylim=c(-25,35),xlab='t',ylab='x(t)',main=paste('Random curve',i),cex.main=2,cex.axis=1.5,cex.lab=1.5)
      }
      dev.off()
    }
  }

  
  
  

  sim=sim_simulation[10]
  dir.create(paste0(folder,sim))
  load(paste0(folder_sobolev,sim,'/base_curves_no_noise.RData'))
  curves_no_noise=generate_curves(N=N,coeff_motifs=curves_no_noise$coeff_motifs,motifs_in_curves=curves_no_noise$motifs_in_curves,
                                  dist_knots=dist_knots,len=200,norder=norder,
                                  sd_noise=0,distrib='beta',
                                  only_der=TRUE,coeff_min_only_der=-10,coeff_max_only_der=10,
                                  amplitude=TRUE,coeff_min_amplitude=0.1,coeff_max_amplitude=5)
  save(N,norder,dist_knots,nmotifs,len_motifs,freq_motifs,curves_no_noise,file=paste0(folder,sim,'/base_curves_no_noise.RData'))
  
  # different levels of noise
  for(jj in seq_along(sd_simulation)){
    sd_noise=sd_simulation[jj]
    load(paste0(folder_sobolev,sim,'/base_curves_sd',sd_noise,'.RData'))
    curves_noise_old=curves_noise
    curves_noise=curves_no_noise
    curves_noise$fd_curves=mapply(function(fd_curves_i,motifs_in_curves_i,fd_curves_old_i){
      len_motifs=unlist(lapply(curves_noise$coeff_motifs,length))[motifs_in_curves_i$id_motifs]
      pos_coeff_motifs=unlist(mapply(function(a,b) seq(a)+b,len_motifs,motifs_in_curves_i$pos_motifs-1,SIMPLIFY=FALSE))
      fd_curves_i$coefs[pos_coeff_motifs,]=fd_curves_i$coefs[pos_coeff_motifs,]+fd_curves_old_i$coefs[pos_coeff_motifs,]-unlist(curves_noise$coeff_motifs[motifs_in_curves_i$id_motifs])
      return(fd_curves_i)
    },curves_noise$fd_curves,curves_noise$motifs_in_curves,curves_noise_old$fd_curves,SIMPLIFY=FALSE)
    save(N,norder,dist_knots,nmotifs,len_motifs,freq_motifs,sd_noise,curves_noise,file=paste0(folder,sim,'/base_curves_sd',sd_noise,'.RData'))
  }
  
  # different curve lengths and levels of noise
  for(j in seq_along(len_simulation)){
    len=len_simulation[j]
    curves_len=generate_curves(N=N,dist_knots=dist_knots,len=len,norder=norder,motifs_in_curves=curves_no_noise$motifs_in_curves,
                               coeff_motifs=curves_no_noise$coeff_motifs,nmotifs=nmotifs,len_motifs=len_motifs,freq_motifs=freq_motifs,
                               sd_noise=0,distrib='beta')
    # different levels of noise
    for(jj in seq_along(sd_simulation)){
      sd_noise=sd_simulation[jj]
      curves=curves_len
      load(paste0(folder,sim,'/base_curves_sd',sd_noise,'.RData'))
      curves$fd_curves=mapply(function(fd_curves_i,fd_curves_noise_i){
        fd_curves_i$coefs[1:length(fd_curves_noise_i$coefs)]=fd_curves_noise_i$coefs
        return(fd_curves_i)
      },curves$fd_curves,curves_noise$fd_curves,SIMPLIFY=FALSE)
      save(N,norder,dist_knots,len,nmotifs,len_motifs,freq_motifs,sd_noise,curves,file=paste0(folder,sim,'/len',len,'_sd',sd_noise,'.RData'))
      
      # plot curves with motifs, motifs and curves
      pdf(file=paste0(folder,sim,'/len',len,'_sd',sd_noise,'.pdf'),10,5*N)
      par(mfrow=c(N,1))
      for(i in 1:N){
        plot(seq(0,curves$fd_curves[[i]]$basis$rangeval[2],by=0.5),eval.fd(seq(0,curves$fd_curves[[i]]$basis$rangeval[2],by=0.5),curves$fd_curves[[i]],Lfdobj=0),type='l',ylim=c(-25,35),xlab='t',ylab='x(t)',main=paste('Random curve',i),cex.main=2,cex.axis=1.5,cex.lab=1.5)
        mapply(function(id_motif,pos_motif){
          lines(seq((pos_motif-1)*dist_knots,(pos_motif-1+length(curves$coeff_motifs[[id_motif]])-norder+1)*dist_knots,by=0.5),eval.fd(seq((pos_motif-1)*dist_knots,(pos_motif-1+length(curves$coeff_motifs[[id_motif]])-norder+1)*dist_knots,by=0.5),curves$fd_curves[[i]],Lfdobj=0),col=1+id_motif,lwd=3)
        },curves$motifs_in_curves[[i]]$id_motifs,curves$motifs_in_curves[[i]]$pos_motifs)
      }
      dev.off()
      pdf(file=paste0(folder,sim,'/len',len,'_sd',sd_noise,'_motifs.pdf'),10,10)
      for(motif in 1:nmotifs){
        xy=mapply(function(fd_curve,motifs_in_curves_i){
          x=lapply(motifs_in_curves_i$pos_motifs[motifs_in_curves_i$id_motifs==motif]-1,function(start) seq(start*dist_knots,(start+length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=0.5))
          if(is.null(x))
            return(NULL)
          xy=lapply(x,function(x) cbind(x,eval.fd(x,fd_curve,Lfdobj=0)))
          return(xy)
        },curves$fd_curves,curves$motifs_in_curves,SIMPLIFY=FALSE)
        plot(1,type='n',xlab='t',ylab='x(t)',xlim=c(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots),ylim=c(-25,35),
             main=paste('Motif',motif),cex.main=2,cex.axis=1.5,cex.lab=1.5)
        lapply(xy,function(xy){
          if(!is.null(xy))
            lapply(xy,function(xy) lines(xy[,1]-xy[1,1],xy[,2]))
        })
        breaks=seq(0,(length(curves$coeff_motifs[[motif]])-norder+1+2*(norder-1))*dist_knots,by=dist_knots)-(norder-1)*dist_knots
        basis=create.bspline.basis(norder=norder,breaks=breaks)
        curve=fd(coef=c(rep(0,norder-1),curves$coeff_motifs[[motif]],rep(0,norder-1)),basisobj=basis)
        lines(seq(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=0.5),eval.fd(seq(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=0.5),curve,Lfdobj=0),col=1+motif,lwd=3)
        curve_average=colMeans(Reduce(rbind,lapply(xy,
                                                   function(xy){
                                                     if(!is.null(xy))
                                                       Reduce(rbind,lapply(xy,function(xy) xy[,2]))
                                                   })))
        lines(seq(0,(length(curves$coeff_motifs[[motif]])-norder+1)*dist_knots,by=0.5),curve_average,col=1+motif,lwd=3,lty=2)
        legend('topleft',legend=c('Simulated','Sample mean','Real'),col=c(1,1+motif,1+motif),lwd=c(1,3,3),lty=c(1,2,1))
      }
      dev.off()
      pdf(file=paste0(folder,sim,'/len',len,'_sd',sd_noise,'_curves.pdf'),10,5*N)
      par(mfrow=c(N,1))
      for(i in 1:N){
        plot(seq(0,curves$fd_curves[[i]]$basis$rangeval[2],by=0.5),eval.fd(seq(0,curves$fd_curves[[i]]$basis$rangeval[2],by=0.5),curves$fd_curves[[i]],Lfdobj=0),type='l',ylim=c(-25,35),xlab='t',ylab='x(t)',main=paste('Random curve',i),cex.main=2,cex.axis=1.5,cex.lab=1.5)
      }
      dev.off()
    }
  }
  