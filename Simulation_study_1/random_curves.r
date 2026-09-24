resample <- function(x, ...) x[sample.int(length(x), ...)]

generate_curves <- function(N,dist_knots=10,len=300,norder=4,
                            coeff_motifs=NULL,nmotifs=2,len_motifs=30,
                            min_dist_motifs=norder*dist_knots,freq_motifs=5,
                            motifs_in_curves=NULL,
                            coeff_min=-15,coeff_max=15,distrib='unif',
                            only_der=FALSE,coeff_min_only_der=-10,coeff_max_only_der=10,
                            amplitude=FALSE,coeff_min_amplitude=0.5,coeff_max_amplitude=1.5,
                            sd_noise=1){
  # Function that generates N (univariate) curves with motifs, using b-spline basis. 
  # If the motifs are not provided by the user they are randomly generated, i.e. the b-spline coefficients that define them 
  # are randomly sampled from [coeff_min,coeff_max] with uniform distribution. 
  # Each of the N curves is then constructed by populating it with the motifs (adding gaussian noise to the coefficients that 
  # define them) and randomly sampling the non-motif coefficients from [coeff_min,coeff_max] with uniform/beta distribution.
  # 
  # Returns list with four elements:
  # curves: a list of N fd objects with the curves.
  # motifs_in_curves: a list with the motif ids and starting position (starting coefficient).
  # coeff_motifs: a list with the coefficients defining the motifs.
  # freq_motifs: a list with the frequencies for the different motifs in the different curves.
  #
  # N: number of curves. 
  #
  # dist_knots: distance between b-spline knots. 
  # len: an integer specifying the curve length, or a vector of N integers with the lengths of the different curves. 
  #      Each length must be a multiple of dist_knots. 
  # norder: the order of b-splines, which is one higher than their degree. 
  #
  # coeff_motifs: list of vectors, each specifying the b-spline coefficients that define a motif (every motif corresponds 
  #               to at least norder coefficients).
  # nmotifs: number of different motifs to be inserted in the curves. If coeff_motifs is not NULL, nmotifs=length(coeff_motifs).
  # len_motifs: an integer specifying the length of the motifs, or a vector of nmotifs integers with the lengths of the 
  #             different motifs. Each length must be a multiple of dist_knots. If coeff_motifs is not NULL, it is derived from it.
  #
  # min_dist_motifs: minimum distance between two motifs in the same curve. Must be at least norder*dist_knots, and a multiple 
  #                  of dist_knots. 
  # freq_motifs: an integer specifying the total frequency of the motifs in the N curves, or a vector of nmotifs integers with 
  #              the frequencies for the different motifs, or a list of N vectors of nmotifs integers with the frequencies for 
  #              the different motifs in the different curves.
  #
  # motifs_in_curves: list of N list, each with id_motifs and pos_motifs vectors giving the ids and positions of the motifs in the curve.
  #                   NOTE: no check is done on the ammissibility of it (provide it only together with coeff_motifs and freq_motifs
  #                   returned from a previous run of the function).
  #
  # coeff_min, coeff_max: minumim and maximum values for the b-spline coefficients.
  # distrib: 'unif' for generating the coefficients from a uniform distribution between coeff_min and coeff_max, 
  #          'beta' for generating them from a Beta(0.45,0.45) rescaled to be between coeff_min and coeff_max.
  #
  # only_der: if TRUE, motifs in curves have the same derivative, but they are not the same (random constant b-spline coefficients added).
  # coeff_min_only_der, coeff_max_only_der: minimum and maximum values to be added to b-spline coefficients when only_der=TRUE.
  #
  # sd_noise: an integer specifying the stardard deviation of the gaussian noise to be added to the coefficients defining 
  #           the motifs, or a vector of nmotifs integers with the standard deviation of the gaussian noise for the different motifs.
  
  require(fda)
  
  #### check N, dist_knots and len ####
  if( (N%%1 != 0) | (N < 1) )
    stop('Invalid \'N\'.')
  if( (dist_knots%%1 != 0) | (dist_knots < 1))
    stop('Invalid \'dist_knots\'.')
  if( TRUE %in% ((len%%1 != 0) | (len < 1) | (sum(len%%dist_knots) > 0)) )
    stop('Invalid \'len\'.')
  
  #### create knots and generate the corresponding b-spline basis for every curve ####
  breaks = lapply(rep(len, length.out = N), seq, from = 0, by = dist_knots)
  basis = lapply(breaks, function(breaks_i) create.bspline.basis(norder = norder, breaks = breaks_i))
  
  #### if motifs are provided, compute nmotifs and len_motifs ####
  if( !is.null(coeff_motifs) ){
    nmotifs=length(coeff_motifs)
    len_motifs=(unlist(lapply(coeff_motifs,length))-norder+1)*dist_knots
  }
  #### check nmotifs and len_motifs ####
  if( (nmotifs%%1 != 0) | (nmotifs < 1) )
    stop('Invalid \'nmotifs\'.')
  if( TRUE %in% ((len_motifs%%1 != 0) | (len_motifs < 1) | (sum(len_motifs%%dist_knots) > 0) | (TRUE %in% (max(len_motifs) > len))) )
    stop('Invalid \'len_motifs\'.')
  #### if motifs are not provided, randomly generate them ####
  if( is.null(coeff_motifs) ){
    if( distrib == 'unif' ){
      coeff_motifs = lapply(rep(len_motifs/dist_knots+norder-1, length.out = nmotifs), function(l) runif(l, min = coeff_min, max = coeff_max))
    }else if( distrib == 'beta' ){
      coeff_motifs = lapply(rep(len_motifs/dist_knots+norder-1, length.out = nmotifs), function(l) coeff_min+rbeta(l, 0.45, 0.45)*(coeff_max-coeff_min))
    }else{
      stop('Wrong \'distrib\'')
    }
  }
  
  #### check min_dist_motifs ####
  if( (min_dist_motifs < norder*dist_knots) | (min_dist_motifs%%dist_knots > 0) )
    stop('Invalid \'min_dist_motifs\'.')
  #### check freq_motifs and convert it into a list with the frequencies for the different motifs in the different curves ####
  if( !is.list(freq_motifs) ){
    if( TRUE %in% ((freq_motifs%%1 != 0) | (freq_motifs < 1) | 
                   (sum(norder-1+rep(len/dist_knots, length.out = N)-2*(norder-1)) < 
                    (sum(rep(len_motifs/dist_knots+norder-1, length.out = nmotifs)*rep(freq_motifs, length.out = nmotifs))
                     +max(0, sum(rep(freq_motifs, length.out = nmotifs))-N)*(min_dist_motifs/dist_knots-norder+1)))) )
      stop('Invalid \'freq_motifs\'.')
    freq_check = 1
    it = 0
    while( (sum(freq_check) > 0) & (it < 10000) ){
      it = it+1
      freq_motifs_new = matrix(unlist(lapply(rep(freq_motifs, length.out = nmotifs),
                                             function(freq_motif) 
                                               as.vector(table(factor(sample(N, freq_motif, replace = TRUE), levels = 1:N))))), nrow=N
                               )
      freq_motifs_new = split(freq_motifs_new, rep(1:N, nmotifs))
      freq_check = mapply(function(freq_motifs_i, len_i) 
                            ((len_i-2*(norder-1)) < 
                               (sum(rep(len_motifs/dist_knots+norder-1, length.out = nmotifs)*rep(freq_motifs_i, length.out = nmotifs))
                                +(sum(freq_motifs_i)-1)*(min_dist_motifs/dist_knots-norder+1))),
                          freq_motifs_new, norder-1+rep(len/dist_knots, length.out = N))
    }
    if( it == 10000 )
      stop('Tried 10.000 random configurations of motif frequencies in the different curves, unable to find a valid configuration. Please select lower \'freq_motifs\'.')
    freq_motifs = freq_motifs_new
  }else{
    freq_motifs = rep(freq_motifs, length.out = N)
    freq_check = mapply(function(freq_motifs_i, len_i) TRUE %in% ((freq_motifs_i%%1 != 0) | (freq_motifs_i < 0) | 
                                                                    ((len_i-2*(norder-1)) < (sum(rep(len_motifs/dist_knots+norder-1, length.out = nmotifs)*rep(freq_motifs_i, length.out = nmotifs))+(sum(freq_motifs_i)-1)*(min_dist_motifs/dist_knots-norder+1)))),
                      freq_motifs, norder-1+rep(len/dist_knots, length.out = N))
    if( sum(freq_check) )
      stop('Invalid \'freq_motifs\'.')
  }
  # randomly assign motif positions in curves
  if( is.null(motifs_in_curves) ){
    motifs_in_curves = mapply(function(freq_motifs_i, len_i){
                                if( sum(freq_motifs_i) == 0 )
                                  return(NULL)
                                id_motifs = resample(rep(seq_along(freq_motifs_i), freq_motifs_i))
                                len_elements = c(norder-1, 
                                                 rep(len_motifs/dist_knots+norder-1, length.out = nmotifs)[id_motifs]
                                                 +c(rep((min_dist_motifs/dist_knots-norder+1), sum(freq_motifs_i)-1),0), 
                                                 norder-1)
                                gaps_tot = len_i-sum(len_elements)
                                gaps = diff(c(0, sort(sample(gaps_tot+sum(freq_motifs_i), sum(freq_motifs_i)))))-1
                                pos_motifs = cumsum(len_elements[seq_along(id_motifs)])+1+cumsum(gaps)
                                return(list(id_motifs = id_motifs, pos_motifs = pos_motifs))
                              },freq_motifs, norder-1+rep(len/dist_knots, length.out=N), SIMPLIFY = FALSE)
  }
  # generate curves
  fd_curves = mapply(function(motifs_in_curves_i, len_i, basis_i){
                     coeff=rep(NA,len_i)
                     if( distrib == 'unif' ){
                       if( is.null(motifs_in_curves_i) ){
                         coeff = runif(len_i, min = coeff_min, max = coeff_max)
                       }else{
                         pos_coeff_motifs = unlist(mapply(function(a,b) seq(a)+b, 
                                                          rep(len_motifs, length.out = nmotifs)[motifs_in_curves_i$id_motifs]/dist_knots+norder-1, motifs_in_curves_i$pos_motif-1, SIMPLIFY = FALSE))
                         if(only_der){
                           if(amplitude){
                             coeff[pos_coeff_motifs] = unlist(coeff_motifs[motifs_in_curves_i$id_motifs])*rep(runif(length(motifs_in_curves_i$id_motifs), min = coeff_min_amplitude, max = coeff_max_amplitude), rep(len_motifs, length.out = nmotifs)[motifs_in_curves_i$id_motifs]/dist_knots+norder-1)+
                               rep(runif(length(motifs_in_curves_i$id_motifs), min = coeff_min_only_der, max = coeff_max_only_der), rep(len_motifs, length.out = nmotifs)[motifs_in_curves_i$id_motifs]/dist_knots+norder-1)+
                               rnorm(length(pos_coeff_motifs), sd = rep(rep(sd_noise, length.out = nmotifs)[motifs_in_curves_i$id_motifs], rep(len_motifs, length.out = nmotifs)[motifs_in_curves_i$id_motifs]/dist_knots+norder-1))
                           }else{
                             coeff[pos_coeff_motifs] = unlist(coeff_motifs[motifs_in_curves_i$id_motifs])+
                             rep(runif(length(motifs_in_curves_i$id_motifs), min = coeff_min_only_der, max = coeff_max_only_der), rep(len_motifs, length.out = nmotifs)[motifs_in_curves_i$id_motifs]/dist_knots+norder-1)+
                             rnorm(length(pos_coeff_motifs), sd = rep(rep(sd_noise, length.out = nmotifs)[motifs_in_curves_i$id_motifs], rep(len_motifs, length.out = nmotifs)[motifs_in_curves_i$id_motifs]/dist_knots+norder-1))
                           }
                         }else{
                           if(amplitude){
                             coeff[pos_coeff_motifs] = unlist(coeff_motifs[motifs_in_curves_i$id_motifs])*rep(runif(length(motifs_in_curves_i$id_motifs), min = coeff_min_amplitude, max = coeff_max_amplitude), rep(len_motifs, length.out = nmotifs)[motifs_in_curves_i$id_motifs]/dist_knots+norder-1)+
                               rnorm(length(pos_coeff_motifs), sd = rep(rep(sd_noise, length.out = nmotifs)[motifs_in_curves_i$id_motifs], rep(len_motifs, length.out = nmotifs)[motifs_in_curves_i$id_motifs]/dist_knots+norder-1))
                           }else{
                             coeff[pos_coeff_motifs] = unlist(coeff_motifs[motifs_in_curves_i$id_motifs])+
                               rnorm(length(pos_coeff_motifs), sd = rep(rep(sd_noise, length.out = nmotifs)[motifs_in_curves_i$id_motifs], rep(len_motifs, length.out = nmotifs)[motifs_in_curves_i$id_motifs]/dist_knots+norder-1))
                             
                           }                         
                         }
                         coeff[-pos_coeff_motifs] = runif(len_i-length(pos_coeff_motifs), min = coeff_min, max = coeff_max)
                       }
                     }else if( distrib == 'beta' ){
                       if( is.null(motifs_in_curves_i) ){
                         coeff = coeff_min+rbeta(len_i, 0.45, 0.45)*(coeff_max-coeff_min)
                       }else{
                         pos_coeff_motifs = unlist(mapply(function(a,b) seq(a)+b,
                                                          rep(len_motifs, length.out = nmotifs)[motifs_in_curves_i$id_motifs]/dist_knots+norder-1, motifs_in_curves_i$pos_motif-1, SIMPLIFY = FALSE))
                         if(only_der){
                           if(amplitude){
                             coeff[pos_coeff_motifs] = unlist(coeff_motifs[motifs_in_curves_i$id_motifs])*rep(runif(length(motifs_in_curves_i$id_motifs), min = coeff_min_amplitude, max = coeff_max_amplitude), rep(len_motifs, length.out = nmotifs)[motifs_in_curves_i$id_motifs]/dist_knots+norder-1)+
                               rep(runif(length(motifs_in_curves_i$id_motifs), min = coeff_min_only_der, max = coeff_max_only_der), rep(len_motifs, length.out = nmotifs)[motifs_in_curves_i$id_motifs]/dist_knots+norder-1)+
                               rnorm(length(pos_coeff_motifs), sd = rep(rep(sd_noise, length.out = nmotifs)[motifs_in_curves_i$id_motifs], rep(len_motifs, length.out = nmotifs)[motifs_in_curves_i$id_motifs]/dist_knots+norder-1))
                           }else{
                             coeff[pos_coeff_motifs] = unlist(coeff_motifs[motifs_in_curves_i$id_motifs])+
                               rep(runif(length(motifs_in_curves_i$id_motifs), min = coeff_min_only_der, max = coeff_max_only_der), rep(len_motifs, length.out = nmotifs)[motifs_in_curves_i$id_motifs]/dist_knots+norder-1)+
                               rnorm(length(pos_coeff_motifs), sd = rep(rep(sd_noise, length.out = nmotifs)[motifs_in_curves_i$id_motifs], rep(len_motifs, length.out = nmotifs)[motifs_in_curves_i$id_motifs]/dist_knots+norder-1))
                           }
                         }else{
                           if(amplitude){
                             coeff[pos_coeff_motifs] = unlist(coeff_motifs[motifs_in_curves_i$id_motifs])*rep(runif(length(motifs_in_curves_i$id_motifs), min = coeff_min_amplitude, max = coeff_max_amplitude), rep(len_motifs, length.out = nmotifs)[motifs_in_curves_i$id_motifs]/dist_knots+norder-1)+
                               rnorm(length(pos_coeff_motifs), sd = rep(rep(sd_noise, length.out = nmotifs)[motifs_in_curves_i$id_motifs], rep(len_motifs, length.out = nmotifs)[motifs_in_curves_i$id_motifs]/dist_knots+norder-1))
                           }else{
                              coeff[pos_coeff_motifs] = unlist(coeff_motifs[motifs_in_curves_i$id_motifs])+
                                rnorm(length(pos_coeff_motifs), sd = rep(rep(sd_noise, length.out = nmotifs)[motifs_in_curves_i$id_motifs], rep(len_motifs, length.out = nmotifs)[motifs_in_curves_i$id_motifs]/dist_knots+norder-1))
                           }
                         }
                         coeff[-pos_coeff_motifs] = coeff_min+rbeta(len_i-length(pos_coeff_motifs), 0.45, 0.45)*(coeff_max-coeff_min)
                       }
                     }else{
                       stop('Wrong \'distrib\'')
                     }
                     
                     curve = fd(coef = coeff, basisobj = basis_i)
                     return(curve)
                   }, motifs_in_curves, norder-1+rep(len/dist_knots, length.out = N), basis, SIMPLIFY = FALSE)
  return(list(fd_curves = fd_curves, motifs_in_curves = motifs_in_curves, coeff_motifs = coeff_motifs, freq_motifs = freq_motifs))
}

evaluate_curves <- function(fd_curves,by=1,Lfdobj=0){
  # Function that evaluates the curves generated with generate_curves. 
  #
  # Returns a list of N vectors.
  #
  # fd_curves: list of functional data objects with the curves to be evaluated.
  # by: distance of (equispaced) grid points at which the functional data objects is to be evaluated. 
  #     It must be a divisor of dist_knots. 
  # Lfdobj: nonnegative integer indicating the degree of derivative for evaluating the curves.
  
  require(fda)
  
  # check by
  if(fd_curves[[1]]$basis$params[1]%%by)
    stop('Invalid \'by\'.')
  # evaluate curves
  Y=lapply(fd_curves,
           function(fd_curve){
             x=seq(0,fd_curve$basis$rangeval[2],by=by)
             y=eval.fd(x,fd_curve,Lfdobj=Lfdobj)
             if(ncol(y)==1)
               return(as.vector(y))
             return(y)
           })
  return(Y)
}
