working.variog <- function(obs,nb,vth,vlen,llim,verbose,Large,compute){
  ## obs   : coordinates and time stamps and data of observed points
  ## nb    : neighborhood list.
  ## compute: whether to compute variogram or just determine the flags
  ## additional arguments documented in ltsk
  ## value : fitted product sum variogram, and a flag denoting any other conditions
  ##       : flag documented in ltsk
  
  ## check data availability
  flag <- rep(0,length(nb))
  value <- vector("list",length(nb))
  for(i in 1:length(value)){
    #q0 <- query[i,]
    ii <- nb[[i]]
    if( length(ii)<=5 ){
      ## Not enough neighbors
      if(verbose) cat('k= ',length(ii),'\n')
      flag[i] <- 3
    }else{
      ## alternative check based upon discussion with Jin Aug 04
      ssout <- dsubsample2(obs[ii,],Large=Large)
      ## update neighbor based on sub-sampling results
      if(!is.null(ssout$sub)){
        nb[[i]] <- nb[[i]][ssout$sub]
      }
      nbr <- ssout$nbr
      if(verbose)
      {
        with(ssout,cat('k= ',nrow(nbr),'ns=',ns,'nt=',nt,'\n'))
      }
      if( (ssout$ns > llim[1]) && (ssout$nt > llim[2]) )
      {
        if(compute){
          ## compute variogram
          vout <- dvariogram(nbr,vth,vlen)
          vout <- dsmooth.variogram(vout)
          fout <- dfitvariogram(vout,nbr)
          if(fout$ret){
            value[[i]] <- fout
          }else{
            flag[i] <- 4
          }
        }else{
          flag[i] <- NA
        }
      }else if(ssout$nt <= llim[2]){
        if (verbose) cat('insufficient time points.\n')
        flag[i]<- 1
      }else if(ssout$ns <= llim[1]){
        if (verbose) cat('insufficient space points.\n')		
        flag[i] <- 2
      }else{
        if (verbose) cat('insufficient space & time points.\n')
        flag[i] <- 3
      }
    }
  }
  list(value=value,flag=flag,nb=nb)
}

# dev_ <- function(){
#   rm(list=ls())
#   library(ltsk)
#   source("dnb_2.R")
#   source("dnb_index_1.R")
#   library(RANN)
#   data(epa_cl)
#   obs_ <- obs[sample(251102,2000),c("x","y","t","pr_pm25")]
#   nb2_ <- dnb2(
#     query[,c("x","y","t")],
#     obs_[,1:3],
#     th=c(0.10,10),
#     by=c(0.01,0.01))
#   source("~/../Downloads/ltsk/R/dvariogram.R")
#   source("~/../Downloads/ltsk/R/dfitvariogram.R")
#   source("~/../Downloads/ltsk/R/dsample.pps.R")
#   source("~/../Downloads/ltsk/R/dsample.strata.R")
#   source("~/../Downloads/ltsk/R/dsmooth.variogram.R")
#   source("~/../Downloads/ltsk/R/working.smoothvariogram.R")
#   source("~/../Downloads/ltsk/R/working.fitvariog1.R")
#   source("~/../Downloads/ltsk/R/firstpeak.R")
#   source("~/../Downloads/ltsk/R/working.compvariogmodels1.R")
#   source("~/../Downloads/ltsk/R/vexpn.R")
#   source("~/../Downloads/ltsk/R/vmten.R")
#   source("~/../Downloads/ltsk/R/vsphn.R")
#   source("~/../Downloads/ltsk/R/vgaun.R")
#   source("~/../Downloads/ltsk/R/dadjustsills.R")
#   
#   source("dsubsample_2.R")
#   source("working_variog_1.R")
#   variog2_ <- working.variog2(query = as.matrix(query[,c("x","y","t")]),
#                              obs = as.matrix(obs_),nb=nb2_,vth = NULL,vlen = NULL,
#                              llim = c(3,3),verbose = T,Large = 2000
#                             )
#   plot(sapply(nb2_,length),sapply(variog2_$nb,length))
#   abline(0,1)
# }