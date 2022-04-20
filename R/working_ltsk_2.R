working.ltsk <- function(query,obs,vgm,nmax,subset,verbose){
  ## query : coordinate and time stamp of the query point
  ## obs   : coordinates and time stamps and data of observed points
  ## vgm   : fitted product sum variogram object
  ## nmax : for local kriging: the number of nearest observations that
  ##        should be used for prediction where nearest is defined
  ##        in terms of semi-variance of local model. By default all used.
  ## subset : subset neighbors using estimated sills 
  ## value : ordinary kriging and standard error for querys with available neighbors
  
  fit <- matrix(NA,nrow(query),2)
  ## extract the variogram
  if(is.null(vgm$group)){
    ## variogram is not grouped
    group_ <- seq(1,length(vgm$value))
  }else{
    ## variogram is grouped
    group_ <- vgm$group
  }
  for(i in which(vgm$flag==0)){
    #browser()
    if(verbose){
      cat(query[i,],"\n")
    }
    j <- group_[i]
    gout <- cal.gamma(query[i,],obs[vgm$nb[[i]],],vgm$value[[j]],nmax,subset)
    fit[i,1:2] <- with(gout,work.kriging(Gamma,gamma,dat[,4]))
  }
  fit
}
# dev_ <- function(){
#   rm(list=ls())
#   library(ltsk)
#   library(fields)
#   source("dnb_2.R")
#   source("working_dnb_2.R")
#   source("dnb_index_1.R")
#   library(RANN)
#   data(epa_cl)
#   obs_ <- obs[sample(251102,2000),c("x","y","t","pr_pm25")]
#   nb2_ <- working.dnb2(
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
#   source("dsubsample_2.R")
#   source("working_variog_1.R")
#   variog2_ <- working.variog2(query = as.matrix(query)[,c("x","y","t")],
#                               obs = as.matrix(obs_),nb=nb2_,vth = NULL,vlen = NULL,
#                               llim = c(3,3),verbose = T,Large = 2000
#   )
#   plot(sapply(nb2_,length),sapply(variog2_$nb,length))
#   abline(0,1)
#   
#   source("~/../Downloads/ltsk/R/tritomat.R")
#   source("~/../Downloads/ltsk/R/work.kriging.R")
#   source("working_ltsk_1.R")
#   source("cal.gamma.R")
#   source("work.calgamma.R")
#   krig2_ <- working.ltsk2(query = as.matrix(query)[,c("x","y","t")],
#                           obs = as.matrix(obs_),vgm=variog2_,
#                           nmax=NULL,subset=T)
#   
#   krig_ <- ltsk(query = query[,c("x","y","t")],
#                 obs = obs_,th = c(0.10,10),zcoord = "pr_pm25")
#   
#   identical(krig_$flag,variog2_$flag)
#   #[1] TRUE
#   with(krig_,plot(fit[flag==0],krig2_[flag==0,1],xlab="ltsk",
#                   ylab="ltsk2"))
#   abline(0,1)
#   with(krig_,plot(se[flag==0],krig2_[flag==0,2],xlab="ltsk",
#                   ylab="ltsk2"))
#   abline(0,1)
#   
# }