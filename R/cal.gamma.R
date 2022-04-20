cal.gamma <-
function(query,obs,fout,nmax=NULL,subset= T)
{
  ## query : query point
  ## obs : neighbor points with a given thresholds
  ## fout : output fitted variogram
  ## nmax : for local kriging: the number of nearest observations that
  ##        should be used for prediction where nearest is defined
  ##        in terms of semi-variance of local model. By default all used.
  ## subset : subset neighbors using estimated sills 
  ## value : estimated gamma vector & matrix in ordinary kriging system
  ## product sum model
  loc0 <- matrix(query[1:2], ncol=2)
  tstamp <- matrix(obs[,3],ncol=1)
  coords <- matrix(obs[,1:2],ncol=2)
  dvec <- as.vector(rdist(loc0,coords))
  tvec <- abs( query[3] - obs[,3])
  if (subset) {
	  ii <- which( (dvec < fout$scoef[2]) | (tvec < fout$tcoef[2]) )
	  nn <- nrow(obs)
	  if((length(ii))<=1){
		  #cat('[cal.gamma: subset to null neighbors, switch to sample 10%]\n')
		  if(nn >= 50) obs <- obs[sort(sample(nn,nn/10)),]
	  } else{
		  obs <- matrix(obs[ii,],ncol=4)
	  }
	  tstamp <- matrix(obs[,3],ncol=1)
  	coords <- matrix(obs[,1:2],ncol=2)
  	dvec <- as.vector(rdist(loc0,coords))
 	  tvec <- abs( query[3] - obs[,3])
  }
  dmat <- dist(coords)
  tmat <- dist(tstamp)
  gout <- work.calgamma(dmat,dvec,tmat,tvec,fout,nmax)
  if(!is.null(nmax)){
    obs <- obs[gout$sub,]
  }
  n <- nrow(obs)
  Gamma <- tritomat(gout$Gamma,n)
  colnames(obs) <- c('x','y','t','z')
  list(Gamma = Gamma, gamma = gout$gamma, dat=obs)
}

# Debug_ <- function(){
#   rm(list=ls())
#   setwd("../")
#   library(ltsk)
#   library(fields)
#   library(RANN)
#   dump_ <- sapply(list.files("./R",full=T),source)
#   source("ltsk_1.R")
#   data(epa_cl)
#   obs$id <- seq(1,nrow(obs))
#   ## load cross-validation sets
#   load(file="test3.RData")
#   i <- 1
#   cat("processing fold ",i,"\n")
#   query_ <- obs[part==i,c("id","x","y","t","pr_pm25")]
#   obs_ <- obs[part!=i,c("x","y","t","pr_pm25")]
#   
#   ii <- with(query_,which(x==-82.04&y==41.03&t==1952))
#   debug(cal.gamma)
#   tmp <- ltsk2(query_[ii,],obs_,th=c(0.10,10),
#         zcoord = "pr_pm25",verbose = T,
#         nmax=20)
# }