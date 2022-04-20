work.calgamma <-
function(dmat,dvec,tmat,tvec,fout,nmax)
{
  ## dmat : spatial distance matrix
  ## dvec : spatial distance to query point
  ## tmat : temporal distance matrix
  ## tvec : temporal distance to query point
  ## fout : fitted variogram object
  ## nmax : for local kriging: the number of nearest observations that
  ##        should be used for prediction where nearest is defined
  ##        in terms of semi-variance of local model. By default all used.
  ## value : estimated Gamma and gamma for ordinary kriging
  vfs0 <- get(paste('v',fout$smodel,sep=''))
  vf0t <- get(paste('v',fout$tmodel,sep=''))
  
  gd0 <- with(fout,as.vector(vfs0(dvec,scoef[1],scoef[2],scoef[3])))
  g0t <- with(fout,as.vector(vf0t(tvec,tcoef[1],tcoef[2],tcoef[3])))
  fitgamma <- gd0 + g0t - fout$k * gd0 * g0t
  
  ## distance matrices
  n_ <- length(fitgamma)
  
  #browser()
  sub_ <- NULL
  if(!is.null(nmax)){
    ## Local Kriging 
    ## following gstat convention to scale time
    ## and space using fitted variogram 
    r_ <- rank(fitgamma,ties.method = "random")
    sub_ <- which(r_<=nmax)
    fitgamma <- fitgamma[sub_]
    helper <- function(x,n,sub){
      m <- tritomat(x,n)
      m2 <- m[sub,sub]
      as.dist(m2)
    }
    dmat <- helper(dmat,n_,sub_)
    tmat <- helper(tmat,n_,sub_)
  }
  
  Gd0 <- with(fout,as.vector(vfs0(dmat,scoef[1],scoef[2],scoef[3])))
  G0t <- with(fout,as.vector(vf0t(tmat,tcoef[1],tcoef[2],tcoef[3])))
  fitGamma <- Gd0 + G0t - fout$k * Gd0 * G0t
  
  chk <- with(fout,(scoef[3]==0) & (tcoef[3]==0))
  if(chk){
    ## zero nuggets add ad-hoc small nuggets effect
    fitGamma <- fitGamma + 1e-5
  }	
  ## fitted variogram and subset for local Kriging.
  list(Gamma = fitGamma, 
       gamma = fitgamma,
       sub=sub_)
}
