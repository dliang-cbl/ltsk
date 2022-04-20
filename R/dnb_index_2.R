## Index data according to x/y/t
## value: original data with the coordinates
##       discretized at the desired resolution
## and meta data containing the time-stamp and locations
dnb.index <- function(obs_,nbin=NULL,by=NULL){
  ## nbin  : number of lon,lat,day to bin data
  
  tmp <- if( is.null(nrow(obs_))) matrix(obs_,1,3) else obs_
  colnames(tmp) <- c('x','y','t')
  
  lonr <- range(tmp[,1])+c(-1e-5,.001)
  latr <- range(tmp[,2])+c(-1e-5,.001)
  tstampr <- range(tmp[,3])+c(-.001,.001)

  if(!is.null(by)){
    lonb <- seq(lonr[1],lonr[2],by=by[1])
    latb <- seq(latr[1],latr[2],by=by[2])
    tstampb <- seq(tstampr[1],tstampr[2],by=by[3])
  }else{
    if(is.null(nbin)) nbin <- c(200,200,floor(diff(tstampr))+1)
    lonb <- seq(lonr[1],lonr[2],len=nbin[1]+1)
    latb <- seq(latr[1],latr[2],len=nbin[2]+1)
    tstampb <- seq(tstampr[1],tstampr[2],len=nbin[3]+1)
  }  
  
  lonid <- as.integer(findInterval(tmp[,1],lonb))
  latid <- as.integer(findInterval(tmp[,2],latb))
  tsid <- as.integer(findInterval(tmp[,3],tstampb))
  
  # ## surrogate location and time
  # helper_ <- function(x){
  #   ## find middle point of an increasing sequence
  #   x[-length(x)]+diff(x)/2
  # }
  # lon_ <- helper_(lonb)
  # lat_ <- helper_(latb)
  # tstamp_ <- helper_(tstampb)
  
  ## surrogates of locations
  if(FALSE){
    selx_ <- !duplicated(cbind(lonid,latid))
    xy_ <- cbind(lonid,latid,tmp[,c("x","y"),drop=F])[selx_,]
    if(is.null(dim(xy_))){
      xy_ <- matrix(xy_,1,length(xy_))
    }
  }else{
    xy__ <- as.matrix(tmp[,c("x","y"),drop=F])
    xy_ <- aggregate(xy__~lonid+latid,FUN=mean)
    # selx_ <- !duplicated(cbind(lonid,latid))
    # xy_ <- cbind(lonid,latid,x=lonb[lonid],y=latb[latid])[selx_,]
    # if(is.null(dim(xy_))){
    #   xy_ <- matrix(xy_,1,length(xy_))
    # }
  }
  colnames(xy_) <- c("xid","yid","x","y")
  
  ## create location id
  line_ <- as.data.frame(cbind(sid=1:nrow(xy_),xy_))
  sid_ <- data.frame(dummy=seq(1,length(lonid)),xid=lonid,yid=latid)
  line2_ <- merge(sid_,line_)
  
  #browser()
  
  ## surrogates of time
  if(FALSE){
    selt_ <- !duplicated(tsid)
    ts_ <- cbind(tsid,tmp[,"t",drop=T])[selt_,]
    if(is.null(dim(ts_))){
      ts_ <- matrix(ts_,1,length(ts_))
    }
  }else{
    ts_ <- aggregate(tmp[,"t",drop=T]~tsid,FUN=mean)
    # selt_ <- !duplicated(tsid)
    # ts_ <- cbind(tsid,tstampb[tsid])[selt_,]
    # if(is.null(dim(ts_))){
    #   ts_ <- matrix(ts_,1,length(ts_))
    # }
  }
  colnames(ts_) <- c("tsid","t")
  
  ## create time id
  line3_ <- as.data.frame(cbind(tid=1:nrow(ts_),ts_))
  tid_ <- data.frame(dummy=seq(1,length(tsid)),tsid=tsid)
  line4_ <- merge(line3_,tid_)  
  
  ## combine
  value_ <- merge(line2_,line4_,by="dummy")
  value_$dummy <- NULL
  
  uid0 <- with(value_,as.integer(factor(interaction(sid,tid))))
  value_$uid <- match(uid0,unique(uid0))
  
  list(xy=xy_,ts=ts_,
       value=value_,
       xbin=lonb,ybin=latb,
       tbin=tstampb)
}

# dev_ <- function(){
#   rm(list=ls())
#   library(ltsk)
#   data(epa_cl)
#   obs_ <- obs[sample(251102,2000),c("x","y","t")]
#   source("dnb_index_2.R")
#   debug(dnb.index)
#   r_ <- dnb.index(obs_,by=c(0.05,0.05,20))
#   with(r_,hist(value$x-xbin[value$xid]))
#   with(r_,summary(value$x-obs_[,"x"]))
#   with(r_,hist(value$y-ybin[value$yid]))
#   with(r_,summary(value$y-obs_[,"y"]))
#   with(r_,hist(value$t-tbin[value$tsid]))
#   with(r_,summary(value$t-obs_[,"t"]))
#   with(r_,range(value$tid))
#   dim(r_$ts)
#   with(r_,range(value$sid))
#   dim(r_$xy)
# }