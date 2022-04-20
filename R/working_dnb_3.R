## (1) allow clustering according to certain spatial threshold
## (2) allow time and space by argument when clustering.
working.dnb <- function(query,obs,th,future=TRUE,by=NULL,nbin=NULL,
                         Large=2000,cluster=FALSE){
  ## query : coordinate and time stamp of the query point
  ## obs   : coordinates and time stamps of observed points
  ## th    : space and time thresholds
  ## future : whether to include future temporal neighbors
  ## Large   : Number of unique locations to search for NN
  ##           AND the maximum number of neighbors to keep
  ## by    : vector long/lat/time resolutions to discretize the domain
  ## nbin  : number of long/lat/day to bin data 
  ##         (if by and nbin is NULL, default 200,200 and unique days)
  ## cluster: Whether to cluster query points
  ## value : a list of of potential neighbors and
  ##         cluster: the unique "clustered" query points (if requested)
  
  ## index query data according to x/y/t
  query.index <- dnb.index(query,by=by,nbin = nbin)
  
  if(cluster){
    ## cluster query points at given resolutions and the corresponding
    ## observed neighbors, in preparation for parallel processing
    
    ## Create a time space id
    uid0 <- with(query.index$value,as.integer(factor(interaction(sid,tid))))
    uid <- match(uid0,unique(uid0))
    
    ## make a list of site and time id for each time/space index point
    coord0_ <- cbind(query.index$value,uid=uid)
    coord_ <- subset(coord0_,!duplicated(uid))
    
  }else{
    ## for approximate neighbor search instead
    uid <- NULL
    coord_ <- query.index$value
  }
  
  ## make a list of neighbor id for each time/space index point
  r_ <- vector("list",nrow(coord_))
  
  ## index observed data according to chosen resolution per unit time
  obs.index <- dnb.index(obs,by=by,nbin = nbin)
  
  ## approximate spatial neighbor searches
  ns_ <- nrow(obs.index$xy)
  ks_ <- min(ns_,Large)
  th_s_ <- ifelse(cluster,2*max(th[1],by[1:2]),th[1])
  nns_ <- nn2(obs.index$xy[,3:4,drop=FALSE],
              query.index$xy[,3:4,drop=FALSE],
              k=ks_,searchtype="radius",
              radius=th_s_-sqrt(.Machine$double.eps))
  
  if(!any(nns_$nn.idx>0)){
    ## return NULL if no spatial neighbor
    return(list(value=r_,cluster=uid))
  }
  
  ## Create a nb list at spatial index level
  s0_ <- split(nns_$nn.idx[nns_$nn.idx>0],
               row(nns_$nn.idx)[nns_$nn.idx>0])
  
  ## approximate temporal neighbor search
  if(cluster){
    nnt_ <- nn2(obs.index$ts[,"t",drop=FALSE],
                query.index$ts[,"t",drop=FALSE],
                k=nrow(obs.index$ts),
                searchtype = "radius",radius=2*max(th[2],by[3])-0.001)
    ## minus 0.001 to be consistent with dnb of not including H in search
    
    if(!any(nnt_$nn.idx>0)){
      ## return NULL if no temporal neighbor
      return(list(value=r_,cluster=uid))
    }
    
    ## Create a nb list at temporal index level
    t0_ <- split(nnt_$nn.idx[nnt_$nn.idx>0],
                 row(nnt_$nn.idx)[nnt_$nn.idx>0])
  }else{
    tmat_ <- outer(obs.index$ts[,"t",drop=TRUE],
                   query.index$ts[,"t",drop=TRUE],
                   FUN="-")
    if(future){
      iit_ <- abs(tmat_)<th[2]
    }else{
      iit_ <- (tmat_ > - th[2]) & (tmat_ <= 0)
    }
    if(!any(iit_)){
      ## return NULL if no temporal neighbor
      return(list(value=r_,cluster=uid))
    }
    ## Create a nb list at temporal index level
    t0_ <- split(iit_,col(iit_))
    t0_ <- lapply(t0_,which)
  }
  
  ## Expand the nb list to spatial observation level
  s1_ <- lapply(s0_,function(elmt){
    which(obs.index$value$sid %in% elmt)
    # if(FALSE){
    #   nb_ <- which(obs.index$value$id %in% elmt)
    #   summary(obs[nb_,"x"]-query.index$xy[1,"x"])
    # }
  })
  
  ## Expand the nb list to temporal observation level
  t1_ <- lapply(t0_,function(elmt){
    which(obs.index$value$tid %in% elmt)
  })
  
  ## Create a query level temporal observation neighbor
  tindex_ <- as.integer(names(t1_))
  for(i in 1:length(t1_)){
    w_ <- which(coord_$tid==tindex_[i])
    r_[w_] <- rep(list(t1_[[i]]),length(w_))
  }
  
  ## subset the query level to spatial observation 
  sindex_ <- as.integer(names(s1_))
  for(i in 1:length(s1_)){
    w_ <- which(coord_$sid==sindex_[i])
    for(j in w_){
      r_[[j]] <- intersect(r_[[j]],s1_[[i]])
    }
  }
  
  list(value=r_,cluster=uid)
  
}

# dev_ <- function(){
#   rm(list=ls())
#   library(ltsk)
#   source("dnb_index_2.R")
#   library(RANN)
#   data(epa_cl)
#   
#   obs_ <- obs[sample(251102,2000),c("x","y","t")]
#   
#   hold_ <- vector("list",nrow(query))
#   for(j in 1:nrow(query)){
#     hold_[[j]] <- dnb(query[j,2:4],obs_,th=c(0.10,10))
#   }
#   nb1 <- sapply(hold_,length)
#   
#   source("working_dnb_3.R")
#   tmp_ <- working.dnb(
#     query[,c("x","y","t")],
#     obs_,
#     th=c(0.10,10),
#     by=c(0.01,0.01,1))
#   
#   d_ <- matrix(NA,length(tmp_),3)
#   for(i in 1:length(tmp_)){
#     d_[i,1] <- max(query[i,"x"]-obs_[tmp_$value[[i]],"x"])
#     d_[i,2] <- max(query[i,"y"]-obs_[tmp_$value[[i]],"y"])
#     d_[i,3] <- max(query[i,"t"]-obs_[tmp_$value[[i]],"t"])
#   }
#   nb2 <- sapply(tmp_$value,length)
#   plot(nb1,nb2,xlim=range(nb1,nb2),ylim=range(nb1,nb2))
#   abline(0,1)
#   mean((nb1-nb2)^2)
#   identical(lapply(nb1,sort),lapply(nb2,sort))
#   Check_ <- function(x){
#     xy_ <- as.matrix(query[x,c("x","y"),drop=FALSE])
#     oxy_ <- as.matrix(obs_[tmp_$value[[x]],c("x","y")])
#     dx_ <- xy_[1,1] - oxy_[,1]
#     dy_ <- xy_[1,2] - oxy_[,2]
#     print(sqrt(sum(dx_^2+dy_^2)))
#   }
#   Check_(2281)
#   
#   ## check Cluster = TRUE
#   source("working_dnb_3.R")
#   #debug(working.dnb)
#   tmp2_ <- working.dnb(
#     query[,c("x","y","t")],
#     obs_,
#     th=c(0.10,10),
#     by=c(0.15,0.15,30),
#     cluster = TRUE)
#   
#   table(tmp2_$cluster)
#   
#   d2_ <- matrix(NA,length(tmp2_$value),3)
#   for(i in 1:length(tmp2_$value)){
#     if(length(tmp2_$value[[i]])>0){
#       d2_[i,1] <- max(outer(query[tmp2_$cluster==i,"x"],
#                             obs_[tmp2_$value[[i]],"x"],FUN="-"))
#       d2_[i,2] <-  max(outer(query[tmp2_$cluster==i,"y"],
#                              obs_[tmp2_$value[[i]],"y"],FUN="-"))
#       d2_[i,3] <-  max(outer(query[tmp2_$cluster==i,"t"],
#                              obs_[tmp2_$value[[i]],"t"],FUN="-"))
#     }
#   }
#   summary(d2_)
# }