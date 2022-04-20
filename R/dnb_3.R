## cl: a parallel object or NULL, Allow parallel computing eventually
## cluster: if parallel, the time and space domain to run in parallel
##          in terms of th
dnb <- function(query,obs,th,xcoord='x',ycoord='y',tcoord='t',
                 future=TRUE,by=NULL,nbin=NULL,Large=2000,
                 cl=NULL,cluster=c(1,3)){
  ## extract input
  l.query <- check_input(query,xcoord,ycoord,tcoord,"zcoord")
  l.obs <- check_input(obs,xcoord,ycoord,tcoord,"zcoord")
  
  ## check missing for query values
  l.query <- check_na(l.query,'query')
  
  ## replace missing values in observed
  replaceNA <- function(x,delta){
    x[is.na(x)] <- min(x,na.rm=T)-2*delta
    x
  }
  if(any(is.na(l.obs))){
    l.obs <- cbind(x=replaceNA(l.obs[,1],th[1]),
                   y=replaceNA(l.obs[,2],th[1]),
                   t=replaceNA(l.obs[,3],th[2]),
                   z=-99999)
  }

  ## neighbor search 
  if(is.null(cl)){ # when not doing parallel computing
    res_ <- working.dnb(query=l.query[,1:3,drop=FALSE],
                         obs=l.obs[,1:3,drop=FALSE],
                        th=th,future = future,by=by,
                        nbin = nbin,Large = Large, 
                        cluster = FALSE)
    nb_ <- res_$value
  }
  else if("cluster" %in% class(cl)){  ## doing parallel computing
    # first cluster points
    part_ <- working.dnb(query=l.query[,1:3,drop=FALSE],
                         obs=l.obs[,1:3,drop=FALSE],
                         th=th,future = TRUE, 
                         by=cluster[c(1,1,2)]*th[c(1,1,2)],
                         nbin=NULL,Large=Large, cluster=TRUE)
    ## return neighbor list
    nb_ <- vector("list",nrow(l.query))
    
    ## number of neighbors to each cluster
    nb_cluster_ <- sapply(part_$value,length)
    if(!any(nb_cluster_>0)){
      ## no neighbor, return NULL
      return(nb_)
    }
    
    ## helper function to run neighbor search
    run_ <- function(x,y,z){
      ## x: matrix of query points
      ## y: matrix of observed points (could be zero rows)
      ## z: index of observed points (could be null)
      ## return: neighbor size
      output <- vector("list",nrow(x))
      if(length(z)>0){
        ## for each cluster find the neighbors
        lnb_ <- working.dnb(
          query=x,obs=y,th=th,future=future,by=by,nbin=nbin,Large=Large,
          cluster = FALSE)
        ## translate the id to original observed
        for(j in 1:length(lnb_$value)){
          output[[j]] <- z[lnb_$value[[j]]]
        }
      }
      return(output)
    }
    
    ## make a list of clusters and corresponding observed
    nc_ <- length(part_$value)
    ll.query <- vector('list',nc_)
    ll.obs <- vector('list',nc_)
    ll.order <- vector("list",nc_)
    for(i in 1:nc_){
      ll.order[[i]] <- which(part_$cluster==i)
      ll.query[[i]] <- l.query[part_$cluster==i,1:3,drop=FALSE]
      ll.obs[[i]] <- l.obs[part_$value[[i]],1:3,drop=FALSE]
    }
    
    ## neighbor query
    if(FALSE){
      for(i in 1:nc_){
        nb_[part_$cluster==i] <- run_(
          x=ll.query[[i]],
          y=ll.obs[[i]],
          z=part_$value[[i]]
        )
      }
    }else{
      ## run in parallel
      out1 <- clusterMap(cl=cl,fun=run_,ll.query,ll.obs,part_$value)
      #browser()
      for(i in 1:nc_){
        nb_[ll.order[[i]]] <- out1[[i]]
      }
    }
  }else{
    stop("Please check cl argument.\n")
  }
  
  ## return
  l.query.na <- attr(l.query,"na.action")
  if(!is.null(l.query.na)){
    nb2_ <- vector("list",nrow(query))
    nb2_[-l.query.na] <- nb_
    return(nb2_)
  }else{
    return(nb_)
  }
}
# dev_ <- function(){
#   rm(list=ls())
#   library(ltsk)
#   library(RANN)
#   data(epa_cl)
#   set.seed(1)
#   obs_ <- obs[sample(251102,2000),c("x","y","t")]
#   source("dnb_index_2.R")
#   source("./R/check_input.R")
#   source("dnb_3.R")
#   source("working_dnb_3.R")
#   #source("working_part_1.R")
#   nb0_ <- dnb(query,obs_,th=c(0.1,10))
#   nb0_[[2]]
#   
#   library(parallel)
#   cl <- makeCluster(6)
#   junk_ <- clusterEvalQ(cl,{
#     source("dnb_index_2.R")
#     source("working_dnb_3.R")
#     library(RANN)
#   })
#   nb1_ <- dnb(query,obs_,th=c(0.1,10),cl=cl,cluster = c(0.5,0.5))
#   nb1_[[2]]
#   n0_ <- sapply(nb0_,length)
#   n1_ <- sapply(nb1_,length)
#   table(n1_-n0_)
#   query[1,"x"] <- NA; query[2,"y"] <- NA
#   obs_[199,"x"] <- NA; obs_[434,"y"] <- NA
#   #debug(dnb2)
#   nb2_ <- dnb(query,obs_,th=c(0.1,10))
#   nb3_ <- dnb(query,obs_,th=c(0.1,10),cluster = NULL)
#   summary(sapply(nb2_,length)-sapply(nb3_,length))
#   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   # 0       0       0       0       0       0 
#   nb2_[[1]]
#   nb3_[[1]]
#   any(sapply(nb2_,function(x) any(x %in% c(199,434))))
#   #[1] FALSE
#   stopCluster(cl)
# }
# 
# debug_ <- function(){
#   rm(list=ls())
#   library(ltsk)
#   library(RANN)
#   source("dnb_index_2.R")
#   source("./R/check_input.R")
#   source("dnb_2.R")
#   source("working_dnb_3.R")
#   load(file="ltsk_test_debug_1.rData")
#   tmp_ <- dnb(query_,obs_,th=c(0.1,10),cluster = c(0.5,0.5))
# }