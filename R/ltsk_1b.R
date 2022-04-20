ltsk <- function(query,obs,th,xcoord='x',ycoord='y',tcoord='t',zcoord='z',
                  by=c(0.001,0.001,1),nbin=NULL,byvariog=c(0.5,0.5,5),
                  subset=T,nmax=NULL,
                  vth=NULL,vlen=NULL,llim=c(3,3),
                  verbose=T,Large=2000,future=T,cl=NULL,cluster=c(2,10)){
  ## cluster: proportion of th[1] to cluster query points into
  # if(FALSE){
  #   rm(list=ls())
  #   library(ltsk)
  #   library(fields)
  #   library(RANN)
  #   dump_ <- sapply(list.files("./R",full=T),source)
  #   data(epa_cl)
  #   set.seed(1)
  #   #query <- query[1,]
  #   obs <- obs[sample(251102,2000),c("x","y","t","pr_pm25")]
  #   xcoord='x';ycoord='y';tcoord='t';zcoord='pr_pm25'
  #   th=c(0.1,10);by=c(0.001,0.001,1);nbin=NULL;
  #   byvariog=c(0.5,0.5,5);
  #   subset=T;nmax=NULL
  #   vth=NULL;vlen=NULL;llim=c(3,3)
  #   verbose=T;Large=2000;future=T  
  # }
  ## maintain interface as the package
  seed <- round(runif(1) * 1000000)
  ## extract input
  l.query <- check_input(query,xcoord,ycoord,tcoord,zcoord)
  l.obs <- check_input(obs,xcoord,ycoord,tcoord,zcoord)
  
  ## check missing for query values
  l.query <- check_na(l.query,'query')
  l.obs <- check_na(l.obs,'observed')
  
  if(is.null(cl)){ # No parallel computing
    if(is.null(byvariog)){ ## sequential runs
      ## Neighbor search
      l.nb <- working.dnb(l.query[,1:3],l.obs[,1:3],th,future,by,nbin,Large)

      ## Variogram estimation
      #source("working_variog_2.R")
      l.vgm <- working.variog(l.obs,l.nb$value,
                               vth,vlen,llim,verbose,Large,T)
      l.vgm$group <- NULL ## no grouping
      flag_ <- l.vgm$flag
      
      ## Kriging output
      #source("working_ltsk_2.R")
      l.krig <- working.ltsk(l.query[,1:3,drop=FALSE],l.obs,l.vgm,
                              nmax = nmax,subset = subset,verbose)
    }else{
      ## Neighbor search
      l.nb <- working.dnb(l.query[,1:3],l.obs[,1:3],th,future,by,nbin,Large)
      
      ## sub-sampling and time/space neighbor around each query
      #source("working_variog_2.R")
      l.sub <- working.variog(obs=l.obs,nb = l.nb$value,vth = vth,
                              vlen = vlen,llim = llim,verbose = verbose,
                              Large = Large,compute = FALSE)
      
      ## Group query points to calculate variogram
      ## index query points according to the given time-space resolution
      lq.index <- dnb.index(l.query[,1:3],by=byvariog)
      
      ## Extract neighbors from each group of query points
      lnb.group.0 <- split(l.nb$value,lq.index$value$uid)
      lnb.group <- lapply(lnb.group.0,function(x){
        unique(do.call(c,x))
      })
      
      ## Grouped variogram estimation
      #source("working_variog_2.R")
      lvgm.group <- working.variog(
        obs=l.obs,nb = lnb.group,vth=vth,vlen = vlen,llim=llim,
        verbose = verbose,Large = Large,compute = T)
      
      ## Impute flag for individual query with sufficient neighbor
      ## with the corresponding group flag:
      ## ASSUMPTION: the group estimates will success
      flag_ <- l.sub$flag
      uid_ <- lq.index$value$uid
      flagg_ <- lvgm.group$flag[uid_]
      flag_[is.na(flag_)] <- flagg_[is.na(flag_)]
      
      ## Group variogram estimates
      l.vgm <- list(value=lvgm.group$value,group=uid_,
                    flag=flag_,nb=l.sub$nb)
      
      ## Kriging output
      #source("working_ltsk_2.R")
      l.krig <- working.ltsk(
        l.query[,1:3,drop=FALSE],l.obs,l.vgm,nmax = nmax,
        subset = subset,verbose)
    }
  }else{ # first cluster points
   ## Parallel computing TBD
  }
  
  ## return
  out <- cbind(fit=l.krig[,1],se=l.krig[,2],flag=flag_)
  
  if(!is.null(attr(l.query,"na.action"))){
    ## remove missing query locations.
    query <- query[-attr(l.query,"na.action"),]
  }
  cbind(query,out)
}

# dev_ <- function(){
#   rm(list=ls())
#   library(ltsk)
#   library(fields)
#   library(RANN)
#   data(epa_cl)
#   set.seed(1)
#   obs_ <- obs[sample(251102,2000),c("x","y","t","pr_pm25")]
#   r0 <- ltsk(query,obs_,th=c(0.1,10),zcoord="pr_pm25")
#   dump_ <- sapply(list.files("./R",full=T),source)
#   source("working_variog_2.R")
#   source("working_ltsk_2.R")
#   source("ltsk_1b.R")
#   r1 <- ltsk2(query,obs_,th=c(0.1,10),zcoord="pr_pm25",byvariog = NULL)
#   table(r0$flag,r1$flag)
#   #     0    1    3
#   # 0  292    0    0
#   # 1    0   79    0
#   # 3    0    0 1910
#   summary(with(r0,fit[flag==0])-with(r1,fit[flag==0]))
#   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   # 0       0       0       0       0       0 
#   summary(with(r0,se[flag==0])-with(r1,se[flag==0]))
#   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   # 0       0       0       0       0       0 
#   
#   r2 <- ltsk2(query,obs_,th=c(0.1,10),zcoord="pr_pm25",byvariog = c(0.5,0.5,1))
#   table(r0$flag,r2$flag)
#   #     0    1    3
#   # 0  292    0    0
#   # 1    0   79    0
#   # 3    0    0 1910
#   plot(with(r0,fit[flag==0]),with(r2,fit[flag==0]),xlab="ltsk",ylab="ltsk2")
#   abline(0,1)
#   summary(1-with(r2,fit[flag==0])/with(r0,fit[flag==0]))
#   # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#   # -1.02229  0.00000  0.00000  0.08595  0.14291  1.19239 
#   plot(with(r0,se[flag==0]),with(r2,se[flag==0]),xlab="ltsk",ylab="ltsk2")
#   abline(0,1)
#   summary(1-with(r2,se[flag==0])/with(r0,se[flag==0]))
#   # Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#   # -0.512870 -0.067256  0.000000 -0.003705  0.107478  0.722037 
#   
# }