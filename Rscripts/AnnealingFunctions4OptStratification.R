getCriterion <- function(bh, b.rand, x, Esh, costs, sim, Zsim){
  strata<-findInterval(x,bh)+1
  Nh <- tapply(x,INDEX=strata,FUN=length)
  ah <- Nh/sum(Nh)
  if(sum(Nh<2)>0 | length(unique(strata))<length(bh)+1){
   criterion <- 1E20 
  }else{
    if(sim==TRUE) {
      s2h <- matrix(nrow=length(bh)+1,ncol=ncol(Zsim))
      for(i in (b.rand:(b.rand+1))){
        ids <- which(strata==i)
        s2h[i,] <- apply(Zsim[ids,],MARGIN=2,FUN=var)
      }
      
    }else{
      s2h <- tapply(x,INDEX=strata,FUN=var)
    }
  sh <- sqrt(s2h)
  Esh.new <- apply(sh,MARGIN=1,FUN=mean)
  Esh[!is.na(Esh.new)] <- Esh.new[!is.na(Esh.new)]
  ch <- tapply(costs,INDEX=strata,FUN=mean)
  sum.ini <- sum(ah*Esh/sqrt(ch))
  criterion <- sum(ah*Esh*sqrt(ch)*sum.ini)
  }
  list(Esh=Esh,crit=criterion)
}


permute.b<-function(bh, x)  {
  id <- sample(x = seq_len(length(bh)), size = 1)
  if (id==1) {
    bh[id] <- runif(n=1,min=min(x),max=bh[id+1])
  }  else if (id==length(bh)) {
    bh[id] <- runif(n=1,min=bh[id-1],max=max(x))
  } else {
    bh[id] <- runif(n=1,min=bh[id-1],max=bh[id+1])
  }
  list(b.rand=id,bh=bh)
}

Ospatsim<-function(bh_ini,
                   x,
                   costs = rep(1,length(x)),
                   sim=FALSE,
                   Zsim=NULL,
                   T_ini = 1,
                   coolingRate = 0.95,
                   maxPermuted=20*length(bh_ini),
                   maxNoChange=20*length(bh_ini),
                   verbose = getOption("verbose")) {
  
  # set initial temperature
  T <- T_ini
  
  # compute the minimisation criterion for initial stratum boundaries
  strata<-findInterval(x,bh_ini)+1
  Nh <- tapply(x,INDEX=strata,FUN=length)
  ah <- Nh/sum(Nh)
  ch <- tapply(costs,INDEX=strata,FUN=mean)
  if(sim==TRUE){
    s2h <- matrix(nrow=length(bh_ini)+1,ncol=ncol(Zsim))
    for(i in (1:(length(bh_ini)+1))) {
      ids <- which(strata==i)
      s2h[i,] <- apply(Zsim[ids,],MARGIN=2,FUN=var)
    }
  }else{
    s2h <- tapply(x,INDEX=strata,FUN=var)
  }
  sh <- sqrt(s2h)
  #compute average across simulations of stratum standard deviations
  Esh_cur <- apply(sh,MARGIN=1,FUN=mean)
  sum.ini <- sum(ah*Esh_cur/sqrt(ch))
  criterion_cur <- sum(ah*Esh_cur*sqrt(ch)*sum.ini)

  # Define structure for storing trace of criterion
  trace<-NULL
  
  # initialize number of zero changes of objective function
  nNoChange <-0
  
  bh_cur <- bh_ini
  
  # start cooling loop
  repeat{
    
    # initialize number of accepted configurations
    nAccepted <- 0
    
    # initialize number of permuted configurations
    nPermuted <- 0
    
    # initialize number of improved configurations
    nImproved <- 0
    
    # start permutation loop
    repeat {
      
      # increase the number of permutations
      nPermuted <- nPermuted + 1
      
      # propose new stratum boundaries
      res <- permute.b(bh_cur, x)
      bh_p <- res$bh
      b.rand <- res$b.rand #randomly selected id of strarum boundary that is changed
      
      # compute the criterion of this new stratification
      res <- getCriterion(bh_p,b.rand,x,Esh_cur,costs,sim,Zsim)
      criterion_p <- res$crit
      Esh_p <- res$Esh

      # accept/reject proposal by means of Metropolis criterion
      delta <- criterion_p - criterion_cur
      if (delta < 0) {
        nImproved <- nImproved + 1
        p <- 1 # always accept improvements
      } else {
        p <- exp(-delta / T) # use Boltzmann to judge if deteriorations should be accepted
      }
      u <- runif(n = 1) # draw uniform deviate
      if (u < p) { # accept proposal
        nAccepted <- nAccepted + 1
        bh_cur <- bh_p
        criterion_cur <- criterion_p
        Esh_cur <- Esh_p
      }
      # are conditions met to lower temperature?
      lowerTemperature <- (nPermuted == maxPermuted)
      if (lowerTemperature) {
        if (nImproved==0)
        {nNoChange<-nNoChange+1}
        else
        {nNoChange<-0}
        trace<-rbind(trace,criterion_cur)
        break  
      }
    }
    
    if (verbose) {
      cat(
        format(Sys.time()), "|",
        sprintf("T = %e  C = %e  permuted = %d  accepted = %d  improved = %d  acceptance rate = %f  \n",
                T, criterion_cur, nPermuted, nAccepted, nImproved, nAccepted / nPermuted)
      )
    }
    
    # check on convergence
    if (nNoChange == maxNoChange) {
      break
    }
    
    # lower temperature
    T <- coolingRate * T
  }
  
  # return result
  list(optbh=bh_cur,criterion = criterion_cur,trace=trace)
}
