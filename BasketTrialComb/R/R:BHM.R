#' Get posterior samples from BHM Method
#' 
#' @import rjags
#' @param cohortsize vector format, sample size for each cohort per disease type 
#' @param ntype number of disease types/subgroups     
#' @param p.null the vector of null hypothesis rate in each type
#' @param p.target the vector of alternative hypothesis rate in each type
#' @param ntrial number of simulated trials needs to be 1 here, repeated simulations are done outside    
#' @param n.burnin control prameters of MCMC sampling using jags
#' @param n.chain control prameters of MCMC sampling using jags
#' @param n.iter control prameters of MCMC sampling using jags
#' @param mu.par the mean parameter of mu's prior distribution    
#' @param eff.ref the efficacy evaluation criteria, calibrated through simulations 
#' @param r.all number of responses   
#' @return Probability to stop
#' @export
berrymethod<-function(cohortsize, ntype = 4, p.true, p.null, p.target, ntrial = 1, n.burnin = 2000, 
                      n.chain = 1, n.iter = 30000, mu.par, eff.ref, r.all){
  
  cat('model
  {
      for (j in 1:ntype){                
      y[j]~ dbin(p[j],n[j])
      p[j]<-exp(theta[j]+log(p.target[j]/(1-p.target[j])))/(1+exp(theta[j]+log(p.target[j]/(1-p.target[j]))))
      
      }
      for (j in 1:ntype){
      theta[j]~ dnorm(mu,tau)      
      }
      sigma2<-1/tau                    #variance=sigma2
      mu ~ dnorm(mu.par, 0.01)         #prior on mu
      tau ~ dgamma(0.0005,0.000005)    #prior on tau
}', file={fbhm <- tempfile()})
  
  p.est<-matrix(0, nrow=ntrial, ncol=ntype)                   # store the estimated response rate for each disease type
  sample.size<-matrix(0, nrow=ntrial,ncol=ntype)              # store the maximum sample size used for each disease type
  efficacy<-matrix(0,nrow=ntrial,ncol=ntype)                  # store the decision go or not go for each disease type
  eff.prob.store<-matrix(0,nrow=ntrial,ncol=ntype)            # store the efficacy probability to go for each disease type
  
  ncohort<-length(cohortsize)
  futstop<-0.05                                              # futility stopping cutoff
  
  
  for (trial in 1:ntrial){
    
    #set.seed(100+trial)
    n<-numeric(ntype)
    y<-numeric(ntype)
    stopping<-numeric(ntype)
    presponse<-p.true
    csize<-matrix(cohortsize, nrow=ncohort, ncol=ntype)
    
    ## added for futility control
    cont<-1:ntype
    stop<-NULL
    
    for (i in 1:ncohort){
      
      # generate data for a new cohort of patients
      y = y + r.all[i,];
      n = n + csize[i,];
      
      if (i != ncohort){ # interim analysis
        jags.data<-list("y"=y,"n"=n,"ntype"=ntype,"p.target"=p.target,"mu.par"=mu.par)
        jags.fit<- jags.model(fbhm, data=jags.data, n.adapt = n.burnin, n.chains = n.chain, quiet = T)
        post.samples <- coda.samples(jags.fit, variable.names = c("p"), n.iter = n.iter/n.chain, progress.bar = "none")
        pres.est <- do.call(rbind, post.samples)
        fut.prob<-sapply(seq(1,ntype,1), function(x) mean(pres.est[,x]>(p.null[x]+p.target[x])/2)) # calculate the futility probability
        
        stop_temp <- cont[which(fut.prob[cont]<futstop)] ##index of arms stopped at 'interim' analysis, early stop for futility.
        stop <- c(stop, stop_temp)  ##store all stopped arms
        cont <- cont[!cont%in%stop] ##only the continued arms remain
        efficacy[trial, stop_temp] <- 0  ##set the futile arms not go decision
        #presponse[stop_temp]<-0 ##stop
        r.all[(i+1),stop_temp]<-0
        csize[(i+1),stop_temp]<-0
        
        
        if (length(cont)==0){ # when all groups have stopped for futility
          efficacy[trial,]<-0                                         
          p.est[trial,]<-apply(pres.est,2,mean)                                    # store the rate estimate for each disease type
          sample.size[trial,]<-n                                                       # store the maximum sample size used for each disease subgroup
          break                                                                        # stop simulating patients if all subgroups have stopped
        }  
      } else{ # final analysis
        jags.data<-list("y"=y,"n"=n,"ntype"=ntype,"p.target"=p.target,"mu.par"=mu.par)
        jags.fit<- jags.model(fbhm, data=jags.data,n.adapt = n.burnin, n.chains = n.chain, quiet = T)
        post.samples <- coda.samples(jags.fit, variable.names = c("p"), n.iter = n.iter/n.chain, progress.bar = "none")
        pres.est <- do.call(rbind, post.samples)
        
        eff.prob<-sapply(seq(1,ntype,1), function(x) mean(pres.est[,x]>p.null[x]))
        eff.prob.store[trial,]<-eff.prob
        
        go <- cont[which((eff.prob > eff.ref)[cont])]  ##index of arm satisfies efficacy boundary
        nogo <- cont[which((eff.prob <= eff.ref)[cont])]
        efficacy[trial,go]<-1
        efficacy[trial,nogo]<-0
        
        p.est[trial,]<-apply(pres.est,2,mean)                    # store the rate estimate for each disease type      
        sample.size[trial,]<-n                                   # store the maximum sample size used for each disease type      
      }  
    }
  }
  
  
  
  results<-cbind(sample.size,efficacy,eff.prob.store,p.est) 
  # sample size is the sample size for each cohort under all trial repetitions. 
  # efficacy stores 0(not efficacious) or 1(efficacious) for each cohort under all repetitions. 
  # eff.prob.store is the proportion of posterior probability of being greater than p.null. p.est is the posterior estimate of response rate. 
  return(results)
  
}


#' Simulation for BHM Method
#'
#' @param cohortsize vector format, sample size for each cohort per disease type 
#' @param ntype number of disease types/subgroups     
#' @param p.true the vector of true response rate in each disease type  
#' @param p.null the vector of null hypothesis rate in each type
#' @param p.target the vector of alternative hypothesis rate in each type
#' @param ntrial the number of simulated trials      
#' @param n.burnin control prameters of MCMC sampling using jags
#' @param n.chain control prameters of MCMC sampling using jags
#' @param n.iter control prameters of MCMC sampling using jags
#' @param mu.par the mean parameter of mu's prior distribution    
#' @param FWER: the pre-specified family wise error rate for efficacy cutoff calculation#
#' @param  precision: the grid space controlled in the grid search, can be refined       #
# if higher precision is needed                                                 #
#' @param  fwer.threshold: the acceptable threshold for the difference between           #
# estimated FWER and proposed FWER 
#' @return probability of reject
#' @export
berrymethod.calibrate<-function(cohortsize,ntype=4,p.true, p.null, p.target,
                                ntrial=10000, n.burnin=2000, n.chain=1,n.iter=30000, mu.par, 
                                FWER, precision=0.0001, fwer.threshold=0.001){
  cat('model
  {
      for (j in 1:ntype){                
      y[j]~ dbin(p[j],n[j])
      p[j]<-exp(theta[j]+log(p.target[j]/(1-p.target[j])))/(1+exp(theta[j]+log(p.target[j]/(1-p.target[j]))))
      
      }
      for (j in 1:ntype){
      theta[j]~ dnorm(mu,tau)      
      }
      sigma2<-1/tau                    #variance=sigma2
      mu ~ dnorm(mu.par, 0.01)         #prior on mu
      tau ~ dgamma(0.0005,0.000005)    #prior on tau
}', file={fbhm <- tempfile()})
  
  eff.prob.store<-matrix(0,nrow=ntrial,ncol=ntype)
  ncohort<-length(cohortsize)
  futstop<-0.05                                              # futulity stopping cutoff
  
  for (trial in 1:ntrial){
    set.seed(10100+trial)
    n<-numeric(ntype)
    y<-numeric(ntype)
    stopping<-numeric(ntype)
    presponse<-p.true
    csize<-matrix(cohortsize, nrow=ncohort, ncol=ntype)
    ## added for futility control
    cont<-1:ntype
    stop<-NULL
    for (i in 1:ncohort){
      # generate data for a new cohort of patients
      y = y + rbinom(rep(1,ntype), cohortsize[i], presponse);
      n = n + csize[i,];
      if (i != ncohort){ # interim analysis
        jags.data<-list("y"=y,"n"=n,"ntype"=ntype,"p.target"=p.target,"mu.par"=mu.par)
        jags.fit<- jags.model(fbhm, data = jags.data,
                              n.adapt = n.burnin, n.chains = n.chain, quiet = T)
        post.samples <- coda.samples(jags.fit, variable.names = c("p"), n.iter = n.iter/n.chain, progress.bar = "none")
        pres.est <- do.call(rbind, post.samples)
        fut.prob<-sapply(seq(1,ntype,1), function(x) mean(pres.est[,x]>(p.null[x]+p.target[x])/2)) # calculate the futility probability
        
        stop_temp <- cont[which(fut.prob[cont]<futstop)] ##index of arms stopped at 'interim' analysis, early stop for futility.
        stop <- c(stop, stop_temp)
        cont <- cont[!cont%in%stop] ##only the continued arms remain
        eff.prob.store[trial, stop_temp] <- 0
        presponse[stop_temp]<-0
        csize[(i+1),stop_temp]<-0
        if (length(cont)==0){ # when all groups have stopped for futility
          break                                                                        # stop simulating patients if all subgroups have stopped
        }
      } else{ # final analysis
        jags.data<-list("y"=y,"n"=n,"ntype"=ntype,"p.target"=p.target,"mu.par"=mu.par)
        jags.fit<- jags.model(fbhm, data = jags.data,
                              n.adapt = n.burnin, n.chains = n.chain, quiet = T)
        post.samples <- coda.samples(jags.fit, variable.names = c("p"), n.iter = n.iter/n.chain, progress.bar = "none")
        pres.est <- do.call(rbind, post.samples)
        
        eff.prob<-sapply(seq(1,ntype,1), function(x) mean(pres.est[,x]>p.null[x]))
        eff.prob.store[trial,cont]<-eff.prob[cont]
      }
    }
  }
  ###grid search for efficacy cutoff
  w0 = 1-(1-FWER)^(1/4)
  Grid = seq(w0-0.01, FWER,by = precision) # one can modify this for different methods
  fwer.pred = sapply(Grid, function(o) fwer(eff.prob.store,o))
  fwer.pred.left = fwer.pred[fwer.pred < FWER]
  anchor.grid = which.min(abs(fwer.pred.left - FWER))
  #anchor.grid = which.min(abs(fwer.pred-FWER))
  if(abs(fwer.pred.left[anchor.grid]-FWER)>fwer.threshold){
    print("Precision is not small enough or simulation number is not big enough, please reput again")
    eff.cut = sapply(1:ntype, FUN = function(x){quantile(eff.prob.store[,x], 1-Grid[anchor.grid])})
    results = list(prob.success = eff.prob.store, eff.ref=eff.cut)
    #return(results)
    #return(eff.prob.store)
  }else{
    eff.cut = sapply(1:ntype, FUN = function(x){quantile(eff.prob.store[,x], 1-Grid[anchor.grid])})
    results = list(prob.success = eff.prob.store, eff.ref=eff.cut)
    #return(results)
  }
  plot_search_type1_bhm = ggplot2::ggplot(data=data.frame(Grid, fwer.pred), ggplot2::aes(x = Grid, y = fwer.pred, group=1)) +
    ggplot2::geom_line()+
    ggplot2::geom_point() + ggplot2::ggtitle("Family wise error rate with grid search")+
    ggplot2::labs(x="Fixed type one error for each cohort", y="FWER")
  
  return(list(results, plot_search_type1_bhm))
}


