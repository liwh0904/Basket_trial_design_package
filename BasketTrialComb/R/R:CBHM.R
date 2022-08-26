#' Function to decide tuning parameter a and b
#' 
#' @import rjags
#' @param cohortsize vector format, sample size for each cohort per disease type 
#' @param ntype number of disease types/subgroups  
#' @param ntrial the number of simulated trials   
#' @param p0 the value of null hypothesis rate   
#' @param p1 the value of alternative hypothesis rate  
#' @param var.small the prespecified small value of shrinkage parameter sigma^2 for calibrating a and b  
#' @param var.big the prespecified big value of shrinkage parameter sigma^2 for calibrating a and b                                       
#' @return tuning parameters a and b
#' @export
decidePar<-function(cohortsize, ntype, ntrial, p0, p1, var.small, var.big){
  
  
  
  presponse <- matrix(NA, nrow= ntype, ncol = ntype)
  for (i in 1:(ntype-1)){
    presponse[i,]<- c(rep(p0, i),rep(p1, ntype-i))
  }
  presponse[ntype,] <-  rep(p1, ntype)
  ncohort <- length(cohortsize)
  medianT<-NULL
  
  #----------------------------#
  for (j in 1:ntype){ ##loop for all possible combinations of response rates
    test.stat<-matrix(NA,nrow=ntrial,ncol=ncohort)
    for (sim in 1:ntrial){
      
      set.seed(20200+ntrial)
      n<-numeric(ntype); y<-numeric(ntype)
      for (i in 1:ncohort){
        
        y = y + rbinom(rep(1,ntype), cohortsize[i], presponse[j,]);
        n = n + cohortsize[i];
        
        p<-sum(y)/sum(n)  ##observed response rate using pooled data; assume all cohorts follow the same distribution.
        x<-cbind(y, n-y)
        E <- cbind(n * p, n * (1 - p)) ##expected responder/non-responder depend on the common response rate p. 
        T <- sum((abs(x - E))^2/E)
        if (is.nan(T)){T<-0}
        test.stat[sim,i]<-T       # store the test statistic  
        
      }
    }
    medianT<-c(medianT, apply(test.stat, 2, median)[ncohort])
    
    
  }
  heteroT <- min(medianT[1:(ntype-1)]); homoT <- medianT[ntype]
  b<-(log(var.big)-log(var.small))/(log(heteroT)-log(homoT)) 
  a<-log(var.small)-b*log(homoT)
  
  results<-list(a = a, b=b)
  results
}



#' CBHMdesign() Function
#'
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
#' @param a the tuning parameter that characterize the relationship between test       #
# statistic T and the shrinkage parameter sigma^2   
#' @param b the tuning parameter that characterize the relationship between test       #
# statistic T and the shrinkage parameter sigma^2    
#' @param r.all number of responses 
#' @return reponses and sample sizes
#' @export
CBHMdesign<-function(cohortsize,ntype = 4, p.true, p.null, p.target, ntrial = 1, n.burnin, n.chain, 
                     n.iter, mu.par, eff.ref, a, b, r.all){
  cat('model
      {
      for (j in 1:ntype){                #J=4, the number of arms
      y[j]~ dbin(p[j],n[j])
      p[j]<-exp(theta[j])/(1+exp(theta[j]))
      
      }
      for (j in 1:ntype){
      theta[j] ~ dnorm(mu,tau)      #hierarchical model for theta
      }
      mu ~ dnorm(mu.par, 0.01)         #prior on mu
}', file={fcbhm <- tempfile()})
  
  ncohort<-length(cohortsize)
  futstop<-0.05                                              # futility stopping cutoff
  
  p.est<-matrix(0, nrow=ntrial, ncol=ntype)                   # store the estimated response rate for each disease type
  sample.size<-matrix(0, nrow=ntrial,ncol=ntype)              # store the maximum sample size used for each disease type
  test.stat<-matrix(0, nrow=ntrial, ncol=ncohort)             # store the test statistic value at each interim analysis
  efficacy<-matrix(0,nrow=ntrial,ncol=ntype)                  # store the decision go or not go for each disease type
  eff.prob.store<-matrix(0,nrow=ntrial,ncol=ntype)            # store the efficacy probability to go for each disease type
  
  
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
        
        phat<-sum(y)/sum(n)
        obs<-cbind(y, n-y); E <- cbind(n * phat, n * (1 - phat))
        T <- sum((abs(obs - E))^2/E)
        if (is.nan(T) | (T<1)) {T<-1}
        test.stat[trial, i]<-T       # store the test statistic
        sigma2<-exp(a+b*log(T))
        if (sigma2==Inf) {sigma2<-10^4}
        tau<-1/sigma2
        jags.data <- list("y" = y, "n" = n, "ntype" = ntype, "tau" = tau, "mu.par" = mu.par)
        jags.fit <- jags.model(fcbhm, data = jags.data,
                               n.adapt = n.burnin, n.chains = n.chain, quiet = T)
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
        phat<-sum(y)/sum(n)
        obs<-cbind(y, n-y); E <- cbind(n * phat, n * (1 - phat))
        T <- sum((abs(obs - E))^2/E)
        if (is.nan(T) | (T<1)) {T<-1}
        test.stat[trial, i]<-T       # store the test statistic
        sigma2<-exp(a+b*log(T))
        if (sigma2==Inf) {sigma2<-10^4}
        tau<-1/sigma2
        jags.data<-list("y"=y,"n"=n,"ntype"=ntype,"tau"=tau,"mu.par"=mu.par)
        jags.fit<- jags.model(fcbhm,data=jags.data,
                              n.adapt = n.burnin, n.chains = n.chain, quiet = T)
        post.samples <- coda.samples(jags.fit, variable.names = c("p"), n.iter = n.iter, progress.bar = "none")
        pres.est <- do.call(rbind, post.samples)
        eff.prob<-sapply(seq(1,ntype,1), function(x) mean(pres.est[,x]>p.null[x]))
        eff.prob.store[trial,]<-eff.prob
        
        go <- cont[which((eff.prob > eff.ref)[cont])]  ##index of arm satisfies efficacy boundary
        nogo <- cont[which((eff.prob <= eff.ref)[cont])]
        efficacy[trial,go]<-1
        efficacy[trial,nogo]<-0
        
        p.est[trial,]<-apply(pres.est,2,mean)                    # store the rate estimate for each disease type      
        sample.size[trial,]<-n   
        
      }  
    }
    
  }
  
  
  results<-cbind(test.stat, sample.size,efficacy,eff.prob.store, p.est)
  return(results)
  
}







#' Get efficacy cutoff
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
#' @param FWER the pre-specified family wise error rate for efficacy cutoff calculation#
#' @param precision: the grid space controlled in the grid search, can be refined       #
# if higher precision is needed                                                 #
#' @param fwer.threshold: the acceptable threshold for the difference between           #
# estimated FWER and proposed FWER    
#' @param a the tuning parameter that characterize the relationship between test       #
# statistic T and the shrinkage parameter sigma^2   
#' @param b the tuning parameter that characterize the relationship between test       #
# statistic T and the shrinkage parameter sigma^2      
#' @return efficacy cutoff
#' @export
CBHMdesign.calibrate <-function(cohortsize,ntype = 4, p.true, p.null, p.target, ntrial=5000, n.burnin, n.chain, 
                                n.iter, mu.par,FWER, precision=0.001, fwer.threshold=0.001, a, b){
  
  cat('model
  {
    for (j in 1:ntype){                #J=4, the number of arms
      y[j]~ dbin(p[j],n[j])
      p[j]<-exp(theta[j])/(1+exp(theta[j]))
      
    }
    for (j in 1:ntype){
      theta[j] ~ dnorm(mu,tau)      #hierarchical model for theta
    }
    mu ~ dnorm(mu.par, 0.01)         #prior on mu
  }', file={fcbhm <- tempfile()})
  
  ncohort<-length(cohortsize)
  test.stat<-matrix(0, nrow=ntrial, ncol=ncohort)             # store the test statistic value at each interim analysis
  eff.prob.store<-matrix(0,nrow=ntrial,ncol=ntype)
  futstop<-0.05  
  
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
        
        phat<-sum(y)/sum(n)
        obs<-cbind(y, n-y); E <- cbind(n * phat, n * (1 - phat))
        T <- sum((abs(obs - E))^2/E)
        if (is.nan(T) | (T<1)) {T<-1}
        test.stat[trial, i]<-T       # store the test statistic
        sigma2<-exp(a+b*log(T))
        if (sigma2==Inf) {sigma2<-10^4}
        tau<-1/sigma2
        jags.data<-list("y"=y,"n"=n,"ntype"=ntype,"tau"=tau,"mu.par"=mu.par)
        jags.fit<- jags.model(fcbhm,data=jags.data,
                              n.adapt = n.burnin, n.chains = n.chain, quiet = T)
        post.samples <- coda.samples(jags.fit, variable.names = c("p"), n.iter = n.iter, progress.bar = "none")
        pres.est <- do.call(rbind, post.samples)                                               # extract the mcmc samples of p
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
        phat<-sum(y)/sum(n)
        obs<-cbind(y, n-y); E <- cbind(n * phat, n * (1 - phat))
        T <- sum((abs(obs - E))^2/E)
        if (is.nan(T) | (T<1)) {T<-1}
        test.stat[trial, i]<-T       # store the test statistic
        sigma2<-exp(a+b*log(T))
        if (sigma2==Inf) {sigma2<-10^4}
        tau<-1/sigma2
        jags.data<-list("y"=y,"n"=n,"ntype"=ntype,"tau"=tau,"mu.par"=mu.par)
        jags.fit<- jags.model(fcbhm,data=jags.data,
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
  Grid = seq(w0-0.02, FWER,by = precision) # one can modify this for different methods
  fwer.pred = sapply(Grid, function(o) fwer(eff.prob.store,o))
  #fwer.pred.left = fwer.pred[fwer.pred < FWER]
  #anchor.grid = which.min(abs(fwer.pred.left - FWER))
  anchor.grid = which.min(abs(fwer.pred-FWER))
  if(abs(fwer.pred[anchor.grid] - FWER) > fwer.threshold){
    print("Precision is not small enough or simulation number is not big enough, please reput again")
    eff.cut = sapply(1:ntype, FUN = function(x){quantile(eff.prob.store[,x], 1 - Grid[anchor.grid])})
    results = list(prob.success = eff.prob.store, eff.ref = eff.cut)
    #return(NULL)
    #return(results)
  }else{
    eff.cut = sapply(1:ntype, FUN = function(x){quantile(eff.prob.store[,x], 1 - Grid[anchor.grid])})
    results = list(prob.success = eff.prob.store, eff.ref = eff.cut)
    #return(results)
  }
  
  plot_search_type1_cbhm = ggplot2::ggplot(data=data.frame(Grid, fwer.pred), ggplot2::aes(x = Grid, y = fwer.pred, group=1)) +
    ggplot2::geom_line()+
    ggplot2::geom_point() + ggplot2::ggtitle("Family wise error rate with grid search")+
    ggplot2::labs(x="Fixed type one error for each cohort", y="FWER")
  
  return(list(results, plot_search_type1_cbhm))
  
}
