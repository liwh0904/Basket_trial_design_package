#' Get posterior samples from MUCE Method
#' 
#' @import rjags
#' @param nk sample size for each cohort at each analysis (can be different). row is for analysis#
#' @param rk ?
#' @param k no. of cohorts/indications
#' @param Nexch the no. of exchangeability distributions  
#' @param Nmix  ==Nexch+1. no. of total candidate prior dists for theta_j. 
#' @param pMix  the probs of each strata following each dist. 
#' @param dirichprior parameter for prior distribution specified in exnex.dirich.txt. but we do not have this file.
#' @param dirichprior.check if TRUE, will use exnex.dirich.txt to specify bayesian model
#' @param mu.mean  mean of prior dist of mu in exchangeable part 
#' @param mu.prec std of prior dist of mu in exchangeable part
#' @param tau.HN.scale parameter to derive prior dist of 
#' @param nex.mean to derive prior dist of theta_j under non-exchangeability
#' @param nex.prec   to derive prior dist of theta_j under non-exchangeability
#' @param n.burnin control prameters of MCMC sampling using jags
#' @param n.chain control prameters of MCMC sampling using jags
#' @param n.iter control prameters of MCMC sampling using jags
#' @return Probability to stop
#' @export
jags.exnex.post.samples <- function(nk = rep(24,6), 
                                    rk = rep(4,6), 
                                    k = 6, 
                                    Nexch = 2,
                                    Nmix = 3,
                                    pMix = c(0.5,0,0.5),
                                    dirichprior = c(1,1,1),
                                    dirichprior.check = TRUE,
                                    mu.mean = c(-1.735,0.847),
                                    mu.prec = c(0.146,0.266),
                                    tau.HN.scale = c(1,1),
                                    nex.mean = rep(-1.734,6),
                                    nex.prec = rep(0.128,6),
                                    n.burnin = 1000, 
                                    n.chain = 1, 
                                    n.iter = 8000){

  cat('model{
    # prior distributions for EX-parameters
    for (jj in 1:Nexch) {
      mu[jj] ~dnorm(mu.mean[jj],mu.prec[jj])
      prior.tau.prec[jj] <- pow(tau.HN.scale[jj],-2)
      tau[jj] ~ dnorm(0,prior.tau.prec[jj])I(0.001,)
      prec.tau[jj] <- pow(tau[jj],-2)
    }
    # log-odds parameters under EX
    for (jj in 1:Nexch) {
      for (j in 1:Nstrata) {
        re[jj,j] ~ dnorm(0,prec.tau[jj])
        LogOdds[jj,j] <- mu[jj]+re[jj,j]
      }
    }
    
    # log-odds parameters under NEX
    for (j in 1:Nstrata) {	
      LogOdds[Nmix,j] ~ dnorm(nex.mean[j],nex.prec[j])
    }
    
    # latent mixture indicators:
    # exch.index: categorial 1,...,Nmix=Nexch+1
    # exch: Nstrata x Nmix matrix of 0/1 elements
    
    for (j in 1:Nstrata) {
      exch.index[j] ~ dcat(pMix[1:Nmix])
      for (jj in 1:Nmix) {
        exch[j,jj] <- equals(exch.index[j],jj)
      }
    }
    # pick theta
    for (j in 1:Nstrata) {
      theta[j] <- LogOdds[exch.index[j],j]
    }
    # likelihood part
    for (i in 1:Nstrata) {
      logit( p[i] ) <- theta[i]
      r[i] ~ dbin(p[i],n[i])
    }
  }', file={fexnex <- tempfile()})
  
  if (dirichprior.check == TRUE){
    jags.data <- list(
      "Nexch" = Nexch, "Nmix" = Nmix, "dirichprior" = dirichprior,
      "Nstrata" = k, "n" = nk, "r" = rk,
      # prior means and precisions for EX parameter mu
      "mu.mean" = mu.mean, "mu.prec" = mu.prec,
      # scale parameter of Half-Normal prior for tau
      "tau.HN.scale" = tau.HN.scale,
      # NEX priors; make them strata-specific if needed
      "nex.mean" = nex.mean, "nex.prec" = nex.prec
    )
    
    jags.fit <- jags.model(fexnex, data = jags.data,
                           n.adapt = n.burnin, n.chains = n.chain, quiet = T)
    #jags.fit <- jags.model(file = "./src/functions/exnex.dirich.txt", data = jags.data,
    #                       n.adapt = n.burnin, n.chains = n.chain, quiet = T)
    ##note: we do not have exnex.dirich.txt
  }else if (dirichprior.check == FALSE){
    jags.data <- list(
      "Nexch" = Nexch, "Nmix" = Nmix, "pMix" = pMix,
      "Nstrata" = k, "n" = nk, "r" = rk,
      # prior means and precisions for EX parameter mu
      "mu.mean" = mu.mean, "mu.prec" = mu.prec,
      # scale parameter of Half-Normal prior for tau
      "tau.HN.scale" = tau.HN.scale,
      # NEX priors; make them strata-specific if needed
      "nex.mean" = nex.mean, "nex.prec" = nex.prec
    )
    
    jags.fit <- jags.model(fexnex, data = jags.data,
                           n.adapt = n.burnin, n.chains = n.chain, quiet = T)
    #jags.fit <- jags.model(file = "src_functions_exnex.nondirich.txt", data = jags.data,
    #                      n.adapt = n.burnin, n.chains = n.chain, quiet = T)
  }
  ##file = "./src/functions/exnex.nondirich.txt"
  post.samples <- coda.samples(jags.fit, variable.names = c("p"), n.iter = n.iter/n.chain, progress.bar = "none")
  exnex.out <- do.call(rbind, post.samples)
  
  return(exnex.out) ##the output is the posterior samples of response rates
  
}


#' Main function for EXNEX Method
#' @param num.sim number of simulated trials needs to be 1 here, repeated simulations are done outside 
#' @param Q.exnex efficacy cutoff for each cohort  
#' @param p0 the vector of true response rate in each disease type      
#' @param H0 the vector of null hypothesis rate in each type 
#' @param H1 the vector of alternative hypothesis rate in each type   
#' @param nk sample size for each cohort at each analysis (can be different). row is for analysis#
#' @param k no. of cohorts/indications
#' @param v1 efficacy cutoff
#' @param Nexch the no. of exchangeability distributions  
#' @param Nmix  ==Nexch+1. no. of total candidate prior dists for theta_j. 
#' @param pMix  the probs of each strata following each dist. 
#' @param dirichprior parameter for prior distribution specified in exnex.dirich.txt. but we do not have this file.
#' @param dirichprior.check if TRUE, will use exnex.dirich.txt to specify bayesian model
#' @param mu.mean  mean of prior dist of mu in exchangeable part 
#' @param mu.prec std of prior dist of mu in exchangeable part
#' @param tau.HN.scale parameter to derive prior dist of 
#' @param nex.mean to derive prior dist of theta_j under non-exchangeability
#' @param nex.prec   to derive prior dist of theta_j under non-exchangeability
#' @param n.burnin control prameters of MCMC sampling using jags
#' @param n.chain control prameters of MCMC sampling using jags
#' @param n.iter control prameters of MCMC sampling using jags
#' @param r.all number of responses 
#' @return probability of reject
#' @export
exnex <- function(num.sim = 1,
                  Q.exnex = rep(0.8,6),
                  p0 = rep(0.2,6),
                  H0 = rep(0.2,6),
                  H1 = rep(0.4,6),
                  nk = matrix(rep(c(12,24),6),nrow=2), 
                  K = 6, 
                  v1, 
                  Nexch = 2,
                  Nmix = 3,
                  pMix = c(0.5,0,0.5),
                  dirichprior = c(1,1,1),
                  dirichprior.check = TRUE,
                  mu.mean = c(-1.735,0.847),
                  mu.prec = c(0.146,0.266),
                  tau.HN.scale = c(1,1),
                  nex.mean = rep(-1.734,6),
                  nex.prec = rep(0.128,6),
                  n.burnin = 1000, 
                  n.chain = 1, 
                  n.iter = 8000,
                  r.all){
  L <- dim(nk)[1]-1 #no. of interim analysis
  Nik <- rbind(nk[1,], apply(nk, 2, diff)) ##the no. of pts newly enrolled before each analysis. row is for analysis
  
  tp <- which(p0 >= H1)  ##index for TP
  tn <- which(p0 < H1)
  Decision <- matrix(NA, num.sim, K)
  TP <- NULL
  FP <- NULL
  TN <- NULL
  FN <- NULL
  perfect <- NULL
  FWER <- NULL
  decision.exnex <- matrix(NA, num.sim, K)
  N.SS <- matrix(NA, num.sim, K) 
  
  for (sim in 1:num.sim){
    #set.seed(100+sim)
    
    stop <- NULL
    cont <- 1:K  ##index of arms continued enrollment before analysis
    interim <- 1
    Nk <- matrix(NA, L+1, K) ##no. of pts newly enrolled before the specific analysis.
    rik <- matrix(NA, L+1, K) 
    
    while ((interim <= L) & (length(cont)>0)){##interim analysis and at least one arm continued enrollment before analysis
      ##Nk[interim, ] is the no. of pts newly enrolled before 'interim' analysis. if the arm has stopped, no. is 0.
      ##rik[interim, ] is the simulated no. of response for 'interim' analysis.
      Nk[interim, ] <- sapply(1:K, FUN = function(x){ifelse(is.element(x, cont), Nik[interim, x], 0)})
      #rik[interim, ] <- sapply(1:K, FUN = function(x){rbinom(n = 1, size = Nk[interim, x], prob = p0[x])})
      rik[interim, ] <- r.all[interim,]
      rik[interim, ][Nk[interim, ]==0] <- 0
      
      post.sample <- jags.exnex.post.samples(colSums(Nk,na.rm=T),colSums(rik,na.rm=T),K,Nexch,Nmix,pMix,dirichprior, dirichprior.check, mu.mean,mu.prec,tau.HN.scale,nex.mean, nex.prec, n.burnin, n.chain, n.iter)
      posterior <- rowMeans(t(post.sample)>H0)  ##the prop for each indication's response rate exceeding null hypothesis 
      
      stop_temp <- cont[which(posterior[cont] < v1)] ##index of arms stopped at 'interim' analysis, early stop for futility.
      stop <- c(stop, stop_temp) 
      cont <- cont[!cont%in%stop] ##only the continued arms remain
      decision.exnex[sim, stop_temp] <- -interim  
      interim <- interim+1
    }	
    
    if (length(cont)>0){##if there is at least one arm proceed to final analysis.
      Nk[L+1, ] <- sapply(1:K, FUN = function(x){ifelse(is.element(x, cont), Nik[L+1, x], 0)})
      #rik[L+1, ] <- sapply(1:K, FUN = function(x){rbinom(n = 1, size = Nk[L+1, x], prob = p0[x])})
      rik[L+1, ] <- r.all[interim,]
      rik[L+1, ][Nk[L+1, ]==0] <- 0
      
      post.sample <- jags.exnex.post.samples(colSums(Nk,na.rm=T),colSums(rik,na.rm=T),K,Nexch,Nmix,pMix,dirichprior, dirichprior.check, mu.mean,mu.prec,tau.HN.scale,nex.mean, nex.prec, n.burnin, n.chain, n.iter)
      posterior <- rowMeans(t(post.sample)>H0)  ##the prop for each indication's response rate exceeding null hypothesis 
      go <- cont[which((posterior > Q.exnex)[cont])]  ##index of arm satisfies efficacy boundary
      nogo <- cont[which((posterior <= Q.exnex)[cont])]
      decision.exnex[sim, go] <- 1  ##==1 calim efficacy
      decision.exnex[sim, nogo] <- 0       
    }
    
    N.SS[sim, ] <- colSums(Nk, na.rm = T) ##total sample size for each arm
    TP[sim] <- sum(decision.exnex[sim,tp]==1)
    FP[sim] <- sum(decision.exnex[sim,tn]==1)
    TN[sim] <- sum(decision.exnex[sim,tn]<=0)
    FN[sim] <- sum(decision.exnex[sim,tp]<=0)
    perfect[sim] <- (TP[sim]+TN[sim])==K   ##perfect==1 means all arms are correctly identified.
    FWER[sim] <- FP[sim]>0  ##if there is at least one FP, then this simulation will be counted in FWER calculation.
    Decision[sim,tp] <- 0*(decision.exnex[sim,tp]==1)+2*(decision.exnex[sim,tp]<=0) 
    Decision[sim,tn] <- 1*(decision.exnex[sim,tn]<=0)+3*(decision.exnex[sim,tn]==1)
  }
  
  N.SS <- colMeans(N.SS) ##average sample size for each arm.  
  out <- list(decision.exnex = decision.exnex, N.SS=N.SS,TP=mean(TP), FP=mean(FP),TN=mean(TN), FN=mean(FN), perfect=mean(perfect), FWER=mean(FWER), Decision=Decision)
  
  return(out)  
}


#' Get efficacy cutoff
#' 
#' @param k no. of cohorts/indications
#' @param v1 efficacy cutoff
#' @param H0 the vector of null hypothesis rate in each type   
#' @param p0 the vector of true response rate in each disease type  
#' @param FWER: the pre-specified family wise error rate for efficacy cutoff calculation
#' @param precision the grid space controlled in the grid search, can be refined       #
# if higher precision is needed                                                 #
#' @param fwer.threshold the acceptable thresholf for the difference between           #
# estimated FWER and proposed FWER  
#' @param nk sample size for each cohort at each analysis (can be different). row is for analysis#
#' @param num.sim number of MC simulations
#' @param Nexch the no. of exchangeability distributions  
#' @param Nmix  ==Nexch+1. no. of total candidate prior dists for theta_j. 
#' @param pMix  the probs of each strata following each dist. 
#' @param dirichprior parameter for prior distribution specified in exnex.dirich.txt. but we do not have this file.
#' @param dirichprior.check if TRUE, will use exnex.dirich.txt to specify bayesian model
#' @param mu.mean  mean of prior dist of mu in exchangeable part 
#' @param mu.prec std of prior dist of mu in exchangeable part
#' @param tau.HN.scale parameter to derive prior dist of 
#' @param nex.mean to derive prior dist of theta_j under non-exchangeability
#' @param nex.prec   to derive prior dist of theta_j under non-exchangeability
#' @param n.burnin control prameters of MCMC sampling using jags
#' @param n.chain control prameters of MCMC sampling using jags
#' @param n.iter control prameters of MCMC sampling using jags
#' @return parameters
#' @export
exnex.calibration <- function(K = 6,v1,
                              H0 = rep(0.2,6),p0,
                              FWER = 0.15,
                              precision = 0.001, 
                              fwer.threshold = 0.001,
                              nk = matrix(rep(c(12,24),6),nrow=2),
                              num.sim = 1000, 
                              Nexch = 2,
                              Nmix = 3,
                              pMix = c(0.5,0,0.5),
                              dirichprior = c(1,1,1),
                              dirichprior.check = TRUE,
                              mu.mean = c(-1.735,0.847),
                              mu.prec = c(0.146,0.266),
                              tau.HN.scale = c(1,1),
                              nex.mean = rep(-1.734,6),
                              nex.prec = rep(0.128,6),
                              n.burnin = 500, 
                              n.chain = 1, 
                              n.iter = 2000){
  
  L <- dim(nk)[1]-1 #no. of interim analysis
  Nik <- rbind(nk[1,], apply(nk, 2, diff)) ##the no. of pts newly enrolled before each analysis. row is for analysis
  posterior.exnex <- matrix(NA, num.sim, K)
  
  for (sim in 1:num.sim){
    set.seed(10100+sim) 
    
    stop <- NULL
    cont <- 1:K
    interim <- 1
    Nk <- matrix(NA, L+1, K) 
    rik <- matrix(NA, L+1, K) 
    
    while ((interim <= L) & (length(cont)>0)){##interim analysis and at least one arm continued enrollment before analysis
      ##Nk[interim, ] is the no. of pts newly enrolled before 'interim' analysis. if the arm has stopped, no. is 0.
      ##rik[interim, ] is the simulated no. of response for 'interim' analysis.
      Nk[interim, ] <- sapply(1:K, FUN = function(x){ifelse(is.element(x, cont), Nik[interim, x], 0)})
      rik[interim, ] <- sapply(1:K, FUN = function(x){rbinom(n = 1, size = Nk[interim, x], prob = p0[x])})
      
      post.sample <- jags.exnex.post.samples(colSums(Nk,na.rm=T),colSums(rik,na.rm=T),K,Nexch,Nmix,pMix,dirichprior, dirichprior.check, mu.mean,mu.prec,tau.HN.scale,nex.mean, nex.prec, n.burnin, n.chain, n.iter)
      posterior <- rowMeans(t(post.sample)>H0)  ##the prop for each indication's response rate exceeding null hypothesis 
      stop_temp <- cont[which(posterior[cont] < v1)] ##index of arms stopped at 'interim' analysis, early stop for futility.
      stop <- c(stop, stop_temp) 
      cont <- cont[!cont%in%stop] ##only the continued arms remain
      posterior.exnex[sim, stop_temp] <- 0  
      interim <- interim+1
    }	
    
    if (length(cont)>0){##if there is at least one arm proceed to final analysis.
      Nk[L+1, ] <- sapply(1:K, FUN = function(x){ifelse(is.element(x, cont), Nik[L+1, x], 0)})
      rik[L+1, ] <- sapply(1:K, FUN = function(x){rbinom(n = 1, size = Nk[L+1, x], prob = p0[x])})
      
      post.sample <- jags.exnex.post.samples(colSums(Nk,na.rm=T),colSums(rik,na.rm=T),K,Nexch,Nmix,pMix,dirichprior, dirichprior.check, mu.mean,mu.prec,tau.HN.scale,nex.mean, nex.prec, n.burnin, n.chain, n.iter)
      posterior <- rowMeans(t(post.sample)>H0)  ##the prop for each indication's response rate exceeding null hypothesis 
      posterior.exnex[sim, cont] <- posterior[cont]       
    }
  }
  
  ###grid search for efficacy cutoff
  w0 = 1 - (1 - FWER)^(1/4)
  Grid = seq(w0 - 0.01, FWER, by = precision) # one can modify this for different methods
  fwer.pred = sapply(Grid, function(o) fwer(posterior.exnex, o))
  fwer.pred.left = fwer.pred[fwer.pred < FWER]
  anchor.grid = which.min(abs(fwer.pred.left - FWER))
  #anchor.grid = which.min(abs(fwer.pred - FWER))
  if(abs(fwer.pred.left[anchor.grid] - FWER) > fwer.threshold){
    print("Precision is not small enough or simulation number is not big enough, please reput again")
    eff.cut = sapply(1:K, FUN = function(x){quantile(posterior.exnex[,x], 1 - Grid[anchor.grid])})
    #return(NULL)
    #return(eff.cut)
  }else{
    eff.cut = sapply(1:K, FUN = function(x){quantile(posterior.exnex[,x], 1 - Grid[anchor.grid])})
    #return(eff.cut)
  }
  plot_search_type1_exnex = ggplot2::ggplot(data=data.frame(Grid, fwer.pred), ggplot2::aes(x = Grid, y = fwer.pred, group=1)) +
    ggplot2::geom_line()+
    ggplot2::geom_point() + ggplot2::ggtitle("Family wise error rate with grid search")+
    ggplot2::labs(x="Fixed type one error for each cohort", y="FWER")
  
  return(list(eff.cut, plot_search_type1_exnex))
}

