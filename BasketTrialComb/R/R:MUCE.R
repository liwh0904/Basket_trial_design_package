#' Get posterior samples from MUCE Method
#' 
#' @import rjags
#' @param k no. of indications
#' @param J no. of dose levels
#' @param Nk ? sample size at each analysis
#' @param rik ?
#' @param H0 largest response rate under null, of length(K)
#' @param gamma a parameter for prior dist of theta_ij
#' @param prec0 std of Z_ij, latent probit regression
#' @param prec.eta std of dose specific effects
#' @param prec.ksi std of indication specific effects 
#' @param mu.eta0 the mean of hyperprior eta_0
#' @param prec.eta0 the std of hyperprior eta_0 
#' @param mu.ksi0 the mean of hyperprior ksi_0 
#' @param prec.ksi0 the std of hyperprior ksi_0
#' @param n.burnin control prameters of MCMC sampling using jags
#' @param n.chain control prameters of MCMC sampling using jags
#' @param n.iter control prameters of MCMC sampling using jags
#' @return Probability to stop
#' @export
jags.muce.post.samples <- function(k, 
                                   J, 
                                   Nk, 
                                   rik, 
                                   H0, 
                                   gamma, 
                                   prec0, 
                                   prec.eta, 
                                   prec.ksi,
                                   mu.eta0,
                                   prec.eta0,
                                   mu.ksi0,
                                   prec.ksi0,
                                   n.burnin,
                                   n.chain,
                                   n.iter){
  jags.data <- list(
    "K" = k, "J" = J,
    "n" = matrix(colSums(Nk,na.rm=T), k, J), "r" = matrix(colSums(rik,na.rm=T), k, J),
    "H0" = H0,
    "gamma" = gamma, "prec0" = prec0,
    "prec.eta" = prec.eta, "prec.ksi" = prec.ksi,
    "mu.eta0" = mu.eta0, "prec.eta0" = prec.eta0,
    "mu.ksi0" = mu.ksi0, "prec.ksi0" = prec.ksi0
  )
  
  
  cat('model
  {
    tau <- 1/pow(gamma,2)
    
    for (i in 1:K)
    {
      for (j in 1:J)
      { 
        r[i,j] ~ dbin(p[i,j],n[i,j])
        p[i,j] <- exp(theta[i,j])/(1+exp(theta[i,j]))
        muz[i,j] <- ksi[i] + eta[j]
        Z[i,j] ~ dnorm(muz[i,j], prec0)
        #      a[i,j] ~ dt(hc[i], tau, 1)T(,hc[i])
        #      b[i,j] ~ dt(hc[i], tau, 1)T(hc[i],)
        a[i,j] ~ dnorm(hc[i], tau)T(,hc[i])
        b[i,j] ~ dnorm(hc[i], tau)T(hc[i],)
        theta[i,j] <- ifelse (Z[i,j]<=0, a[i,j], b[i,j])
        lambda[i,j] <- ifelse(Z[i,j]<=0, 0, 1)
        
      }
      
    }
    
    for (i in 1:K){
      hc[i] <- logit(H0[i])
      ksi[i] ~ dnorm(ksi0, prec.ksi)
      
    }
    
    for (j in 1:J){
      eta[j] ~ dnorm(eta0, prec.eta)
    }
    eta0 ~ dnorm(mu.eta0, prec.eta0)
    ksi0 ~ dnorm(mu.ksi0, prec.ksi0)
  }', file={fmuce <- tempfile()})
  
  jags.fit <- jags.model(fmuce, data = jags.data,
                         n.adapt = 5000, n.chains = n.chain, quiet = T)
  #jags.fit <- jags.model(file = "fmuce.txt", data = jags.data,
  #                       n.adapt = 5000, n.chains = n.chain, quiet = T)
  update(jags.fit,n.burnin)
  post.samples <- coda.samples(jags.fit, variable.names = c("lambda"), n.iter = n.iter/n.chain, progress.bar = "none")
  muce.out <- do.call(rbind, post.samples) ##lambda==1 means reject the null hypothesis
  ##posterior is the posterior estimate of rejection rate of null hypothesis.
  posterior <- colMeans(muce.out)  
  return(posterior)
}


#' MUCE() design Function
#'
#' @param k no. of indications
#' @param J no. of dose levels
#' @param L no. of interim analysis; 
#' @param nik sample size at each analysis 
#' @param H0 largest response rate under null, of length(K)
#' @param H1 lowest response rate under alternative, of length(K)
#' @param p0 true response rate for each arm
#' @param v1 phi_1, the futility boundary at interim 
#' @param v2 phi_2, the efficacy boundary at final
#' @param gamma a parameter for prior dist of theta_ij
#' @param prec0 std of Z_ij, latent probit regression
#' @param prec.eta std of dose specific effects
#' @param prec.ksi std of indication specific effects 
#' @param mu.eta0 the mean of hyperprior eta_0
#' @param prec.eta0 the std of hyperprior eta_0 
#' @param mu.ksi0 the mean of hyperprior ksi_0 
#' @param prec.ksi0 the std of hyperprior ksi_0
#' @param n.burnin control prameters of MCMC sampling using jags
#' @param n.chain control prameters of MCMC sampling using jags
#' @param n.iter control prameters of MCMC sampling using jags
#' @param r.all number of responses  
#' @param num.sim number of simulated trials needs to be 1 here, repeated simulations are done outside 
#' @return probability of reject
#' @export
muce <- function(K, 
                 J,
                 L,
                 nik,
                 H0, 
                 H1,
                 p0,
                 v1,
                 v2,
                 gamma, 
                 prec0, 
                 prec.eta, 
                 prec.ksi,
                 mu.eta0,
                 prec.eta0,
                 mu.ksi0,
                 prec.ksi0,
                 n.burnin,
                 n.chain,
                 n.iter,
                 r.all,
                 num.sim = 1){
  
  Nik <- rbind(nik[1,], apply(nik, 2, diff)) ##the no. of pts newly enrolled before each analysis. row is for analysis
  ##decision.muce indicates the analysis at which the arm stopped. '-interim' refers to stopping at 'interim' analysis. 1 means claim efficacy at final analysis, while 0 means not efficious at final analysis.
  decision.muce <- matrix(NA, num.sim, K*J)#each row is for each repetition. each col is for one analysis of one indication.
  N.SS <- matrix(NA, num.sim, K*J)
  tp <- which(p0 >= rep(H1,J))  ##index of true indications. should be K*J because H1 should be of length K.
  tn <- which(p0 < rep(H1,J))
  Decision <- matrix(NA, num.sim, K*J)  ##0 means tp, 2 means fn, 1 means tn, 3 means fp
  TP <- NULL  ##total no. of TP in each repetition
  FP <- NULL
  TN <- NULL
  FN <- NULL
  perfect <- NULL
  FWER <- NULL
  for (sim in 1:1){##loop for repetitions
    #set.seed(100+sim)  
    
    stop <- NULL
    cont <- 1:(K*J)  ##index of arms continued enrollment before analysis
    interim <- 1
    Nk <- matrix(NA, L+1, (K*J)) ##no. of pts enrolled before the specific analysis.
    rik <- matrix(NA, L+1, (K*J)) 
    
    while ((interim <= L) & (length(cont)>0)){##interim analysis and at least one arm continued enrollment before analysis
      ##Nk[interim, ] is the no. of pts newly enrolled before 'interim' analysis. if the arm has stopped, no. is 0.
      ##rik[interim, ] is the simulated no. of response for 'interim' analysis.
      Nk[interim, ] <- sapply(1:(K*J), FUN = function(x){ifelse(is.element(x, cont), Nik[interim, x], 0)})
      rik[interim, ] <- r.all[interim,]
      rik[interim, ][Nk[interim, ]==0] <- 0
      
      posterior <- jags.muce.post.samples(K,J,Nk,rik,H0,gamma,prec0,prec.eta,prec.ksi,mu.eta0,prec.eta0,mu.ksi0,prec.ksi0,n.burnin,n.chain,n.iter)
      stop_temp <- cont[which(posterior[cont] < v1)] ##index of arms stopped at 'interim' analysis, early stop for futility.
      stop <- c(stop, stop_temp) 
      cont <- cont[!cont%in%stop] ##only the continued arms remain
      decision.muce[sim, stop_temp] <- -interim  
      interim <- interim+1
    }
    if (length(cont)>0){##if there is at least one arm proceed to final analysis.
      Nk[L+1, ] <- sapply(1:(K*J), FUN = function(x){ifelse(is.element(x, cont), Nik[L+1, x], 0)})
      rik[L+1, ] <- r.all[interim,]
      rik[L+1, ][Nk[L+1, ]==0] <- 0
      
      posterior <-  jags.muce.post.samples(K,J,Nk,rik,H0,gamma,prec0,prec.eta,prec.ksi,mu.eta0,prec.eta0,mu.ksi0,prec.ksi0,n.burnin,n.chain,n.iter)
      go <- cont[which(posterior[cont] > v2[cont])]  ##index of arm satisfies efficacy boundary
      nogo <- cont[which(posterior[cont] <= v2[cont])]
      decision.muce[sim, go] <- 1
      decision.muce[sim, nogo] <- 0
      
    }
    N.SS[sim, ] <- colSums(Nk, na.rm = T) ##total sample size for each arm
    TP[sim] <- sum(decision.muce[sim,tp]==1)
    FP[sim] <- sum(decision.muce[sim,tn]==1)
    TN[sim] <- sum(decision.muce[sim,tn]<=0)
    FN[sim] <- sum(decision.muce[sim,tp]<=0)
    perfect[sim] <- (TP[sim]+TN[sim])==K*J   ##perfect==1 means all arms are correctly identified.
    FWER[sim] <- FP[sim]>0  ##if there is at least one FP, then this simulation will be counted in FWER calculation.
    Decision[sim,tp] <- 0*(decision.muce[sim,tp]==1)+2*(decision.muce[sim,tp]<=0) 
    Decision[sim,tn] <- 1*(decision.muce[sim,tn]<=0)+3*(decision.muce[sim,tn]==1)
  }
  N.SS <- colMeans(N.SS) ##average sample size for each arm.
  out <- list(decision.muce=decision.muce, N.SS=N.SS,
              TP=mean(TP), FP=mean(FP),TN=mean(TN), FN=mean(FN), perfect=mean(perfect), FWER=mean(FWER), Decision=Decision)
  return(out)
}


#' Get efficacy cutoff
#' 
#' @param k no. of indications
#' @param J no. of dose levels
#' @param L no. of interim analysis
#' @param num.sim  no. of repetitions 
#' @param FWER the pre-specified family wise error rate for efficacy cutoff calculatio
#' @param precisio the grid space controlled in the grid search, can be refined if higher precision is needed 
#' @param fwer.threshold the acceptable thresholf for the difference between estimated FWER and proposed FWER
#' @param nik sample size at each analysis 
#' @param H0 largest response rate under null, of length(K)
#' @param H1 lowest response rate under alternative, of length(K)
#' @param p0 true response rate for each arm
#' @param v1 phi_1, the futility boundary at interim 
#' @param gamma a parameter for prior dist of theta_ij
#' @param prec0 std of Z_ij, latent probit regression
#' @param prec.eta std of dose specific effects
#' @param prec.ksi std of indication specific effects 
#' @param mu.eta0 the mean of hyperprior eta_0
#' @param prec.eta0 the std of hyperprior eta_0 
#' @param mu.ksi0 the mean of hyperprior ksi_0 
#' @param prec.ksi0 the std of hyperprior ksi_0
#' @param n.burnin control prameters of MCMC sampling using jags
#' @param n.chain control prameters of MCMC sampling using jags
#' @param n.iter control prameters of MCMC sampling using jags
#' @return efficacy cutoff
#' @export
muce.calibrate <- function(K, 
                           J, 
                           L,
                           num.sim,
                           FWER,
                           precision, 
                           fwer.threshold,
                           nik,
                           H0, 
                           H1,
                           p0,
                           v1,
                           gamma, 
                           prec0, 
                           prec.eta, 
                           prec.ksi,
                           mu.eta0,
                           prec.eta0,
                           mu.ksi0,
                           prec.ksi0,
                           n.burnin,
                           n.chain,
                           n.iter){
  Nik <- rbind(nik[1,], apply(nik, 2, diff)) ##the no. of pts newly enrolled before each analysis. row is for analysis
  posterior.muce <- matrix(NA, num.sim, K*J)
  
  for (sim in 1:num.sim){
    set.seed(10100+sim) 
    
    stop <- NULL
    cont <- 1:(K*J)
    interim <- 1
    Nk <- matrix(NA, L+1, (K*J)) 
    rik <- matrix(NA, L+1, (K*J)) 
    
    while ((interim <= L) & (length(cont)>0)){
      
      Nk[interim, ] <- sapply(1:(K*J), FUN = function(x){ifelse(is.element(x, cont), Nik[interim, x], 0)})
      rik[interim, ] <- sapply(1:(K*J), FUN = function(x){rbinom(n = 1, size = Nk[interim, x], prob = p0[x])})
      
      posterior <- jags.muce.post.samples(K,J,  Nk, rik, H0,gamma,prec0,prec.eta,prec.ksi,mu.eta0,
                                          prec.eta0,mu.ksi0,prec.ksi0,n.burnin,n.chain,n.iter)
      stop_temp <- cont[which(posterior[cont] < v1)]
      stop <- c(stop, stop_temp)
      cont <- cont[!cont%in%stop]
      posterior.muce[sim, stop_temp] <- 0
      interim <- interim+1
    }
    if (length(cont)>0){
      Nk[L+1, ] <- sapply(1:(K*J), FUN = function(x){ifelse(is.element(x, cont), Nik[L+1, x], 0)})
      rik[L+1, ] <- sapply(1:(K*J), FUN = function(x){rbinom(n = 1, size = Nk[L+1, x], prob = p0[x])})
      
      posterior <-  jags.muce.post.samples(K, J, Nk, rik, H0, gamma, prec0, prec.eta, prec.ksi, mu.eta0,
                                           prec.eta0,mu.ksi0,prec.ksi0,n.burnin,n.chain,n.iter)
      posterior.muce[sim, cont] <- posterior[cont]
      
    }
  }
  
  ###grid search for efficacy cutoff
  w0 = 1-(1-FWER)^(1/4)
  Grid = seq(w0-0.01, FWER,by = precision) # one can modify this for different methods
  fwer.pred = sapply(Grid, function(o) fwer(posterior.muce,o))
  fwer.pred.left = fwer.pred[fwer.pred < FWER]
  anchor.grid = which.min(abs(fwer.pred.left - FWER))
  #anchor.grid = which.min(abs(fwer.pred-FWER))
  if(abs(fwer.pred.left[anchor.grid]-FWER)>fwer.threshold){
    print("Precision is not small enough or simulation number is not big enough, please reput again")
    eff.cut = sapply(1:K, FUN = function(x){quantile(posterior.muce[,x], 1-Grid[anchor.grid])})
    #return(NULL)
    #return(eff.cut)
  }else{
    eff.cut = sapply(1:K, FUN = function(x){quantile(posterior.muce[,x], 1-Grid[anchor.grid])})
    #return(eff.cut)
  }
  plot_search_type1_muce = ggplot2::ggplot(data=data.frame(Grid, fwer.pred), ggplot2::aes(x = Grid, y = fwer.pred, group=1)) +
    ggplot2::geom_line()+
    ggplot2::geom_point() + ggplot2::ggtitle("Family wise error rate with grid search")+
    ggplot2::labs(x="Fixed type one error for each cohort", y="FWER")
  
  return(list(eff.cut, plot_search_type1_muce))
  
}

