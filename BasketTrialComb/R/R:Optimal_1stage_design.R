#' This function finds alphastar
#' 
#' @import dplyr
#' @import poibin
#' @import gtools
#' @param K number of indications
#' @param g the prior probability of H0: that the strata are completely correlated
#' @param r number of responses
#' @param n number of patients
#' @param p0 probability vector of the null hypothesis
#' @param p1 probability vector of the alternative hypothesis
#' @param alphastar adjusted alpha
#' @param nsim number of simulation used for type I and type II calculation
#' @return power
#' @export
PowerEmpirical.1stg <- function(K, g, r, n, p0, p1, alphastar, nsim){
  # r is the final pooling criteria for individual tumor type
  accept = rej = 0
  set.seed(10)
  pval <- c()
  # ?
  probs <- p1*g+p0*(1-g)
  for (i in 1:nsim){
    # number of responses
    resp = rbinom(K,n,probs)
    #r <- qbinom(1-alpha0, n, p1) # alpha0 is the pruning bar;
    index = which (resp > r)
    resp_pool <- resp[index]
    #wt.p0 <- sum(p0[index]*n[index])/sum(n[index])
    # CDF of Poisson Binomial Distribution
    if(length(index)==0){pval[i]=1} else {pval[i] <- 1- poibin::ppoibin(sum(resp_pool)-1, pp=p0[index],wts=n[index])}
    if(pval[i]<alphastar) {rej = rej + 1}  # if pval less than alphastar, we claim effective drug
    else {accept = accept + 1}
  }    
  return(rej/nsim) # prob of claim effective at alphastar level, corresponding to global type I error under H0 and power at H1;
}


#' This function finds alphastar
#'
#' @param K number of indications
#' @param g the prior probability of H0: that the strata are completely correlated
#' @param alphastar adjusted alpha
#' @param r response probability for inactive drug
#' @param n number of patients
#' @param p probability vector, p[1] is the null hypothesis probability, p[2] is the alternative hypothesis probability
#' @param alpha  global type I error
#' @param nsim number of simulation used for type I and type II calculation
#' @return return alphastar
#' @export
diff0.emp.1stg <- function(K, g, alphastar, r, n, p, alpha, nsim){
  # calculate alphastar
  PowerEmpirical.1stg(K=K,g=g,r=r,n=n, p0=p, p1=p, alphastar, nsim=nsim)-alpha
}

#' Optim paramerters setting
#'
#' @param K number of indications
#' @param n number of patients
#' @param p0 probability vector of the null hypothesis
#' @param p1 probability vector of the alternative hypothesis
#' @param Alpha global type I error
#' @param nsim number of simulation used for type I and type II calculation
#' @return optim paramerters setting for optimal one stage design
#' @export
Optimal_pooling_hetero <- function (K, n, p0, p1, Alpha, nsim=100000){
  res <- NULL
  # grid search for alpha \talpha and \tbeta
  for (alpha0seq in seq(0.5, 0.01, by=-0.01)){
    # a sequence of r_k, using bionamial quantile to calculate r_k
    r.seq <- qbinom(1-alpha0seq, n, p0) + 1
    # \t beta
    beta0 <- pbinom(r.seq-1, n, p1)
    # \t alpha
    alpha0 <- 1-pbinom(r.seq-1,n,p0)
    #diff <- diff(typeIIerror[!duplicated(typeIIerror)])
    diff <- diff(range(alpha0))+diff(range(beta0))
    
    if(PowerEmpirical.1stg(K=K,g=0,r=r.seq,n=n, p0=p0, p1=p0, alphastar=1, nsim=nsim) < Alpha) break
    
    alphastar <- uniroot(diff0.emp.1stg,c(0,1),p=p0,K=K,g=0,r=r.seq, n=n, alpha=Alpha, nsim=nsim)$root
    
    typeIerror <-  round(PowerEmpirical.1stg(K=K,g=0,r=r.seq,n=n, p0=p0, p1=p0, alphastar=alphastar, nsim=nsim), 3)
    
    while(typeIerror > Alpha +0.001 ) { # if typeIerror still exceeds 0.05, then fine tune alphastar;
      # ?
      alphastar.new = alphastar - 0.001
      alphastar <- alphastar.new
      typeIerror <- round(PowerEmpirical.1stg(K=K,g=0,r=r.seq,n=n, p0=p0, p1=p0, alphastar=alphastar, nsim=nsim),3)
    }
    # all combinations of 0 and 1 in a 4 elements vector
    g.list <- permutations(2,K,0:1,repeats.allowed = TRUE)[-1,] # scenarios with at least one positive
    power.seq <- n.active <-  rep(0,nrow(g.list))
    for (j in 1:nrow(g.list)){
      # K-length vector of indicators
      g <- g.list[j,]
      # number of truly active
      n.active[j] <- sum(g)
      power.seq[j] <- PowerEmpirical.1stg(K=K,g=g,r=r.seq,n=n, p0=p0, p1=p1, alphastar=alphastar, nsim=nsim)
    } 
    temp.list <- cbind(n.active, power.seq)
    # ?
    power.unif <- mean(sapply(split(power.seq,n.active),mean))
    tmp1 <- c(n=n, r=r.seq, alpha0seq=alpha0seq, alpha0=round(alpha0,3), beta0=round(beta0,3), alphastar=round(alphastar,4), typeIerror=typeIerror,power=round(power.unif,4), diff=round(diff,3))
    res <- rbind(res, tmp1)
  }   
  # with max power
  opt_power <- res[res[,"power"]==max(res[,"power"]),]
  # with multiple max power
  if (is.null(nrow(opt_power))==FALSE) {opt_power <- opt_power[!duplicated(opt_power[,"power"]),]}
  # ?
  min_diff <- res[res[,"diff"]==min(res[,"diff"]),]
  if (nrow(min_diff)>1){  min_diff <- min_diff[!duplicated(min_diff[,"power"]),]}
  res <- res[!duplicated(res[,"power"]),]
  return(list(res=res, opt_power=opt_power, min_diff=min_diff))
}

#' main function decision making
#' 
#' @param n number of patients
#' @param k number of indications
#' @param r.1stg number of responses for interium 1st stage
#' @param alpha.adj.1stg adjusted alpha for 1st stage
#' @param resp.cntl probability under the null
#' @return return drug active ones, total number of subjects
#' @export
optimal.1stage <- function(n, k, r.1stg, alpha.adj.1stg,r,resp.cntl){
  fut.id.1stg <- which(r<r.1stg)
  pool <- c(1:k)[!(1:k) %in% fut.id.1stg]
  resp_pool <- sum(r[!(1:k) %in% fut.id.1stg])
  
  if(length(pool)==0){
    pval <- 1
    rej <- rep(0,k)
    
  } else {
    rej <- rep(0,k)
    pval <- 1- poibin::ppoibin(resp_pool-1, pp=resp.cntl[pool],wts=n[pool])
    if(pval<alpha.adj.1stg){
      rej[pool] <- 1
    }
    
  }
  n.tot <- n
  res <- c(rej,n.tot)
  
  return(  res )
}
