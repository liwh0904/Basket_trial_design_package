#'
#' This function find the post prob
#'
#' @param lambda the prior probability of H0: that the strata are completely correlated
#' @param gamma the prior probability that the drug is active in any specific stratum
#' @param r number of responses
#' @param n number of patients
#' @param plo response probability for inactive drug
#' @param phi response probability for active drug
#' @return post prob of activity in each strata
#' @export
find_postk <- function(lambda, gamma, r, n, plo, phi){
  ## P0: probability of drug is active in all strata assuming strata are homogeneous
  P0.num <- (1-gamma)*prod(dbinom(r, n, plo))
  P0.deno <- gamma*prod(dbinom(r, n, phi))
  P0 <- 1/(1 + (P0.num/P0.deno))
  
  ## Pk_ind: probability of drug is active in all strata assuming strata are independent
  Pk.num <- gamma*dbinom(r, n, phi)
  # ?? 1-gamma typo in the paper
  Pk.deno <- gamma*dbinom(r, n, phi) + (1-gamma)*dbinom(r, n, plo)
  Pk <- Pk.num/Pk.deno
  
  ## P[pi = 1|data]: the posterior probability for H_0(posterior probaability treatment effects are equal among strata)
  pi1.num <- (1-lambda)*prod(gamma*dbinom(r, n, phi) + (1-gamma)*dbinom(r, n, plo))
  pi1.deno <- lambda*(gamma*prod(dbinom(r, n, phi)) + (1-gamma)*prod(dbinom(r, n, plo)))
  pi1 <- 1/(1+(pi1.num/pi1.deno))
  
  ## Finally, post prob of activity in each strata
  postk <- P0*pi1 + Pk*(1 - pi1)
  postk
}


#' main function decision making
#'
#' @param r.ia number of responses in interium analysis
#' @param r  number of responses
#' @param lambda the prior probability of H0: that the strata are completely correlated
#' @param gamma the prior probability that the drug is active in any specific stratum
#' @param plo response probability for inactive drug
#' @param phi response probability for active drug
#' @param T threshold for conclusive posterior probability
#' @param n.ia  number of patients for interium analysis
#' @param k number of cohorts
#' @param n Path to the input file
#' @return return drug active ones, total number of subjects
#' @export
SimonBayes.2stg <- function(r.ia,r,lambda, gamma, plo, phi,T=0.8,n.ia,k,n){
  
  #?ia: interim analysis
  pos.ia <- find_postk(lambda, gamma, r.ia, n.ia, plo, phi)
  # declare drug inactive in stratum k
  fut.ia.id <- which(pos.ia<(1-T))
  # declare drug active in stratum k
  rej.ia.id <- which(pos.ia>T)
  # terminate accrual to stratum k for both cases
  stop.ia.id <- c(fut.ia.id,rej.ia.id)
  if(length(stop.ia.id)==0){
    # strat are still open
    go.ia.id <- 1:k
  } else {
    go.ia.id <- c(1:k)[-stop.ia.id]
  }
  
  # K cohorts
  if(length(fut.ia.id)==k){
    rej <- rep(0,4)
    n.tot <- n.ia
  } else if(length(rej.ia.id)==k){
    rej <- rep(1,4)
    n.tot <- n.ia
  } else if(length(go.ia.id)==0){
    
    rej <- rep(NA,k)
    rej[fut.ia.id] <- 0
    rej[rej.ia.id] <- 1
    
    # n.ia?: number of subjects in interium analysis
    n.tot <- n.ia
    
  } else {
    r.go <- r[go.ia.id]
    n.go <- n[go.ia.id]
    # ?fa: fa means final analysis
    pos.fa <- find_postk(lambda, gamma, r.go, n.go, plo, phi)
    
    ### If Pk>T declare drug active in stratum k
    ## T=threshold for conclusive posterior probability
    rej.fa.id <- go.ia.id[pos.fa>T]
    #?1-T
    fut.fa.id <- go.ia.id[pos.fa<=T]
    
    fut.id <- c(fut.ia.id,fut.fa.id)
    rej.id <- c(rej.ia.id,rej.fa.id)
    
    n.tot <- n
    #?: number of patients not go
    n.tot[-go.ia.id] <- n.ia[-go.ia.id]
    if(length(fut.id)==k){
      rej <- rep(0,k)
      
    } else if(length(rej.id)==k){
      #?
      rej <- rep(1,k)
      
    } else {
      rej <- rep(1,k)
      rej[fut.id] <- 0
    }
    
  }
  ## return number of drug active, total number of subjects
  res <- c(rej,n.tot)
  return(res)
}

#' Simulate data and run main function decision making
#'
#' @param rep number of repetition
#' @param k number of cohorts
#' @param n number of patients
#' @param n.ia number of patients for interium analysis
#' @param resp.true probability under the alternative
#' @param resp.cntl probability under the null
#' @param resp.trt probability under the alternative
#' @param lambda the prior probability of H0: that the strata are completely correlated
#' @param gamma the prior probability that the drug is active in any specific stratum
#' @param plo response probability for inactive drug
#' @param phi response probability for active drug
#' @param T threshold for conclusive posterior probability
#' @return return drug active ones, total number of subjects
#' @export
SimonBayes.type1 <- function(rep, k,n,n.ia,resp.true,resp.cntl,resp.trt, lambda,gamma,plo, phi,T){
  
  ## Simulate data
  dat <- data.frame(id=NA,cohort=NA,orr=NA)
  for(i in 1:k){
    cohort <- c(rep(i,n[i]))
    dat.tmp <- data.frame(cohort)
    id <- 1:n[i]
    # Simulate ORR
    # ? overall response rate
    orr <- rbinom(n=n[i], size=1, prob=resp.true[i])
    
    
    dat.tmp <- data.frame(dat.tmp,id,orr)
    dat <- rbind(dat,dat.tmp)
    
  }
  
  dat <- dat[-1,]
  
  r.ia <- rep(NA,k)
  r <- rep(NA,k)
  for(i in 1: k){
    dat.tmp <- subset(dat,cohort==i)
    #?n.ia
    r.ia[i] <- nrow(subset(dat.tmp[1:n.ia[i],], orr==1))
    r[i] <- nrow(subset(dat.tmp, orr==1))
    
  }
  
  ## Apply methods
  # (1) Simon's Bayesian
  Simon.Bayes.res <- SimonBayes.2stg(r.ia,r,lambda, gamma, plo, phi,T,n.ia,k,n)
  
  #res <- data.frame(SimBayes=Simon.Bayes.res,Opt2Stg=opt.2stg.res,Opt1Stg=opt.1stg.res,Simon2Stg=simon.2stg.res,Pool= simple.pool.res  )
  return(Simon.Bayes.res)
  
}

#' search for typeI error under differen T
#'
#' @param k number of cohorts
#' @param n number of patients
#' @param n.ia number of patients for interium analysis
#' @param resp.true probability under the alternative
#' @param resp.cntl probability under the null
#' @param resp.trt probability under the alternative
#' @param lambda the prior probability of H0: that the strata are completely correlated
#' @param gamma the prior probability that the drug is active in any specific stratum
#' @param plo response probability for inactive drug
#' @param phi response probability for active drug
#' @param T threshold for conclusive posterior probability
#' @return return typeI error under differen T
#' @export
SimonBayes.type1.search <- function(k,n,n.ia,resp.true,resp.cntl,resp.trt, lambda,gamma,plo, phi,T){
  typ1.lst <-  rep(NA,length(T.lst))
  for(j in 1:length(T.lst)){
    T <- T.lst[j]
    res <- sapply(1:100,SimonBayes.type1,k,n,n.ia,resp.true,resp.cntl,resp.trt, lambda,gamma,plo, phi,T)
    typ1 <- res[1:k,]
    typ1.lst[j] <- mean(colSums( typ1)>0) 
  }
  
  data.frame(T.lst,typ1.lst)
}
