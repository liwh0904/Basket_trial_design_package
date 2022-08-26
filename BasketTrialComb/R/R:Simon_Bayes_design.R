#'
#' This function find the post prob
#'
#' @import ggplot2
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
SimonBayes.2stg <- function(r.ia = NULL, r, lambda, gamma, plo, phi, T=0.8, n.ia = NULL, k, n){
  
  #ia: interim analysis
  if(is.null(n.ia)){
    fut.ia.id <- NULL
    rej.ia.id <- NULL
  } else {
    pos.ia <- find_postk(lambda, gamma, r.ia, n.ia, plo, phi)
    fut.ia.id <- which(pos.ia<(1-T))
    rej.ia.id <- which(pos.ia>T)
  }
  
  if(length(fut.ia.id)==0){
    go.ia.id <- 1:k
  } else {
    go.ia.id <- c(1:k)[-fut.ia.id]
  }
  
  if(length(go.ia.id)==0){
    rej <- rep(0,4)
    n.tot <- n.ia
  }  else {
    r.go <- r[go.ia.id]
    n.go <- n[go.ia.id]
    pos.fa <- find_postk(lambda, gamma, r.go, n.go, plo, phi)
    
    
    rej.fa.id <- go.ia.id[pos.fa>T]
    fut.fa.id <- go.ia.id[pos.fa<=T]
    
    fut.id <- c(fut.ia.id,fut.fa.id)
    rej.id <- c(rej.fa.id)
    
    n.tot <- n
    n.tot[-go.ia.id] <- n.ia[-go.ia.id]
    if(length(fut.id)==k){
      rej <- rep(0,k)
      
      
    } else if(length(rej.id)==k){
      rej <- rep(1,k)
      
    } else {
      rej <- rep(1,k)
      rej[fut.id] <- 0
    }
    
  }
  
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
#' @param lambda the prior probability of H0: that the strata are completely correlated
#' @param gamma the prior probability that the drug is active in any specific stratum
#' @param plo response probability for inactive drug
#' @param phi response probability for active drug
#' @param T threshold for conclusive posterior probability
#' @return return drug active ones, total number of subjects
#' @export
SimonBayes.type1 <- function(rep, k, n, n.ia = NULL, resp.true, lambda, gamma, plo, phi, T){
  
  ## Simulate data
  dat <- data.frame(id=NA,cohort=NA,orr=NA)
  for(i in 1:k){
    cohort <- c(rep(i,n[i]))
    dat.tmp <- data.frame(cohort)
    id <- 1:n[i]
    # Simulate ORR： overall response rate
    orr <- rbinom(n=n[i], size=1, prob=resp.true[i])
    
    
    dat.tmp <- data.frame(dat.tmp,id,orr)
    dat <- rbind(dat,dat.tmp)
    
  }
  
  dat <- dat[-1,]
  
  if(is.null(n.ia)){
    
    r <- rep(NA,k)
    for(i in 1: k){
      dat.tmp <- subset(dat,cohort==i)
      r[i] <- nrow(subset(dat.tmp,orr==1))
      
    }
    
    Simon.Bayes.res <- SimonBayes.2stg(r.ia=NULL,r,lambda, gamma, plo, phi,T,n.ia=NULL,k,n)
    
  } else {
    r.ia <- rep(NA,k)
    r <- rep(NA,k)
    for(i in 1: k){
      dat.tmp <- subset(dat,cohort==i)
      r.ia[i] <- nrow(subset(dat.tmp[1:n.ia[i],],orr==1))
      r[i] <- nrow(subset(dat.tmp,orr==1))
      
    }
    
    Simon.Bayes.res <- SimonBayes.2stg(r.ia,r,lambda, gamma, plo, phi,T,n.ia,k,n)
    
  }
  
  return(Simon.Bayes.res)
  
}


#' search for typeI error under differen T
#'
#' @param k number of cohorts
#' @param n number of patients
#' @param n.ia number of patients for interium analysis
#' @param resp.true probability under the null
#' @param lambda the prior probability of H0: that the strata are completely correlated
#' @param gamma the prior probability that the drug is active in any specific stratum
#' @param plo response probability for inactive drug
#' @param phi response probability for active drug
#' @param T threshold for conclusive posterior probability
#' @param nsim Number of simulations for each T
#' @return return typeI error under differen T
#' @export
SimonBayes.type1.search <- function(k, n, n.ia, resp.true, lambda, gamma, plo, phi, T, nsim){
  typ1.lst <-  rep(NA,length(T.lst))
  for(j in 1:length(T.lst)){
    T <- T.lst[j]
    res <- sapply(1:nsim, SimonBayes.type1, k, n, n.ia, resp.true, lambda, gamma, plo, phi, T)
    typ1 <- res[1:k,]
    typ1.lst[j] <- mean(colSums(typ1)>0) 
  }
  plot_search_type1 = ggplot2::ggplot(data=data.frame(T.lst, typ1.lst), ggplot2::aes(x=T.lst, y=typ1.lst, group=1)) +
    ggplot2::geom_line()+
    ggplot2::geom_point() + ggplot2::ggtitle("Global type I error with different threshold")+
    ggplot2::labs(x="Posterior cutoff at final analysis", y="Global type I error")
  
  #setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  #png(paste0("Global type I error with different threshold",".png"),800,800)
  #par(mar=c(5,6,6,2)+0.1, oma=c(1, 1, 1, 1)) ##c(5,6,4,2)+0.1
  #plot_search_type1
  #dev.off()
  
  return(list(data.frame(T.lst, typ1.lst), plot_search_type1))
}
