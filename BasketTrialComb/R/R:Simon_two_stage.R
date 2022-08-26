#' Probability to stop for futility at IA
#' 
#' @import dplyr
#' @import poibin
#' @import gtools
#' @param r1 number of responses at IA
#' @param size sample size
#' @param prob probability vector
#' @return Probability to stop
#' @export
PET <- function(r1,n1,p){
  pet <- pbinom(r1, size=n1, prob=p)
  return(pet)
}

#' Prob for stop futility at end of the trial
#'
#' @param r1 number of responses at IA
#' @param n1 number of patients at IA
#' @param p probability vector
#' @param r number of responses
#' @param n number of patients
#' @return probability of reject
#' @export
ProbRej <- function(r1, n1, p, r, n){
  
  cp<-0
  
  for(j in (r1+1):min(n1, r) ){
    cp <- cp + dbinom(j, size=n1, prob=p)*pbinom(r-j, size=n-n1, prob=p)
  }
  
  
  prej <- PET(r1,n1,p) + cp
  return(prej)
}

#' Expected N
#'
#' @param r1 number of responses at IA
#' @param n1 number of patients at IA
#' @param p probability vector
#' @param n number of patients
#' @return Expected sample size
#' @export
EN <- function(r1,n1,p,n){
  en <- n1 + (1-PET(r1,n1,p))*(n-n1)
  return(en)
}

#' Optim parameter
#' 
#' @param n.s2stg number of patients at 2stgs
#' @param n1.s2stg number of patients
#' @param p0.s2stg number of responses for interium 1st stage
#' @param p1.s2stg adjusted alpha for 1st stage
#' @param r1.s2stg number of responses
#' @param r.s2stg number of responses
#' @return parameters
#' @export
Simon.2stg.search <- function(n.s2stg, n1.s2stg, p0.s2stg, p1.s2stg, r1.s2stg, r.s2stg){
  
  sum.2stage <- data.frame(n=NA,n1=NA,r1=NA,r=NA,Type1=NA, Power=NA, EN=NA, PET=NA )
  
  for(i in 1:length(r1.s2stg)){
    r1.tmp <- r1.s2stg [i]
    for(j in (r1.tmp+1):r.s2stg){
      # Type 1 error
      typ1 <- 1-ProbRej(r1.tmp,n1.s2stg,p=p0.s2stg, j,n.s2stg)
      
      # Power
      power <- 1-ProbRej(r1.tmp,n1.s2stg,p=p1.s2stg, j,n.s2stg)
      
      # Early stopping for futility
      early.stop <- PET(r1.tmp,n1.s2stg,p=p0.s2stg)
      
      # EN
      
      en <- EN(r1.tmp,n1.s2stg,p=p0.s2stg,n.s2stg)
      
      sum.2stage.tmp <- data.frame(n=n.s2stg,n1=n1.s2stg,r1=r1.tmp,r=j,Type1=round(typ1,6), Power=round(power,6), EN=round(en,6), PET=round(early.stop,6) )
      sum.2stage <- rbind(sum.2stage,sum.2stage.tmp)
    }
  }
  
  return(sum.2stage)
  
}

#' main function decision making
#' 
#' @param n1.s2stg number of patients
#' @param n.s2stg number of patients at 2stgs
#' @param r1.s2stg number of responses
#' @param r.s2stg number of responses
#' @param r.ia number of responses at IA
#' @param r number of responses
#' @return return drug active ones, total number of subjects
#' @export
Simon.2stage <- function(n1.s2stg, n.s2stg, r1.s2stg, r.s2stg, r.ia, r){
  if(r.ia <=r1.s2stg ){
    rej <- 0
    n.tot <- n1.s2stg
  } else {
    n.tot <- n.s2stg
    if(r <=r.s2stg){
      rej <- 0
      
    } else {
      rej <- 1
    }
  }
  res <- c(rej,n.tot)
  return(res)
}
