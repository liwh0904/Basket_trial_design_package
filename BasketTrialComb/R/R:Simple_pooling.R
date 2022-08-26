#' main function decision making
#' 
#' @import dplyr
#' @import poibin
#' @import gtools
#' @param n number of patients
#' @param K number of indications
#' @param r number of responses
#' @param resp.cntl probability under the null
#' @param alpha global type I error
#' @return  return drug active ones, total number of subjects
#' @export
Simple.Pool <- function(n,k,r,resp.cntl,alpha){
  pval <- 1- poibin::ppoibin(sum(r)-1, pp=resp.cntl,wts=n)
  n.tot <- n
  if(pval<alpha){
    rej<-rep(1,k)
  } else {
    rej<-rep(0,k)
  }
  res <- c(rej,n.tot)
  return(res)
}

