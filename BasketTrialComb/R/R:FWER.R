#' Get FWER
#' 
#' @param post.p posterior probability of efficacy
#' @param type_one fixed type one error for each cohort 
#' @return family wise error rate for fixed type one error
#' @export
fwer = function(post.p,type_one){
  ntype <- dim(post.p)[2]
  eff.cut <- sapply(1:ntype, FUN = function(x){quantile(post.p[,x], 1-type_one)})
  Decision <- post.p > eff.cut
  fwer <- mean(apply(Decision, 1, function(o) sum(o)>0))
  return(fwer)
}
