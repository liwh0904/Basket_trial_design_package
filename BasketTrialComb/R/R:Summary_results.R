#'
#' Summary Results 
#' 
#' @param res results from Monte Carlo Simulation
#' @param p is the number of methods
#' @param q is the number of cohorts
#' @param mnames names of methods
#' @param tp True positive indicator
#' @return summary results 
#' @export
summary_sim <- function(res, p = 5, q = 4, mnames = c('SimonBayes','Opt2stg','Opt1stg','Simon2stg','Pool'), 
                        tp){
  
  outdat <- data.frame(plch=rep(NA,8))
  for(i in 1:p){
    row <- 1 + (i - 1) * (q * 2)
    tmp.mthd.rej <- res[(row:(row+q-1)),]
    tmp.mthd.n <- res[( (row+q):(row+2*q-1)),]
    tmp.rej.mean <- rowMeans(tmp.mthd.rej)
    tmp.rej.sd <- apply(tmp.mthd.rej, 1, sd)
    tmp.n.mean <- mean(colSums(tmp.mthd.n))
    tmp.n.sd <- sd(colSums(tmp.mthd.n))
    
    tmp.grej.mean <- mean(colSums( tmp.mthd.rej)>0)
    tmp.grej.sd <- sd(colSums( tmp.mthd.rej)>0)
    if(is.null(tp)){
      tmp.tp.mean <- NA
      tmp.tp.sd <- NA
      tmp.fp.mean <- mean(apply(tmp.mthd.rej, 2, function(x) sum(x==1) ))
      tmp.fp.sd <- sd(apply(tmp.mthd.rej, 2, function(x) sum(x==1) ))
    } else {
      tmp.tp.mean <- mean(apply(tmp.mthd.rej, 2, function(x) sum(x[tp]==1) ))
      tmp.fp.mean <- mean(apply(tmp.mthd.rej, 2, function(x) sum(x[-tp]==1) ))
      
      tmp.tp.sd <- sd(apply(tmp.mthd.rej, 2, function(x) sum(x[tp]==1) ))
      tmp.fp.sd <- sd(apply(tmp.mthd.rej, 2, function(x) sum(x[-tp]==1) ))
      
    }
    
    m.lst <- c(tmp.grej.mean,tmp.rej.mean,tmp.tp.mean,tmp.fp.mean,tmp.n.mean)
    sd.lst <- c(tmp.rej.sd,tmp.grej.sd,tmp.tp.sd,tmp.fp.sd,tmp.n.sd)
    
    outdat <- cbind(outdat, data.frame(m.lst))
    
  }
  outdat <- outdat[,-1]
  if(length(mnames) == 1){
     names(outdat) =  c('Global', sapply(c(1:k),function(x) paste0('Ind',x,sep='') ),'TP','FP','EN')
  }else{
    colnames(outdat) <- mnames
    row.names(outdat) <- c('Global', sapply(c(1:k),function(x) paste0('Ind',x,sep='') ),'TP','FP','EN')
  }
 return(outdat)
}

