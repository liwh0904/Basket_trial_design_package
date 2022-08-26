#'
#' This function use binomial exact test to calculate the sample size. 
#' 
#' @param p probability vector, p[1] is the null hypothesis probability, p[2] is the alternative hypothesis probability
#' @param alpha type I error
#' @param beta type II error
#' @return sample size
#' @export
cal.N.2stg <- function(p, alpha, beta){
  
  for(i in 1:1000000000){
    # quantile of 1 - alpha
    r<-qbinom(p=1-alpha,size=i,prob=p[1])
    
    #reject region is great than r
    #control type I error then get power greater than 1- beta with minimum sample size
    power<-(1-pbinom(r,i,p[2]))
    if(power>(1-beta)){
      return(i)
    }
  }
  
}


#' This helper function to calculate the power
#' 
#' @param K number of indications
#' @param g the prior probability of H0: that the strata are completely correlated
#' @param r1 number of responses for stage I 
#' @param n1 number of patients for stage I 
#' @param n2 number of patients for stage II
#' @param p0 probability vector of the null hypothesis
#' @param p1 probability vector of the alternative hypothesis
#' @param alphastar adjusted alpha
#' @param nsim number of simulation used for type I and type II calculation
#' @return power
#' @export
PowerEmpirical.helper.2stg<-function(K,g,r1,n1,n2,p0,p1,alphastar,nsim=10000){
  # g is vector
  # n1 is stage I sample size
  # n2 is stage II sample size
  
  accept = rej = 0
  n0pts=n1
  n1pts=n2
  n12pts=n1+n2
  
  for (sim in 1:nsim)
    
  {
    #?
    probs <- p1*g+p0*(1-g)
    
    resp0 = rbinom(K,n0pts,probs)
    
    # stage I type I error control thru r1 for each indication
    
    # r1 = qbinom(1-alpha0,n0pts,p0)
    
    if(all(resp0 <= r1)) rej = rej + 1
    
    else {
      
      index = which(resp0 > r1)
      # allocate extra from inactive to remaining active tumors
      #extra=sum(n1pts[-index])
      extra <- 0
      #modulo
      # no extra, no reallocation
      common_add<-extra%/%length(index)
      #numerical division
      remainder<-extra%%length(index)
      #N2 is the sample size of stage II
      N2<-n1pts[index]+common_add
      
      #random assign remainder to indications
      if(remainder!=0){
        remainder_add_index<-sample(1:length(index),remainder,replace=FALSE)
        
        
        N2<-N2[remainder_add_index]+1
      }
      resp1=rbinom(length(index),N2,probs[index])
      
      # control the combined two stage type I error for each indication
      
      #r2<-qbinom(1-alpha1,n1pts,p0)
      
      resp_combine<-resp0[index]+resp1
      r<- qpoibin(qq=1-alphastar,pp=p0[index],wts=N2+n0pts[index])
      
      resp.pool = sum(resp_combine)
      
      if(resp.pool <= r) {rej = rej + 1}
      
      else{ accept = accept + 1}
      
      
    }
    
  }
  
  return(accept/nsim)
}

#' This function uses simulation to calculate the power
#' 
#' @param K number of indications
#' @param g the prior probability of H0: that the strata are completely correlated
#' @param alpha0 the type I error parameter for stage1
#' @param beta0 the type II error parameter for stage1
#' @param n12pts total sample size for each indication
#' @param p0 probability vector of the null hypothesis
#' @param p1 probability vector of the alternative hypothesis
#' @param alphastar adjusted alpha
#' @param nsim number of simulation used for type I and type II calculation
#' @return power
#' @export
PowerEmpirical.2stg <- function(K,g,alpha0,beta0,n12pts,alphastar,p0,p1,nsim=10000){
  
  accept = rej = 0
  
  ### grid search for sample size meet type I and type II error requirments
  n0pts = apply(cbind(p0,p1),1,cal.N.2stg,alpha=alpha0,beta=beta0)
  
  #n12pts = apply(cbind(p0,p1),1,cal.N,alpha=alpha1,beta=beta1)
  # control type 1 error 
  r1 = qbinom(1-alpha0,n0pts,p0)
  
  
  return(PowerEmpirical.helper.2stg(K,g,r1,n1=n0pts,n2=n12pts-n0pts,p0,p1,alphastar,nsim))
  
  
  
}


#' This function finds alphastar
#'
#' @param K number of indications
#' @param alpha0 the type I error parameter for stage1
#' @param beta0 the type II error parameter for stage1
#' @param n12pts total sample size for each indication
#' @param p0 probability vector of the null hypothesis
#' @param p1 probability vector of the alternative hypothesis
#' @param alphastar adjusted alpha
#' @param nsim number of simulation used for type I and type II calculation
#' @param alpha type I error

#' @return return alphastar
#' @export
diff0.emp.2stg <- function(alphastar,K,alpha0,beta0,n12pts,p0,p1,nsim,alpha){
  
  PowerEmpirical.2stg(K=K,g=0,alpha0=alpha0,beta0=beta0,n12pts=n12pts, alphastar=alphastar,p0=p0, p1=p1, nsim=nsim)-alpha
  
}

#' This function finds alphastar
#' @param alphastar adjusted alpha
#' @param alpbeta0 alpha
#' @param n12pts total sample size for each indication
#' @param K number of indications
#' @param p0 probability vector of the null hypothesis
#' @param p1 probability vector of the alternative hypothesis
#' @param nsim number of simulation used for type I and type II calculation
#' @param alpha type I error
#' @return return alphastar
#' @export
power.cal.2stg <- function(alphastar, alpbeta0,n12pts,K,p0,p1,nsim,alpha)
{
  
  # scenarios with at least one positive
  g.list <- permutations(2,K,0:1,repeats.allowed = TRUE)[-1,] 
  
  power.seq <- n.active <-  rep(0,nrow(g.list))
  
  for (j in 1:nrow(g.list)){
    
    g <- g.list[j,]
    
    n.active[j] <- sum(g)
    
    power.seq[j] <- PowerEmpirical.2stg(K=K,g=g,alpha0=alpbeta0[1],beta0=alpbeta0[2],n12pts=n12pts, alphastar=alphastar, p0=p0, p1=p1, nsim=nsim)
    
  }
  
  temp.list <- cbind(n.active, power.seq)#useless
  #? each cohort mean, every cohorts mean
  power.unif <- mean(sapply(split(power.seq,n.active),mean))
  
  
  
  return(power.unif)
  
}

#' Optim paramerters setting
#'
#' @param K number of indications
#' @param Alpha global type I error
#' @param p0 probability vector of the null hypothesis
#' @param p1 probability vector of the alternative hypothesis
#' @param n12pts number of patients
#' @param a0range range for grid search 
#' @param b0range range for grid search 
#' @param nsim number of simulation used for type I and type II calculation
#' @param track tracking
#' @return optim paramerters setting for optimal one stage design
#' @export
fixed_twostage <- function(K, Alpha, p0, p1, a0range=seq(alpha,0.3,by=0.1), b0range=seq(0.05,0.3,by=0.1), 
                           n12pts ,nsim=10000, track=TRUE){
  
  
  pow.cur = 0
  
  this <- NULL
  
  ntotal<-length(a0range)*length(b0range)
  
  count=0
  
  validcount=0
  
  for (a0 in a0range){
    
    for(b0 in b0range){
      
      
      
      ## check whether the a0,b0,a1,b1 satifies contraint that total sample size should be larger than stage I sample size for each indication
      
      
      
      # n0pts is the sample size of stage I
      
      # n1pts is the sample size of stage I and II combined
      # for loop for alpha0 and beta0
      n0pts = apply(cbind(p0,p1),1,cal.N.2stg,alpha=a0,beta=b0)
      if(!all(n12pts>n0pts)){
        cat("n12pts less than  n0pts","\n")
        next
      }else{
        
        
        stage_II_patients=n12pts-n0pts
        
        alphastar <- tryCatch(uniroot(diff0.emp.2stg,c(0,1),K=K,alpha0=a0,beta0=b0,n12pts=n12pts,p0=p0,
                                      p1=p1,nsim=nsim,alpha=Alpha)$root, error=function(e) return("alphastar Uniroot unsolvable"))
        if(is.character(alphastar)){
          cat("alphastar not solved here","\n")
          next
        }
        #? uniform power
        pow.temp <- tryCatch(power.cal.2stg(alphastar, alpbeta0=c(a0,b0),n12pts=n12pts,K,p0,p1,nsim,
                                            alpha=Alpha), error=function(e) return("Uniroot unsolvable"))
        
        
        if (is.character(pow.temp)){
          print("power not solvable")
          next
          
        }else { # record the power, do not check value
          
          pow.cur <- pow.temp;
          
          
        }
        
        
        
        typeIerror = PowerEmpirical.2stg(K=K,g=0,alpha0=a0,beta0=b0,n12pts=n12pts,alphastar=alphastar,p0=p0,p1=p1,nsim=nsim)
        
        while(typeIerror > Alpha +0.001 && alphastar >0.001) { # if typeIerror still exceeds 0.05, then fine tune alphastar;
          #?
          alphastar.new = alphastar - 0.001
          
          alphastar <- alphastar.new
          
          typeIerror <- round(PowerEmpirical.2stg(K=K,g=0,alpha0=a0,beta0=b0,n12pts=n12pts,alphastar=alphastar,p0=p0,p1=p1,nsim=nsim),3)
          
        }
        r <- qbinom(1-a0, n0pts, p0)  # r for the stage I
        
        
        this <- rbind(this, c(alpha0=a0,beta0=b0,stage_I_patients=n0pts, stage_I_pooling_bar_r=r,stage_II_patients=stage_II_patients, stage_I_II_combine=n12pts,alphastar=alphastar, typeIerror=typeIerror, power=pow.cur,Npatients=sum(n12pts)))
        
        count=count+1
        
        validcount=validcount+1
        
        if(track==TRUE) print(paste(" currently calculate case#=", count,  "valid case #=", validcount,"out of total=",ntotal, "alpha0= ", a0,"beta0= ",b0, "power=",pow.cur, sep=" "))
        
      } 
      
    }
  }
  
  
  
  
  
  
  
  print(paste("finally legetimate case numbers is ", validcount, "out of ", ntotal,sep = " "))
  
  
  
  #generate design table using the searched parameter.
  
  result <- NULL
  
  
  result <- as.data.frame(this)
  
  typeIerror_list<-result$typeIerror< Alpha
  candi_result<-result[typeIerror_list,]
  best_id<-which(result[typeIerror_list,"power"]==max(result[typeIerror_list,"power"]))
  
  return(list(result=result, largest_pow=candi_result[best_id,] , large_pow=candi_result ,parameter=this))
  
}

#' main function decision making
#' 
#' @param k number of indications
#' @param n number of patients
#' @param r.stg1 number of responses for interium 1st stage
#' @param n.stg1 number of responses for interium 1st stage
#' @param alphastar.2stg adjusted alpha for 2stg 
#' @param r.ia.2stg number of responses for interium 2stg
#' @param r number of responses 
#' @param resp.cntl probability under the null
#' @return return drug active ones, total number of subjects
#' @export
Optimal.2stg<- function(k, n, r.stg1, n.stg1, alphastar.2stg, r.ia.2stg, r,resp.cntl){
  
  
  fut.id.stg1 <- which(r.ia.2stg<r.stg1)
  pool <- c(1:k)[!(1:k) %in% fut.id.stg1]
  resp_pool <- sum(r[!(1:k) %in% fut.id.stg1])
  
  if(length(pool)==0){
    pval <- 1
    
    rej <- rep(0,k)
    n.tot <- n.stg1
    
  } else {
    rej <- rep(0,k)
    n.tot <- n.stg1
    n.tot[pool] <- n[pool]
    pval <- 1- poibin::ppoibin(resp_pool-1, pp=resp.cntl[pool],wts=n[pool])
    
    if(pval<alphastar.2stg){
      rej[pool] <- 1
      
    } 
  }
  
  res <- c(rej,n.tot)
  return(res)
  
  
}
