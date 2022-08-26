devtools::document()
devtools::install()
library(BasketTrialComb)

##################################################
##### Simulation performance metrics #####
##### 
##################################################

### Monte Carlo Simulation for posterior mean
Bask.Sim.Test <- function(rep, k = 4, n = rep(40, 4), n.ia = rep(15, 4), resp.true =  rep(0.65, 4), 
                          resp.cntl =  rep(0.5, 4), resp.trt =  rep(0.65, 4), 
                          
                          lambda = 0.33, gamma = 0.5, plo = 0.2, phi = 0.4, T =  seq(0.4, 0.8, 0.025),
                          r.stg1 = NULL, n.stg1 = NULL, alphastar.2stg = NULL, r.1stg = NULL, 
                          alphastar.1stg = NULL, n1.s2stg = 10, n.s2stg = 20, 
                          r1.s2stg = NULL, r.s2stg = NULL, alpha = 0.1, 
                          # muce
                          J = 1, L = 1, nik = matrix(rep(c(15, 40), 2), nrow = 2, ncol = 4), v1 = 0.3, 
                          eff.ref.muce = NULL, 
                          gamma_muce = 2.5, prec0 = 1, prec.eta = 1, prec.ksi = 1, mu.eta0 = 0,
                          prec.eta0 = 1, mu.ksi0 = 0, prec.ksi0 = 1, n.burnin = 1000, n.chain = 1, n.iter = 8000,
                          # exnex
                          Q.exnex = NULL, nk = matrix(rep(c(15, 40), 2), nrow = 2, ncol = 4), Nexch = 1, 
                          Nmix = 2, pMix = c(0.5, 0.5), 
                          dirichprior = c(1,1,1), dirichprior.check = FALSE, mu.mean = NULL, 
                          mu.prec = mean(1/(1/resp.cntl+1/(1-resp.cntl)-1)),
                          tau.HN.scale = 1, nex.mean = NULL, nex.prec = NULL,
                          # BHM
                          cohortsize = c(15, 25), mu.par.BHM = NULL, eff.ref.bhm = NULL,
                          # CBHM
                          mu.par.CBHM = NULL, eff.ref.CBHM = NULL, a = NULL, b = NULL,
                          # Indications
                          Simon.Bayes.Ind = 0, Optimal.2stg.Ind = 0, Optimal.1stg.Ind = 0, 
                          Simon.2stg.Ind = 0, Simple.pooling.Ind = 0, MUCE.Ind = 0, EXNEX.Ind = 0, 
                          BHM.Ind = 0, CBHM.Ind = 0, Binom.Test.Ind = 0){
  
  ## Simulate data
  dat <- data.frame(id = NA, cohort = NA, orr = NA)
  for(i in 1:k){
    cohort <- c(rep(i, n[i]))
    dat.tmp <- data.frame(cohort)
    id <- 1:n[i]
    # Simulate ORR
    orr <- rbinom(n = n[i], size = 1, prob = resp.true[i])
    
    dat.tmp <- data.frame(dat.tmp, id, orr)
    dat <- rbind(dat, dat.tmp)
    
  }
  
  dat <- dat[-1,]
  
  r.ia <- rep(NA,k)
  r <- rep(NA,k)
  r.ia.2stg <- rep(NA,k)
  for(i in 1: k){
    dat.tmp <- subset(dat, cohort == i)
    if(is.null(n.ia)){
      r.ia = NULL
      
      r[i] <- nrow(subset(dat.tmp, orr == 1))
      r.all = matrix(r, nrow = length(cohortsize), ncol = k)
    }else{
      r.ia[i] <- nrow(subset(dat.tmp[1:n.ia[i],], orr == 1))
      
      r[i] <- nrow(subset(dat.tmp, orr == 1))
      r.all = rbind(r.ia, r - r.ia)
    }
    if(Optimal.2stg.Ind == 1) {
      r.ia.2stg[i] <- nrow(subset(dat.tmp[1:n.stg1[i],], orr == 1))
    }
  }
  
  ## Binomimal test
  binom.test.res = NULL
  if(Binom.Test.Ind == 1){
    binom.test.1 = sum(binom.test(x = r[1], n = n[1], p = resp.cntl[i])$p.value < FWER/k)
    binom.test.2 = sum(binom.test(x = r[2], n = n[2], p = resp.cntl[i])$p.value < FWER/k)
    binom.test.3 = sum(binom.test(x = r[3], n = n[3], p = resp.cntl[i])$p.value < FWER/k)
    binom.test.4 = sum(binom.test(x = r[4], n = n[4], p = resp.cntl[i])$p.value < FWER/k)
    binom.test.res = c(binom.test.1, binom.test.2, binom.test.3, binom.test.4, n)
  } 
  
  ## Apply methods Simon's Bayesian
  Simon.Bayes.res = NULL
  if(Simon.Bayes.Ind == 1){
    Simon.Bayes.res <- SimonBayes.2stg(r.ia, r, lambda, gamma, plo, phi, T, n.ia, k, n)
  } 
  
  ## Apply methods Optim1stg 
  opt.1stg.res = NULL
  if(Optimal.1stg.Ind == 1){
    opt.1stg.res <- optimal.1stage(n, k, r.1stg, alphastar.1stg, r, resp.cntl)
  } 
  
  ## Apply methods Optim2stg
  opt.2stg.res = NULL
  if(Optimal.2stg.Ind == 1){
    opt.2stg.res <- Optimal.2stg(k, n, r.stg1, n.stg1, alphastar.2stg, r.ia.2stg, r, resp.cntl)
  } 
    
  ## Apply methods Simon's 2-stage
  simon.2stg.res = NULL
  if(Simon.2stg.Ind == 1){
    alpha.simon2stg <- 1 - (1 - alpha)^(1/k)
    simon2stage.rej.lst <- rep(NA, k)
    simon2stage.n.lst <- rep(NA, k)
    for(f in 1:k){
      simon2stage.res.lst<-Simon.2stage(n1.s2stg, n.s2stg, r1.s2stg, r.s2stg, r.ia[f], r[f])
      simon2stage.rej.lst[f] <- simon2stage.res.lst[1]
      simon2stage.n.lst[f] <- simon2stage.res.lst[2]
    }
    simon.2stg.res <- c(simon2stage.rej.lst, simon2stage.n.lst)
  } 
  
  ## Apply methods Simple pooling
  simple.pool.res = NULL
  if(Simple.pooling.Ind == 1){
    simple.pool.res <- Simple.Pool(n, k, r, resp.cntl, alpha = alpha)
  } 
  
  ## Apply methods MUCE
  MUCE.res = NULL
  if(MUCE.Ind == 1){
    MUCE.res <- muce(K = k, J = J, L = L, nik = nik, H0 = resp.cntl, H1 = resp.trt, p0 = resp.true,
                     v1 = v1, v2 = eff.ref.muce, gamma = gamma_muce, prec0 = prec0, prec.eta = prec.eta, 
                     prec.ksi = prec.ksi, mu.eta0 = mu.eta0, prec.eta0 = prec.eta0, mu.ksi0 = mu.ksi0,
                     prec.ksi0 = prec.ksi0, n.burnin = n.burnin, n.chain = n.chain, n.iter = n.iter, 
                     r.all, num.sim = 1)
    MUCE.res$decision.muce[MUCE.res$decision.muce == -1] = 0
    MUCE.res = c(MUCE.res$decision.muce, MUCE.res$N.SS)
  } 
  
  ## Apply methods EXNEX
  EXNEX.res = NULL
  if(EXNEX.Ind == 1){
    EXNEX.res <- exnex(num.sim = 1, Q.exnex = eff.ref.exnex, p0 = resp.true, H0 = resp.cntl, H1 = resp.trt, 
                           nk = nik, k, v1, Nexch = Nexch, Nmix = Nmix, pMix = pMix, dirichprior = c(1,1,1),
                           dirichprior.check = FALSE, mu.mean = mu.mean, mu.prec = mu.prec,
                           tau.HN.scale = tau.HN.scale, nex.mean = nex.mean, nex.prec = nex.prec,
                           n.burnin = n.burnin, n.chain = n.chain, n.iter = n.iter, r.all)
    EXNEX.res$decision.exnex[EXNEX.res$decision.exnex == -1] = 0
    EXNEX.res = c(EXNEX.res$decision.exnex, EXNEX.res$N.SS)
  } 
  
  ## Apply methods BHM
  BHM.res = NULL
  if(BHM.Ind == 1){
    BHM.res <- berrymethod(cohortsize = cohortsize, ntype = k, p.true = resp.true, p.null = resp.cntl, 
                           p.target = resp.trt, ntrial = 1, n.burnin = n.burnin, n.chain = n.chain, 
                           n.iter = n.iter, mu.par = mu.par.BHM, eff.ref.bhm, r.all)
    efficacy <- BHM.res[(k + 1):(2 * k )]
    sample.size <- BHM.res[1:k]
    BHM.res = c(efficacy, sample.size)
  } 
  
  ## Apply methods CBHM
  CBHM.res = NULL
  if(CBHM.Ind == 1){
    CBHM.res <- CBHMdesign(cohortsize = cohortsize, ntype = k, r.all, p.null = resp.cntl, 
                           p.target = resp.trt, ntrial = 1, n.burnin = n.burnin, n.chain = n.chain,
                           n.iter = n.iter, mu.par = mu.par.CBHM, eff.ref = eff.ref.CBHM, 
                           a = a.par, b = b.par, r.all)
    if(L == 1){
      efficacy <- CBHM.res[7:10]
      sample.size <- CBHM.res[3:6]
    }else{
      efficacy <- CBHM.res[6:9]
      sample.size <- CBHM.res[2:5]
    }

    #efficacy <- CBHM.res[(2 * k  + 1):(3 * k)]
    #sample.size <- CBHM.res[(k+1):(2*k)]
    CBHM.res = c(efficacy, sample.size)
  } 
  
  res <- c(Simon.Bayes.res, opt.2stg.res, opt.1stg.res, simon.2stg.res, simple.pool.res, 
           MUCE.res, EXNEX.res, BHM.res, CBHM.res, binom.test.res) 
  return(res)
  
}

#### Simulation Settings  --------------
set.seed(10) 
# Indications for which method to use
Simon.Bayes.Ind = 1
Optimal.2stg.Ind = 1
Optimal.1stg.Ind = 1
Simon.2stg.Ind = 1
Simple.pooling.Ind = 1
MUCE.Ind = 1
EXNEX.Ind = 1
BHM.Ind = 1
CBHM.Ind = 1
Binom.Test.Ind = 0


# K cohorts/indications
k <- 4
nsim = 1000 # number of MC simulations
# total samples
n <- rep(40, 4) #20
# probability under the alternative
resp.trt <- rep(0.65, 4) #0.4
# probability under the null
resp.cntl <- rep(0.5, 4) #0.2
resp.true <- resp.cntl
alpha <- 0.1

#### settings for interim analysis
## for Simon's Bayes samples for interium analysis
#without IA
#n.ia <- NULL
#with IA
n.ia <- rep(15, 4) #

## for Simon 2-stage ####
n.s2stg <- 40
n1.s2stg  <- 15

## for muce and exnex
## without ia
#nik <- matrix(rep(c(20), 2), nrow =1 , ncol = 4) 
## with ia
nik <-  matrix(rep(c(15, 40), 2), nrow = 2, ncol = 4) #rbind(n.ia, n) ##the sample size matrix. each row is the sample size at each analysis, each col is cohort.
## for MUCE
##no. of interim analysis
L <- 1  

## for BHM and CBHM
# with ia
cohortsize <- c(15, 25) #
#cohortsize = c(n.ia[1], n[1])

#### (1) Simon's Bayes ####
plo <- resp.cntl[1]
phi <- resp.trt[1]
lambda <- 0.33
gamma <- 0.5
# T=threshold for conclusive posterior probability
T.lst <- seq(0.4, 0.8, 0.025)
T.search.result = SimonBayes.type1.search(k, n, n.ia, resp.true, lambda, gamma, plo, phi, 
                                          T = T.lst, nsim)

#### (2) Optimal 2-stage ####
optimal.2stg.params <- fixed_twostage(K = k, Alpha = alpha, p0 = resp.cntl, p1 = resp.trt, 
                                      a0range = seq(0.05, 0.5, by = 0.05), b0range = seq(0.05, 0.5, by = 0.05), 
                                      n12pts = n, nsim, track = TRUE)
optimal.2stg.parm <- optimal.2stg.params$largest_pow
# n.ia no use
n.stg1 <- as.numeric(optimal.2stg.parm[1,3:6])
r.stg1 <- as.numeric(optimal.2stg.parm[1,7:10])
alphastar.2stg <- as.numeric(optimal.2stg.parm[1,'alphastar'])

#### (3) Optimal 1-stage ####
optimal.1stg.params <- Optimal_pooling_hetero(K = k, n = n, p0 = resp.cntl, p1 = resp.trt, 
                                              Alpha = alpha, nsim)
optimal.1stg.parm <- optimal.1stg.params$opt_power

r.1stg <- as.numeric(optimal.1stg.parm[c('r1','r2','r3','r4')])
alphastar.1stg <- as.numeric(optimal.1stg.parm['alphastar'])

#### (4) Simon 2-stage ####
p0.s2stg  <- resp.cntl[1] #0.2
p1.s2stg  <- resp.trt[1] #0.4
r1.s2stg.lst <- seq(1, 15, 1)
r.s2stg.lst  <- 40

alpha.simon2stg <- 1 - (1 - alpha)^(1/k)
sum.2stage <- Simon.2stg.search(n.s2stg, n1.s2stg, p0.s2stg, p1.s2stg, r1.s2stg.lst, r.s2stg.lst )
sum.2stage <- sum.2stage[-1,]
sum.2stage.val <- subset(sum.2stage, Type1 < alpha.simon2stg)
# plot y sample size x power
# plot y sample size r1 power
sum.2stage.val = sum.2stage.val[abs(sum.2stage.val$Power  - max(sum.2stage.val$Power)) < max(sum.2stage.val$Power)/15,]
sum.2stage.pick <- sum.2stage.val[sum.2stage.val$EN == min(sum.2stage.val$EN),][1,]

#sum.2stage.val <- subset(sum.2stage,Type1<alpha.simon2stg)
#sum.2stage.pick <- sum.2stage.val[sum.2stage.val$Power==max(sum.2stage.val$Power),]
r1.s2stg <- sum.2stage.pick$r1
r.s2stg <- sum.2stage.pick$r

#### (5) Simple pool ####
# no need settings

#### For muce, exnex, BHM, CBHM ####
n.burnin <- 1000
n.chain <- 1
n.iter <- 8000

FWER = 0.1
v1 <- 0.3 ##threshold for futility.

#### (6) MUCE ####
##no. of dose levels
J <- 1 
# prior distribution
gamma_muce <- 2.5
prec0 <- 1
prec.eta <- 1
prec.ksi <- 1
mu.eta0 <- 0
prec.eta0 <- 1
mu.ksi0 <- 0
prec.ksi0 <- 1

##search for efficacy cutoff 
eff.ref.muce.v <- muce.calibrate(K= k, J = J, L = L, num.sim = 10000, FWER, precision = 0.0001, 
                               fwer.threshold = 0.001, nik = nik, H0 = resp.cntl, H1 = resp.trt, 
                               p0 = resp.true, v1 = v1, gamma = gamma_muce, prec0 = prec0, 
                               prec.eta = prec.eta, prec.ksi = prec.ksi, mu.eta0 = mu.eta0, 
                               prec.eta0 = prec.eta0, mu.ksi0 = mu.ksi0, prec.ksi0 = prec.ksi0, 
                               n.burnin = n.burnin, n.chain = n.chain, n.iter = n.iter)     
eff.ref.muce = eff.ref.muce.v[[1]] #mean(eff.ref.muce.v[[1]])

#### (7) EXNEX ####
Nexch <- 1
Nmix <- 2
## weight of exchange and nonexchange
pMix <- c(0.5, 0.5)
dirichprior.check <- FALSE
tau.HN.scale <- 1

prec.eta0 <- 1
mu.ksi0 <- 0
prec.ksi0 <- 1

mu.mean <- mean(log(resp.cntl/(1-resp.cntl)))
mu.prec <- mean(1/(1/resp.cntl+1/(1-resp.cntl)-1))
nex.mean <- log(resp.cntl/(1-resp.cntl))
nex.prec <- 1/(1/resp.cntl+1/(1-resp.cntl)-1)
nk = nik

##search for efficacy cutoff
## main efficacy for different cutoff???? cutoff for each indicators
Q.exnex <- exnex.calibration(k, v1 = v1, H0 = resp.cntl, p0 = resp.true, FWER, 
                             precision = 0.0001, fwer.threshold = 0.001,
                             nk = nik, num.sim = 10000, Nexch = Nexch, Nmix = Nmix, pMix = pMix,
                             dirichprior = c(1,1,1), dirichprior.check = FALSE, mu.mean = mu.mean, 
                             mu.prec = mu.prec, tau.HN.scale = tau.HN.scale, nex.mean = nex.mean,
                             nex.prec = nex.prec, n.burnin = n.burnin, n.chain = n.chain, 
                             n.iter = n.iter)
eff.ref.exnex <- mean(Q.exnex[[1]])

#### (8) BHM ####
mu.par.BHM <- mean(log(resp.cntl/(1-resp.cntl))-log(resp.trt/(1-resp.trt)))

BHM.results <- berrymethod.calibrate(cohortsize = cohortsize, ntype = k, p.true = resp.true, 
                                 p.null = resp.cntl, p.target = resp.trt, ntrial = 10000, 
                                 n.burnin = n.burnin, n.chain = n.chain, n.iter = n.iter, 
                                 mu.par = mu.par.BHM, FWER, precision = 0.0001, fwer.threshold = 0.001)
eff.ref.bhm = mean(BHM.results[[1]]$eff.ref)

#### (9) CBHM ####
var.small <- 1
var.big <- 80
mu.par.CBHM <- mean(log(resp.cntl/(1 - resp.cntl))) ##the testing parameter is different from BHM
 
Tun <- decidePar(cohortsize = cohortsize, ntype = k, ntrial = 10000, p0 = mean(resp.cntl), 
                        p1 = mean(resp.trt), var.small = var.small, var.big = var.big)
a_par_list <- Tun$a
b_par_list <- Tun$b
CBHM.results <- CBHMdesign.calibrate(cohortsize = cohortsize, ntype = k, p.true = resp.cntl, 
                                p.null = resp.cntl, p.target = resp.trt, ntrial = 10000, 
                                n.burnin = n.burnin, n.chain = n.chain, n.iter = n.iter, 
                                mu.par = mu.par.CBHM, FWER, precision = 0.0001, fwer.threshold = 0.001, 
                                a = a_par_list,  b = b_par_list)  
eff.ref.CBHM =  mean(CBHM.results[[1]]$eff.ref)
a.par <- a_par_list
b.par <- b_par_list

#### Saving data ####
save(k, nsim, n, resp.trt, resp.true, resp.cntl, n.ia, alpha,
     #Simon's Bayes
     plo, phi, lambda, gamma, T.lst, T.search.result,
     #Optimal 2-stage
     optimal.2stg.params, optimal.2stg.parm, n.stg1, r.stg1, alphastar.2stg,
     #Optimal 1-stage
     optimal.1stg.params, optimal.1stg.parm, r.1stg, alphastar.1stg,
     #Simon 2-stage
     n.s2stg, n1.s2stg, p0.s2stg, p1.s2stg, r1.s2stg.lst, r.s2stg.lst, alpha.simon2stg, sum.2stage,
     sum.2stage.val, sum.2stage.pick, r1.s2stg, r.s2stg,
     nik, n.burnin, n.chain, n.iter, FWER, v1,
     #MUCE,
     J, L, gamma_muce, prec0, prec.eta, prec.ksi, mu.eta0, prec.eta0, mu.ksi0, prec.ksi0, 
     eff.ref.muce.v, eff.ref.muce,
     #EXNEX
     Nexch, Nmix, pMix, dirichprior.check, tau.HN.scale, prec.eta0, mu.ksi0, prec.ksi0, mu.mean,
     mu.prec, nex.mean, nex.prec, Q.exnex, eff.ref.exnex,
     #BHM
     cohortsize, mu.par.BHM, BHM.results, eff.ref.bhm,
     #CBHM
     var.small, var.big, mu.par.CBHM, a_par_list, b_par_list, CBHM.results, eff.ref.CBHM, a.par, b.par,
     
     file = paste0("Basket_trial_design_setting_parameters_sample_20_with_ia", ".RData"))

load('Basket_trial_design_setting_parameters_sample_20_with_ia.RData')

#### 5 cases ####
set.seed(15) 
mnames = c('SimonBayes','Opt2stg','Opt1stg','Simon2stg','Pool', 'MUCE', 'EXNEX', 'BHM', 
           'CBHM', 'Independent')[c(Simon.Bayes.Ind == 1, Optimal.2stg.Ind == 1, Optimal.1stg.Ind == 1, 
                              Simon.2stg.Ind == 1, Simple.pooling.Ind == 1, MUCE.Ind == 1, EXNEX.Ind == 1, 
                              BHM.Ind == 1, CBHM.Ind == 1, Binom.Test.Ind == 1)]
p = length(mnames)
#### S1: resp.true=c(0.2,0.2,0.2,0.2) ####
resp.true <- resp.cntl
s1.res <- sapply(1:10000, Bask.Sim.Test, k, n, n.ia, resp.true, resp.cntl, resp.trt, 
                 #Tchange
                 lambda, gamma, plo, phi, T = 0.8,
                 r.stg1, n.stg1, alphastar.2stg,r.1stg, alphastar.1stg, n1.s2stg, n.s2stg, r1.s2stg, 
                 r.s2stg, alpha, 
                 
                 J, L, nik, v1, eff.ref.muce, gamma_muce, prec0, prec.eta, prec.ksi, mu.eta0,
                 prec.eta0, mu.ksi0, prec.ksi0, n.burnin, n.chain, n.iter,
                 
                 eff.ref.exnex, nik, Nexch = Nexch, Nmix = Nmix, pMix = pMix, dirichprior = c(1,1,1),
                 dirichprior.check = FALSE, mu.mean = mu.mean, mu.prec = mu.prec,
                 tau.HN.scale = tau.HN.scale, nex.mean = nex.mean, nex.prec = nex.prec,
                 
                 cohortsize = cohortsize, mu.par.BHM, eff.ref.bhm,
                 
                 mu.par.CBHM, eff.ref.CBHM, a = a.par, b = b.par,
                 
                 Simon.Bayes.Ind, Optimal.2stg.Ind, Optimal.1stg.Ind, Simon.2stg.Ind, 
                 Simple.pooling.Ind, MUCE.Ind, EXNEX.Ind, BHM.Ind, CBHM.Ind, Binom.Test.Ind)

s1.res.output = round(summary_sim(s1.res, p, q = 4, mnames, tp = NULL), 3)


#### S2: resp.true=c(0.4,0.4,0.4,0.4) ####
# True probability under the alternative

resp.true <- resp.trt
s2.res <- sapply(1:nsim, Bask.Sim.Test, k, n, n.ia, resp.true, resp.cntl, resp.trt, 
                 
                 lambda, gamma, plo, phi, T = 0.8,
                 r.stg1, n.stg1, alphastar.2stg,r.1stg, alphastar.1stg, n1.s2stg, n.s2stg, r1.s2stg, 
                 r.s2stg,alpha, 
                 
                 J, L, nik, v1, eff.ref.muce, gamma_muce, prec0, prec.eta, prec.ksi, mu.eta0,
                 prec.eta0, mu.ksi0, prec.ksi0, n.burnin, n.chain, n.iter,
                 
                 eff.ref.exnex, nik, Nexch = Nexch, Nmix = Nmix, pMix = pMix, dirichprior = c(1,1,1),
                 dirichprior.check = FALSE, mu.mean = mu.mean, mu.prec = mu.prec,
                 tau.HN.scale = tau.HN.scale, nex.mean = nex.mean, nex.prec = nex.prec,
                 
                 cohortsize = cohortsize, mu.par.BHM, eff.ref.bhm,
                 
                 mu.par.CBHM, eff.ref.CBHM, a = a.par, b = b.par,
                 
                 Simon.Bayes.Ind, Optimal.2stg.Ind, Optimal.1stg.Ind, Simon.2stg.Ind, 
                 Simple.pooling.Ind, MUCE.Ind, EXNEX.Ind, BHM.Ind, CBHM.Ind, Binom.Test.Ind)

s2.res.output = round(summary_sim(s2.res, p, q = 4,mnames, tp = c(1:4)), 3)

#### S3: resp.true=c(0.2,0.2,0.2,0.4) ####
# True probability under the alternative

resp.true <- c(0.5,0.5,0.5,0.65) #c(0.2,0.2,0.2,0.4)
s3.res <- sapply(1:nsim, Bask.Sim.Test, k, n, n.ia, resp.true, resp.cntl, resp.trt, 
                 
                 lambda, gamma, plo, phi, T = 0.8,
                 r.stg1, n.stg1, alphastar.2stg,r.1stg, alphastar.1stg, n1.s2stg, n.s2stg, r1.s2stg, 
                 r.s2stg,alpha, 
                 
                 J, L, nik, v1, eff.ref.muce, gamma_muce, prec0, prec.eta, prec.ksi, mu.eta0,
                 prec.eta0, mu.ksi0, prec.ksi0, n.burnin, n.chain, n.iter,
                 
                 eff.ref.exnex, nik, Nexch = Nexch, Nmix = Nmix, pMix = pMix, dirichprior = c(1,1,1),
                 dirichprior.check = FALSE, mu.mean = mu.mean, mu.prec = mu.prec,
                 tau.HN.scale = tau.HN.scale, nex.mean = nex.mean, nex.prec = nex.prec,
                 
                 cohortsize = cohortsize, mu.par.BHM, eff.ref.bhm,
                 
                 mu.par.CBHM, eff.ref.CBHM, a = a.par, b = b.par,
                 
                 Simon.Bayes.Ind, Optimal.2stg.Ind, Optimal.1stg.Ind, Simon.2stg.Ind, 
                 Simple.pooling.Ind, MUCE.Ind, EXNEX.Ind, BHM.Ind, CBHM.Ind, Binom.Test.Ind)

s3.res.output = round(summary_sim(s3.res, p, q = 4,mnames, tp = 4), 3)

#### S4: resp.true=c(0.2,0.2,0.4,0.4) ####
# True probability under the alternative

resp.true <- c(0.5,0.5,0.65,0.65)
s4.res <- sapply(1:nsim, Bask.Sim.Test, k, n, n.ia, resp.true, resp.cntl, resp.trt, 
                 
                 lambda, gamma, plo, phi, T = 0.8,
                 r.stg1, n.stg1, alphastar.2stg,r.1stg, alphastar.1stg, n1.s2stg, n.s2stg, r1.s2stg, 
                 r.s2stg,alpha, 
                 
                 J, L, nik, v1, eff.ref.muce, gamma_muce, prec0, prec.eta, prec.ksi, mu.eta0,
                 prec.eta0, mu.ksi0, prec.ksi0, n.burnin, n.chain, n.iter,
                 
                 eff.ref.exnex, nik, Nexch = Nexch, Nmix = Nmix, pMix = pMix, dirichprior = c(1,1,1),
                 dirichprior.check = FALSE, mu.mean = mu.mean, mu.prec = mu.prec,
                 tau.HN.scale = tau.HN.scale, nex.mean = nex.mean, nex.prec = nex.prec,
                 
                 cohortsize = cohortsize, mu.par.BHM, eff.ref.bhm,
                 
                 mu.par.CBHM, eff.ref.CBHM, a = a.par, b = b.par,
                 
                 Simon.Bayes.Ind, Optimal.2stg.Ind, Optimal.1stg.Ind, Simon.2stg.Ind, 
                 Simple.pooling.Ind, MUCE.Ind, EXNEX.Ind, BHM.Ind, CBHM.Ind, Binom.Test.Ind)

s4.res.output = round(summary_sim(s4.res, p, q = 4, mnames, tp = c(3, 4)), 3)

#### S5: resp.true=c(0.2,0.4,0.4,0.4) ####
# True probability under the alternative

resp.true <- c(0.5,0.65,0.65,0.65)
s5.res <- sapply(1:nsim, Bask.Sim.Test, k, n, n.ia, resp.true, resp.cntl, resp.trt, 
                 
                 lambda, gamma, plo, phi, T = 0.8,
                 r.stg1, n.stg1, alphastar.2stg,r.1stg, alphastar.1stg, n1.s2stg, n.s2stg, r1.s2stg, 
                 r.s2stg,alpha, 
                 
                 J, L, nik, v1, eff.ref.muce, gamma_muce, prec0, prec.eta, prec.ksi, mu.eta0,
                 prec.eta0, mu.ksi0, prec.ksi0, n.burnin, n.chain, n.iter,
                 
                 eff.ref.exnex, nik, Nexch = Nexch, Nmix = Nmix, pMix = pMix, dirichprior = c(1,1,1),
                 dirichprior.check = FALSE, mu.mean = mu.mean, mu.prec = mu.prec,
                 tau.HN.scale = tau.HN.scale, nex.mean = nex.mean, nex.prec = nex.prec,
                 
                 cohortsize = cohortsize, mu.par.BHM, eff.ref.bhm,
                 
                 mu.par.CBHM, eff.ref.CBHM, a = a.par, b = b.par,
                 
                 Simon.Bayes.Ind, Optimal.2stg.Ind, Optimal.1stg.Ind, Simon.2stg.Ind, 
                 Simple.pooling.Ind, MUCE.Ind, EXNEX.Ind, BHM.Ind, CBHM.Ind, Binom.Test.Ind)

s5.res.output = round(summary_sim(s5.res, p, q = 4,mnames, tp = c(2,3,4)), 3)

#### Saving data ####
save(k, nsim, n, resp.trt, resp.true, n.ia, alpha,
     #Simon's Bayes
     plo, phi, lambda, gamma, T.lst, T.search.result,
     #Optimal 2-stage
     optimal.2stg.params, optimal.2stg.parm, n.stg1, r.stg1, alphastar.2stg,
     #Optimal 1-stage
     optimal.1stg.params, optimal.1stg.parm, r.1stg, alphastar.1stg,
     #Simon 2-stage
     n.s2stg, n1.s2stg, p0.s2stg, p1.s2stg, r1.s2stg.lst, r.s2stg.lst, alpha.simon2stg, sum.2stage,
     sum.2stage.val, sum.2stage.pick, r1.s2stg, r.s2stg,
     nik, n.burnin, n.chain, n.iter, FWER, v1,
     #MUCE,
     J, L, gamma_muce, prec0, prec.eta, prec.ksi, mu.eta0, prec.eta0, mu.ksi0, prec.ksi0, eff.ref.muce,
     #EXNEX
     Nexch, Nmix, pMix, dirichprior.check, tau.HN.scale, prec.eta0, mu.ksi0, prec.ksi0, mu.mean,
     mu.prec, nex.mean, nex.prec, Q.exnex, eff.ref.exnex,
     #BHM
     cohortsize, mu.par.BHM, BHM.results, eff.ref.bhm,
     #CBHM
     var.small, var.big, mu.par.CBHM, a_par_list, b_par_list, eff.ref.CBHM, a.par, b.par,
     # results
     s1.res, s1.res.output, s2.res, s2.res.output, s3.res, s3.res.output, s4.res, s4.res.output,
     s5.res, s5.res.output,
     
     file = paste0("Basket_trial_design_setting_parameters_and_results_set_1_5_cases_20samples_with_ia", ".RData"))

load('Basket_trial_design_setting_parameters_and_results_set_1_5_cases_20samples_with_ia.RData')

write.csv(s1.res.output,"s1_res_output_ia.csv", row.names = TRUE)
write.csv(s2.res.output,"s2_res_output_ia.csv", row.names = TRUE)
write.csv(s3.res.output,"s3_res_output_ia.csv", row.names = TRUE)
write.csv(s4.res.output,"s4_res_output_ia.csv", row.names = TRUE)
write.csv(s5.res.output,"s5_res_output_ia.csv", row.names = TRUE)

#### Plot ####
library(ggplot2)
p = 9 #8
type1.collection <- data.frame(Scenarios = rep(c("0 active: global null"), each = p),
                                Methods = rep(colnames(s1.res.output), 1),
                               type1 = c(unlist(s1.res.output[1,])))
type1.collection$Methods = factor(type1.collection$Methods, levels = type1.collection$Methods)
ggplot(type1.collection, aes(x = Scenarios, y = type1, fill = Methods)) + 
  geom_bar(stat = "identity", position = position_dodge(width=0.8), width = 0.7) + 
  ggtitle("Global type 1 error rate")  + geom_hline(yintercept = 0.1, show.legend = TRUE)+
  xlab("") + ylab("Type 1 error rate")


global.collection <- data.frame(Scenarios = rep(c("1 active", "2 active", 
                                                  "3 active", "4 active: global alternative"), each = p),
                     Methods = rep(colnames(s1.res.output), 4),
                     Global = c(unlist(s3.res.output[1,]), unlist(s4.res.output[1,]),
                                unlist(s5.res.output[1,]), unlist(s2.res.output[1,])))
global.collection$Methods = factor(unique(global.collection$Methods), levels = unique(global.collection$Methods))
ggplot(global.collection, aes(x = Scenarios, y = Global, fill = Methods)) + 
  geom_bar(stat = "identity", position = position_dodge(width=0.8), width = 0.7) + 
  ggtitle("Global power") +
  xlab("4 scenarios") + ylab("Global power")

TP.collection <- data.frame(Scenarios = rep(c("1 active", "2 active", 
                                              "3 active", "4 active: global alternative"), each = p),
                                Methods = rep(colnames(s1.res.output), 4),
                                TP = c(unlist(s3.res.output[6,]), unlist(s4.res.output[6,]),
                                           unlist(s5.res.output[6,]), unlist(s2.res.output[6,])))
TP.collection$Methods = factor(unique(TP.collection$Methods), levels = unique(TP.collection$Methods))
ggplot(TP.collection, aes(x = Scenarios, y = TP, fill = Methods)) + 
  geom_bar(stat = "identity", position = position_dodge(width=0.8), width = 0.7) + 
  ggtitle("True positive")  +
  xlab("5 scenarios") + ylab("True positive")

FP.collection <- data.frame(Scenarios = rep(c("0 active: global null", "1 active", "2 active", 
                                              "3 active"), each = p),
                            Methods = rep(colnames(s1.res.output), 4),
                            FP = c(unlist(s1.res.output[7,]),
                                   unlist(s3.res.output[7,]), unlist(s4.res.output[7,]),
                                   unlist(s5.res.output[7,])))
FP.collection$Methods = factor(unique(FP.collection$Methods), levels = unique(FP.collection$Methods))
ggplot(FP.collection, aes(x = Scenarios, y = FP, fill = Methods)) + 
  geom_bar(stat = "identity", position = position_dodge(width=0.8), width = 0.7) + 
  ggtitle("False positive") +
  xlab("5 scenarios") + ylab("False positive")

EN.collection <- data.frame(Scenarios = rep(c("0 active: global null", "4 active: global alternative", 
                                              "1 active", "2 active", 
                                              "3 active"), each = p),
                            Methods = rep(colnames(s1.res.output), 5),
                            EN = c(unlist(s1.res.output[8,]), unlist(s2.res.output[8,]),
                                   unlist(s3.res.output[8,]), unlist(s4.res.output[8,]),
                                   unlist(s5.res.output[8,])))
EN.collection$Methods = factor(unique(EN.collection$Methods), levels = unique(EN.collection$Methods))
ggplot(EN.collection, aes(x = Scenarios, y = EN, fill = Methods)) + 
  geom_bar(stat = "identity", position = position_dodge(width=0.8), width = 0.7) + 
  ggtitle("Expected sample size") 
