#####Load packages#####
#-------------------------#

library(dplyr)
library(jagsUI)
library(MCMCvis)
library(ggdag)
library(stringr)
library(matrixStats)

#####Load data#####
#-------------------------#
load("squirrelObservations.RData")
load("urbanCovariates.RData")

#####Bundle data for JAGS#####
#-------------------------#

all.data <- list(dis = dis.s, pdn = pdn.log.s, bld = bld.s,
                 tre = tre.s, rod = rod.s, 
                 frg = frg.log.s, prd = prd.s,
                 yg.occ = y.occ.g, ym.occ =y.occ.m, 
                 yg.ct = y.ct.g, ym.ct =y.ct.m, 
                 nsites = nsites,
                 nsurveys.occ = nsurveys.occ, nsurveys.ct = nsurveys.ct, 
                 tmp.occ = temp.avg.occ, tmp.ct = temp.avg.ct)
str(all.data)

#####Fit urbanization and abundance submodels#####
#-------------------------#

sink("fullModel.txt")
cat("

model {

#####URBANIZATION SUBMODEL#####

###population density###
pdn.tau <- pow(pdn.sd, -2)                        #precision for observed variable
pdn.sd ~ dt(0, 1/pow(25, 2), 1)I(0,)              #prior for standard deviation

b.pdn.0 ~ dnorm(0, 0.01)                          #prior for intercept
b.pdn.dis ~ dnorm(0, 0.01)                        #prior for slope

for(i in 1:nsites){                               #regression model
pdn[i] ~ dnorm(pdn.mu[i], pdn.tau)                #stochastic component
pdn.mu[i] <- b.pdn.0 +                            #deterministic component
             b.pdn.dis*dis[i]
pdn.res[i] <- pdn[i] - pdn.mu[i]                  #residuals             
}

###building cover###
bld.tau <- pow(bld.sd, -2)                        #precision for observed variable
bld.sd ~ dt(0, 1/pow(25, 2), 1)I(0,)              #prior for standard deviation

b.bld.0 ~ dnorm(0, 0.01)                          #prior for intercept
b.bld.pdn ~ dnorm(0, 0.01)                        #prior for slope

for(i in 1:nsites){                               #regression model
bld[i] ~ dnorm(bld.mu[i], bld.tau)                #stochastic component
bld.mu[i] <- b.bld.0 +                            #deterministic component
             b.bld.pdn*pdn[i] 
bld.res[i] <- bld[i] - bld.mu[i]                  #residuals
}


###tree area###
tre.tau <- pow(tre.sd, -2)                        #precision for observed variable
tre.sd ~ dt(0, 1/pow(25, 2), 1)I(0,)              #prior for standard deviation

b.tre.0 ~ dnorm(0, 0.01)                          #prior for intercept
b.tre.bld ~ dnorm(0, 0.01)

for(i in 1:nsites){                               #regression model
tre[i] ~ dnorm(tre.mu[i], tre.tau)                #stochastic component
tre.mu[i] <- b.tre.0 +                            #deterministic component
             b.tre.bld*bld[i]
tre.res[i] <- tre[i] - tre.mu[i]                  #residuals
}

###fragmentation###
frg.tau <- pow(frg.sd, -2)                        #precision for observed variable
frg.sd ~ dt(0, 1/pow(25, 2), 1)I(0,)              #prior for standard deviation

b.frg.0 ~ dnorm(0, 0.01)                          #prior for intercept
b.frg.dis ~ dnorm(0, 0.01)                        #prior for slopes
b.frg.pdn ~ dnorm(0, 0.01)
b.frg.tre ~ dnorm(0, 0.01)                        

for(i in 1:nsites){                               #regression model
frg[i] ~ dnorm(frg.mu[i], frg.tau)                #stochastic component
frg.mu[i] <- b.frg.0 +                            #deterministic component
             b.frg.dis*dis[i] +
             b.frg.pdn*pdn[i] +
             b.frg.tre*tre[i] 
frg.res[i] <- frg[i] - frg.mu[i]                  #residuals
}

###predators###

prd.tau <- pow(prd.sd, -2)                        #precision for observed variable
prd.sd ~ dt(0, 1/pow(25, 2), 1)I(0,)              #prior for standard deviation
  
b.prd.0 ~ dnorm(0, 0.01)                          #prior for intercept
b.prd.pdn ~ dnorm(0, 0.01)                        #prior for slope

for(i in 1:nsites){                               #regression model
prd[i] ~ dnorm(prd.mu[i], prd.tau)                #stochastic component
prd.mu[i] <- b.prd.0 +                            #deterministic component
             b.prd.pdn*pdn[i]
prd.res[i] <- prd[i] - prd.mu[i]                  #residuals
}

#####ABUNDANCE SUBMODEL INTEGRATING COUNT AND DETECTION DATA#####
    
    #Priors
    
    for(m in 1:2) {                 #m = morphs, 1 = melanic 2 = gray
    a.0[m] <- logit(mean.p[m])      #detection intercept
    mean.p[m] ~ dunif(0.01, 0.30)   #mean detection probabilty at value = 0 for covariates (mean temperature)
    a.tmp[m] ~ dnorm(0, 0.01)       #linear temp slope
    a.tmp2[m] ~ dnorm(0, 0.01)      #squadratic temp slope
    }
    
    for(m in 1:2) {
    b.0[m] ~ dunif(0,4)             #abundance intercept
    b.dis[m] ~ dnorm(0, 0.01)       #slopes
    b.pdn[m] ~ dnorm(0, 0.01)       
    b.rod[m] ~ dnorm(0, 0.01)       
    b.bld[m] ~ dnorm(0, 0.01)
    b.tre[m] ~ dnorm(0, 0.01)       
    b.frg[m] ~ dnorm(0, 0.01)       
    b.prd[m] ~ dnorm(0, 0.01)
    mean.abu[m] <- exp(b.0[m])      #mean abundance at value 0 for covariates (mean distance)
    }

    #Likelihood
    
    #Ecological process model for abundance
    for(i in 1:nsites) {
    N1[i] ~ dpois(lambda.abu1[i])                     #melanic abundance
    log(lambda.abu1[i]) <-  b.0[1] +                  #expected abundance as a function of covariates
                            b.dis[1]*dis[i] + 
                            b.pdn[1]*pdn[i] + 
                            b.rod[1]*rod[i] +
                            b.bld[1]*bld[i] +
                            b.tre[1]*tre[i] +
                            b.frg[1]*frg[i] +
                            b.prd[1]*prd[i]
    N1.res[i] <- N1[i] - lambda.abu1[i]               #residual abundance
    
    N2[i] ~ dpois(lambda.abu2[i])                     #gray abundance
    log(lambda.abu2[i]) <-  b.0[2] +                  #expected abundance as a function of covariates
                            b.dis[2]*dis[i] + 
                            b.pdn[2]*pdn[i] + 
                            b.rod[2]*rod[i] +
                            b.bld[2]*bld[i] +
                            b.tre[2]*tre[i] +
                            b.frg[2]*frg[i] +
                            b.prd[2]*prd[i]
    N2.res[i] <- N2[i] - lambda.abu2[i]                #residual abundance
    }

    #Observation model for detection probability - occupancy data
    for(i in 1:nsites) {
    for(j in 1:nsurveys.occ) {
    ym.occ[i,j] ~ dbern(pstar.m.occ[i,j])                       #observed detections for melanic morph              
    pstar.m.occ[i,j] <- 1-(1-pdet.m.occ[i,j])^N1[i]             #Pstar = P(detect melanic morph), pdet = ind. detection prob.
    logit(pdet.m.occ[i,j]) <- a.0[1] +                          #ind. detection prob. as a function of survey temp
                              a.tmp[1]*tmp.occ[i,j] +
                              a.tmp2[1]*pow(tmp.occ[i,j], 2)
                                
    yg.occ[i,j] ~ dbern(pstar.g.occ[i,j])                       #observed detections for gray morph         
    pstar.g.occ[i,j] <- 1-(1-pdet.g.occ[i,j])^N2[i]             #Pstar = P(detect gray morph), pdet = ind. detection prob.     
    logit(pdet.g.occ[i,j]) <- a.0[2] +                          #ind. detection prob. as a function of survey temp
                              a.tmp[2]*tmp.occ[i,j] +
                              a.tmp2[2]*pow(tmp.occ[i,j], 2)
    }
    }
    
    #Observation model for detection probability - count data
    for(i in 1:nsites) {
    for(j in 1:nsurveys.ct) {
    ym.ct[i,j] ~ dbinom(pdet.m.ct[i,j], N1[i])                  #observed melanic count; pdet = ind detection prob.
    logit(pdet.m.ct[i,j]) <- a.0[1] +                           #ind. detection prob. as a function of survey temp
                             a.tmp[1]*tmp.ct[i,j] +
                             a.tmp2[1]*pow(tmp.ct[i,j], 2)
    
    yg.ct[i,j] ~ dbinom(pdet.g.ct[i,j], N2[i])                  #observed gray count; pdet = ind detection prob.
    logit(pdet.g.ct[i,j]) <- a.0[2] +                           #ind. detection prob. as a function of survey temp
                             a.tmp[2]*tmp.ct[i,j] +
                             a.tmp2[2]*pow(tmp.ct[i,j], 2)
    }
    }
    
    #Derived quantities
    for(i in 1:nsites) {
    N[i] <- N1[i] + N2[i]                          #derived total abundance
    pm.obs[i] <- N1[i]/ifelse(N[i]==0, 1e6, N[i])  #derived proportion melanic at each site
    }
    
    d.a.0    <- a.0[1] - a.0[2]                    #delta alphas for morphs (melanic - gray)
    d.a.tmp  <- a.tmp[1] - a.tmp[2]
    d.a.tmp2 <- a.tmp2[1] - a.tmp2[2]

    d.b.0 <- b.0[1] - b.0[2]                       #delta betas for morphs (melanic - gray)
    d.b.dis <- b.dis[1] - b.dis[2]
    d.b.pdn <- b.pdn[1] - b.pdn[2]
    d.b.rod <- b.rod[1] - b.rod[2]
    d.b.bld <- b.bld[1] - b.bld[2]
    d.b.tre <- b.tre[1] - b.tre[2]
    d.b.frg <- b.frg[1] - b.frg[2]
    d.b.prd <- b.prd[1] - b.prd[2]

}

    ", fill=TRUE)
sink()

##Initial values----
y.g.max.ct <- apply(y.ct.g, 1, max, na.rm=TRUE, warn=F)
y.g.max.ct[!is.finite(y.g.max.ct)] <- 0
y.m.max.ct <- apply(y.ct.m, 1, max, na.rm=TRUE, warn=F)
y.m.max.ct[!is.finite(y.m.max.ct)] <- 0
y.g.max.occ <- apply(y.occ.g, 1, max, na.rm=TRUE, warn=F)
y.g.max.occ[!is.finite(y.g.max.occ)] <- 0
y.m.max.occ <- apply(y.occ.m, 1, max, na.rm=TRUE, warn=F)
y.m.max.occ[!is.finite(y.m.max.occ)] <- 0

inits <- function(){list(pdn.sd = 1, bld.sd = 1, 
                         imp.sd = 1, tre.sd = 1, frg.sd =1, prd.sd = 1,
                         N1 = apply(cbind(y.m.max.occ,y.m.max.ct), 1, max), 
                         N2 = apply(cbind(y.g.max.occ,y.g.max.ct), 1, max))}

##Parameters monitored
params <- c("pdn.sd", "pdn.mu", "pdn.res", 
            "bld.sd", "bld.mu", "bld.res", 
            "tre.sd", "tre.mu", "tre.res", 
            "frg.sd", "frg.mu", "frg.res", 
            "prd.sd", "prd.mu", "prd.res", 
            "b.pdn.0", "b.pdn.dis",
            "b.bld.0", "b.bld.pdn",
            "b.tre.0", "b.tre.bld",
            "b.frg.0", "b.frg.dis", "b.frg.pdn", "b.frg.tre",
            "b.prd.0", "b.prd.pdn",
            "mean.p", "a.0", "a.tmp", "a.tmp2", 
            "d.a.0", "d.a.tmp", "d.a.tmp2", 
            "mean.abu", "b.0", "b.dis", "b.pdn", "b.bld", "b.rod", "b.tre", "b.frg", "b.prd",
            "d.b.0", "d.b.dis", "d.b.pdn", "d.b.bld", "d.b.rod", "d.b.tre", "d.b.frg", "d.b.prd", 
            "N", "N1", "N2", "pm.obs", "lambda.abu1", "lambda.abu2",
            "N1.res", "N2.res")

##MCMC
ni <- 4000 ; nt <- 3 ;  nb <- 1000 ; nc <- 3 #6 min

m1 <- jags(all.data, inits, params, "fullModel.txt", 
           n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=20000,
           parallel=TRUE, n.cores=3)

#Parameter estimates with Rhat
ModSummary <- MCMCsummary(m1, 
            params = c("b.pdn.0", "b.pdn.dis",
                       "b.bld.0", "b.bld.pdn",
                       "b.tre.0", "b.tre.bld", 
                       "b.frg.0", "b.frg.dis", "b.frg.pdn", "b.frg.tre",
                       "b.prd.0", "b.prd.pdn",
                       "mean.p", "a.0", "a.tmp", "a.tmp2", 
                       "d.a.0", "d.a.tmp", "d.a.tmp2", 
                       "mean.abu", "b.0", "b.dis", "b.pdn", "b.bld", "b.rod", "b.tre", "b.frg", "b.prd",
                       "d.b.0", "d.b.dis", "d.b.pdn", "d.b.bld", "d.b.rod", "d.b.tre", "d.b.frg", "d.b.prd"),
            Rhat = TRUE,
            n.eff = TRUE,
            probs = c(0.025, 0.5, 0.975),
            round = 4)

#####R2#####
#-------------------------#

r2.end  <- matrix(NA, nrow=7, ncol=5)
r2.end[,1] <- c("pdn", "bld", "tre", "frg", "prd", "mel", "gry")
colnames(r2.end) <- c("variable", "mean", "sd", "lcl.95", "ucl.95")

#population density
r2.sims.pdn <- apply(m1$sims.list$pdn.mu, 1, var)/(apply(m1$sims.list$pdn.mu, 1, var) + apply(m1$sims.list$pdn.res, 1, var))
r2.end[1,2] <- mean(r2.sims.pdn)
r2.end[1,3] <- sd(r2.sims.pdn)
r2.end[1,4] <- quantile(r2.sims.pdn, prob=c(0.025))
r2.end[1,5] <- quantile(r2.sims.pdn, prob=0.975)

#building density
r2.sims.bld <- apply(m1$sims.list$bld.mu, 1, var)/(apply(m1$sims.list$bld.mu, 1, var) + apply(m1$sims.list$bld.res, 1, var))
r2.end[2,2] <- mean(r2.sims.bld)
r2.end[2,3] <- sd(r2.sims.bld)
r2.end[2,4] <- quantile(r2.sims.bld, prob=c(0.025))
r2.end[2,5] <- quantile(r2.sims.bld, prob=0.975)

#tree
r2.sims.tre <- apply(m1$sims.list$tre.mu, 1, var)/(apply(m1$sims.list$tre.mu, 1, var) + apply(m1$sims.list$tre.res, 1, var))
r2.end[3,2] <- mean(r2.sims.tre)
r2.end[3,3] <- sd(r2.sims.tre)
r2.end[3,4] <- quantile(r2.sims.tre, prob=c(0.025))
r2.end[3,5] <- quantile(r2.sims.tre, prob=0.975)

#frag
r2.sims.frg <- apply(m1$sims.list$frg.mu, 1, var)/(apply(m1$sims.list$frg.mu, 1, var) + apply(m1$sims.list$frg.res, 1, var))
r2.end[4,2] <- mean(r2.sims.frg)
r2.end[4,3] <- sd(r2.sims.frg)
r2.end[4,4] <- quantile(r2.sims.frg, prob=c(0.025))
r2.end[4,5] <- quantile(r2.sims.frg, prob=0.975)

#predators
r2.sims.prd <- apply(m1$sims.list$prd.mu, 1, var)/(apply(m1$sims.list$prd.mu, 1, var) + apply(m1$sims.list$prd.res, 1, var))
r2.end[5,2] <- mean(r2.sims.prd)
r2.end[5,3] <- sd(r2.sims.prd)
r2.end[5,4] <- quantile(r2.sims.prd, prob=c(0.025))
r2.end[5,5] <- quantile(r2.sims.prd, prob=0.975)

#mel
r2.sims.n1 <- apply(m1$sims.list$lambda.abu1, 1, var)/(apply(m1$sims.list$lambda.abu1, 1, var) + apply(m1$sims.list$N1.res, 1, var))
r2.end[6,2] <- mean(r2.sims.n1)
r2.end[6,3] <- sd(r2.sims.n1)
r2.end[6,4] <- quantile(r2.sims.n1, prob=c(0.025))
r2.end[6,5] <- quantile(r2.sims.n1, prob=0.975)

#gry
r2.sims.n2 <- apply(m1$sims.list$lambda.abu2, 1, var)/(apply(m1$sims.list$lambda.abu2, 1, var) + apply(m1$sims.list$N2.res, 1, var))
r2.end[7,2] <- mean(r2.sims.n2)
r2.end[7,3] <- sd(r2.sims.n2)
r2.end[7,4] <- quantile(r2.sims.n2, prob=c(0.025))
r2.end[7,5] <- quantile(r2.sims.n2, prob=0.975)

r2.end <- as.data.frame(r2.end)
r2.end[2:5] <- sapply(r2.end[2:5],as.numeric)

#####Standardized coefficients for melanic abundance#####
#-------------------------#

std.beta.n1 <- matrix(NA, nrow=7, ncol=nrow(m1$sims.list$N))
rownames(std.beta.n1) <- c("dis", "pdn", "bld", "rod", "tre", "frg", "prd")

for (i in 1:nrow(m1$sims.list$N1)) {

  yhat <- m1$sims.list$lambda.abu1[i,] 
  r2 <- r2.sims.n1[i]
  sd.y <- sqrt(var(log(yhat))/r2) 
  
  #dist effect
  std.beta.n1[1,i] <- m1$sims.list$b.dis[i,1] * (sd(all.data$dis)/sd.y)  
  
  #pdn effect
  std.beta.n1[2,i] <- m1$sims.list$b.pdn[i,1] * (sd(all.data$pdn)/sd.y)  
  
  #building effect
  std.beta.n1[3,i] <- m1$sims.list$b.bld[i,1] * (sd(all.data$bld)/sd.y)  
  
  #road effect
  std.beta.n1[4,i] <- m1$sims.list$b.rod[i,1] * (sd(all.data$rod)/sd.y)  

  #tree effect
  std.beta.n1[5,i] <- m1$sims.list$b.tre[i,1] * (sd(all.data$tre)/sd.y)  
  
  #frag effect
  std.beta.n1[6,i] <- m1$sims.list$b.frg[i,1] * (sd(all.data$frg)/sd.y)  
  
  #predator effect
  std.beta.n1[7,i] <- m1$sims.list$b.prd[i,1] * (sd(all.data$prd)/sd.y)  
  
}

std.beta.n1.sum <- matrix(NA, nrow=7, ncol=4)
std.beta.n1.sum[,1] <- apply(std.beta.n1, 1, mean)
std.beta.n1.sum[,2] <- apply(std.beta.n1, 1, sd)
std.beta.n1.sum[,3] <- apply(std.beta.n1, 1, function(x) quantile(x, prob=c(0.025)))
std.beta.n1.sum[,4] <- apply(std.beta.n1, 1, function(x) quantile(x, prob=c(0.975)))
colnames(std.beta.n1.sum) <- c("mean", "psd", "lcl.95", "ucl.95")
rownames(std.beta.n1.sum) <- c("dis", "pdn", "bld", "rod", "tre", "frg", "prd")

#####Standardized coefficients for gray morph#####
#-------------------------#

std.beta.n2 <- matrix(NA, nrow=7, ncol=nrow(m1$sims.list$N))
rownames(std.beta.n2) <- c("dis", "pdn", "bld", "rod", "tre", "frg", "prd")

for (i in 1:nrow(m1$sims.list$N2)) {
  
  yhat <- m1$sims.list$lambda.abu2[i,] 
  r2 <- r2.sims.n2[i]
  sd.y <- sqrt(var(log(yhat))/r2) 
  
  #dist effect
  std.beta.n2[1,i] <- m1$sims.list$b.dis[i,2] * (sd(all.data$dis)/sd.y)  
  
  #pdn effect
  std.beta.n2[2,i] <- m1$sims.list$b.pdn[i,2] * (sd(all.data$pdn)/sd.y) 
  
  #bldtlement effect
  std.beta.n2[3,i] <- m1$sims.list$b.bld[i,2] * (sd(all.data$bld)/sd.y)  
  
  #road effect
  std.beta.n2[4,i] <- m1$sims.list$b.rod[i,2] * (sd(all.data$rod)/sd.y) 
  
  #tree effect
  std.beta.n2[5,i] <- m1$sims.list$b.tre[i,2] * (sd(all.data$tre)/sd.y) 
  
  #frag effect
  std.beta.n2[6,i] <- m1$sims.list$b.frg[i,2] * (sd(all.data$frg)/sd.y) 
  
  #predator effect
  std.beta.n2[7,i] <- m1$sims.list$b.prd[i,2] * (sd(all.data$prd)/sd.y) 
  
}

std.beta.n2.sum <- matrix(NA, nrow=7, ncol=4)
std.beta.n2.sum[,1] <- apply(std.beta.n2, 1, mean)
std.beta.n2.sum[,2] <- apply(std.beta.n2, 1, sd)
std.beta.n2.sum[,3] <- apply(std.beta.n2, 1, function(x) quantile(x, prob=c(0.025)))
std.beta.n2.sum[,4] <- apply(std.beta.n2, 1, function(x) quantile(x, prob=c(0.975)))
colnames(std.beta.n2.sum) <- c("mean", "psd", "lcl.95", "ucl.95")
rownames(std.beta.n2.sum) <- c("dis", "pdn", "bld", "rod", "tre", "frg", "prd")

#####Direct, indirect, and total effects#####
#-------------------------#

#specify DAG
dagified <- dagify(
  pdn ~ dis,
  bld ~ pdn,
  tre ~ bld,
  frg ~ dis + pdn + tre,
  prd ~ pdn,
  mel ~ dis + pdn + bld + rod + tre + frg + prd,
  gry ~ dis + pdn + bld + rod + tre + frg + prd
)

#path effects from urbanization submodel
dis.pdn <- m1$sims.list$b.pdn.dis 
dis.frg <- m1$sims.list$b.frg.dis
pdn.bld <-   m1$sims.list$b.bld.pdn
pdn.frg <-   m1$sims.list$b.frg.pdn
pdn.prd <-   m1$sims.list$b.prd.pdn
bld.tre <-  m1$sims.list$b.tre.bld
tre.frg <- m1$sims.list$b.frg.tre

#direct effects on melanic abundance
dis.mel <- std.beta.n1["dis",]
pdn.mel <- std.beta.n1["pdn",]
bld.mel <- std.beta.n1["bld",]
rod.mel <- std.beta.n1["rod",]
tre.mel <- std.beta.n1["tre",]
frg.mel <- std.beta.n1["frg",]
prd.mel <- std.beta.n1["prd",] 

#direct effects on gray abundance
dis.gry <- std.beta.n2["dis",]
pdn.gry <- std.beta.n2["pdn",]
bld.gry <- std.beta.n2["bld",]
rod.gry <- std.beta.n2["rod",]
tre.gry <- std.beta.n2["tre",]
frg.gry <- std.beta.n2["frg",]
prd.gry <- std.beta.n2["prd",] 

#####create a table for direct, indirect, and total effects - melanic abundance####
effects.n1 <- matrix(NA, nrow=0, ncol = 7)
colnames(effects.n1) <- c("response", "predictor", "effect", "mean", "sd", "95% lcl", "95% ucl")

predictors <- c("dis", "pdn", "bld", "rod", "tre", "frg", "prd")
response <- "mel"

dir.effect.n1 <- matrix(NA, nrow=3000, ncol=length(predictors))
colnames(dir.effect.n1) <- c("dis", "pdn", "bld", "rod", "tre", "frg", "prd")

ind.effect.n1 <- matrix(NA, nrow=3000, ncol=length(predictors))
colnames(ind.effect.n1) <- c("dis", "pdn", "bld", "rod", "tre", "frg", "prd")

tot.effect.n1 <- matrix(NA, nrow=3000, ncol=length(predictors))
colnames(tot.effect.n1) <- c("dis", "pdn", "bld", "rod", "tre", "frg", "prd")

for(k in 1:length(predictors)) {
  
  #list of all paths
  paths.x <- paths(dagified, from = predictors[k], to = response, directed=TRUE)["paths"] 
  paths.x <- as.character(paths.x$paths)
  paths.x <- paths.x[order(nchar(paths.x), paths.x)]
  
  #is there a direct effect?
  is.direct <- nrow(str_locate_all(paths.x[1], ' -> ')[[1]])==1 
  
  #number of indirect paths
  if (is.direct == FALSE) 
  {n.ind.paths <- length(paths.x)
  start <- 1
  } else 
  {n.ind.paths <- length(paths.x) - 1
  start <- 2
  }          
  
  #are there indirect paths?
  is.indirect <- n.ind.paths>0
  
  #quantify all indirect paths
  if (is.indirect ==TRUE) {
    ind.effects <- matrix(NA, nrow=3000, ncol=n.ind.paths)        
    for(i in start:length(paths.x)) {
      path <- paths.x[i]                                        
      a <- str_locate_all(path, ' -> ')                                     
      connections <- nrow(a[[1]])
      betas <- matrix(NA, nrow=3000, ncol=connections)
      for(j in 1:connections){
        x <- substring(path, a[[1]][j,]["start"]-3, a[[1]][j,]["start"]-1)
        y <- substring(path, a[[1]][j,]["end"]+1, a[[1]][j,]["end"]+3)
        betas[,j] <- eval(parse(text=paste0(x, ".", y))) 
      }
      if (is.direct == FALSE) {ind.effects[,(i)] <- rowProds(betas)
      } else {ind.effects[,(i-1)] <- rowProds(betas)
      }
    }
  } else {ind.effects <- rep(0, 3000)}
  
  #quantify direct effect
  dir.effect.n1[,k] <- if (is.direct == TRUE) {eval(parse(text = str_replace_all(paths.x[1], " -> ", "."))) 
  } else {rep(0, 3000)}
  
  #quantify indirect effect
  ind.effect.n1[,k] <- if (is.indirect == TRUE) {rowSums(ind.effects)} else {ind.effects}
  
  #quantify total effect
  tot.effect.n1[,k] <- dir.effect.n1[,k] + ind.effect.n1[,k]
  
  #summarize results and add to table
  effects.temp <- matrix(NA, nrow=3, ncol = 7)
  
  effects.temp[1,] <-  cbind(response, predictors[k], "direct", 
                             mean(dir.effect.n1[,k]), 
                             sd(dir.effect.n1[,k]), 
                             quantile(dir.effect.n1[,k], prob=c(0.025)), 
                             quantile(dir.effect.n1[,k], prob=c(0.975)))
  effects.temp[2,] <- cbind(response, predictors[k], "indirect", 
                            mean(ind.effect.n1[,k]), 
                            sd(ind.effect.n1[,k]), 
                            quantile(ind.effect.n1[,k], prob=c(0.025)), 
                            quantile(ind.effect.n1[,k], prob=c(0.975)))
  effects.temp[3,] <- cbind(response, predictors[k], "total",
                            mean(tot.effect.n1[,k]), 
                            sd(tot.effect.n1[,k]), 
                            quantile(tot.effect.n1[,k], prob=c(0.025)), 
                            quantile(tot.effect.n1[,k], prob=c(0.975)))
  
  effects.n1 <- rbind(effects.n1, effects.temp)
  
}

effects.n1 <- as.data.frame(effects.n1)
effects.n1[4:7] <- sapply(effects.n1[4:7],as.numeric)
print(effects.n1)


#####create a table for direct, indirect, and total effects - gray abundance####

effects.n2 <- matrix(NA, nrow=0, ncol = 7)
colnames(effects.n2) <- c("response", "predictor", "effect", "mean", "sd", "95% lcl", "95% ucl")

predictors <- c("dis", "pdn", "bld", "rod", "tre", "frg", "prd")
response <- "gry"

dir.effect.n2 <- matrix(NA, nrow=3000, ncol=length(predictors))
colnames(dir.effect.n2) <- c("dis", "pdn", "bld", "rod", "tre", "frg", "prd")

ind.effect.n2 <- matrix(NA, nrow=3000, ncol=length(predictors))
colnames(ind.effect.n2) <- c("dis", "pdn", "bld", "rod", "tre", "frg", "prd")

tot.effect.n2 <- matrix(NA, nrow=3000, ncol=length(predictors))
colnames(tot.effect.n2) <- c("dis", "pdn", "bld", "rod", "tre", "frg", "prd")

for(k in 1:length(predictors)) {
  
  #list of all paths
  paths.x <- paths(dagified, from = predictors[k], to = response, directed=TRUE)["paths"] 
  paths.x <- as.character(paths.x$paths)
  paths.x <- paths.x[order(nchar(paths.x), paths.x)]
  
  #is there a direct effect?
  is.direct <- nrow(str_locate_all(paths.x[1], ' -> ')[[1]])==1 
  
  #number of indirect paths
  if (is.direct == FALSE) 
  {n.ind.paths <- length(paths.x)
  start <- 1
  } else 
  {n.ind.paths <- length(paths.x) - 1
  start <- 2
  }          
  
  #are there indirect paths?
  is.indirect <- n.ind.paths>0
  
  #quantify all indirect paths
  if (is.indirect ==TRUE) {
    ind.effects <- matrix(NA, nrow=3000, ncol=n.ind.paths)        
    for(i in start:length(paths.x)) {
      path <- paths.x[i]                                        
      a <- str_locate_all(path, ' -> ')                                     
      connections <- nrow(a[[1]])
      betas <- matrix(NA, nrow=3000, ncol=connections)
      for(j in 1:connections){
        x <- substring(path, a[[1]][j,]["start"]-3, a[[1]][j,]["start"]-1)
        y <- substring(path, a[[1]][j,]["end"]+1, a[[1]][j,]["end"]+3)
        betas[,j] <- eval(parse(text=paste0(x, ".", y))) 
      }
      if (is.direct == FALSE) {ind.effects[,(i)] <- rowProds(betas)
      } else {ind.effects[,(i-1)] <- rowProds(betas)
      }
    }
  } else {ind.effects <- rep(0, 3000)}
  
  #quantify direct effect
  dir.effect.n2[,k] <- if (is.direct == TRUE) {eval(parse(text = str_replace_all(paths.x[1], " -> ", "."))) 
  } else {rep(0, 3000)}
  
  #quantify indirect effect
  ind.effect.n2[,k] <- if (is.indirect == TRUE) {rowSums(ind.effects)} else {ind.effects}
  
  #quantify total effect
  tot.effect.n2[,k] <- dir.effect.n2[,k] + ind.effect.n2[,k]
  
  #summarize results and add to table
  effects.temp <- matrix(NA, nrow=3, ncol = 7)
  
  effects.temp[1,] <-  cbind(response, predictors[k], "direct", 
                             mean(dir.effect.n2[,k]), 
                             sd(dir.effect.n2[,k]), 
                             quantile(dir.effect.n2[,k], prob=c(0.025)), 
                             quantile(dir.effect.n2[,k], prob=c(0.975)))
  effects.temp[2,] <- cbind(response, predictors[k], "indirect", 
                            mean(ind.effect.n2[,k]), 
                            sd(ind.effect.n2[,k]), 
                            quantile(ind.effect.n2[,k], prob=c(0.025)), 
                            quantile(ind.effect.n2[,k], prob=c(0.975)))
  effects.temp[3,] <- cbind(response, predictors[k], "total",
                            mean(tot.effect.n2[,k]), 
                            sd(tot.effect.n2[,k]), 
                            quantile(tot.effect.n2[,k], prob=c(0.025)), 
                            quantile(tot.effect.n2[,k], prob=c(0.975)))
  
  effects.n2 <- rbind(effects.n2, effects.temp)
  
}

effects.n2 <- as.data.frame(effects.n2)
effects.n2[4:7] <- sapply(effects.n2[4:7],as.numeric)
print(effects.n2)

#####create a table for direct, indirect, and total effects - predators####
effects.prd <- matrix(NA, nrow=0, ncol = 7)
colnames(effects.prd) <- c("response", "predictor", "effect", "mean", "sd", "95% lcl", "95% ucl")

predictors <- c("dis", "pdn")
response <- "prd"

for(k in 1:length(predictors)) {
  
  #list of all paths
  paths.x <- paths(dagified, from = predictors[k], to = response, directed=TRUE)["paths"] 
  paths.x <- as.character(paths.x$paths)
  paths.x <- paths.x[order(nchar(paths.x), paths.x)]
  
  #is there a direct effect?
  is.direct <- nrow(str_locate_all(paths.x[1], ' -> ')[[1]])==1 
  
  #number of indirect paths
  if (is.direct == FALSE) 
  {n.ind.paths <- length(paths.x)
  start <- 1
  } else 
  {n.ind.paths <- length(paths.x) - 1
  start <- 2
  }          
  
  #are there indirect paths?
  is.indirect <- n.ind.paths>0
  
  #quantify all indirect paths
  if (is.indirect ==TRUE) {
    ind.effects <- matrix(NA, nrow=3000, ncol=n.ind.paths)        
    for(i in start:length(paths.x)) {
      path <- paths.x[i]                                        
      a <- str_locate_all(path, ' -> ')                                     
      connections <- nrow(a[[1]])
      betas <- matrix(NA, nrow=3000, ncol=connections)
      for(j in 1:connections){
        x <- substring(path, a[[1]][j,]["start"]-3, a[[1]][j,]["start"]-1)
        y <- substring(path, a[[1]][j,]["end"]+1, a[[1]][j,]["end"]+3)
        betas[,j] <- eval(parse(text=paste0(x, ".", y))) 
      }
      if (is.direct == FALSE) {ind.effects[,(i)] <- rowProds(betas)
      } else {ind.effects[,(i-1)] <- rowProds(betas)
      }
    }
  } else {ind.effects <- rep(0, 3000)}
  
  #quantify direct effect
  dir.effect <- if (is.direct == TRUE) {eval(parse(text = str_replace_all(paths.x[1], " -> ", "."))) 
  } else {rep(0, 3000)}
  
  #quantify indirect effect
  ind.effect <- if (is.indirect == TRUE) {rowSums(ind.effects)} else {ind.effects}
  
  #quantify total effect
  tot.effect <- dir.effect + ind.effect
  
  #summarize results and add to table
  effects.temp <- matrix(NA, nrow=3, ncol = 7)
  
  effects.temp[1,] <-  cbind(response, predictors[k], "direct", 
                             mean(dir.effect), 
                             sd(dir.effect), 
                             quantile(dir.effect, prob=c(0.025)), 
                             quantile(dir.effect, prob=c(0.975)))
  effects.temp[2,] <- cbind(response, predictors[k], "indirect", 
                            mean(ind.effect), 
                            sd(ind.effect), 
                            quantile(ind.effect, prob=c(0.025)), 
                            quantile(ind.effect, prob=c(0.975)))
  effects.temp[3,] <- cbind(response, predictors[k], "total",
                            mean(tot.effect), 
                            sd(tot.effect), 
                            quantile(tot.effect, prob=c(0.025)), 
                            quantile(tot.effect, prob=c(0.975)))
  
  effects.prd <- rbind(effects.prd, effects.temp)
  
}

effects.prd <- as.data.frame(effects.prd)
effects.prd[4:7] <- sapply(effects.prd[4:7],as.numeric)
print(effects.prd)

#####create a table for direct, indirect, and total effects - fragmentation####
effects.frg <- matrix(NA, nrow=0, ncol = 7)
colnames(effects.frg) <- c("response", "predictor", "effect", "mean", "sd", "95% lcl", "95% ucl")

predictors <- c("dis", "pdn", "bld", "tre")
response <- "frg"

for(k in 1:length(predictors)) {
  
  #list of all paths
  paths.x <- paths(dagified, from = predictors[k], to = response, directed=TRUE)["paths"] 
  paths.x <- as.character(paths.x$paths)
  paths.x <- paths.x[order(nchar(paths.x), paths.x)]
  
  #is there a direct effect?
  is.direct <- nrow(str_locate_all(paths.x[1], ' -> ')[[1]])==1 
  
  #number of indirect paths
  if (is.direct == FALSE) 
  {n.ind.paths <- length(paths.x)
  start <- 1
  } else 
  {n.ind.paths <- length(paths.x) - 1
  start <- 2
  }          
  
  #are there indirect paths?
  is.indirect <- n.ind.paths>0
  
  #quantify all indirect paths
  if (is.indirect ==TRUE) {
    ind.effects <- matrix(NA, nrow=3000, ncol=n.ind.paths)        
    for(i in start:length(paths.x)) {
      path <- paths.x[i]                                        
      a <- str_locate_all(path, ' -> ')                                     
      connections <- nrow(a[[1]])
      betas <- matrix(NA, nrow=3000, ncol=connections)
      for(j in 1:connections){
        x <- substring(path, a[[1]][j,]["start"]-3, a[[1]][j,]["start"]-1)
        y <- substring(path, a[[1]][j,]["end"]+1, a[[1]][j,]["end"]+3)
        betas[,j] <- eval(parse(text=paste0(x, ".", y))) 
      }
      if (is.direct == FALSE) {ind.effects[,(i)] <- rowProds(betas)
      } else {ind.effects[,(i-1)] <- rowProds(betas)
      }
    }
  } else {ind.effects <- rep(0, 3000)}
  
  #quantify direct effect
  dir.effect <- if (is.direct == TRUE) {eval(parse(text = str_replace_all(paths.x[1], " -> ", "."))) 
  } else {rep(0, 3000)}
  
  #quantify indirect effect
  ind.effect <- if (is.indirect == TRUE) {rowSums(ind.effects)} else {ind.effects}
  
  #quantify total effect
  tot.effect <- dir.effect + ind.effect
  
  #summarize results and add to table
  effects.temp <- matrix(NA, nrow=3, ncol = 7)
  
  effects.temp[1,] <-  cbind(response, predictors[k], "direct", 
                             mean(dir.effect), 
                             sd(dir.effect), 
                             quantile(dir.effect, prob=c(0.025)), 
                             quantile(dir.effect, prob=c(0.975)))
  effects.temp[2,] <- cbind(response, predictors[k], "indirect", 
                            mean(ind.effect), 
                            sd(ind.effect), 
                            quantile(ind.effect, prob=c(0.025)), 
                            quantile(ind.effect, prob=c(0.975)))
  effects.temp[3,] <- cbind(response, predictors[k], "total",
                            mean(tot.effect), 
                            sd(tot.effect), 
                            quantile(tot.effect, prob=c(0.025)), 
                            quantile(tot.effect, prob=c(0.975)))
  
  effects.frg <- rbind(effects.frg, effects.temp)
  
}

effects.frg <- as.data.frame(effects.frg)
effects.frg[4:7] <- sapply(effects.frg[4:7],as.numeric)
print(effects.frg)

#####create a table for direct, indirect, and total effects - tree cover####
effects.tre <- matrix(NA, nrow=0, ncol = 7)
colnames(effects.tre) <- c("response", "predictor", "effect", "mean", "sd", "95% lcl", "95% ucl")

predictors <- c("dis", "pdn", "bld")
response <- "tre"

for(k in 1:length(predictors)) {
  
  #list of all paths
  paths.x <- paths(dagified, from = predictors[k], to = response, directed=TRUE)["paths"] 
  paths.x <- as.character(paths.x$paths)
  paths.x <- paths.x[order(nchar(paths.x), paths.x)]
  
  #is there a direct effect?
  is.direct <- nrow(str_locate_all(paths.x[1], ' -> ')[[1]])==1 
  
  #number of indirect paths
  if (is.direct == FALSE) 
  {n.ind.paths <- length(paths.x)
  start <- 1
  } else 
  {n.ind.paths <- length(paths.x) - 1
  start <- 2
  }          
  
  #are there indirect paths?
  is.indirect <- n.ind.paths>0
  
  #quantify all indirect paths
  if (is.indirect ==TRUE) {
    ind.effects <- matrix(NA, nrow=3000, ncol=n.ind.paths)        
    for(i in start:length(paths.x)) {
      path <- paths.x[i]                                        
      a <- str_locate_all(path, ' -> ')                                     
      connections <- nrow(a[[1]])
      betas <- matrix(NA, nrow=3000, ncol=connections)
      for(j in 1:connections){
        x <- substring(path, a[[1]][j,]["start"]-3, a[[1]][j,]["start"]-1)
        y <- substring(path, a[[1]][j,]["end"]+1, a[[1]][j,]["end"]+3)
        betas[,j] <- eval(parse(text=paste0(x, ".", y))) 
      }
      if (is.direct == FALSE) {ind.effects[,(i)] <- rowProds(betas)
      } else {ind.effects[,(i-1)] <- rowProds(betas)
      }
    }
  } else {ind.effects <- rep(0, 3000)}
  
  #quantify direct effect
  dir.effect <- if (is.direct == TRUE) {eval(parse(text = str_replace_all(paths.x[1], " -> ", "."))) 
  } else {rep(0, 3000)}
  
  #quantify indirect effect
  ind.effect <- if (is.indirect == TRUE) {rowSums(ind.effects)} else {ind.effects}
  
  #quantify total effect
  tot.effect <- dir.effect + ind.effect
  
  #summarize results and add to table
  effects.temp <- matrix(NA, nrow=3, ncol = 7)
  
  effects.temp[1,] <-  cbind(response, predictors[k], "direct", 
                             mean(dir.effect), 
                             sd(dir.effect), 
                             quantile(dir.effect, prob=c(0.025)), 
                             quantile(dir.effect, prob=c(0.975)))
  effects.temp[2,] <- cbind(response, predictors[k], "indirect", 
                            mean(ind.effect), 
                            sd(ind.effect), 
                            quantile(ind.effect, prob=c(0.025)), 
                            quantile(ind.effect, prob=c(0.975)))
  effects.temp[3,] <- cbind(response, predictors[k], "total",
                            mean(tot.effect), 
                            sd(tot.effect), 
                            quantile(tot.effect, prob=c(0.025)), 
                            quantile(tot.effect, prob=c(0.975)))
  
  effects.tre <- rbind(effects.tre, effects.temp)
  
}

effects.tre <- as.data.frame(effects.tre)
effects.tre[4:7] <- sapply(effects.tre[4:7],as.numeric)
print(effects.tre)


#####create a table for direct, indirect, and total effects - building density####
effects.bld <- matrix(NA, nrow=0, ncol = 7)
colnames(effects.bld) <- c("response", "predictor", "effect", "mean", "sd", "95% lcl", "95% ucl")

predictors <- c("dis", "pdn")
response <- "bld"

for(k in 1:length(predictors)) {
  
  #list of all paths
  paths.x <- paths(dagified, from = predictors[k], to = response, directed=TRUE)["paths"] 
  paths.x <- as.character(paths.x$paths)
  paths.x <- paths.x[order(nchar(paths.x), paths.x)]
  
  #is there a direct effect?
  is.direct <- nrow(str_locate_all(paths.x[1], ' -> ')[[1]])==1 
  
  #number of indirect paths
  if (is.direct == FALSE) 
  {n.ind.paths <- length(paths.x)
  start <- 1
  } else 
  {n.ind.paths <- length(paths.x) - 1
  start <- 2
  }          
  
  #are there indirect paths?
  is.indirect <- n.ind.paths>0
  
  #quantify all indirect paths
  if (is.indirect ==TRUE) {
    ind.effects <- matrix(NA, nrow=3000, ncol=n.ind.paths)        
    for(i in start:length(paths.x)) {
      path <- paths.x[i]                                        
      a <- str_locate_all(path, ' -> ')                                     
      connections <- nrow(a[[1]])
      betas <- matrix(NA, nrow=3000, ncol=connections)
      for(j in 1:connections){
        x <- substring(path, a[[1]][j,]["start"]-3, a[[1]][j,]["start"]-1)
        y <- substring(path, a[[1]][j,]["end"]+1, a[[1]][j,]["end"]+3)
        betas[,j] <- eval(parse(text=paste0(x, ".", y))) 
      }
      if (is.direct == FALSE) {ind.effects[,(i)] <- rowProds(betas)
      } else {ind.effects[,(i-1)] <- rowProds(betas)
      }
    }
  } else {ind.effects <- rep(0, 3000)}
  
  #quantify direct effect
  dir.effect <- if (is.direct == TRUE) {eval(parse(text = str_replace_all(paths.x[1], " -> ", "."))) 
  } else {rep(0, 3000)}
  
  #quantify indirect effect
  ind.effect <- if (is.indirect == TRUE) {rowSums(ind.effects)} else {ind.effects}
  
  #quantify total effect
  tot.effect <- dir.effect + ind.effect
  
  #summarize results and add to table
  effects.temp <- matrix(NA, nrow=3, ncol = 7)
  
  effects.temp[1,] <-  cbind(response, predictors[k], "direct", 
                             mean(dir.effect), 
                             sd(dir.effect), 
                             quantile(dir.effect, prob=c(0.025)), 
                             quantile(dir.effect, prob=c(0.975)))
  effects.temp[2,] <- cbind(response, predictors[k], "indirect", 
                            mean(ind.effect), 
                            sd(ind.effect), 
                            quantile(ind.effect, prob=c(0.025)), 
                            quantile(ind.effect, prob=c(0.975)))
  effects.temp[3,] <- cbind(response, predictors[k], "total",
                            mean(tot.effect), 
                            sd(tot.effect), 
                            quantile(tot.effect, prob=c(0.025)), 
                            quantile(tot.effect, prob=c(0.975)))
  
  effects.bld <- rbind(effects.bld, effects.temp)
  
}

effects.bld <- as.data.frame(effects.bld)
effects.bld[4:7] <- sapply(effects.bld[4:7],as.numeric)
print(effects.bld, digits = 2)

#####create a table for direct, indirect, and total effects - pdn density####
effects.pdn <- matrix(NA, nrow=0, ncol = 7)
colnames(effects.pdn) <- c("response", "predictor", "effect", "mean", "sd", "95% lcl", "95% ucl")

predictors <- c("dis")
response <- "pdn"

for(k in 1:length(predictors)) {
  
  #list of all paths
  paths.x <- paths(dagified, from = predictors[k], to = response, directed=TRUE)["paths"] 
  paths.x <- as.character(paths.x$paths)
  paths.x <- paths.x[order(nchar(paths.x), paths.x)]
  
  #is there a direct effect?
  is.direct <- nrow(str_locate_all(paths.x[1], ' -> ')[[1]])==1 
  
  #number of indirect paths
  if (is.direct == FALSE) 
  {n.ind.paths <- length(paths.x)
  start <- 1
  } else 
  {n.ind.paths <- length(paths.x) - 1
  start <- 2
  }          
  
  #are there indirect paths?
  is.indirect <- n.ind.paths>0
  
  #quantify all indirect paths
  if (is.indirect ==TRUE) {
    ind.effects <- matrix(NA, nrow=3000, ncol=n.ind.paths)        
    for(i in start:length(paths.x)) {
      path <- paths.x[i]                                        
      a <- str_locate_all(path, ' -> ')                                     
      connections <- nrow(a[[1]])
      betas <- matrix(NA, nrow=3000, ncol=connections)
      for(j in 1:connections){
        x <- substring(path, a[[1]][j,]["start"]-3, a[[1]][j,]["start"]-1)
        y <- substring(path, a[[1]][j,]["end"]+1, a[[1]][j,]["end"]+3)
        betas[,j] <- eval(parse(text=paste0(x, ".", y))) 
      }
      if (is.direct == FALSE) {ind.effects[,(i)] <- rowProds(betas)
      } else {ind.effects[,(i-1)] <- rowProds(betas)
      }
    }
  } else {ind.effects <- rep(0, 3000)}
  
  #quantify direct effect
  dir.effect <- if (is.direct == TRUE) {eval(parse(text = str_replace_all(paths.x[1], " -> ", "."))) 
  } else {rep(0, 3000)}
  
  #quantify indirect effect
  ind.effect <- if (is.indirect == TRUE) {rowSums(ind.effects)} else {ind.effects}
  
  #quantify total effect
  tot.effect <- dir.effect + ind.effect
  
  #summarize results and add to table
  effects.temp <- matrix(NA, nrow=3, ncol = 7)
  
  effects.temp[1,] <-  cbind(response, predictors[k], "direct", 
                             mean(dir.effect), 
                             sd(dir.effect), 
                             quantile(dir.effect, prob=c(0.025)), 
                             quantile(dir.effect, prob=c(0.975)))
  effects.temp[2,] <- cbind(response, predictors[k], "indirect", 
                            mean(ind.effect), 
                            sd(ind.effect), 
                            quantile(ind.effect, prob=c(0.025)), 
                            quantile(ind.effect, prob=c(0.975)))
  effects.temp[3,] <- cbind(response, predictors[k], "total",
                            mean(tot.effect), 
                            sd(tot.effect), 
                            quantile(tot.effect, prob=c(0.025)), 
                            quantile(tot.effect, prob=c(0.975)))
  
  effects.pdn <- rbind(effects.pdn, effects.temp)
  
}

effects.pdn <- as.data.frame(effects.pdn)
effects.pdn[4:7] <- sapply(effects.pdn[4:7],as.numeric)
print(effects.pdn, digits = 2)

#table of all effects

effects <- rbind(effects.n1, effects.n2, effects.pdn, effects.bld,
                 effects.tre, effects.frg, effects.prd)
effects$sig <- ifelse((effects$`95% lcl` <0 & effects$`95% ucl` <0) | (effects$`95% lcl` >0 & effects$`95% ucl` >0), "*", "")

#####Difference in effects between morphs#####
#-------------------------#
#need to be careful here as variances differ between morphs - more variance among sites for mel than gry
#https://data.library.virginia.edu/the-shortcomings-of-standardized-regression-coefficients/
predictors <- c("dis", "pdn", "bld", "rod", "tre", "frg", "prd")
dir.effect.dif <- matrix(NA, nrow=3000, ncol=length(predictors))
ind.effect.dif <- matrix(NA, nrow=3000, ncol=length(predictors))
tot.effect.dif <- matrix(NA, nrow=3000, ncol=length(predictors))

for(k in 1:length(predictors)){
  dir.effect.dif[,k] <- dir.effect.n1[,k] - dir.effect.n2[,k]
  ind.effect.dif[,k] <- ind.effect.n1[,k] - ind.effect.n2[,k]
  tot.effect.dif[,k] <- tot.effect.n1[,k] - tot.effect.n2[,k]
}

effect.dif <- matrix(NA, nrow=0, ncol=6)
colnames(effect.dif) <- c("predictor", "effect", "mean", "sd", "95% lcl", "95% ucl")
for(k in 1:length(predictors)){
  
  dir.effect.temp <- cbind(predictors[k],
                           "direct",
                           mean(dir.effect.dif[,k]),
                           sd(dir.effect.dif[,k]),
                           as.numeric(quantile(dir.effect.dif[,k], prob=c(0.025))),
                           as.numeric(quantile(dir.effect.dif[,k], prob=c(0.975))))
  effect.dif <- rbind(effect.dif, dir.effect.temp)
  
  ind.effect.temp <- cbind(predictors[k],
                           "indirect",
                           mean(ind.effect.dif[,k]),
                           sd(ind.effect.dif[,k]),
                           as.numeric(quantile(ind.effect.dif[,k], prob=c(0.025))),
                           as.numeric(quantile(ind.effect.dif[,k], prob=c(0.975))))
  effect.dif <- rbind(effect.dif, ind.effect.temp)
  
  tot.effect.temp <- cbind(predictors[k],
                           "total",
                           mean(tot.effect.dif[,k]),
                           sd(tot.effect.dif[,k]),
                           as.numeric(quantile(tot.effect.dif[,k], prob=c(0.025))),
                           as.numeric(quantile(tot.effect.dif[,k], prob=c(0.975))))
  effect.dif <- rbind(effect.dif, tot.effect.temp)
  
}

effect.dif <- as.data.frame(effect.dif)
effect.dif[3:6] <- sapply(effect.dif[3:6],as.numeric)


