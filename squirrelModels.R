# Load packages ------------------------------------------------

library(dplyr)
library(jagsUI)
library(MCMCvis)
library(dagitty)

# Load data ------------------------------------------------

load("squirrelObservations.RData")
load("measuredVariables.RData")

# Bundle data for JAGS ------------------------------------------------

all.data <- list(dis = dis.s, pdn = pdn.log.s, bld = bld.logit.s,
                 rdl = rdl.s, imp = imp.logit.s,
                 tre = tre.s, frg = frg.log.s, 
                 detections = rai.sum$detections, days = rai.sum$active.days,
                 yg.occ = y.occ.g, ym.occ =y.occ.m, 
                 yg.ct = y.ct.g, ym.ct =y.ct.m, 
                 nsites = nsites,
                 nsurveys.occ = nsurveys.occ, nsurveys.ct = nsurveys.ct, 
                 tmp.occ = temp.avg.occ, tmp.ct = temp.avg.ct)
str(all.data)


# Initial abundance values ------------------------------------------------

##Initial values----
y.g.max.ct <- apply(y.ct.g, 1, max, na.rm=TRUE, warn=F)
y.g.max.ct[!is.finite(y.g.max.ct)] <- 0
y.m.max.ct <- apply(y.ct.m, 1, max, na.rm=TRUE, warn=F)
y.m.max.ct[!is.finite(y.m.max.ct)] <- 0
y.g.max.occ <- apply(y.occ.g, 1, max, na.rm=TRUE, warn=F)
y.g.max.occ[!is.finite(y.g.max.occ)] <- 0
y.m.max.occ <- apply(y.occ.m, 1, max, na.rm=TRUE, warn=F)
y.m.max.occ[!is.finite(y.m.max.occ)] <- 0

inits <- function(){list(N1 = apply(cbind(y.m.max.occ,y.m.max.ct), 1, max), 
                         N = apply(cbind(y.m.max.occ,y.m.max.ct), 1, max) + apply(cbind(y.g.max.occ,y.g.max.ct), 1, max))}

# MCMC settings --------------------------------------------------------------
na <- 5000
ni <- 30000
nt <- 5
nb <- 25000
nc <- 4 

# dag ------------------------------------------------
dag <- dagitty("dag {pdn -> bld; pdn -> imp; pdn -> frg; pdn -> U; pdn -> prd;
                     bld -> rdl; bld -> imp; bld -> mel;
                     rdl -> imp; rdl-> mel;
                     imp -> tre; imp -> frg;
                     tre -> frg; tre -> mel;
                     frg -> prd; frg -> mel; 
                     prd -> mel;
                     U -> mel; 
                     U [unobserved];
}")

# HUMAN POP DENSITY ------------------------------------------------------------

## total effect  -------------------------------------------

adjustmentSets(dag, exposure = "pdn", outcome = "mel", effect="total")

model.name <- "pdn.t"
sink(paste0(model.name, ".txt"))
cat("

model {

    ###Priors###
    
    #detection submodel parameters
    for(m in 1:2) {                 #m = morphs, 1 = melanic 2 = gray
    a.0[m] <- logit(mean.p[m])      #detection intercept
    mean.p[m] ~ dunif(0.01, 0.30)   #mean detection probabilty at mean temp
    a.tmp[m] ~ dnorm(0, 0.01)       #linear temp slope
    a.tmp2[m] ~ dnorm(0, 0.01)      #quadratic temp slope
    }
    
    #abundance submodel parameters
    mean.abu ~ dunif(1, 50)     #abundance at mean of covariates
    b.a.0 <- log(mean.abu)      #intercept
    b.a.dis ~ dnorm(0, 0.01)    #slopes

    #melanism submodel parameters
    mean.pm ~ dunif(0.01, 0.5)  #proportion melanic at value = 0 for covariates (mean distance)
    b.m.0 <- logit(mean.pm)     #mean p(melanic) at mean of covariates 
    b.m.pdn ~ dnorm(0, 0.01)    #slopes
    
    ###Likelihood###
    
    #Process model for abundance & coat color
    for(i in 1:nsites) {
    N[i] ~ dpois(lambda[i])               #total squirrel abundance
    log(lambda[i]) <- b.a.0 +
                      b.a.dis*dis[i] 

    N1[i] ~ dbinom(pm[i], N[i])           #coat color process 
    logit(pm[i]) <- b.m.0 +
                    b.m.pdn*pdn[i] 
    
    #derived:
    N2[i] <- N[i] - N1[i]                 #abundance of gray morph
    pm.obs[i] <- N1[i]/ifelse(N[i]==0, 1e6, N[i])  #derived proportion melanic at each site
    N.res[i] <- N[i] - lambda[i]          #abundance residual
    pm.res[i] <- pm.obs[i] - pm[i]       #pm residual
    
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
}

    ", fill=TRUE)
sink()

##Parameters monitored
params <- c("mean.p", "a.0", "a.tmp", "a.tmp2", 
            "mean.abu", "b.a.0", "b.a.dis",
            "mean.pm", "b.m.0", "b.m.pdn",
            "N", "N1", "N2", "pm.obs", "lambda", "N.res", "pm.res")

##MCMC
m <- jags(all.data, inits, params, paste0(model.name, ".txt"), 
          n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na,
          parallel=TRUE, n.cores=4)
assign(paste0("m.", model.name), m)

#Parameter estimates with Rhat
MCMCsummary(get(paste0("m.", model.name)), 
            params = c("mean.abu", "b.a.0", "b.a.dis", 
                       "mean.pm", "b.m.0", "b.m.pdn"),
            Rhat = TRUE,
            n.eff = TRUE,
            probs = c(0.025, 0.5, 0.975),
            round = 4)

# BUILDING COVER ------------------------------------------------------------

## direct effect  -------------------------------------------

adjustmentSets(dag, exposure = "bld", outcome = "mel", effect="direct")

model.name <- "bld.d"
sink(paste0(model.name, ".txt"))
cat("

model {

    ###Priors###
    
    #detection submodel parameters
    for(m in 1:2) {                 #m = morphs, 1 = melanic 2 = gray
    a.0[m] <- logit(mean.p[m])      #detection intercept
    mean.p[m] ~ dunif(0.01, 0.30)   #mean detection probabilty at mean temp
    a.tmp[m] ~ dnorm(0, 0.01)       #linear temp slope
    a.tmp2[m] ~ dnorm(0, 0.01)      #quadratic temp slope
    }
    
    #abundance submodel parameters
    mean.abu ~ dunif(1, 50)     #abundance at mean of covariates
    b.a.0 <- log(mean.abu)      #intercept
    b.a.dis ~ dnorm(0, 0.01)    #slopes

    #melanism submodel parameters
    mean.pm ~ dunif(0.01, 0.5)  #proportion melanic at value = 0 for covariates (mean distance)
    b.m.0 <- logit(mean.pm)     #mean p(melanic) at mean of covariates 
    b.m.bld ~ dnorm(0, 0.01)    #slopes
    b.m.pdn ~ dnorm(0, 0.01)
    b.m.rdl ~ dnorm(0, 0.01)     
    b.m.tre ~ dnorm(0, 0.01)    
    b.m.frg ~ dnorm(0, 0.01)    

    ###Likelihood###
    
    #Process model for abundance & coat color
    for(i in 1:nsites) {
    N[i] ~ dpois(lambda[i])               #total squirrel abundance
    log(lambda[i]) <- b.a.0 +
                      b.a.dis*dis[i]

    N1[i] ~ dbinom(pm[i], N[i])           #coat color process 
    logit(pm[i]) <- b.m.0 +
                    b.m.bld*bld[i] +      
                    b.m.pdn*pdn[i] +
                    b.m.rdl*rdl[i] +
                    b.m.tre*tre[i] + 
                    b.m.frg*frg[i]
    
    #derived:
    N2[i] <- N[i] - N1[i]                 #abundance of gray morph
    pm.obs[i] <- N1[i]/ifelse(N[i]==0, 1e6, N[i])  #derived proportion melanic at each site
    N.res[i] <- N[i] - lambda[i]          #abundance residual
    pm.res[i] <- pm.obs[i] - pm[i]       #pm residual
    
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
}

    ", fill=TRUE)
sink()


##Parameters monitored
params <- c("mean.p", "a.0", "a.tmp", "a.tmp2", 
            "mean.abu", "b.a.0", "b.a.dis",
            "mean.pm", "b.m.0", "b.m.bld", "b.m.pdn", "b.m.rdl", "b.m.tre", "b.m.frg",
            "N", "N1", "N2", "pm.obs", "lambda", "N.res", "pm.res")

##MCMC
m <- jags(all.data, inits, params, paste0(model.name, ".txt"), 
          n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na,
          parallel=TRUE, n.cores=4)
assign(paste0("m.", model.name), m)

m.bld.d <- update(m.bld.d, n.iter=10000, n.thin=10, parallel=TRUE)

#Parameter estimates with Rhat
MCMCsummary(get(paste0("m.", model.name)), 
            params = c("mean.abu", "b.a.0", "b.a.dis",
                       "mean.pm", "b.m.0", "b.m.bld", "b.m.pdn", "b.m.rdl", "b.m.tre", "b.m.frg"),
            Rhat = TRUE,
            n.eff = TRUE,
            probs = c(0.025, 0.5, 0.975),
            round = 4)

## total effect  -------------------------------------------

adjustmentSets(dag, exposure = "bld", outcome = "mel", effect="total")

model.name <- "bld.t"
sink(paste0(model.name, ".txt"))
cat("

model {

    ###Priors###
    
    #detection submodel parameters
    for(m in 1:2) {                 #m = morphs, 1 = melanic 2 = gray
    a.0[m] <- logit(mean.p[m])      #detection intercept
    mean.p[m] ~ dunif(0.01, 0.30)   #mean detection probabilty at mean temp
    a.tmp[m] ~ dnorm(0, 0.01)       #linear temp slope
    a.tmp2[m] ~ dnorm(0, 0.01)      #quadratic temp slope
    }
    
    #abundance submodel parameters
    mean.abu ~ dunif(1, 50)     #abundance at mean of covariates
    b.a.0 <- log(mean.abu)      #intercept
    b.a.dis ~ dnorm(0, 0.01)    #slopes

    #melanism submodel parameters
    mean.pm ~ dunif(0.01, 0.5)  #proportion melanic at value = 0 for covariates (mean distance)
    b.m.0 <- logit(mean.pm)     #mean p(melanic) at mean of covariates 
    b.m.bld ~ dnorm(0, 0.01)    #slopes
    b.m.pdn ~ dnorm(0, 0.01)

    ###Likelihood###
    
    #Process model for abundance & coat color
    for(i in 1:nsites) {
    N[i] ~ dpois(lambda[i])               #total squirrel abundance
    log(lambda[i]) <- b.a.0 +
                      b.a.dis*dis[i]
                      
    N1[i] ~ dbinom(pm[i], N[i])           #coat color process 
    logit(pm[i]) <- b.m.0 +
                    b.m.bld*bld[i] +      
                    b.m.pdn*pdn[i] 
    
    #derived:
    N2[i] <- N[i] - N1[i]                 #abundance of gray morph
    pm.obs[i] <- N1[i]/ifelse(N[i]==0, 1e6, N[i])  #derived proportion melanic at each site
    N.res[i] <- N[i] - lambda[i]          #abundance residual
    pm.res[i] <- pm.obs[i] - pm[i]       #pm residual
    
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
}

    ", fill=TRUE)
sink()

##Parameters monitored
params <- c("mean.p", "a.0", "a.tmp", "a.tmp2", 
            "mean.abu", "b.a.0", "b.a.dis", 
            "mean.pm", "b.m.0", "b.m.bld", "b.m.pdn", 
            "N", "N1", "N2", "pm.obs", "lambda", "N.res", "pm.res")

##MCMC
m <- jags(all.data, inits, params, paste0(model.name, ".txt"), 
          n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na,
          parallel=TRUE, n.cores=4)
assign(paste0("m.", model.name), m)

#Parameter estimates with Rhat
MCMCsummary(get(paste0("m.", model.name)), 
            params = c("mean.abu", "b.a.0", "b.a.dis", 
                       "mean.pm", "b.m.0", "b.m.bld", "b.m.pdn"),
            Rhat = TRUE,
            n.eff = TRUE,
            probs = c(0.025, 0.5, 0.975),
            round = 4)

## univariate model  -------------------------------------------

model.name <- "bld.u"
sink(paste0(model.name, ".txt"))
cat("

model {

    ###Priors###
    
    #detection submodel parameters
    for(m in 1:2) {                 #m = morphs, 1 = melanic 2 = gray
    a.0[m] <- logit(mean.p[m])      #detection intercept
    mean.p[m] ~ dunif(0.01, 0.30)   #mean detection probabilty at mean temp
    a.tmp[m] ~ dnorm(0, 0.01)       #linear temp slope
    a.tmp2[m] ~ dnorm(0, 0.01)      #quadratic temp slope
    }
    
    #abundance submodel parameters
    mean.abu ~ dunif(1, 50)     #abundance at mean of covariates
    b.a.0 <- log(mean.abu)      #intercept
    b.a.dis ~ dnorm(0, 0.01)    #slopes

    #melanism submodel parameters
    mean.pm ~ dunif(0.01, 0.5)  #proportion melanic at value = 0 for covariates (mean distance)
    b.m.0 <- logit(mean.pm)     #mean p(melanic) at mean of covariates 
    b.m.bld ~ dnorm(0, 0.01)    #slopes

    ###Likelihood###
    
    #Process model for abundance & coat color
    for(i in 1:nsites) {
    N[i] ~ dpois(lambda[i])               #total squirrel abundance
    log(lambda[i]) <- b.a.0 +
                      b.a.dis*dis[i] 
                      
    N1[i] ~ dbinom(pm[i], N[i])           #coat color process 
    logit(pm[i]) <- b.m.0 +
                    b.m.bld*bld[i] 
                    
    #derived:
    N2[i] <- N[i] - N1[i]                 #abundance of gray morph
    pm.obs[i] <- N1[i]/ifelse(N[i]==0, 1e6, N[i])  #derived proportion melanic at each site
    N.res[i] <- N[i] - lambda[i]          #abundance residual
    pm.res[i] <- pm.obs[i] - pm[i]       #pm residual
    
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
}

    ", fill=TRUE)
sink()

##Parameters monitored
params <- c("mean.p", "a.0", "a.tmp", "a.tmp2", 
            "mean.abu", "b.a.0", "b.a.dis", 
            "mean.pm", "b.m.0", "b.m.bld",
            "N", "N1", "N2", "pm.obs", "lambda", "N.res", "pm.res")

##MCMC
m <- jags(all.data, inits, params, paste0(model.name, ".txt"), 
          n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na,
          parallel=TRUE, n.cores=4)
assign(paste0("m.", model.name), m)

m.bld.u <- update(m.bld.u, n.iter=10000, n.thin=10, parallel=TRUE)

#Parameter estimates with Rhat
MCMCsummary(get(paste0("m.", model.name)), 
            params = c("mean.abu", "b.a.0", "b.a.dis", 
                       "mean.pm", "b.m.0", "b.m.bld"),
            Rhat = TRUE,
            n.eff = TRUE,
            probs = c(0.025, 0.5, 0.975),
            round = 4)

# ROAD COVER ------------------------------------------------------------

## direct effect  -------------------------------------------

adjustmentSets(dag, exposure = "rdl", outcome = "mel", effect="direct")

m.rdl.d <- m.bld.d #identical models

## total effect   -------------------------------------------

adjustmentSets(dag, exposure = "rdl", outcome = "mel", effect="total")

model.name <- "rdl.t"
sink(paste0(model.name, ".txt"))
cat("

model {

    ###Priors###
    
    #detection submodel parameters
    for(m in 1:2) {                 #m = morphs, 1 = melanic 2 = gray
    a.0[m] <- logit(mean.p[m])      #detection intercept
    mean.p[m] ~ dunif(0.01, 0.30)   #mean detection probabilty at mean temp
    a.tmp[m] ~ dnorm(0, 0.01)       #linear temp slope
    a.tmp2[m] ~ dnorm(0, 0.01)      #quadratic temp slope
    }
    
    #abundance submodel parameters
    mean.abu ~ dunif(1, 50)     #abundance at mean of covariates
    b.a.0 <- log(mean.abu)      #intercept
    b.a.dis ~ dnorm(0, 0.01)    #slopes

    #melanism submodel parameters
    mean.pm ~ dunif(0.01, 0.5)  #proportion melanic at value = 0 for covariates (mean distance)
    b.m.0 <- logit(mean.pm)     #mean p(melanic) at mean of covariates 
    b.m.bld ~ dnorm(0, 0.01)    #slopes
    b.m.rdl ~ dnorm(0, 0.01)

    ###Likelihood###
    
    #Process model for abundance & coat color
    for(i in 1:nsites) {
    N[i] ~ dpois(lambda[i])               #total squirrel abundance
    log(lambda[i]) <- b.a.0 +
                      b.a.dis*dis[i] 
                      
    N1[i] ~ dbinom(pm[i], N[i])           #coat color process 
    logit(pm[i]) <- b.m.0 +
                    b.m.bld*bld[i] +      
                    b.m.rdl*rdl[i]
    
    #derived:
    N2[i] <- N[i] - N1[i]                 #abundance of gray morph
    pm.obs[i] <- N1[i]/ifelse(N[i]==0, 1e6, N[i])  #derived proportion melanic at each site
    N.res[i] <- N[i] - lambda[i]          #abundance residual
    pm.res[i] <- pm.obs[i] - pm[i]       #pm residual
    
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
}

    ", fill=TRUE)
sink()

##Parameters monitored
params <- c("mean.p", "a.0", "a.tmp", "a.tmp2", 
            "mean.abu", "b.a.0", "b.a.dis", 
            "mean.pm", "b.m.0", "b.m.bld","b.m.rdl", 
            "N", "N1", "N2", "pm.obs", "lambda", "N.res", "pm.res")

##MCMC
m <- jags(all.data, inits, params, paste0(model.name, ".txt"), 
          n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na,
          parallel=TRUE, n.cores=4)
assign(paste0("m.", model.name), m)

m.rdl.t <- update(m.rdl.t, n.iter=10000, n.thin=10, parallel=TRUE)

#Parameter estimates with Rhat
MCMCsummary(get(paste0("m.", model.name)), 
            params = c("mean.abu", "b.a.0", "b.a.dis", 
                       "mean.pm", "b.m.0", "b.m.bld", "b.m.rdl"),
            Rhat = TRUE,
            n.eff = TRUE,
            probs = c(0.025, 0.5, 0.975),
            round = 4)


## univariate model -------------------------------------------
model.name <- "rdl.u"
sink(paste0(model.name, ".txt"))
cat("

model {

    ###Priors###
    
    #detection submodel parameters
    for(m in 1:2) {                 #m = morphs, 1 = melanic 2 = gray
    a.0[m] <- logit(mean.p[m])      #detection intercept
    mean.p[m] ~ dunif(0.01, 0.30)   #mean detection probabilty at mean temp
    a.tmp[m] ~ dnorm(0, 0.01)       #linear temp slope
    a.tmp2[m] ~ dnorm(0, 0.01)      #quadratic temp slope
    }
    
    #abundance submodel parameters
    mean.abu ~ dunif(1, 50)     #abundance at mean of covariates
    b.a.0 <- log(mean.abu)      #intercept
    b.a.dis ~ dnorm(0, 0.01)    #slopes

    #melanism submodel parameters
    mean.pm ~ dunif(0.01, 0.5)  #proportion melanic at value = 0 for covariates (mean distance)
    b.m.0 <- logit(mean.pm)     #mean p(melanic) at mean of covariates 
    b.m.rdl ~ dnorm(0, 0.01)    #slopes

    ###Likelihood###
    
    #Process model for abundance & coat color
    for(i in 1:nsites) {
    N[i] ~ dpois(lambda[i])               #total squirrel abundance
    log(lambda[i]) <- b.a.0 +
                      b.a.dis*dis[i]

    N1[i] ~ dbinom(pm[i], N[i])           #coat color process 
    logit(pm[i]) <- b.m.0 +
                    b.m.rdl*rdl[i]
    
    #derived:
    N2[i] <- N[i] - N1[i]                 #abundance of gray morph
    pm.obs[i] <- N1[i]/ifelse(N[i]==0, 1e6, N[i])  #derived proportion melanic at each site
    N.res[i] <- N[i] - lambda[i]          #abundance residual
    pm.res[i] <- pm.obs[i] - pm[i]       #pm residual
    
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
}

    ", fill=TRUE)
sink()

##Parameters monitored
params <- c("mean.p", "a.0", "a.tmp", "a.tmp2", 
            "mean.abu", "b.a.0", "b.a.dis",
            "mean.pm", "b.m.0", "b.m.rdl",
            "N", "N1", "N2", "pm.obs", "lambda", "N.res", "pm.res")

##MCMC
m <- jags(all.data, inits, params, paste0(model.name, ".txt"), 
          n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na,
          parallel=TRUE, n.cores=4)
assign(paste0("m.", model.name), m)

m.rdl.u <- update(m.rdl.u, n.iter=10000, n.thin=10, parallel=TRUE)

#Parameter estimates with Rhat
MCMCsummary(get(paste0("m.", model.name)), 
            params = c("mean.abu", "b.a.0", "b.a.dis", 
                       "mean.pm", "b.m.0", "b.m.rdl"),
            Rhat = TRUE,
            n.eff = TRUE,
            probs = c(0.025, 0.5, 0.975),
            round = 4)

# Impervious cover ------------------------------------------------------------

## total effect  -------------------------------------------

adjustmentSets(dag, exposure = "imp", outcome = "mel", effect="total")

model.name <- "imp.t"
sink(paste0(model.name, ".txt"))
cat("

model {

    ###Priors###
    
    #detection submodel parameters
    for(m in 1:2) {                 #m = morphs, 1 = melanic 2 = gray
    a.0[m] <- logit(mean.p[m])      #detection intercept
    mean.p[m] ~ dunif(0.01, 0.30)   #mean detection probabilty at mean temp
    a.tmp[m] ~ dnorm(0, 0.01)       #linear temp slope
    a.tmp2[m] ~ dnorm(0, 0.01)      #quadratic temp slope
    }
    
    #abundance submodel parameters
    mean.abu ~ dunif(1, 50)     #abundance at mean of covariates
    b.a.0 <- log(mean.abu)      #intercept
    b.a.dis ~ dnorm(0, 0.01)    #slopes

    #melanism submodel parameters
    mean.pm ~ dunif(0.01, 0.5)  #proportion melanic at value = 0 for covariates (mean distance)
    b.m.0 <- logit(mean.pm)     #mean p(melanic) at mean of covariates 
    b.m.imp ~ dnorm(0, 0.01)    #slopes
    b.m.rdl ~ dnorm(0, 0.01)
    b.m.pdn ~ dnorm(0, 0.01)
    b.m.bld ~ dnorm(0, 0.01)


    ###Likelihood###
    
    #Process model for abundance & coat color
    for(i in 1:nsites) {
    N[i] ~ dpois(lambda[i])               #total squirrel abundance
    log(lambda[i]) <- b.a.0 +
                      b.a.dis*dis[i]

    N1[i] ~ dbinom(pm[i], N[i])           #coat color process 
    logit(pm[i]) <- b.m.0 +
                    b.m.imp*imp[i] +      
                    b.m.rdl*rdl[i] +
                    b.m.pdn*pdn[i] +
                    b.m.bld*bld[i]
    
    #derived:
    N2[i] <- N[i] - N1[i]                 #abundance of gray morph
    pm.obs[i] <- N1[i]/ifelse(N[i]==0, 1e6, N[i])  #derived proportion melanic at each site
    N.res[i] <- N[i] - lambda[i]          #abundance residual
    pm.res[i] <- pm.obs[i] - pm[i]       #pm residual
    
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
}

    ", fill=TRUE)
sink()

##Parameters monitored
params <- c("mean.p", "a.0", "a.tmp", "a.tmp2", 
            "mean.abu", "b.a.0", "b.a.dis", 
            "mean.pm", "b.m.0", "b.m.imp", "b.m.bld", "b.m.rdl", "b.m.pdn", 
            "N", "N1", "N2", "pm.obs", "lambda", "N.res", "pm.res")

##MCMC
m <- jags(all.data, inits, params, paste0(model.name, ".txt"), 
          n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na,
          parallel=TRUE, n.cores=4)
assign(paste0("m.", model.name), m)

m.imp.t <- update(m.imp.t, n.iter=10000, n.thin=10, parallel=TRUE)

#Parameter estimates with Rhat
MCMCsummary(get(paste0("m.", model.name)), 
            params = c("mean.abu", "b.a.0", "b.a.dis", 
                       "mean.pm", "b.m.0", "b.m.imp", "b.m.rdl", "b.m.pdn", "b.m.bld"),
            Rhat = TRUE,
            n.eff = TRUE,
            probs = c(0.025, 0.5, 0.975),
            round = 4)

## univariate model  -------------------------------------------

model.name <- "imp.u"
sink(paste0(model.name, ".txt"))
cat("

model {

    ###Priors###
    
    #detection submodel parameters
    for(m in 1:2) {                 #m = morphs, 1 = melanic 2 = gray
    a.0[m] <- logit(mean.p[m])      #detection intercept
    mean.p[m] ~ dunif(0.01, 0.30)   #mean detection probabilty at mean temp
    a.tmp[m] ~ dnorm(0, 0.01)       #linear temp slope
    a.tmp2[m] ~ dnorm(0, 0.01)      #quadratic temp slope
    }
    
    #abundance submodel parameters
    mean.abu ~ dunif(1, 50)     #abundance at mean of covariates
    b.a.0 <- log(mean.abu)      #intercept
    b.a.dis ~ dnorm(0, 0.01)    #slopes

    #melanism submodel parameters
    mean.pm ~ dunif(0.01, 0.5)  #proportion melanic at value = 0 for covariates (mean distance)
    b.m.0 <- logit(mean.pm)     #mean p(melanic) at mean of covariates 
    b.m.imp ~ dnorm(0, 0.01)    #slopes

    ###Likelihood###
    
    #Process model for abundance & coat color
    for(i in 1:nsites) {
    N[i] ~ dpois(lambda[i])               #total squirrel abundance
    log(lambda[i]) <- b.a.0 +
                      b.a.dis*dis[i] 
                      
    N1[i] ~ dbinom(pm[i], N[i])           #coat color process 
    logit(pm[i]) <- b.m.0 +
                    b.m.imp*imp[i] 
                    
    #derived:
    N2[i] <- N[i] - N1[i]                 #abundance of gray morph
    pm.obs[i] <- N1[i]/ifelse(N[i]==0, 1e6, N[i])  #derived proportion melanic at each site
    N.res[i] <- N[i] - lambda[i]          #abundance residual
    pm.res[i] <- pm.obs[i] - pm[i]       #pm residual
    
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
}

    ", fill=TRUE)
sink()

##Parameters monitored
params <- c("mean.p", "a.0", "a.tmp", "a.tmp2", 
            "mean.abu", "b.a.0", "b.a.dis", 
            "mean.pm", "b.m.0", "b.m.imp",
            "N", "N1", "N2", "pm.obs", "lambda", "N.res", "pm.res")

##MCMC
m <- jags(all.data, inits, params, paste0(model.name, ".txt"), 
          n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na,
          parallel=TRUE, n.cores=4)
assign(paste0("m.", model.name), m)

#Parameter estimates with Rhat
MCMCsummary(get(paste0("m.", model.name)), 
            params = c("mean.abu", "b.a.0", "b.a.dis", 
                       "mean.pm", "b.m.0", "b.m.imp"),
            Rhat = TRUE,
            n.eff = TRUE,
            probs = c(0.025, 0.5, 0.975),
            round = 4)


# FOREST COVER ------------------------------------------------------------

## direct effect  ---------
adjustmentSets(dag, exposure = "tre", outcome = "mel", effect="direct")

m.tre.d <- m.bld.d #identical models

## total effect -------------------------------------------

adjustmentSets(dag, exposure = "tre", outcome = "mel", effect="total")

model.name <- "tre.t"
sink(paste0(model.name, ".txt"))
cat("

model {

    ###Priors###
    
    #detection submodel parameters
    for(m in 1:2) {                 #m = morphs, 1 = melanic 2 = gray
    a.0[m] <- logit(mean.p[m])      #detection intercept
    mean.p[m] ~ dunif(0.01, 0.30)   #mean detection probabilty at mean temp
    a.tmp[m] ~ dnorm(0, 0.01)       #linear temp slope
    a.tmp2[m] ~ dnorm(0, 0.01)      #quadratic temp slope
    }
    
    #abundance submodel parameters
    mean.abu ~ dunif(1, 50)     #abundance at mean of covariates
    b.a.0 <- log(mean.abu)      #intercept
    b.a.dis ~ dnorm(0, 0.01)    #slopes

    #melanism submodel parameters
    mean.pm ~ dunif(0.01, 0.5)  #proportion melanic at value = 0 for covariates (mean distance)
    b.m.0 <- logit(mean.pm)     #mean p(melanic) at mean of covariates 
    b.m.tre ~ dnorm(0, 0.01)    #slopes
    b.m.imp ~ dnorm(0, 0.01)

    ###Likelihood###
    
    #Process model for abundance & coat color
    for(i in 1:nsites) {
    N[i] ~ dpois(lambda[i])               #total squirrel abundance
    log(lambda[i]) <- b.a.0 +
                      b.a.dis*dis[i]
                      
    N1[i] ~ dbinom(pm[i], N[i])           #coat color process 
    logit(pm[i]) <- b.m.0 +
                    b.m.tre*tre[i] +      
                    b.m.imp*imp[i] 
    
    #derived:
    N2[i] <- N[i] - N1[i]                 #abundance of gray morph
    pm.obs[i] <- N1[i]/ifelse(N[i]==0, 1e6, N[i])  #derived proportion melanic at each site
    N.res[i] <- N[i] - lambda[i]          #abundance residual
    pm.res[i] <- pm.obs[i] - pm[i]       #pm residual
    
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
}

    ", fill=TRUE)
sink()

##Parameters monitored
params <- c("mean.p", "a.0", "a.tmp", "a.tmp2", 
            "mean.abu", "b.a.0", "b.a.dis", 
            "mean.pm", "b.m.0", "b.m.tre", "b.m.imp", 
            "N", "N1", "N2", "pm.obs", "lambda", "N.res", "pm.res")

##MCMC
m <- jags(all.data, inits, params, paste0(model.name, ".txt"), 
          n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na,
          parallel=TRUE, n.cores=4)
assign(paste0("m.", model.name), m)

m.tre.t <- update(m.tre.t, n.iter=10000, n.thin=10, parallel=TRUE)

#Parameter estimates with Rhat
MCMCsummary(get(paste0("m.", model.name)), 
            params = c("mean.abu", "b.a.0", "b.a.dis",
                       "mean.pm", "b.m.0", "b.m.tre", "b.m.imp"),
            Rhat = TRUE,
            n.eff = TRUE,
            probs = c(0.025, 0.5, 0.975),
            round = 4)

save.image("~/Library/CloudStorage/Dropbox/projects/sccaColor/nsf2018249/manuscripts/causalInference/analysis/squirrelPaths_out.RData")


## univariate model  -------------------------------------------

model.name <- "tre.u"
sink(paste0(model.name, ".txt"))
cat("

model {

    ###Priors###
    
    #detection submodel parameters
    for(m in 1:2) {                 #m = morphs, 1 = melanic 2 = gray
    a.0[m] <- logit(mean.p[m])      #detection intercept
    mean.p[m] ~ dunif(0.01, 0.30)   #mean detection probabilty at mean temp
    a.tmp[m] ~ dnorm(0, 0.01)       #linear temp slope
    a.tmp2[m] ~ dnorm(0, 0.01)      #quadratic temp slope
    }
    
    #abundance submodel parameters
    mean.abu ~ dunif(1, 50)     #abundance at mean of covariates
    b.a.0 <- log(mean.abu)      #intercept
    b.a.dis ~ dnorm(0, 0.01)    #slopes

    #melanism submodel parameters
    mean.pm ~ dunif(0.01, 0.5)  #proportion melanic at value = 0 for covariates (mean distance)
    b.m.0 <- logit(mean.pm)     #mean p(melanic) at mean of covariates 
    b.m.tre ~ dnorm(0, 0.01)    #slopes
    
    ###Likelihood###
    
    #Process model for abundance & coat color
    for(i in 1:nsites) {
    N[i] ~ dpois(lambda[i])               #total squirrel abundance
    log(lambda[i]) <- b.a.0 +
                      b.a.dis*dis[i] 

    N1[i] ~ dbinom(pm[i], N[i])           #coat color process 
    logit(pm[i]) <- b.m.0 +
                    b.m.tre*tre[i] 
    
    #derived:
    N2[i] <- N[i] - N1[i]                 #abundance of gray morph
    pm.obs[i] <- N1[i]/ifelse(N[i]==0, 1e6, N[i])  #derived proportion melanic at each site
    N.res[i] <- N[i] - lambda[i]          #abundance residual
    pm.res[i] <- pm.obs[i] - pm[i]       #pm residual
    
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
}

    ", fill=TRUE)
sink()

##Parameters monitored
params <- c("mean.p", "a.0", "a.tmp", "a.tmp2", 
            "mean.abu", "b.a.0", "b.a.dis",
            "mean.pm", "b.m.0", "b.m.tre",
            "N", "N1", "N2", "pm.obs", "lambda", "N.res", "pm.res")

##MCMC
m <- jags(all.data, inits, params, paste0(model.name, ".txt"), 
          n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na,
          parallel=TRUE, n.cores=4)
assign(paste0("m.", model.name), m)

m.tre.u <- update(m.tre.u, n.iter=10000, n.thin=10, parallel=TRUE)

#Parameter estimates with Rhat
MCMCsummary(get(paste0("m.", model.name)), 
            params = c("mean.abu", "b.a.0", "b.a.dis", 
                       "mean.pm", "b.m.0", "b.m.tre"),
            Rhat = TRUE,
            n.eff = TRUE,
            probs = c(0.025, 0.5, 0.975),
            round = 4)

# FRAGMENTATION ------------------------------------------------------------

#more iterations
na <- 20000
ni <- 200000
nt <- 25
nb <- 175000
nc <- 4 

## direct effect  ---------
adjustmentSets(dag, exposure = "frg", outcome = "mel", effect="direct")

model.name <- "frg.d"
sink(paste0(model.name, ".txt"))
cat("
    
    model {
      
      ###Priors###
      
      #detection submodel parameters
      for(m in 1:2) {                 #m = morphs, 1 = melanic 2 = gray
        a.0[m] <- logit(mean.p[m])      #detection intercept
        mean.p[m] ~ dunif(0.01, 0.30)   #mean detection probabilty at mean temp
        a.tmp[m] ~ dnorm(0, 0.01)       #linear temp slope
        a.tmp2[m] ~ dnorm(0, 0.01)      #quadratic temp slope
      }
      
      #abundance submodel parameters
      mean.abu ~ dunif(1, 50)     #abundance at mean of covariates
      b.a.0 <- log(mean.abu)      #intercept
      b.a.dis ~ dnorm(0, 0.01)    #slopes
      
      #melanism submodel parameters
      mean.pm ~ dunif(0.01, 0.5)  #proportion melanic at value = 0 for covariates (mean distance)
      b.m.0 <- logit(mean.pm)     #mean p(melanic) at mean of covariates 
      b.m.prd ~ dnorm(0, 0.01)    #slopes
      b.m.pdn ~ dnorm(0, 0.01)    #slopes
      b.m.tre ~ dnorm(0, 0.01)    #slopes
      b.m.frg ~ dnorm(0, 0.01)    #slopes
      b.m.imp ~ dnorm(0, 0.01)
      
      #daily carnivore detection probability
      for (i in 1:nsites) {
        prd[i] ~ dunif(0, 1)  # uniform prior for parameter
        detections[i] ~ dbinom(prd[i], days[i])    # binomial likelihood for detections
      }
      
      ###Likelihood###
      
      #Process model for abundance & coat color
      for(i in 1:nsites) {
        N[i] ~ dpois(lambda[i])               #total squirrel abundance
        log(lambda[i]) <- b.a.0 +
          b.a.dis*dis[i]
        
        N1[i] ~ dbinom(pm[i], N[i])           #coat color process 
        logit(pm[i]) <- b.m.0 +
          b.m.prd*prd[i] + 
          b.m.pdn*pdn[i] + 
          b.m.tre*tre[i] + 
          b.m.frg*frg[i] +
          b.m.imp*imp[i] 
        
        #derived:
        N2[i] <- N[i] - N1[i]                 #abundance of gray morph
        pm.obs[i] <- N1[i]/ifelse(N[i]==0, 1e6, N[i])  #derived proportion melanic at each site
        N.res[i] <- N[i] - lambda[i]          #abundance residual
        pm.res[i] <- pm.obs[i] - pm[i]       #pm residual
        
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
    }
    
    ", fill=TRUE)
sink()

##Parameters monitored
params <- c("mean.p", "a.0", "a.tmp", "a.tmp2", 
            "mean.abu", "b.a.0", "b.a.dis",
            "mean.pm", "b.m.0", "b.m.prd", "b.m.pdn", "b.m.tre", "b.m.frg", "b.m.bld", "b.m.imp",
            "N", "N1", "N2", "pm.obs", "lambda", "N.res", "pm.res", "prd")

##MCMC
m <- jags(all.data, inits, params, paste0(model.name, ".txt"), 
          n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na,
          parallel=TRUE, n.cores=4)
assign(paste0("m.", model.name), m)

m.frg.d <- update(m.frg.d, n.iter=10000, n.thin=10, parallel=TRUE)

#Parameter estimates with Rhat
MCMCsummary(get(paste0("m.", model.name)), 
            params = c("mean.abu", "b.a.0", "b.a.dis",
                       "mean.pm", "b.m.0", "b.m.prd", "b.m.pdn", "b.m.tre", "b.m.frg", "b.m.imp"),
            Rhat = TRUE,
            n.eff = TRUE,
            probs = c(0.025, 0.5, 0.975),
            round = 4)


## total effect  ---------
adjustmentSets(dag, exposure = "frg", outcome = "mel", effect="total")

model.name <- "frg.t"
sink(paste0(model.name, ".txt"))
cat("
    
    model {
      
      ###Priors###
      
      #detection submodel parameters
      for(m in 1:2) {                 #m = morphs, 1 = melanic 2 = gray
        a.0[m] <- logit(mean.p[m])      #detection intercept
        mean.p[m] ~ dunif(0.01, 0.30)   #mean detection probabilty at mean temp
        a.tmp[m] ~ dnorm(0, 0.01)       #linear temp slope
        a.tmp2[m] ~ dnorm(0, 0.01)      #quadratic temp slope
      }
      
      #abundance submodel parameters
      mean.abu ~ dunif(1, 50)     #abundance at mean of covariates
      b.a.0 <- log(mean.abu)      #intercept
      b.a.dis ~ dnorm(0, 0.01)    #slopes
      
      #melanism submodel parameters
      mean.pm ~ dunif(0.01, 0.5)  #proportion melanic at value = 0 for covariates (mean distance)
      b.m.0 <- logit(mean.pm)     #mean p(melanic) at mean of covariates 
      b.m.imp ~ dnorm(0, 0.01)    #slopes
      b.m.pdn ~ dnorm(0, 0.01)    #slopes
      b.m.tre ~ dnorm(0, 0.01)    #slopes
      b.m.frg ~ dnorm(0, 0.01)    #slopes
      
      
      ###Likelihood###
      
      #Process model for abundance & coat color
      for(i in 1:nsites) {
        N[i] ~ dpois(lambda[i])               #total squirrel abundance
        log(lambda[i]) <- b.a.0 +
          b.a.dis*dis[i] 
        
        N1[i] ~ dbinom(pm[i], N[i])           #coat color process 
        logit(pm[i]) <- b.m.0 +
          b.m.imp*imp[i] + 
          b.m.pdn*pdn[i] + 
          b.m.tre*tre[i] + 
          b.m.frg*frg[i]
        
        #derived:
        N2[i] <- N[i] - N1[i]                 #abundance of gray morph
        pm.obs[i] <- N1[i]/ifelse(N[i]==0, 1e6, N[i])  #derived proportion melanic at each site
        N.res[i] <- N[i] - lambda[i]          #abundance residual
        pm.res[i] <- pm.obs[i] - pm[i]       #pm residual
        
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
    }
    
    ", fill=TRUE)
sink()

##Parameters monitored
params <- c("mean.p", "a.0", "a.tmp", "a.tmp2", 
            "mean.abu", "b.a.0","b.a.dis",
            "mean.pm", "b.m.0", "b.m.imp", "b.m.pdn", "b.m.tre", "b.m.frg", 
            "N", "N1", "N2", "pm.obs", "lambda", "N.res", "pm.res", "imp")

##MCMC
m <- jags(all.data, inits, params, paste0(model.name, ".txt"), 
          n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na,
          parallel=TRUE, n.cores=4)
assign(paste0("m.", model.name), m)

m.frg.t <- update(m.frg.t, n.iter=10000, n.thin=10, parallel=TRUE)

#Parameter estimates with Rhat
MCMCsummary(get(paste0("m.", model.name)), 
            params = c("mean.abu", "b.a.0", "b.a.dis",
                       "mean.pm", "b.m.0", "b.m.imp", "b.m.pdn", "b.m.tre", "b.m.frg" ),
            Rhat = TRUE,
            n.eff = TRUE,
            probs = c(0.025, 0.5, 0.975),
            round = 4)

## univariate model  -------------------------------------------

model.name <- "frg.u"
sink(paste0(model.name, ".txt"))
cat("

model {

    ###Priors###
    
    #detection submodel parameters
    for(m in 1:2) {                 #m = morphs, 1 = melanic 2 = gray
    a.0[m] <- logit(mean.p[m])      #detection intercept
    mean.p[m] ~ dunif(0.01, 0.30)   #mean detection probabilty at mean temp
    a.tmp[m] ~ dnorm(0, 0.01)       #linear temp slope
    a.tmp2[m] ~ dnorm(0, 0.01)      #quadratic temp slope
    }
    
    #abundance submodel parameters
    mean.abu ~ dunif(1, 50)     #abundance at mean of covariates
    b.a.0 <- log(mean.abu)      #intercept
    b.a.dis ~ dnorm(0, 0.01)    #slopes

    #melanism submodel parameters
    mean.pm ~ dunif(0.01, 0.5)  #proportion melanic at value = 0 for covariates (mean distance)
    b.m.0 <- logit(mean.pm)     #mean p(melanic) at mean of covariates 
    b.m.frg ~ dnorm(0, 0.01)    #slopes
    
    ###Likelihood###
    
    #Process model for abundance & coat color
    for(i in 1:nsites) {
    N[i] ~ dpois(lambda[i])               #total squirrel abundance
    log(lambda[i]) <- b.a.0 +
                      b.a.dis*dis[i] 

    N1[i] ~ dbinom(pm[i], N[i])           #coat color process 
    logit(pm[i]) <- b.m.0 +
                    b.m.frg*frg[i] 
    
    #derived:
    N2[i] <- N[i] - N1[i]                 #abundance of gray morph
    pm.obs[i] <- N1[i]/ifelse(N[i]==0, 1e6, N[i])  #derived proportion melanic at each site
    N.res[i] <- N[i] - lambda[i]          #abundance residual
    pm.res[i] <- pm.obs[i] - pm[i]       #pm residual
    
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
}

    ", fill=TRUE)
sink()

##Parameters monitored
params <- c("mean.p", "a.0", "a.tmp", "a.tmp2", 
            "mean.abu", "b.a.0", "b.a.dis",
            "mean.pm", "b.m.0", "b.m.frg",
            "N", "N1", "N2", "pm.obs", "lambda", "N.res", "pm.res")

##MCMC
m <- jags(all.data, inits, params, paste0(model.name, ".txt"), 
          n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na,
          parallel=TRUE, n.cores=4)
assign(paste0("m.", model.name), m)

#Parameter estimates with Rhat
MCMCsummary(get(paste0("m.", model.name)), 
            params = c("mean.abu", "b.a.0", "b.a.dis", 
                       "mean.pm", "b.m.0", "b.m.frg"),
            Rhat = TRUE,
            n.eff = TRUE,
            probs = c(0.025, 0.5, 0.975),
            round = 4)

# CARNIVORE ACTIVITY ------------------------------------------------------------

## direct effect  -------------------------------------------
adjustmentSets(dag, exposure = "prd", outcome = "mel", effect="direct")

model.name <- "prd.d"
sink(paste0(model.name, ".txt"))
cat("
    
    model {
      
      ###Priors###
      
      #detection submodel parameters
      for(m in 1:2) {                 #m = morphs, 1 = melanic 2 = gray
        a.0[m] <- logit(mean.p[m])      #detection intercept
        mean.p[m] ~ dunif(0.01, 0.30)   #mean detection probabilty at mean temp
        a.tmp[m] ~ dnorm(0, 0.01)       #linear temp slope
        a.tmp2[m] ~ dnorm(0, 0.01)      #quadratic temp slope
      }
      
      #abundance submodel parameters
      mean.abu ~ dunif(1, 50)     #abundance at mean of covariates
      b.a.0 <- log(mean.abu)      #intercept
      b.a.dis ~ dnorm(0, 0.01)    #slopes
      
      #melanism submodel parameters
      mean.pm ~ dunif(0.01, 0.5)  #proportion melanic at value = 0 for covariates (mean distance)
      b.m.0 <- logit(mean.pm)     #mean p(melanic) at mean of covariates 
      b.m.prd ~ dnorm(0, 0.01)    #slopes
      b.m.pdn ~ dnorm(0, 0.01)    #slopes
      b.m.frg ~ dnorm(0, 0.01)    #slopes
      
      #daily carnivore detection probability
      for (i in 1:nsites) {
        prd[i] ~ dunif(0, 1)  # uniform prior for parameter
        detections[i] ~ dbinom(prd[i], days[i])    # binomial likelihood for detections
      }
      
      ###Likelihood###
      
      #Process model for abundance & coat color
      for(i in 1:nsites) {
        N[i] ~ dpois(lambda[i])               #total squirrel abundance
        log(lambda[i]) <- b.a.0 +
          b.a.dis*dis[i]
        
        N1[i] ~ dbinom(pm[i], N[i])           #coat color process 
        logit(pm[i]) <- b.m.0 +
          b.m.prd*prd[i] + 
          b.m.pdn*pdn[i] + 
          b.m.frg*frg[i] 
        
        #derived:
        N2[i] <- N[i] - N1[i]                 #abundance of gray morph
        pm.obs[i] <- N1[i]/ifelse(N[i]==0, 1e6, N[i])  #derived proportion melanic at each site
        N.res[i] <- N[i] - lambda[i]          #abundance residual
        pm.res[i] <- pm.obs[i] - pm[i]       #pm residual
        
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
    }
    
    ", fill=TRUE)
sink()

##Parameters monitored
params <- c("mean.p", "a.0", "a.tmp", "a.tmp2", 
            "mean.abu", "b.a.0", "b.a.dis",
            "mean.pm", "b.m.0", "b.m.prd", "b.m.pdn",  "b.m.frg",
            "N", "N1", "N2", "pm.obs", "lambda", "N.res", "pm.res", "prd")

##MCMC
m <- jags(all.data, inits, params, paste0(model.name, ".txt"), 
          n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na,
          parallel=TRUE, n.cores=4)
assign(paste0("m.", model.name), m)

#Parameter estimates with Rhat
MCMCsummary(get(paste0("m.", model.name)), 
            params = c("mean.abu", "b.a.0", "b.a.dis",
                       "mean.pm", "b.m.0", "b.m.prd", "b.m.pdn",  "b.m.frg", "prd"),
            Rhat = TRUE,
            n.eff = TRUE,
            probs = c(0.025, 0.5, 0.975),
            round = 4)

#standardize the predator coefficient
m.prd.d$sims.list$prd
prd.sd <- apply(m.prd.d$sims.list$prd, 1, sd)
m.prd.d$sims.list$b.m.prd.original <- m.prd.d$sims.list$b.m.prd
m.prd.d$sims.list$b.m.prd <- m.prd.d$sims.list$b.m.prd * prd.sd


## univariate model -------------------------------------------

na <- 20000
ni <- 150000
nt <- 25
nb <- 125000
nc <- 4 

model.name <- "prd.u"
sink(paste0(model.name, ".txt"))
cat("
    
    model {
      
      ###Priors###
      
      #detection submodel parameters
      for(m in 1:2) {                 #m = morphs, 1 = melanic 2 = gray
        a.0[m] <- logit(mean.p[m])      #detection intercept
        mean.p[m] ~ dunif(0.01, 0.30)   #mean detection probabilty at mean temp
        a.tmp[m] ~ dnorm(0, 0.01)       #linear temp slope
        a.tmp2[m] ~ dnorm(0, 0.01)      #quadratic temp slope
      }
      
      #abundance submodel parameters
      mean.abu ~ dunif(1, 50)     #abundance at mean of covariates
      b.a.0 <- log(mean.abu)      #intercept
      b.a.dis ~ dnorm(0, 0.01)    #slopes
      
      #melanism submodel parameters
      mean.pm ~ dunif(0.01, 0.5)  #proportion melanic at value = 0 for covariates (mean distance)
      b.m.0 <- logit(mean.pm)     #mean p(melanic) at mean of covariates 
      b.m.prd ~ dnorm(0, 0.01)    #slopes
      
      #daily carnivore detection probability
      for (i in 1:nsites) {
        prd[i] ~ dunif(0, 1)  # uniform prior for parameter
        detections[i] ~ dbinom(prd[i], days[i])    # binomial likelihood for detections
      }
      
      ###Likelihood###
      
      #Process model for abundance & coat color
      for(i in 1:nsites) {
        N[i] ~ dpois(lambda[i])               #total squirrel abundance
        log(lambda[i]) <- b.a.0 +
          b.a.dis*dis[i]
        
        N1[i] ~ dbinom(pm[i], N[i])           #coat color process 
        logit(pm[i]) <- b.m.0 +
          b.m.prd*prd[i]
        
        #derived:
        N2[i] <- N[i] - N1[i]                 #abundance of gray morph
        pm.obs[i] <- N1[i]/ifelse(N[i]==0, 1e6, N[i])  #derived proportion melanic at each site
        N.res[i] <- N[i] - lambda[i]          #abundance residual
        pm.res[i] <- pm.obs[i] - pm[i]       #pm residual
        
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
    }
    
    ", fill=TRUE)
sink()

##Parameters monitored
params <- c("mean.p", "a.0", "a.tmp", "a.tmp2", 
            "mean.abu", "b.a.0", "b.a.dis",
            "mean.pm", "b.m.0", "b.m.prd",
            "N", "N1", "N2", "pm.obs", "lambda", "N.res", "pm.res", "prd")

##MCMC
m <- jags(all.data, inits, params, paste0(model.name, ".txt"), 
          n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na,
          parallel=TRUE, n.cores=4)
assign(paste0("m.", model.name), m)

m.prd.u <- update(m.prd.u, n.iter=25000, n.thin=25, parallel=TRUE)

#Parameter estimates with Rhat
MCMCsummary(get(paste0("m.", model.name)), 
            params = c("mean.abu", "b.a.0","b.a.dis",
                       "mean.pm", "b.m.0", "b.m.prd", "prd"),
            Rhat = TRUE,
            n.eff = TRUE,
            probs = c(0.025, 0.5, 0.975),
            round = 4)

#standardize the predator coefficient
m.prd.u$sims.list$prd
prd.sd <- apply(m.prd.u$sims.list$prd, 1, sd)
m.prd.u$sims.list$b.m.prd.original <- m.prd.u$sims.list$b.m.prd
m.prd.u$sims.list$b.m.prd <- m.prd.u$sims.list$b.m.prd * prd.sd