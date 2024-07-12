#####Load packages#####
#-------------------------#

library(dplyr)
library(jagsUI)
library(MCMCvis)

#####Load site covariates#####
#-------------------------#
sites <- read.csv("siteCovariates.csv", header=T)
sites$site <- as.factor(sites$site)

#####Load predator data and create activity index#####
#-------------------------#
rai <- read.csv("predators.csv", header=T)
rai <- rai[,-1]
rai$site <- as.factor(rai$site)

##Combine detection histories for predators of interest
rai.sum <- as.data.frame(rai %>% group_by(site) %>%
                           summarise(across(everything(), sum)))
rai.sum <- rai.sum %>% mutate_if(is.numeric, ~1 * (. > 0))

##Quantify relative activity index
rai.sum$detections <- rowSums(rai.sum[,2:333], na.rm=TRUE)
rai.sum$active.days <- rowSums(!is.na(rai.sum[,2:333]))
rai.sum$rai <- (rai.sum$detections/rai.sum$active.days)*100
rai.sum$site <- factor(rai.sum$site)

#####Prep covariates#####
#-------------------------#

#distance to city center
dis <- sites$DistCityCen_km
dis.s <- (dis - mean(dis))/sd(dis)

#pop density with log transformation, 300-m scale
pdn <- sites$pdn300m
pdn.log <- log(pdn)
pdn.log.s <- (pdn.log - mean(pdn.log))/sd(pdn.log)

#building cover, 300-m scale
bld <- sites$bld300m/(pi*300^2) #building area
bld.logit <- log((bld + 0.0004311796)/(1 - bld + 0.0004311796))
bld.s <- (bld.logit - mean(bld.logit))/sd(bld.logit)

#road length with speed>=30 mph, 300-m scale with log transformation
rod <- log(sites$rod30near)
rod.s <- (rod - mean(rod))/sd(rod)

#tree cover, 300-m scale
tre <- sites$tre300m
tre.s <- (tre - mean(tre))/sd(tre)

#fragmentation with log transformation
frg <- sites$dcores1000m
frg.log <- log(frg)
frg.log.s <- (frg.log - mean(frg.log))/sd(frg.log)

#predators
prd <- rai.sum$rai
prd.s <- (prd - mean(prd))/sd(prd)

#save urban covariate data
save(dis, dis.s, pdn, pdn.log, pdn.log.s, 
     bld, bld.logit, bld.s, rod, rod.s, tre, tre.s,
     frg, frg.log, frg.log.s, prd, prd.s,
     file = "urbanCovariates.RData")

#####Bundle data for JAGS#####
#-------------------------#
urbanCov.data <- list(dis = dis.s, pdn = pdn.log.s, bld = bld.s,
                 tre = tre.s, rod = rod.s, 
                 frg = frg.log.s, prd = prd.s, nsites = nrow(sites))


#####Fit initial urbanization submodel#####
#-------------------------#

sink("urbanSubmodel.txt")
cat("

model {

#####PATHS FOR ENVIRONMENTAL COVARIATES#####

###human population density###
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

###distance to road###
rod.tau <- pow(rod.sd, -2)                        #precision for observed variable
rod.sd ~ dt(0, 1/pow(25, 2), 1)I(0,)              #prior for standard deviation

b.rod.0 ~ dnorm(0, 0.01)                          #prior for intercept
b.rod.pdn ~ dnorm(0, 0.01)                        #prior for slopes
b.rod.bld ~ dnorm(0, 0.01)                        

for(i in 1:nsites){                               #regression model
rod[i] ~ dnorm(rod.mu[i], rod.tau)                #stochastic component
rod.mu[i] <- b.rod.0 +                            #deterministic component
             b.rod.pdn*pdn[i] +
             b.rod.bld*bld[i]
rod.res[i] <- rod[i] - rod.mu[i]                  #residuals
}

###forest cover###
tre.tau <- pow(tre.sd, -2)                        #precision for observed variable
tre.sd ~ dt(0, 1/pow(25, 2), 1)I(0,)              #prior for standard deviation

b.tre.0 ~ dnorm(0, 0.01)                          #prior for intercept
b.tre.bld ~ dnorm(0, 0.01)
b.tre.rod ~ dnorm(0, 0.01)                        
                 
for(i in 1:nsites){                               #regression model
tre[i] ~ dnorm(tre.mu[i], tre.tau)                #stochastic component
tre.mu[i] <- b.tre.0 +                            #deterministic component
             b.tre.bld*bld[i] +
             b.tre.rod*rod[i]
tre.res[i] <- tre[i] - tre.mu[i]                  #residuals
}

###fragmentation###
frg.tau <- pow(frg.sd, -2)                        #precision for observed variable
frg.sd ~ dt(0, 1/pow(25, 2), 1)I(0,)              #prior for standard deviation

b.frg.0 ~ dnorm(0, 0.01)                          #prior for intercept
b.frg.bld ~ dnorm(0, 0.01)                        #prior for slopes
b.frg.rod ~ dnorm(0, 0.01)                        
b.frg.tre ~ dnorm(0, 0.01)                        

for(i in 1:nsites){                               #regression model
frg[i] ~ dnorm(frg.mu[i], frg.tau)                #stochastic component
frg.mu[i] <- b.frg.0 +                            #deterministic component
             b.frg.bld*bld[i] +
             b.frg.rod*rod[i] + 
             b.frg.tre*tre[i]
frg.res[i] <- frg[i] - frg.mu[i]                  #residuals
}

###predators###

prd.tau <- pow(prd.sd, -2)                        #precision for observed variable
prd.sd ~ dt(0, 1/pow(25, 2), 1)I(0,)              #prior for standard deviation
  
b.prd.0 ~ dnorm(0, 0.01)                          #prior for intercept
b.prd.pdn ~ dnorm(0, 0.01)                        #prior for slope
b.prd.bld ~ dnorm(0, 0.01)  
b.prd.rod ~ dnorm(0, 0.01) 
b.prd.tre ~ dnorm(0, 0.01)                        
b.prd.frg ~ dnorm(0, 0.01)    

for(i in 1:nsites){                               #regression model
prd[i] ~ dnorm(prd.mu[i], prd.tau)                #stochastic component
prd.mu[i] <- b.prd.0 +                            #deterministic component
             b.prd.pdn*pdn[i] +
             b.prd.bld*bld[i] +
             b.prd.rod*rod[i] + 
             b.prd.tre*tre[i] +
             b.prd.frg*frg[i]
prd.res[i] <- prd[i] - prd.mu[i]                  #residuals
}

}

    ", fill=TRUE)
sink()

##Initial values----

inits <- function(){list(pdn.sd = 1, bld.sd = 1, rod.sd = 1, tre.sd = 1, frg.sd =1, prd.sd = 1)}

##Parameters monitored
params <- c("pdn.sd", "pdn.mu", "pdn.res", 
            "bld.sd", "bld.mu", "bld.res", 
            "rod.sd", "rod.mu", "rod.res", 
            "tre.sd", "tre.mu", "tre.res", 
            "frg.sd", "frg.mu", "frg.res", 
            "prd.sd", "prd.mu", "prd.res", 
            "b.pdn.0", "b.pdn.dis",
            "b.bld.0", "b.bld.pdn",
            "b.rod.0", "b.rod.pdn", "b.rod.bld",
            "b.tre.0", "b.tre.bld", "b.tre.rod", 
            "b.frg.0", "b.frg.bld", "b.frg.rod", "b.frg.tre",
            "b.prd.0", "b.prd.pdn", "b.prd.bld", "b.prd.rod", "b.prd.tre", "b.prd.frg")

##MCMC
ni <- 4000 ; nt <- 3 ;  nb <- 1000 ; nc <- 3 #1 min

m1 <- jags(urbanCov.data, inits, params, "urbanSubmodel.txt", 
           n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=20000,
           parallel=TRUE, n.cores=3)

#parameter estimates with Rhat
MCMCsummary(m1, 
            params = c("b.pdn.0", "b.pdn.dis",
                       "b.bld.0", "b.bld.pdn",
                       "b.rod.0", "b.rod.pdn", "b.rod.bld",
                       "b.tre.0", "b.tre.bld", "b.tre.rod", 
                       "b.frg.0", "b.frg.bld", "b.frg.rod", "b.frg.tre",
                       "b.prd.0", "b.prd.pdn", "b.prd.bld", "b.prd.rod", "b.prd.tre", "b.prd.frg"),
            Rhat = TRUE,
            n.eff = TRUE,
            probs = c(0.025, 0.5, 0.975),
            round = 2)



#####D-separation tests#####
#-------------------------#

###bld ~ dis

#conditioning set: pdn

sink("urbanSubmodel_dsep.txt")
cat("

model {

bld.tau <- pow(bld.sd, -2)                      #precision for observed variable
bld.sd ~ dt(0, 1/pow(25, 2), 1)I(0,)             #prior for standard deviation
  
b.bld.0 ~ dnorm(0, 0.01)                          #prior for intercept
b.bld.dis ~ dnorm(0, 0.01)
b.bld.pdn ~ dnorm(0, 0.01)                        
  
  for(i in 1:nsites){                               #regression model
    bld[i] ~ dnorm(bld.mu[i], bld.tau)                #stochastic component
    bld.mu[i] <-b.bld.0 +                            #deterministic component
                b.bld.dis*dis[i] +
                b.bld.pdn*pdn[i]
  }
}
    ", fill=TRUE)
sink()

params <- c("b.bld.dis", "b.bld.pdn")

ni <- 4000 ; nt <- 3 ;  nb <- 1000 ; nc <- 3 

ds.bld.dis <- jags(urbanCov.data, inits, params, "urbanSubmodel_dsep.txt", 
                   n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=5000,
                   parallel=TRUE, n.cores=3)

MCMCsummary(ds.bld.dis, params = params, Rhat = TRUE, n.eff = TRUE,
            probs = c(0.05, 0.5, 0.95), round = 2)

###rod ~ dis 

#conditioning set: pdn, bld

sink("urbanSubmodel_dsep.txt")
cat("

model {

rod.tau <- pow(rod.sd, -2)                      #precision for observed variable
rod.sd ~ dt(0, 1/pow(25, 2), 1)I(0,)             #prior for standard deviation
  
b.rod.0 ~ dnorm(0, 0.01)                          #prior for intercept
b.rod.dis ~ dnorm(0, 0.01)
b.rod.pdn ~ dnorm(0, 0.01)                        
b.rod.bld ~ dnorm(0, 0.01)
  
  for(i in 1:nsites){                               #regression model
    rod[i] ~ dnorm(rod.mu[i], rod.tau)                #stochastic component
    rod.mu[i] <-b.rod.0 +                            #deterministic component
                b.rod.dis*dis[i] +
                b.rod.pdn*pdn[i] +
                b.rod.bld*bld[i]
  }
}
    ", fill=TRUE)
sink()

params <- c("b.rod.dis", "b.rod.pdn", "b.rod.bld")

ni <- 4000 ; nt <- 3 ;  nb <- 1000 ; nc <- 3 

ds.rod.dis <- jags(urbanCov.data, inits, params, "urbanSubmodel_dsep.txt", 
                   n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=5000,
                   parallel=TRUE, n.cores=3)

MCMCsummary(ds.rod.dis, params = params, Rhat = TRUE, n.eff = TRUE,
            probs = c(0.05, 0.5, 0.95), round = 2)

###tre ~ dis

#conditioning set: bld, rod

sink("urbanSubmodel_dsep.txt")
cat("

model {

tre.tau <- pow(tre.sd, -2)                      #precision for observed variable
tre.sd ~ dt(0, 1/pow(25, 2), 1)I(0,)             #prior for standard deviation
  
b.tre.0 ~ dnorm(0, 0.01)                          #prior for intercept
b.tre.dis ~ dnorm(0, 0.01)
b.tre.bld ~ dnorm(0, 0.01)
b.tre.rod ~ dnorm(0, 0.01)
  
  for(i in 1:nsites){                               #regression model
    tre[i] ~ dnorm(tre.mu[i], tre.tau)                #stochastic component
    tre.mu[i] <-b.tre.0 +                            #deterministic component
                b.tre.dis*dis[i] +
                b.tre.bld*bld[i] +
                b.tre.rod*rod[i]
  }
}
    ", fill=TRUE)
sink()

params <- c("b.tre.dis", "b.tre.bld", "b.tre.rod")

ni <- 4000 ; nt <- 3 ;  nb <- 1000 ; nc <- 3 

ds.tre.dis <- jags(urbanCov.data, inits, params, "urbanSubmodel_dsep.txt", 
                   n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=5000,
                   parallel=TRUE, n.cores=3)

###tre ~ pdn

#conditioning set: dis, bld, rod

sink("urbanSubmodel_dsep.txt")
cat("

model {

tre.tau <- pow(tre.sd, -2)                      #precision for observed variable
tre.sd ~ dt(0, 1/pow(25, 2), 1)I(0,)             #prior for standard deviation
  
b.tre.0 ~ dnorm(0, 0.01)                          #prior for intercept
b.tre.pdn ~ dnorm(0, 0.01)
b.tre.bld ~ dnorm(0, 0.01)
b.tre.rod ~ dnorm(0, 0.01)
b.tre.dis ~ dnorm(0, 0.01)
  
  for(i in 1:nsites){                               #regression model
    tre[i] ~ dnorm(tre.mu[i], tre.tau)                #stochastic component
    tre.mu[i] <-b.tre.0 +                            #deterministic component
                b.tre.pdn*pdn[i] +
                b.tre.bld*bld[i] +
                b.tre.rod*rod[i] +
                b.tre.dis*dis[i]
  }
}
    ", fill=TRUE)
sink()

params <- c("b.tre.pdn", "b.tre.bld", "b.tre.rod", "b.tre.dis")

ni <- 4000 ; nt <- 3 ;  nb <- 1000 ; nc <- 3 #20 min

ds.tre.pdn <- jags(urbanCov.data, inits, params, "urbanSubmodel_dsep.txt", 
                   n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=5000,
                   parallel=TRUE, n.cores=3)

MCMCsummary(ds.tre.pdn, params = params, Rhat = TRUE, n.eff = TRUE,
            probs = c(0.05, 0.5, 0.95), round = 2)

###frg ~ dis

#conditioning set: bld, rod, tre

sink("urbanSubmodel_dsep.txt")
cat("

model {

frg.tau <- pow(frg.sd, -2)                      #precision for observed variable
frg.sd ~ dt(0, 1/pow(25, 2), 1)I(0,)             #prior for standard deviation
  
b.frg.0 ~ dnorm(0, 0.01)                          #prior for intercept
b.frg.dis ~ dnorm(0, 0.01)
b.frg.tre ~ dnorm(0, 0.01)                        
b.frg.bld ~ dnorm(0, 0.01)
b.frg.rod ~ dnorm(0, 0.01)
  
  for(i in 1:nsites){                               #regression model
    frg[i] ~ dnorm(frg.mu[i], frg.tau)                #stochastic component
    frg.mu[i] <-b.frg.0 +                            #deterministic component
                b.frg.dis*dis[i] +
                b.frg.tre*tre[i] +
                b.frg.bld*bld[i] +
                b.frg.rod*rod[i]
  }
}
    ", fill=TRUE)
sink()

params <- c("b.frg.dis", "b.frg.tre", "b.frg.bld", "b.frg.rod")

ni <- 4000 ; nt <- 3 ;  nb <- 1000 ; nc <- 3 #20 min

ds.frg.dis <- jags(urbanCov.data, inits, params, "urbanSubmodel_dsep.txt", 
                   n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=5000,
                   parallel=TRUE, n.cores=3)

MCMCsummary(ds.frg.dis, params = params, Rhat = TRUE, n.eff = TRUE,
            probs = c(0.05, 0.5, 0.95), round = 2)

###frg ~ pdn

#conditioning set: dis, bld, rod, tre

sink("urbanSubmodel_dsep.txt")
cat("

model {

frg.tau <- pow(frg.sd, -2)                      #precision for observed variable
frg.sd ~ dt(0, 1/pow(25, 2), 1)I(0,)             #prior for standard deviation
  
b.frg.0 ~ dnorm(0, 0.01)                          #prior for intercept
b.frg.bld ~ dnorm(0, 0.01)                        #slopes
b.frg.rod ~ dnorm(0, 0.01)
b.frg.tre ~ dnorm(0, 0.01)
b.frg.pdn ~ dnorm(0, 0.01)
b.frg.dis ~ dnorm(0, 0.01)
  
  for(i in 1:nsites){                               #regression model
    frg[i] ~ dnorm(frg.mu[i], frg.tau)                #stochastic component
    frg.mu[i] <-b.frg.0 +                            #deterministic component
                b.frg.pdn*pdn[i] +
                b.frg.bld*bld[i] +
                b.frg.rod*rod[i] + 
                b.frg.tre*tre[i] + 
                b.frg.dis*dis[i]
  }
}
    ", fill=TRUE)
sink()

params <- c("b.frg.pdn", "b.frg.bld", "b.frg.rod", "b.frg.tre", "b.frg.dis")

#MCMC
ni <- 4000 ; nt <- 3 ;  nb <- 1000 ; nc <- 3

ds.frg.pdn <- jags(urbanCov.data, inits, params, "urbanSubmodel_dsep.txt", 
              n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=5000,
              parallel=TRUE, n.cores=3)

MCMCsummary(ds.frg.pdn, params = params, Rhat = TRUE, n.eff = TRUE,
            probs = c(0.05, 0.5, 0.95), round = 2)


###prd ~ dis

#conditioning set: bld, rod, tre

sink("urbanSubmodel_dsep.txt")
cat("

model {

prd.tau <- pow(prd.sd, -2)                      #precision for observed variable
prd.sd ~ dt(0, 1/pow(25, 2), 1)I(0,)             #prior for standard deviation
  
b.prd.0 ~ dnorm(0, 0.01)                          #prior for intercept
b.prd.dis ~ dnorm(0, 0.01)
b.prd.pdn ~ dnorm(0, 0.01)
b.prd.tre ~ dnorm(0, 0.01)                        
b.prd.bld ~ dnorm(0, 0.01)
b.prd.rod ~ dnorm(0, 0.01)
b.prd.frg ~ dnorm(0, 0.01)
  
  for(i in 1:nsites){                               #regression model
    prd[i] ~ dnorm(prd.mu[i], prd.tau)                #stochastic component
    prd.mu[i] <-b.prd.0 +                            #deterministic component
                b.prd.dis*dis[i] +
                b.prd.pdn*pdn[i] +
                b.prd.tre*tre[i] +
                b.prd.bld*bld[i] +
                b.prd.rod*rod[i] +
                b.prd.frg*frg[i]
  }
}
    ", fill=TRUE)
sink()

params <- c("b.prd.dis", "b.prd.pdn", "b.prd.tre", "b.prd.bld", "b.prd.rod", "b.prd.frg")

ni <- 4000 ; nt <- 3 ;  nb <- 1000 ; nc <- 3 

ds.prd.dis <- jags(urbanCov.data, inits, params, "urbanSubmodel_dsep.txt", 
                   n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=5000,
                   parallel=TRUE, n.cores=3)

MCMCsummary(ds.prd.dis, params = params, Rhat = TRUE, n.eff = TRUE,
            probs = c(0.05, 0.5, 0.95), round = 2)



#####Revised urbanization submodel - add missing links#####
#-------------------------#

##Model in BUGS language----
sink("urbanSubmodel_addedLinks.txt")
cat("

model {

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

###distance to road###
rod.tau <- pow(rod.sd, -2)                        #precision for observed variable
rod.sd ~ dt(0, 1/pow(25, 2), 1)I(0,)              #prior for standard deviation

b.rod.0 ~ dnorm(0, 0.01)                          #prior for intercept
b.rod.pdn ~ dnorm(0, 0.01)                        #prior for slopes
b.rod.bld ~ dnorm(0, 0.01)                        

for(i in 1:nsites){                               #regression model
rod[i] ~ dnorm(rod.mu[i], rod.tau)                #stochastic component
rod.mu[i] <- b.rod.0 +                            #deterministic component
             b.rod.pdn*pdn[i] +
             b.rod.bld*bld[i]
rod.res[i] <- rod[i] - rod.mu[i]                  #residuals
}

###forest cover###
tre.tau <- pow(tre.sd, -2)                        #precision for observed variable
tre.sd ~ dt(0, 1/pow(25, 2), 1)I(0,)              #prior for standard deviation

b.tre.0 ~ dnorm(0, 0.01)                          #prior for intercept
b.tre.bld ~ dnorm(0, 0.01)
b.tre.rod ~ dnorm(0, 0.01)                        
                 
for(i in 1:nsites){                               #regression model
tre[i] ~ dnorm(tre.mu[i], tre.tau)                #stochastic component
tre.mu[i] <- b.tre.0 +                            #deterministic component
             b.tre.bld*bld[i] +
             b.tre.rod*rod[i]
tre.res[i] <- tre[i] - tre.mu[i]                  #residuals
}

###fragmentation###
frg.tau <- pow(frg.sd, -2)                        #precision for observed variable
frg.sd ~ dt(0, 1/pow(25, 2), 1)I(0,)              #prior for standard deviation

b.frg.0 ~ dnorm(0, 0.01)                          #prior for intercept
b.frg.dis ~ dnorm(0, 0.01)                        #prior for slopes
b.frg.pdn ~ dnorm(0, 0.01)
b.frg.bld ~ dnorm(0, 0.01)                        
b.frg.rod ~ dnorm(0, 0.01)                        
b.frg.tre ~ dnorm(0, 0.01)                        

for(i in 1:nsites){                               #regression model
frg[i] ~ dnorm(frg.mu[i], frg.tau)                #stochastic component
frg.mu[i] <- b.frg.0 +                            #deterministic component
             b.frg.dis*dis[i] +
             b.frg.pdn*pdn[i] +
             b.frg.bld*bld[i] +
             b.frg.rod*rod[i] + 
             b.frg.tre*tre[i] 
frg.res[i] <- frg[i] - frg.mu[i]                  #residuals
}

###predators###
prd.tau <- pow(prd.sd, -2)                        #precision for observed variable
prd.sd ~ dt(0, 1/pow(25, 2), 1)I(0,)              #prior for standard deviation
  
b.prd.0 ~ dnorm(0, 0.01)                          #prior for intercept
b.prd.pdn ~ dnorm(0, 0.01)                        #prior for slope
b.prd.bld ~ dnorm(0, 0.01)  
b.prd.rod ~ dnorm(0, 0.01) 
b.prd.tre ~ dnorm(0, 0.01)                        
b.prd.frg ~ dnorm(0, 0.01)    

for(i in 1:nsites){                               #regression model
prd[i] ~ dnorm(prd.mu[i], prd.tau)                #stochastic component
prd.mu[i] <- b.prd.0 +                            #deterministic component
             b.prd.pdn*pdn[i] +
             b.prd.bld*bld[i] +
             b.prd.rod*rod[i] + 
             b.prd.tre*tre[i] +
             b.prd.frg*frg[i]
prd.res[i] <- prd[i] - prd.mu[i]                  #residuals
}

}

    ", fill=TRUE)
sink()

##Initial values----

inits <- function(){list(pdn.sd = 1, bld.sd = 1, rod.sd = 1, tre.sd = 1, frg.sd =1, prd.sd = 1)}

##Parameters monitored
params <- c("pdn.sd", "pdn.mu", "pdn.res", 
            "bld.sd", "bld.mu", "bld.res", 
            "rod.sd", "rod.mu", "rod.res", 
            "tre.sd", "tre.mu", "tre.res", 
            "frg.sd", "frg.mu", "frg.res", 
            "prd.sd", "prd.mu", "prd.res", 
            "b.pdn.0", "b.pdn.dis",
            "b.bld.0", "b.bld.pdn",
            "b.rod.0", "b.rod.pdn", "b.rod.bld",
            "b.tre.0", "b.tre.bld", "b.tre.rod", 
            "b.frg.0", "b.frg.dis", "b.frg.pdn", "b.frg.bld", "b.frg.rod", "b.frg.tre",
            "b.prd.0", "b.prd.pdn", "b.prd.bld", "b.prd.rod", "b.prd.tre", "b.prd.frg")

##MCMC
ni <- 4000 ; nt <- 3 ;  nb <- 1000 ; nc <- 3 

m2 <- jags(urbanCov.data, inits, params, "urbanSubmodel_addedLinks.txt", 
           n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=20000,
           parallel=TRUE, n.cores=3)

#parameter estimates with Rhat
MCMCsummary(m2, 
            params = c("b.pdn.0", "b.pdn.dis",
                       "b.bld.0", "b.bld.pdn",
                       "b.rod.0", "b.rod.pdn", "b.rod.bld",
                       "b.tre.0", "b.tre.bld", "b.tre.rod", 
                       "b.frg.0", "b.frg.dis", "b.frg.pdn", "b.frg.bld", "b.frg.rod", "b.frg.tre",
                       "b.prd.0", "b.prd.pdn", "b.prd.bld", "b.prd.rod", "b.prd.tre", "b.prd.frg"),
            Rhat = TRUE,
            n.eff = TRUE,
            probs = c(0.05, 0.5, 0.95),
            round = 3)


#####Drop terms with 90% CI including 0#####
#-------------------------#

##Model in BUGS language----
sink("urbanSubmodel_pruned.txt")
cat("

model {

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

###forest cover###
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

}

    ", fill=TRUE)
sink()

##Initial values----

inits <- function(){list(pdn.sd = 1, bld.sd = 1, tre.sd = 1, frg.sd =1, prd.sd = 1)}

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
            "b.prd.0", "b.prd.pdn")

##MCMC
ni <- 4000 ; nt <- 3 ;  nb <- 1000 ; nc <- 3 #6 min

m3 <- jags(urbanCov.data, inits, params, "urbanSubmodel_pruned.txt", 
           n.chains=nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=20000,
           parallel=TRUE, n.cores=3)

#parameter estimates with Rhat
MCMCsummary(m3, 
            params = c("b.pdn.0", "b.pdn.dis",
                       "b.bld.0", "b.bld.pdn",
                       "b.tre.0", "b.tre.bld", 
                       "b.frg.0", "b.frg.dis", "b.frg.pdn", "b.frg.tre",
                       "b.prd.0", "b.prd.pdn"),
            Rhat = TRUE,
            n.eff = TRUE,
            probs = c(0.025, 0.5, 0.975),
            round = 2)