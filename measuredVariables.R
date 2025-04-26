# Load packages ----------------------------------------------------------------

library(dplyr)
library(dagitty)
library(brms) 

# Load data --------------------------------------------------------

load("measuredVariables.RData")

# Initial DAG ------------------------------------------------------------------

dag <- dagitty("dag {pdn -> bld; pdn -> imp; pdn -> tre; pdn -> frg; pdn -> U; pdn -> prd;
                     bld -> rdl; bld -> imp; bld -> mel;
                     rdl -> imp; rdl -> prd; rdl-> mel;
                     imp -> tre; imp -> frg;
                     tre -> frg; tre -> prd; tre -> mel;
                     frg -> prd; frg -> mel; 
                     prd -> mel;
                     U -> mel; 
                     U [unobserved];
}")


# Test links among measured variables -------------------------------------------

## pdn -> bld ;  -----------

### adjustment sets
adjustmentSets(dag, exposure = "pdn", outcome = "bld", effect="direct")

### model name
model <- "m.bld"

### bundle data
d <- cbind.data.frame(bld, pdn.log.s)

### specify model formula
assign(paste0(model, ".f"), bf(bld ~ 0 + Intercept + pdn.log.s, family = Beta()))  

### specify priors
assign(paste0(model, ".prior"), 
       c(prior(normal(0, 2.5), class = "b", coef="Intercept"), 
         prior(normal(0, 2), class = "b", coef = "pdn.log.s"),
         prior(gamma(0.01, 0.01), class = "phi")))

## fit model
assign(model, brm(formula = get(paste0(model, ".f")), 
                  data = d, 
                  prior = get(paste0(model, ".prior")),
                  chains = 4, iter = 3000, warmup = 2000, cores = 4))
get(model) #output

## bld -> rdl;-----------

### adjustment sets
adjustmentSets(dag, exposure = "bld", outcome = "rdl", effect="direct")

### model name
model <- "m.rdl.bld"

### bundle data
d <- cbind.data.frame(rdl.c, bld.logit.s)

### specify model formula
assign(paste0(model, ".f"), bf(rdl.c ~ 0 + Intercept + bld.logit.s, family = Beta()))  

### specify priors
assign(paste0(model, ".prior"), 
       c(prior(normal(0, 2.5), class = "b", coef="Intercept"), 
         prior(normal(0, 2), class = "b", coef = "bld.logit.s"),
         prior(gamma(0.01, 0.01), class = "phi")))

## fit model
assign(model, brm(formula = get(paste0(model, ".f")), 
                  data = d, 
                  prior = get(paste0(model, ".prior")),
                  chains = 4, iter = 3000, warmup = 2000, cores = 4))
get(model) #output

## bld -> imp ------

adjustmentSets(dag, exposure = "bld", outcome = "imp", effect="direct")

### model name
model <- "m.imp.bld"

### bundle data
d <- cbind.data.frame(imp, rdl.s, bld.logit.s, pdn.log.s)

### specify model formula
assign(paste0(model, ".f"), bf(imp ~ 0 + Intercept + bld.logit.s + pdn.log.s + rdl.s, family = Beta()))  

### specify priors
assign(paste0(model, ".prior"), 
       c(prior(normal(0, 2), class = "b", coef="Intercept"), 
         prior(normal(0, 2), class = "b", coef = "bld.logit.s"),
         prior(normal(0, 2), class = "b", coef = "rdl.s"),
         prior(normal(0, 2), class = "b", coef = "pdn.log.s"),
         prior(gamma(0.01, 0.01), class = "phi")))

## fit model
assign(model, brm(formula = get(paste0(model, ".f")), 
                  data = d, 
                  prior = get(paste0(model, ".prior")),
                  chains = 4, iter = 3000, warmup = 2000, cores = 4))
get(model) #output


## rdl -> imp; -----
adjustmentSets(dag, exposure = "rdl", outcome = "imp", effect="direct")

### model name
model <- "m.imp.rdl"

### bundle data
d <- cbind.data.frame(imp, rdl.s, bld.logit.s)

### specify model formula
assign(paste0(model, ".f"), bf(imp ~ 0 + Intercept + bld.logit.s + rdl.s, family = Beta()))  

### specify priors
assign(paste0(model, ".prior"), 
       c(prior(normal(0, 2), class = "b", coef="Intercept"), 
         prior(normal(0, 2), class = "b", coef = "bld.logit.s"),
         prior(normal(0, 2), class = "b", coef = "rdl.s"),
         prior(gamma(0.01, 0.01), class = "phi")))

## fit model
assign(model, brm(formula = get(paste0(model, ".f")), 
                  data = d, 
                  prior = get(paste0(model, ".prior")),
                  chains = 4, iter = 3000, warmup = 2000, cores = 4))
get(model) #output

## pdn -> imp; -----
adjustmentSets(dag, exposure = "pdn", outcome = "imp", effect="direct")

### model name
model <- "m.imp.pdn"

### bundle data
d <- cbind.data.frame(imp, pdn.log.s, bld.logit.s)

### specify model formula
assign(paste0(model, ".f"), bf(imp ~ 0 + Intercept + bld.logit.s + pdn.log.s, family = Beta()))  

### specify priors
assign(paste0(model, ".prior"), 
       c(prior(normal(0, 2), class = "b", coef="Intercept"), 
         prior(normal(0, 2), class = "b", coef = "pdn.log.s"),
         prior(normal(0, 2), class = "b", coef = "bld.logit.s"),
         prior(gamma(0.01, 0.01), class = "phi")))

## fit model
assign(model, brm(formula = get(paste0(model, ".f")), 
                  data = d, 
                  prior = get(paste0(model, ".prior")),
                  chains = 4, iter = 3000, warmup = 2000, cores = 4))
get(model) #output


## imp -> tre; pdn -> tre -----------

### adjustment sets
adjustmentSets(dag, exposure = "imp", outcome = "tre", effect="direct")
adjustmentSets(dag, exposure = "pdn", outcome = "tre", effect="direct")

### model name
model <- "m.tre"

### bundle data
d <- cbind.data.frame(tre, imp.logit.s, pdn.log.s)

### specify model formula
assign(paste0(model, ".f"), bf(tre ~ 0 + Intercept + imp.logit.s + pdn.log.s, family = Beta()))  

### specify priors
assign(paste0(model, ".prior"), 
       c(prior(normal(0, 2.5), class = "b", coef="Intercept"), 
         prior(normal(0, 2), class = "b", coef = "imp.logit.s"),
         prior(normal(0, 2), class = "b", coef = "pdn.log.s"),
         prior(gamma(0.01, 0.01), class = "phi")))

## fit model
assign(model, brm(formula = get(paste0(model, ".f")), 
                  data = d, 
                  prior = get(paste0(model, ".prior")),
                  chains = 4, iter = 3000, warmup = 2000, cores = 4))
get(model) #output


## imp -> frg; tre -> frg ; pdn -> frg -----------

### adjustment sets
adjustmentSets(dag, exposure = "imp", outcome = "frg", effect="direct")
adjustmentSets(dag, exposure = "tre", outcome = "frg", effect="direct")
adjustmentSets(dag, exposure = "pdn", outcome = "frg", effect="direct")

### model name
model <- "m.frg"

### bundle data
d <- cbind.data.frame(frg.log.s, imp.logit.s, tre.s, pdn.log.s)

### specify model formula
assign(paste0(model, ".f"), bf(frg.log.s ~ 0 + Intercept + imp.logit.s + tre.s + pdn.log.s, family = gaussian()))  

### specify priors
assign(paste0(model, ".prior"), 
       c(prior(normal(0, 1), class = "b", coef="Intercept"), 
         prior(normal(0, 2), class = "b", coef = "imp.logit.s"),
         prior(normal(0, 2), class = "b", coef = "pdn.log.s"),
         prior(normal(0, 2), class = "b", coef = "tre.s"),
         prior(student_t(3, 0, 2.5), class="sigma")))

## fit model
assign(model, brm(formula = get(paste0(model, ".f")), 
                  data = d, 
                  prior = get(paste0(model, ".prior")),
                  chains = 4, iter = 3000, warmup = 2000, cores = 4))
get(model) #output

## rdl -> prd; tre -> prd; frg -> prd; pdn -> prd -----------

adjustmentSets(dag, exposure = "pdn", outcome = "prd", effect="direct")
adjustmentSets(dag, exposure = "rdl", outcome = "prd", effect="direct")
adjustmentSets(dag, exposure = "tre", outcome = "prd", effect="direct")
adjustmentSets(dag, exposure = "frg", outcome = "prd", effect="direct")

### model name
model <- "m.prd"

### bundle data
d <- cbind.data.frame(pdn.log.s, bld.logit.s, detections = rai.sum$detections, days = rai.sum$active.days,
                      tre.s, rdl.s, frg.log.s, pdn.log.s)

### specify model formula
assign(paste0(model, ".f"), bf(detections | trials(days) ~ 0 + Intercept + 
                                 frg.log.s + tre.s + rdl.s + pdn.log.s , family = binomial()))  

### specify priors
assign(paste0(model, ".prior"), 
       c(prior(normal(0, 2), class = "b", coef="Intercept"), 
         prior(normal(0, 2), class = "b", coef = "tre.s"),
         prior(normal(0, 2), class = "b", coef = "rdl.s"),
         prior(normal(0, 2), class = "b", coef = "pdn.log.s"),
         prior(normal(0, 2), class = "b", coef = "frg.log.s")))

## fit model
assign(model, brm(formula = get(paste0(model, ".f")), 
                  data = d, 
                  prior = get(paste0(model, ".prior")),
                  chains = 4, iter = 3000, warmup = 2000, cores = 4))
get(model) #output


# Revised DAG  ----------------------------------

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
