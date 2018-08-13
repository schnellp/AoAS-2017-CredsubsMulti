# setwd("/home/patrick/Documents/umn/abbvie/credsubs-multi/repo/R/")

library(nimble)
library(matrixStats)

alzCode <- nimbleCode({ 
  # likelihood
  for (i in 1:N){
    for (j in 1:2) {
      # wasn't compiling with linear algebra expression
      # maybe couldn't drop extra dimensions
      eta[i, j] <- X[i, 1] * beta[j, 1] +
        X[i, 2] * beta[j, 2] +
        X[i, 3] * beta[j, 3] +
        X[i, 4] * beta[j, 4] +
        X[i, 5] * beta[j, 5] +
        d[t[i]] * (
          Z[i, 1] * gamma[t[i], j, 1] +
            Z[i, 2] * gamma[t[i], j, 2] +
            Z[i, 3] * gamma[t[i], j, 3] +
            Z[i, 4] * gamma[t[i], j, 4] +
            Z[i, 5] * gamma[t[i], j, 5]
        )
      #eta.mat <- (X[i, 1:P] %*% beta[j, 1:P] + 
      #              d[t[i]] * Z[i, 1:P] %*% gamma[t[i], j, 1:P])
      #eta[i, j] <- eta.mat[1, 1]
    }
    
    # efficacy
    y[i, 1] ~ dnorm(mean=eta[i, 1], tau=sqrt(tau.sq))
    # safety
    y[i, 2] ~ dbern(p=1 / (1 + exp(-eta[i, 2])))
  }
  
  # prior
  tau.sq ~ dgamma(shape=0.001, rate=0.001)
  d[1] <- 0
  d[2] ~ dunif(0, d[3])
  d[3] ~ dunif(0, 1)
  d[4] <- 1
  d[5] <- 1
  
  # regression coefficients
  for (j in 1:2) {
    for (p in 1:P) {
      beta[j, p] ~ dnorm(mean=0, tau=1/10^4)
      gamma[1, j, p] <- 1
      gamma[2, j, p] <- gamma[4, j, p]
      gamma[3, j, p] <- gamma[4, j, p]
      
      gamma[4, j, p] ~ dnorm(mean=0, tau=1/10^4)
      gamma[5, j, p] ~ dnorm(mean=0, tau=1/10^4)
    }
  }
})

set.seed(1)

load("data-simple.RData")

alzConsts <- list(N = nrow(data.simple),
                  P = 5,
                  t = data.simple$treatment + 1)

covariates <- scale(as.matrix(data.simple[, c("severity", "age", "sex", "carrier")]))
X <- cbind(1, covariates)
y = as.matrix(data.simple[, c("improve", "ae")])
y[, "improve"] <- scale(y[, "improve"])

alzData <- list(X = X,
                Z = X,
                y = y)

alzInits <- list(d = c(0, 0.25, 0.5, 1, 1),
                  beta=matrix(0, nrow=2, ncol=5),
                 gamma=array(0, c(5, 2, 5)), # arm, endpoint, predictor
                 tau.sq = 1)

alzModel <- nimbleModel(code = alzCode, name = 'alz', constants = alzConsts,
                   data = alzData,
                   inits = alzInits)

alzSpec <- configureMCMC(alzModel)
alzSpec$addMonitors('gamma')

alzMCMC <- buildMCMC(alzSpec)

CalzModel <- compileNimble(alzModel)
CalzMCMC <- compileNimble(alzMCMC, project = alzModel)

CalzMCMC$run(11000)

MCMCsamples <- as.matrix(CalzMCMC$mvSamples)[1001:11000, ]

save(MCMCsamples, file="nimble-standard-sample-select.RData")

for (k in 1:2) {
  for (p in 2:5) {
    cname1 <- paste0("gamma[4, ", k, ", ", p, "]")
    cname2 <- paste0("gamma[5, ", k, ", ", p, "]")
    print(paste(cname1, "-", cname2))
    col <- MCMCsamples[, cname1] - MCMCsamples[, cname2]
    print(mean(col) - sd(col) * qnorm(0.75))
    print(mean(col) - sd(col) * qnorm(0.25))
  }
}

for (k in 1:2) {
  for (p in 2:5) {
    cname1 <- paste0("gamma[4, ", k, ", ", p, "]")
    print(paste(cname1))
    col <- MCMCsamples[, cname1]
    print(mean(col) - sd(col) * qnorm(0.75))
    print(mean(col) - sd(col) * qnorm(0.25))
  }
}

for (k in 1:2) {
  for (p in 2:5) {
    cname1 <- paste0("gamma[4, ", k, ", ", p, "]")
    cname2 <- paste0("gamma[5, ", k, ", ", p, "]")
    print(paste(cname1, "-", cname2))
    col <- MCMCsamples[, cname1] - MCMCsamples[, cname2]
    print(mean(col) - sd(col) * qnorm(0.95))
    print(mean(col) - sd(col) * qnorm(0.05))
  }
}

for (k in 1:2) {
  for (p in 2:5) {
    cname1 <- paste0("gamma[4, ", k, ", ", p, "]")
    print(paste(cname1))
    col <- MCMCsamples[, cname1]
    print(mean(col) - sd(col) * qnorm(0.95))
    print(mean(col) - sd(col) * qnorm(0.05))
  }
}