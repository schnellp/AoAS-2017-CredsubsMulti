library(nimble)

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
    }
    
    # efficacy
    y[i, 1] ~ dnorm(mean=eta[i, 1], var=sigma.sq)
    # safety
    y[i, 2] ~ dbern(p=1 / (1 + exp(-eta[i, 2])))
  }
  
  # prior
  tau ~ dgamma(shape=0.001, rate=0.001)
  sigma.sq <- 1 / tau
  d[1] <- 0             # placebo
  d[2] ~ dunif(0, d[3]) # low dose test
  d[3] ~ dunif(0, 1)    # medium dose test
  d[4] <- 1             # high dose test
  d[5] <- 1             # active control
  
  # regression coefficients
  for (j in 1:2) {
    beta[j, 1] ~ dnorm(mean=0, tau=1/10^4)
    gamma[1, j, 1] <- 0
    gamma[2, j, 1] <- gamma[4, j, 1]
    gamma[3, j, 1] <- gamma[4, j, 1]
    
    gamma[4, j, 1] ~ dnorm(mean=0, tau=1/10^4)
    gamma[5, j, 1] ~ dnorm(mean=0, tau=1/10^4)
    
    for (p in 2:P) {
      beta[j, p] ~ dnorm(mean=0, tau=1/10^4)
      gamma[1, j, p] <- 0
      gamma[2, j, p] <- gamma[4, j, p]
      gamma[3, j, p] <- gamma[4, j, p]
      
      # Flat priors
      # gamma[4, j, p] ~ dnorm(mean=0, tau=1/10^4)
      # gamma[5, j, p] ~ dnorm(mean=0, tau=1/10^4)
      
      # Shrinkage priors
      gamma[4, j, p] ~ dnorm(mean=0, tau=1)
      gamma[5, j, p] ~ dnorm(mean=0, tau=1)
      
      # Spike-and-slab priors
      # spike[4, j, p] ~ dbern(p=0.1)
      # norm[4, j, p] ~ dnorm(0, 1/10^4)
      # gamma[4, j, p] <- (1 - (p > 1) * spike[4, j, p]) * norm[4, j, p]
      # 
      # spike[5, j, p] ~ dbern(p=0.1)
      # norm[5, j, p] ~ dnorm(0, 1/10^4)
      # gamma[5, j, p] <- (1 - (p > 1) * spike[5, j, p]) * norm[5, j, p]
    }
  }
})

set.seed(1)

data.simple <- read.csv("data-simple.csv")
#load("data-simple.RData")

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
                 spike=array(0, c(5, 2, 5)),
                 norm=array(0, c(5, 2, 5)),
                 tau = 1)

alzModel <- nimbleModel(code = alzCode, name = 'alz', constants = alzConsts,
                   data = alzData,
                   inits = alzInits)

alzSpec <- configureMCMC(alzModel)
alzSpec$addMonitors(c('gamma', 'd'))

alzMCMC <- buildMCMC(alzSpec)

CalzModel <- compileNimble(alzModel)
CalzMCMC <- compileNimble(alzMCMC, project = alzModel)

CalzMCMC$run(110000)

MCMCsamples <- as.matrix(CalzMCMC$mvSamples)

save(MCMCsamples, file="MCMCsamples.RData")
# load("MCMCsamples.RData")

spike.and.slab <- FALSE
if (spike.and.slab) {
  result <- MCMCsamples[, c(166, 11, 11, 11, 
                            c(seq(1, 9, by=2), seq(2, 10, by=2)),
                            (0:4) * 5 * 2 + 16,
                            (0:4) * 5 * 2 + 16 + 5,
                            (0:4) * 5 * 2 + 17,
                            (0:4) * 5 * 2 + 17 + 5,
                            (0:4) * 5 * 2 + 18,
                            (0:4) * 5 * 2 + 18 + 5,
                            (0:4) * 5 * 2 + 19,
                            (0:4) * 5 * 2 + 19 + 5,
                            (0:4) * 5 * 2 + 20,
                            (0:4) * 5 * 2 + 20 + 5,
                            11:15)]
} else {
  result <- MCMCsamples[, c(66, 11, 11, 11,
                            seq(1, 9, by=2), seq(2, 10, by=2),
                            (0:4) * 5 * 2 + 16,
                            (0:4) * 5 * 2 + 16 + 5,
                            (0:4) * 5 * 2 + 17,
                            (0:4) * 5 * 2 + 17 + 5,
                            (0:4) * 5 * 2 + 18,
                            (0:4) * 5 * 2 + 18 + 5,
                            (0:4) * 5 * 2 + 19,
                            (0:4) * 5 * 2 + 19 + 5,
                            (0:4) * 5 * 2 + 20,
                            (0:4) * 5 * 2 + 20 + 5,
                            11:15)]
}

colnames(result) <- c("sigma.sq", "tau1sq", "tau2sq", "rho",
                      paste0("beta", rep(1:2, each=5), 1:5),
                      paste0("gamma", rep(1:5, each=2 * 5), rep(1:2, each=5), 1:5),
                      paste0("d", rep(1:5, each=1)))

result[, 1] <- 1 / result[, 1]
rm(MCMCsamples)


##########################
### Credible Subgroups ###
##########################

keep <- seq(10100, 110000, by=10)

gamma.cols <- array(15:64, dim=c(5, 2, 5))
gamma.cols <- aperm(gamma.cols, 3:1)


grid <- expand.grid(SEVERITY=5:45,
                    AGE=55:90,
                    SEX=0:1,
                    CARRIER=0:1)


grid <- scale(grid,
              center=attr(covariates, "scaled:center"),
              scale=attr(covariates, "scaled:scale"))

X <- as.matrix(expand.grid(
  INTERCEPT=1,
  SEVERITY=5:45,
  AGE=55:90,
  SEX=0:1,
  CARRIER=0:1
))

X <- cbind(1, grid)

grid <- expand.grid(SEVERITY=5:45,
                    AGE=55:90,
                    SEX=factor(c("F", "M")),
                    CARRIER=factor(c("NON-CARRIER", "CARRIER")))

source("sim-cred-bands.R")

library(MCMCpack)
library(vcd)

cross.hatch <- function(m) {
  m2 <- m
  
  for (i in 1:nrow(m)) {
    for (j in 1:ncol(m)) {
      if (m[i, j] < 0 && xor((i + 1) %% 4 < 2, (j + 1) %% 4 < 2)) {
        m2[i, j] <- 0
      }
    }
  }
  
  m2
}

comp.plot.credsubs <- function(coeffs, X, sup.threshold=0, noninf.threshold=0,
                               in.color=TRUE) {
  grid.sample <- coeffs %*% t(X)
  
  band <- quantile.band(grid.sample, 0.5, verbose=TRUE, low.mem=TRUE)
  
  par(mar=c(4.1, 4.1, 4.1, 2.1))
  
  titles.sex <- c("M"="Male", "F"="Female")
  titles.car <- c("NON-CARRIER"="Non-Carriers", "CARRIER"="Carriers")
  
  in.d <- band$lower > sup.threshold
  in.s <- band$upper > sup.threshold
  
  in.d.prime <- band$lower > noninf.threshold
  in.s.prime <- band$upper >= noninf.threshold
  
  layout(matrix(c(1, 2,
                  3, 4), 
                2, 2, byrow=TRUE),
         width=c(1, 1))
  
  for (i in c("F", "M")) {
    for (j in c("NON-CARRIER", "CARRIER")) {
      w <- which(
        grid[, 3] == i &
          grid[, 4] == j)
      imm <- cross.hatch(
        matrix(-2 + in.d[w] + in.s[w] + in.d.prime[w] + in.s.prime[w],
                    nrow=36, ncol=41, byrow=TRUE))
      image(t(imm),
            col=c("gray", "#010101", "white", "gray", "#010101"),
            zlim=c(-2, 2),
            xlab="Severity", ylab="Age",
            xaxt='n', yaxt='n',
            useRaster=FALSE,
            cex.lab=1.5)
      
      title(paste(titles.sex[i], titles.car[j]), cex.main=1.8)
      axis(1, at=(seq(5, 45, by=5) - 5) / (40), labels=seq(5, 45, by=5))
      axis(2, at=(seq(55, 90, by=5) - 55) / (35), labels=seq(55, 90, by=5))
      box()
    }
  }
  
  layout(1)
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}

plot.adm.credsubs <- function(band, in.color=TRUE) {
  par(mar=c(4.1, 4.1, 4.1, 2.1))
  
  titles.sex <- c("M"="Male", "F"="Female")
  titles.car <- c("NON-CARRIER"="Non-Carriers", "CARRIER"="Carriers")
  
  in.d <- band$lower > 0
  in.s <- band$upper > 0
  
  
  layout(matrix(c(1, 2,
                  3, 4), 
                2, 2, byrow=TRUE),
         width=c(1, 1))
  
  for (i in c("F", "M")) {
    for (j in c("NON-CARRIER", "CARRIER")) {
      w <- which(
        grid[, 3] == i &
          grid[, 4] == j)
      imm <- cross.hatch(
        matrix(-1 + in.d[w] + in.s[w],
               nrow=36, ncol=41, byrow=TRUE))
      image(t(imm),
            col=c("#010101", "white", "#010101"),
            zlim=c(-1, 1),
            xlab="Severity", ylab="Age",
            xaxt='n', yaxt='n',
            useRaster=FALSE,
            cex.lab=1.5)
      
      title(paste(titles.sex[i], titles.car[j]), cex.main=1.8)
      axis(1, at=(seq(5, 45, by=5) - 5) / (40), labels=seq(5, 45, by=5))
      axis(2, at=(seq(55, 90, by=5) - 55) / (35), labels=seq(55, 90, by=5))
      box()
    }
  }
  
  layout(1)
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}

sup.eff <- 0
noninf.eff <- -0.5 # sds
sup.saf <- 0
noninf.saf <- -0.18 # 1.20 OR

# efficacy vs placebo
setEPS()
postscript("eff-tst-pcb-gs.eps", width=7, height=7)
comp.plot.credsubs( result[keep, gamma.cols[4, 1, ]], X,
                    sup.threshold=sup.eff, noninf.threshold=noninf.eff)
dev.off()
# safety vs placebo
setEPS()
postscript("saf-tst-pcb-gs.eps", width=7, height=7)
comp.plot.credsubs(-(result[keep, gamma.cols[4, 2, ]]), X,
                   sup.threshold=sup.saf, noninf.threshold=noninf.saf)
dev.off()
# efficacy vs standard
setEPS()
postscript("eff-tst-act-gs.eps", width=7, height=7)
comp.plot.credsubs( (result[keep, gamma.cols[4, 1, ]] -
                       result[keep, gamma.cols[5, 1, ]]), X,
                    sup.threshold=sup.eff, noninf.threshold=noninf.eff)
dev.off()
# safety vs standard
setEPS()
postscript("saf-tst-act-gs.eps", width=7, height=7)
comp.plot.credsubs(-(result[keep, gamma.cols[4, 2, ]] -
                       result[keep, gamma.cols[5, 2, ]]), X,
                   sup.threshold=sup.saf, noninf.threshold=noninf.saf)
dev.off()


########################################
### ADMISSIBILITY CREDIBLE SUBGROUPS ### (Weak)
########################################

gamma.1v0b1 <- result[keep, gamma.cols[4, 1, ]]
grid.sample.1v0b1 <- gamma.1v0b1 %*% t(X)
rm(gamma.1v0b1)

gamma.1v0b2 <- -result[keep, gamma.cols[4, 2, ]]
grid.sample.1v0b2 <- gamma.1v0b2 %*% t(X)
rm(gamma.1v0b2)

grid.sample.1v0bA <- (pmax(grid.sample.1v0b1, grid.sample.1v0b2) > 0) |
  (pmin(grid.sample.1v0b1 - noninf.eff, grid.sample.1v0b2 - noninf.saf) >= 0)
rm(grid.sample.1v0b1, grid.sample.1v0b2)

gamma.1v2b1 <- result[keep, gamma.cols[4, 1, ]] - result[keep, gamma.cols[5, 1, ]]
grid.sample.1v2b1 <- gamma.1v2b1 %*% t(X)
rm(gamma.1v2b1)

gamma.1v2b2 <- -(result[keep, gamma.cols[4, 2, ]] - result[keep, gamma.cols[5, 2, ]])
grid.sample.1v2b2 <- gamma.1v2b2 %*% t(X)
rm(gamma.1v2b2)

grid.sample.1v2bA <- (pmax(grid.sample.1v2b1, grid.sample.1v2b2) > 0) |
  (pmin(grid.sample.1v2b1 - noninf.eff, grid.sample.1v2b2 - noninf.saf) >= 0)
rm(grid.sample.1v2b1, grid.sample.1v2b2)

grid.sample <- grid.sample.1v0bA * grid.sample.1v2bA
rm(grid.sample.1v0bA, grid.sample.1v2bA)

est <- colMeans(grid.sample)

band <- quantile.band(grid.sample, 0.5, verbose=TRUE, low.mem=TRUE)

### PLOTTING ###

setEPS()
postscript("adm-weak-gs.eps", width=7, height=7)
plot.adm.credsubs(band)
dev.off()


########################################
### ADMISSIBILITY CREDIBLE SUBGROUPS ### (Strong)
########################################

gamma.1v0b1 <- result[keep, gamma.cols[4, 1, ]]
grid.sample.1v0b1 <- gamma.1v0b1 %*% t(X)
rm(gamma.1v0b1)

gamma.1v0b2 <- -result[keep, gamma.cols[4, 2, ]]
grid.sample.1v0b2 <- gamma.1v0b2 %*% t(X)
rm(gamma.1v0b2)

grid.sample.1v0bA <- (pmax(grid.sample.1v0b1, grid.sample.1v0b2) > 0) *
  (pmin(grid.sample.1v0b1 - noninf.eff, grid.sample.1v0b2 - noninf.saf) >= 0)
rm(grid.sample.1v0b1, grid.sample.1v0b2)

gamma.1v2b1 <- result[keep, gamma.cols[4, 1, ]] - result[keep, gamma.cols[5, 1, ]]
grid.sample.1v2b1 <- gamma.1v2b1 %*% t(X)
rm(gamma.1v2b1)

gamma.1v2b2 <- -(result[keep, gamma.cols[4, 2, ]] - result[keep, gamma.cols[5, 2, ]])
grid.sample.1v2b2 <- gamma.1v2b2 %*% t(X)
rm(gamma.1v2b2)

grid.sample.1v2bA <- (pmax(grid.sample.1v2b1, grid.sample.1v2b2) > 0) *
  (pmin(grid.sample.1v2b1 - noninf.eff, grid.sample.1v2b2 - noninf.saf) >= 0)
rm(grid.sample.1v2b1, grid.sample.1v2b2)

grid.sample <- grid.sample.1v0bA * grid.sample.1v2bA
rm(grid.sample.1v0bA, grid.sample.1v2bA)

est <- colMeans(grid.sample)

band <- quantile.band(grid.sample, 0.5, verbose=TRUE, low.mem=TRUE)


### PLOTTING ###
setEPS()
postscript("adm-strong-gs.eps", width=7, height=7)

plot.adm.credsubs(band)

dev.off()
