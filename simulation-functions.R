library(nimble)

source("sim-cred-bands.R")

rowAlls <- function(X) {
  apply(as.matrix(X), 1, all)
}

rowAnys <- function(X) {
  apply(as.matrix(X), 1, any)
}

colAlls <- function(X) {
  apply(as.matrix(X), 2, all)
}

colAnys <- function(X) {
  apply(as.matrix(X), 2, any)
}

simulate.covariates <- function(n.arms, n.patients.per.arm,
                                n.endpoints, n.covariates,
                                covariate.range=-2:2) {
  
  covariate.range <- range(covariate.range)[2] - range(covariate.range)[1] + 1
  
  X <- array(NA, dim=c(n.arms, n.patients.per.arm, n.covariates))
  
  intercept.dim <- 1
  X[, , -intercept.dim] <-
    rbinom(n = n.arms * n.patients.per.arm * (n.covariates - 1),
           size = covariate.range - 1,
           prob = 1 / 2) -
    (covariate.range - 1) / 2
  
  X[, , intercept.dim] <- 1
  
  X
}

simulate.outcomes <- function(X, beta, gamma) {
  
  n.arms <- dim(X)[1]
  n.patients.per.arm <- dim(X)[2]
  n.covariates <- dim(X)[3]
  
  n.endpoints <- dim(beta)[1]
  
  y <- array(NA, dim=c(n.arms, n.patients.per.arm, n.endpoints))
  
  if (is.null(n.endpoints)) {
    print("NULL ENDPOINTS")
    traceback()
  }
  
  for (arm in 1:n.arms) {
    for (endpoint in 1:n.endpoints) {
      y[arm, 1:n.patients.per.arm, endpoint] <-
        rnorm(n    = n.patients.per.arm,
              mean = X[arm, 1:n.patients.per.arm, 1:n.covariates] %*%
                beta[endpoint, 1:n.covariates] +
                X[arm, 1:n.patients.per.arm, 1:n.covariates] %*%
                gamma[arm, endpoint, 1:n.covariates],
              sd   = 1)
    }
  }
  
  y
}

simulate.data <- function(seed,
                         n.arms, n.patients.per.arm,
                         n.endpoints, n.covariates,
                         covariate.range,
                         beta, gamma) {
  set.seed(seed)
  
  X <- simulate.covariates(n.arms, n.patients.per.arm,
                           n.endpoints, n.covariates, covariate.range)
  y <- simulate.outcomes(X, beta, gamma)
  
  list(X=X, y=y)
}

compile.nimble <- function(y, X,
                           arms, patients.in.arm,
                           endpoints,
                           covariates, covariate.range) {
  
  simCode <- nimbleCode({ 
    # likelihood
    for (arm in 1:n.arms) {
      for (patient in 1:n.patients.in.arm) {
        for (endpoint in 1:n.endpoints) {
          
          lin.pred[arm, patient, endpoint] <-
            inprod(X[arm, patient, 1:n.covariates],
                   beta[endpoint, 1:n.covariates]) +
            inprod(X[arm, patient, 1:n.covariates],
                   gamma[arm, endpoint, 1:n.covariates])
          
          y[arm, patient, endpoint] ~
            dnorm(mean = lin.pred[arm, patient, endpoint],
                  tau = tau[endpoint])
        }
      }
    }
    
    # prior
    for (endpoint in 1:n.endpoints) {
      tau[endpoint] ~ dgamma(shape=0.001, rate=0.001)
      for (covariate in 1:n.covariates) {
        beta[endpoint, covariate] ~ dnorm(mean=0, sd=100)
      }
    }
    
    for (endpoint in 1:n.endpoints) {
      for (covariate in 1:n.covariates) {
        gamma[1, endpoint, covariate] <- 0
      }
      for (arm in 2:n.arms) {
        gamma[arm, endpoint, 1] ~ dnorm(mean=0, sd=100)
        for (covariate in 2:n.covariates) {
          gamma[arm, endpoint, covariate] ~ dnorm(mean=0, sd=100)
        }
      }
    }
  })
  
  simConsts <- list(n.arms = length(arms),
                    n.patients.in.arm = length(patients.in.arm),
                    n.endpoints = length(endpoints),
                    n.covariates = length(covariates))
  
  beta.init <- array(0, dim=c(length(endpoints),
                              length(covariates)))
  
  gamma.init <- array(0, dim=c(length(arms),
                               length(endpoints),
                               length(covariates)))
  
  simInits <- list(beta = beta.init,
                   gamma = gamma.init,
                   tau = rep(1, length(endpoints)))
  
  init.data <- simulate.data(seed = 1,
                             n.arms = length(arms),
                             n.patients.per.arm = length(patients.in.arm),
                             n.covariates = length(covariates),
                             covariate.range = covariate.range,
                             beta = beta.init, gamma = gamma.init)

  simData <- list(X = init.data$X,
                  y = init.data$y)
  
  simModel <- nimbleModel(code = simCode, name="sim",
                          constants = simConsts,
                          data = simData,
                          inits = simInits)
  
  simSpec <- configureMCMC(simModel)
  
  simMCMC <- buildMCMC(simSpec)
  
  CsimModel <- compileNimble(simModel)
  CsimMCMC <- compileNimble(simMCMC, project = simModel)
  
  environment()
}

fit.nimble <- function(nimble.environment, y, X, n.iters) {
  
  simData <- list(X = X,
                  y = y)
  nimble.environment$CsimModel$setData(simData)
  
  nimble.environment$CsimMCMC$run(n.iters)
  
  param.sample <- as.matrix(nimble.environment$CsimMCMC$mvSamples)
  
  param.sample
}

find.gamma.col.indices <- function(param.sample,
                                   arms, endpoints, covariates) {
  gamma.cols <- array(NA,
                      dim=c(length(arms),
                            length(endpoints),
                            length(covariates)))
  
  for (arm in arms) {
    for (endpoint in endpoints) {
      for (covariate in covariates) {
        gamma.cols[arm, endpoint, covariate] <-
          which(colnames(param.sample) ==
                  paste0("gamma[", arm, ", ", endpoint, ", ", covariate, "]"))
      }
    }
  }
  
  gamma.cols
}

find.points.w.any.superiority <- function(covariate.space, gamma,
                                          test.arm, control.arm,
                                          endpoints, thresholds) {
  any.superior <- rep(FALSE, nrow(covariate.space))
  
  for (endpoint in endpoints) {
    any.superior <- any.superior | covariate.space %*%
      (gamma[test.arm, endpoint, ] - gamma[control.arm, endpoint, ]) >
      thresholds[endpoint]
  }
  
  any.superior
}

find.points.w.all.noninferiority <- function(covariate.space, gamma,
                                            test.arm, control.arm,
                                            endpoints, thresholds) {
  all.noninferior <- rep(TRUE, nrow(covariate.space))
  
  for (endpoint in endpoints) {
    all.noninferior <- all.noninferior & covariate.space %*%
      (gamma[test.arm, endpoint, ] - gamma[control.arm, endpoint, ]) >=
      thresholds[endpoint]
  }
  
  all.noninferior
}

find.weakly.adm.points <- function(covariate.space, gamma,
                                   arms, test.arm,
                                   endpoints,
                                   superiority.thresholds,
                                   noninferiority.thresholds) {
  
  weakly.adm.vs <- matrix(NA, ncol=nrow(covariate.space),
                          nrow=length(arms[-test.arm]))
  
  for (arm in arms[-test.arm]) {
    any.superior <-
      find.points.w.any.superiority(covariate.space, gamma,
                                    test.arm = test.arm,
                                    control.arm = arm,
                                    endpoints,
                                    superiority.thresholds)
    
    all.noninferior <-
      find.points.w.all.noninferiority(covariate.space, gamma,
                                       test.arm = test.arm,
                                       control.arm = arm,
                                       endpoints,
                                       noninferiority.thresholds)
    
    weakly.adm.vs[arm, ] <- any.superior | all.noninferior
  }
  
  colAlls(weakly.adm.vs)
}

find.strongly.adm.points <- function(covariate.space, gamma,
                                     arms, test.arm,
                                     endpoints,
                                     superiority.thresholds,
                                     noninferiority.thresholds) {
  
  strongly.adm.vs <- matrix(NA, ncol=nrow(covariate.space),
                            nrow=length(arms[-test.arm]))
  
  for (arm in arms[-test.arm]) {
    any.superior <-
      find.points.w.any.superiority(covariate.space, gamma,
                                    test.arm = test.arm,
                                    control.arm = arm,
                                    endpoints,
                                    superiority.thresholds)
    
    all.noninferior <-
      find.points.w.all.noninferiority(covariate.space, gamma,
                                       test.arm = test.arm,
                                       control.arm = arm,
                                       endpoints,
                                       noninferiority.thresholds)
    
    strongly.adm.vs[arm, ] <- any.superior & all.noninferior
  }
  
  colAlls(strongly.adm.vs)
}

compute.param.diff.sample <- function(param.sample, gamma.cols,
                                      arms, test.arm) {
  param.diff.sample <- array(NA, dim=dim(param.sample),
                             dimnames=dimnames(param.sample))
  for (arm in arms) {
    param.diff.sample[, gamma.cols[arm, , ]] <-
      param.sample[, gamma.cols[test.arm, , ]] -
      param.sample[, gamma.cols[arm, , ]]
  }
  
  param.diff.sample
}

compute.pte.draws <- function(X, param.sample, gamma.cols,
                              arms, test.arm, endpoints) {
  
  param.diff.sample <- compute.param.diff.sample(
    param.sample = param.sample,
    gamma.cols = gamma.cols,
    arms = arms,
    test.arm = test.arm
  )
  
  pte.draws <- array(NA, dim=c(nrow(param.diff.sample),
                               nrow(X),
                               length(arms[-test.arm]),
                               length(endpoints)))
  
  for (arm in arms[-test.arm]) {
    for (endpoint in endpoints) {
      #     ________       ________      
      #  d |        |   d |        |        ___________________
      #  r | local  |   r |        |       |                   |
      #  a | mean   | = a | params |  X  p | transposed design |
      #  w | draws  |   w |        |       |___________________|
      #  s |________|   s |________|             location
      #     location          p
      pte.draws[, , arm, endpoint] <-
        param.diff.sample[, gamma.cols[arm, endpoint, ]] %*% t(X)
    }
  }
  
  pte.draws
}

bound.adm.adjusted <- function(covariate.space, param.sample, gamma.cols,
                                  arms, test.arm,
                                  endpoints,
                                  superiority.thresholds,
                                  noninferiority.thresholds,
                                  cred=0.95) {
  
  multi.sim.cred.band <- multi.quantile.band(
    sample = compute.pte.draws(covariate.space, param.sample, gamma.cols,
                               arms, test.arm, endpoints),
    alpha= 1 - cred
  )
  
  excl.credsub <- array(NA, dim=c(2, nrow(covariate.space)),
                        dimnames=list(c("weak", "strong"), NULL))
  
  working <- list()
  working$superior.wrt <- working$noninferior.wrt <-
    array(NA, dim=c(length(arms[-test.arm]),
                    length(endpoints),
                    nrow(covariate.space)))
  
  for (arm in arms[-test.arm]) {
    for (endpoint in endpoints) {
      working$superior.wrt[arm, endpoint, ] <-
        multi.sim.cred.band$lower[, arm, endpoint] >
        superiority.thresholds[endpoint]
      working$noninferior.wrt[arm, endpoint, ] <-
        multi.sim.cred.band$lower[, arm, endpoint] >=
        noninferiority.thresholds[endpoint]
    }
  }
  
  working$any.superior <- apply(working$superior.wrt, c(1, 3), any)
  working$all.noninferior <- apply(working$noninferior.wrt, c(1, 3), all)

  working$weakly.adm.vs <- working$any.superior | working$all.noninferior
  working$strongly.adm.vs <- working$any.superior & working$all.noninferior
  
  excl.credsub["weak", ] <- colAlls(working$weakly.adm.vs)
  excl.credsub["strong", ] <- colAlls(working$strongly.adm.vs)
  
  excl.credsub
}

bound.adm.direct <- function(covariate.space, param.sample, gamma.cols,
                             arms, test.arm,
                             endpoints,
                             superiority.thresholds,
                             noninferiority.thresholds,
                             cred=0.95) {
  
  adm.draws.by.arm <- array(NA,
                            dim=c(2,
                                  length(arms[-test.arm]),
                                  nrow(param.sample),
                                  nrow(covariate.space)),
                            dimnames=list(adm.type=c("weak", "strong"),
                                          arm=NULL,
                                          draw=NULL,
                                          location=NULL))
  
  for (arm in arms[-test.arm]) {
    pte.draws.by.endpoint <- array(NA,
                                   dim=c(length(endpoints),
                                         nrow(param.sample),
                                         nrow(covariate.space)),
                                   dimnames=list(endpoint=NULL,
                                                 draw=NULL,
                                                 location=NULL))
    
    for (endpoint in endpoints) {
      hm <- compute.pte.draws(
        X = covariate.space,
        param.sample = param.sample,
        gamma.cols = gamma.cols,
        arms = arms,
        test.arm = test.arm,
        endpoints = endpoints
      )
      pte.draws.by.endpoint[endpoint=endpoint, , ] <- hm[, , arm, endpoint]
    }
    
    any.superior <- apply(
      X = pte.draws.by.endpoint - superiority.thresholds,
      MARGIN = which(names(dimnames(pte.draws.by.endpoint)) %in%
                       c("draw", "location")),
      FUN = max
    ) > 0
    all.noninferior <- apply(
      X = pte.draws.by.endpoint - noninferiority.thresholds,
      MARGIN = which(names(dimnames(pte.draws.by.endpoint)) %in%
                       c("draw", "location")),
      FUN = min
    ) > 0
    
    adm.draws.by.arm["weak", arm, , ] <- any.superior | all.noninferior
    adm.draws.by.arm["strong", arm, , ] <- any.superior & all.noninferior
  }
  
  adm.draws <- apply(
    X = adm.draws.by.arm,
    MARGIN = which(names(dimnames(adm.draws.by.arm)) != "arm"),
    FUN=all
  )
  
  excl.credsub <- array(NA, dim=c(2, nrow(covariate.space)),
                        dimnames=list(c("weak", "strong"), NULL))
  
  for (adm.type in c("weak", "strong")) {
    excl.credsub[adm.type, ] <- quantile.band(
      sample = adm.draws[adm.type, , ],
      alpha = 1 - cred
    )$lower > 0
  }
  
  excl.credsub
}

bound.adm.average <- function(covariate.space, param.sample, gamma.cols,
                              arms, test.arm,
                              endpoints,
                              superiority.thresholds,
                              noninferiority.thresholds,
                              cred=0.95, X=covariate.space) {
  
  X.mat <- X
  dim(X.mat) <- c(dim(X)[1] * dim(X)[2], dim(X)[3])
  
  pte.draws <- compute.pte.draws(X.mat, param.sample, gamma.cols,
                                 arms, test.arm, endpoints)
  
  ate.draws <- apply(
    X = pte.draws,
    MARGIN = c(1, 3, 4),
    FUN = mean
  )
  
  any.superior <- apply(
    X = aperm(aperm(ate.draws) - superiority.thresholds),
    MARGIN = c(1, 2),
    FUN = max
  ) > 0
  all.noninferior <- apply(
    X = aperm(aperm(ate.draws) - noninferiority.thresholds),
    MARGIN = c(1, 2),
    FUN = min
  ) > 0
  
  adm.draws.by.arm <- array(NA, dim=c(2, length(arms[-test.arm]),
                                      nrow(param.sample)),
                            dimnames=list(c("weak", "strong"), NULL, NULL))
  
  adm.draws.by.arm["weak", , ] <- t(any.superior | all.noninferior)
  adm.draws.by.arm["strong", , ] <- t(any.superior & all.noninferior)
  
  adm.draws.single <- apply(
    X = adm.draws.by.arm,
    MARGIN = c(1, 3),
    all
  )
  
  adm.draws <- array(NA, c(2, nrow(param.sample), nrow(covariate.space)))
  dimnames(adm.draws)[[1]] <- c("weak", "strong")
  for (covariate.point in covariate.space) {
    adm.draws[, , covariate.point] <- adm.draws.single
  }
  
  excl.credsub <- array(NA, dim=c(2, nrow(covariate.space)),
                        dimnames=list(c("weak", "strong"), NULL))
  
  for (adm.type in c("weak", "strong")) {
    csl <- quantile.band(
      sample = adm.draws[adm.type, , ],
      alpha = 1 - cred
    )$lower > 0
    excl.credsub[adm.type, ] <- csl
  }
  
  excl.credsub
}
