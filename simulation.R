source("simulation-functions.R")

simulate.one <- function(seed, nimble.environment,
                     arms, patients.in.arm,
                     endpoints, covariates,
                     covariate.range,
                     beta, gamma,
                     n.iters = 1000, n.burn.in = 100) {
  
  set.seed(seed)
  
  adm.types <- c("weak", "strong")
  methods <- c("adjusted", "direct", "average")
  
  superiority.thresholds <- rep(1 / 2, length(endpoints))
  noninferiority.thresholds <- rep(-1 / 2, length(endpoints))
  
  simulated.data <- simulate.data(seed = seed,
                                  n.arms = length(arms),
                                  n.patients.per.arm = length(patients.in.arm),
                                  n.endpoints = length(endpoints),
                                  n.covariates = length(covariates),
                                  covariate.range = covariate.range,
                                  beta, gamma)
  
  param.sample <- fit.nimble(nimble.environment = nimble.environment,
                             y = simulated.data$y,
                             X = simulated.data$X,
                             n.iters = n.iters + n.burn.in)
  
  use <- n.burn.in + 1:n.iters
  
  covariate.ranges <- list()
  covariate.ranges[[1]] <- 1
  for (covariate in covariates[-1]) {
    covariate.ranges[[covariate]] <- covariate.range
  }
  covariate.space <- as.matrix(expand.grid(covariate.ranges))
  
  gamma.cols <- find.gamma.col.indices(param.sample,
                                       arms, endpoints, covariates)
  
  
  ##########################
  ### Credible subgroups ###
  ##########################
  
  excl.credsubs <- array(NA, dim=c(length(adm.types),
                                    length(methods),
                                    nrow(covariate.space)),
                          dimnames=list(adm.types,
                                        methods,
                                        NULL))
  
  excl.credsubs[adm.types, "adjusted", ] <- bound.adm.adjusted(
    covariate.space, param.sample, gamma.cols,
    arms = arms, test.arm = length(arms),
    endpoints = endpoints,
    superiority.thresholds = superiority.thresholds,
    noninferiority.thresholds = noninferiority.thresholds,
    cred=0.50
  )
  
  excl.credsubs[adm.types, "direct", ] <- bound.adm.direct(
    covariate.space, param.sample, gamma.cols,
    arms = arms, test.arm = length(arms),
    endpoints = endpoints,
    superiority.thresholds = superiority.thresholds,
    noninferiority.thresholds = noninferiority.thresholds,
    cred=0.50
  )
  
  excl.credsubs[adm.types, "average", ] <- bound.adm.average(
    covariate.space, param.sample, gamma.cols,
    arms = arms, test.arm = length(arms),
    endpoints = endpoints,
    superiority.thresholds = superiority.thresholds,
    noninferiority.thresholds = noninferiority.thresholds,
    cred=0.50, X=simulated.data$X
  )
  
  excl.credsubs
}

simulate <- function(seeds,
                     arms, patients.in.arm,
                     endpoints, covariates,
                     covariate.range,
                     beta, gamma,
                     n.iters = 1000, n.burn.in = 100) {
  
  adm.types <- c("weak", "strong")
  methods <- c("adjusted", "direct", "average")
  
  summaries <- array(NA, dim=c(2, length(adm.types),
                               length(methods),
                               max(seeds)),
                     dimnames=list(
                       c("sensitivity", "specificity"),
                       adm.types,
                       methods,
                       NULL
                     ))
  
  covariate.ranges <- list()
  covariate.ranges[[1]] <- 1
  for (covariate in covariates[-1]) {
    covariate.ranges[[covariate]] <- covariate.range
  }
  covariate.space <- as.matrix(expand.grid(covariate.ranges))
  
  init.data <- simulate.data(1, length(arms), length(patients.in.arm),
                             length(endpoints), length(covariates),
                             covariate.range,
                             beta, gamma)
  
  nimble.env <- compile.nimble(y=init.data$y,
                               X=init.data$X,
                               arms=arms,
                               patients.in.arm=patients.in.arm,
                               endpoints=endpoints,
                               covariates=covariates,
                               covariate.range=covariate.range)
  
  ### True admissibility
  
  adm <- array(NA, dim=c(length(adm.types), nrow(covariate.space)),
               dimnames=list(adm.types, NULL))
  
  superiority.thresholds <- rep(1 / 2, length(endpoints))
  noninferiority.thresholds <- rep(-1 / 2, length(endpoints))
  
  adm["weak", ] <- find.weakly.adm.points(
    covariate.space = covariate.space,
    gamma = gamma,
    arms = arms,
    test.arm = length(arms),
    endpoints = endpoints,
    superiority.thresholds = superiority.thresholds,
    noninferiority.thresholds = noninferiority.thresholds
  )
  
  adm["strong", ] <- find.strongly.adm.points(
    covariate.space = covariate.space,
    gamma = gamma,
    arms = arms,
    test.arm = length(arms),
    endpoints = endpoints,
    superiority.thresholds = superiority.thresholds,
    noninferiority.thresholds = noninferiority.thresholds
  )
  
  
  excl.credsubs <- array(NA, dim=c(length(seeds),
                                   length(adm.types),
                                   length(methods),
                                   nrow(covariate.space)),
                         dimnames=list(NULL,
                                       adm.types,
                                       methods,
                                       NULL))
  for (seed in seeds) {
    print(paste("run", length(arms), length(endpoints), seed))
    excl.credsubs[seed, , , ] <-
      simulate.one(seed, nimble.env,
                   arms, patients.in.arm,
                   endpoints, covariates,
                   covariate.range = covariate.range,
                   beta=beta, gamma=gamma,
                   n.iters = 1000, n.burn.in = 100)
    
    for (method in methods) {
      for (adm.type in adm.types) {
        summaries["sensitivity", adm.type, method, seed] <-
          mean(excl.credsubs[seed, adm.type, method, ][adm[adm.type, ]])
        
        summaries["specificity", adm.type, method, seed] <-
          mean(!excl.credsubs[seed, adm.type, method, ][!adm[adm.type, ]])
      }
    }
  }
  summaries
}

simulate.config <- function(n.arms.vec, n.endpoints.vec, seeds,
                            n.patients.in.arm=100,
                            methods=c("adjusted", "direct", "average"),
                            adm.types=c("weak", "strong")) {
  
  summaries.all <- array(NA, dim=c(2,
                                   length(adm.types),
                                   length(methods),
                                   max(n.arms.vec),
                                   max(n.endpoints.vec),
                                   max(seeds)),
                         dimnames=list(
                           summary=c("sensitivity", "specificity"),
                           adm.type=adm.types,
                           method=methods,
                           n.arms=NULL,
                           n.endpoints=NULL,
                           seed=NULL
                         ))
  
  for (n.arms in n.arms.vec) {
    for (n.endpoints in n.endpoints.vec) {
      arms <- 1:n.arms
      endpoints <- 1:n.endpoints
      patients.in.arm <- 1:n.patients.in.arm
      covariates <- 1:3
      
      first.non.intercept <- 2
      
      beta <- array(1, dim=c(length(endpoints),
                             length(covariates)))
      
      gamma <- array(0, dim=c(length(arms),
                              length(endpoints),
                              length(covariates)))
      gamma[-arms[1], endpoints[1], covariates[first.non.intercept]] <- 1 / 3
      gamma[arms[length(arms)], endpoints[1], covariates[first.non.intercept]] <- 1
      
      summaries <- simulate(seeds=seeds,
                            arms=arms, patients.in.arm=patients.in.arm,
                            endpoints=endpoints, covariates=covariates,
                            covariate.range=-2:2,
                            beta=beta, gamma=gamma,
                            n.iters = 1000, n.burn.in = 100)
      
      summaries.all[, , , n.arms=n.arms, n.endpoints=n.endpoints, ] <- summaries
    }
  }
  
  summaries.all
}

n.arms.vec <- 2:8
n.endpoints.vec <- 1:8
seeds <- 1:1000

summaries.arms <- simulate.config(n.arms.vec=n.arms.vec,
                                 n.endpoints.vec=n.endpoints.vec[1],
                                 seeds=seeds)

summaries.endpoints <- simulate.config(n.arms.vec=n.arms.vec[1],
                                 n.endpoints.vec=n.endpoints.vec,
                                 seeds=seeds)

print(apply(summaries.arms["sensitivity", "weak", , , 1, ], c(1, 2), mean))
print(apply(summaries.arms["specificity", "weak", , , 1, ], c(1, 2), mean))

print(apply(summaries.arms["sensitivity", "strong", , , 1, ], c(1, 2), mean))
print(apply(summaries.arms["specificity", "strong", , , 1, ], c(1, 2), mean))

print(apply(summaries.endpoints["sensitivity", "weak", , 2, , ], c(1, 2), mean))
print(apply(summaries.endpoints["specificity", "weak", , 2, , ], c(1, 2), mean))

print(apply(summaries.endpoints["sensitivity", "strong", , 2, , ], c(1, 2), mean))
print(apply(summaries.endpoints["specificity", "strong", , 2, , ], c(1, 2), mean))

diagnostic.names <- c("sensitivity", "specificity")
adm.types <- c("weak", "strong")
method.names <- c("adjusted", "direct", "average")

summaries.all <- array(NA, dim=c(2, 2, 3,
                                 max(n.arms.vec), max(n.endpoints.vec),
                                 max(seeds)))

summaries.all[, , , n.arms.vec, n.endpoints.vec[1], seeds] <-
  summaries.arms[, , , n.arms.vec, n.endpoints.vec[1], seeds]

summaries.all[, , , n.arms.vec[1], n.endpoints.vec, seeds] <-
  summaries.endpoints[, , , n.arms.vec[1], n.endpoints.vec, seeds]

means <- apply(summaries.all, 1:5, mean)
dimnames(means) <- list(
  diagnostic.names,
  adm.types,
  method.names,
  NULL,
  NULL
)

# save(summaries.all, file="sim-summaries.RData")

plots <- FALSE

if (plots) {
  setEPS()
  postscript("sim-sens-arms.eps", width=5, height=5)
  plot(0, type='n',
       xlim=c(1, 8), ylim=c(0, 1),
       xlab="Arms", ylab="Avg Sensitivity",
       main="1 Endpoint, Varying Arms")
  for (adm.type in adm.types) {
    for (method in method.names) {
      points(1:max(n.arms.vec), means["sensitivity", adm.type, method, , 1],
             pch=which(adm.types == adm.type))
      lines(1:max(n.arms.vec), means["sensitivity", adm.type, method, , 1],
            lty=which(method.names == method))
    }
  }
  dev.off()
  
  setEPS()
  postscript("sim-spec-arms.eps", width=5, height=5)
  plot(0, type='n',
       xlim=c(1, 8), ylim=c(0, 1),
       xlab="Arms", ylab="Avg Specificity",
       main="1 Endpoint, Varying Arms")
  for (adm.type in adm.types) {
    for (method in method.names) {
      points(1:max(n.arms.vec), means["specificity", adm.type, method, , 1],
             pch=which(adm.types == adm.type))
      lines(1:max(n.arms.vec), means["specificity", adm.type, method, , 1],
            lty=which(method.names == method))
    }
  }
  legend("left", c("Weak adm", "Strong adm"), pch=1:2)
  legend("right", c("Fully adjusted", "Direct", "Naive"), lty=1:3)
  dev.off()
  
  setEPS()
  postscript("sim-sens-ends.eps", width=5, height=5)
  plot(0, type='n',
       xlim=c(1, 8), ylim=c(0, 1),
       xlab="Endpoints", ylab="Avg Sensitivity",
       main="2 Arms, Varying Endpoints")
  for (adm.type in adm.types) {
    for (method in method.names) {
      points(1:max(n.endpoints.vec), means["sensitivity", adm.type, method, 2, ],
             pch=which(adm.types == adm.type))
      lines(1:max(n.endpoints.vec), means["sensitivity", adm.type, method, 2, ],
            lty=which(method.names == method))
    }
  }
  dev.off()
  
  setEPS()
  postscript("sim-spec-ends.eps", width=5, height=5)
  plot(0, type='n',
       xlim=c(1, 8), ylim=c(0, 1),
       xlab="Endpoints", ylab="Avg Specificity",
       main="2 Arms, Varying Endpoints")
  for (adm.type in adm.types) {
    for (method in method.names) {
      points(1:max(n.endpoints.vec), means["specificity", adm.type, method, 2, ],
             pch=which(adm.types == adm.type))
      lines(1:max(n.endpoints.vec), means["specificity", adm.type, method, 2, ],
            lty=which(method.names == method))
    }
  }
  dev.off()
}
