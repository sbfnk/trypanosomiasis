library('truncnorm')

##' Produe a single simuation run
##'
##' @param theta parameter vector
##' @return data frame with trajctor
##' @author seb
##' @export
sim_trajectory <- function(theta, init, village, stoptime, ...)
{
  if (missing(stoptime))
  {
    stoptime <-
      village_screening[village_screening$village.number == village]$stoptime
  }

    chronic_options <- list(params = theta, init = init, times = seq(0, stoptime))
    for (stage in 1:2)
    {
        if (!paste0("p", stage) %in% names(theta))
        {
            chronic_options[[paste0("stage", stage, "_passive")]] <-
                passive[[stage]]
        }
    }

    return(do.call(tryp, c(chronic_options, list(...))))
}

##' Draw a parameter sample from the prior density
##'
##' Default parameters reflect Rebecca's study
##' @param villages Number of villages with different parameters
##' @param passive whether to simulate passive detection
##' @param background whether to simuate backgroud infection rates
##' @param transmitted  whether to simulate transmitted infection rates
##' @param chronic whether chronic carriers transmit
##' @return parameter vector
##' @author Sebastian Funk
rprior <- function(villages = 1, passive = TRUE, background = TRUE,
                   transmitted = FALSE, chronic = TRUE)
{
    param_vector <-
        c(pc = runif(1, 0, 0.5), alpha = runif(1, 0, 0.5))
    param_vector["delta"] <- ifelse(chronic, runif(1, 0, 1), 0)

    if (length(villages) == 1)
    {
        param_vector <- c(param_vector,
                          p1 = ifelse(passive, runif(1, 0, 30), 0),
                          p2 = ifelse(passive, runif(1, 0, 30), 0))
        if (background)
        {
            param_vector["lambda"] <- runif(1, 0, 1e-2)
        }
        if (transmitted)
        {
            param_vector["beta"] <- runif(1, 0, 1e-2)
        }
    } else
    {
        village_vector <- c(ifelse(passive,
                                   runif(length(villages), 0, 30),
                                   rep(0, length(villages))),
                            ifelse(passive,
                                   runif(length(villages), 0, 30),
                                   rep(0, length(villages))))
        name_vector <- c("p1", "p2")
        if (background)
        {
            village_vector <- c(village_vector, 0,  1e-2)
            name_vector <- c(name_vector, "lambda")
        }
        if (transmitted)
        {
            village_vector <- c(village_vector, 0,  1e-2)
            name_vector <- c(name_vector, "beta")
        }
        names(village_vector) <-
            paste(rep(name_vector, each = length(villages)),
                  villages, sep = ".")
        param_vector <- c(param_vector, village_vector)
    }
    param_vector <- c(param_vector,
                      rc = 1/120,
                      r1 = 0.0019 * 30.42,
                      r2 = 0.0040 * 30.42,
                      screen1 = 0.95,
                      screen2 = 0.99,
                      N = length(villages))
    return(param_vector)
}

##' Evaluate prior probability density of a parameter vector for the trypanosomiasis model
##'
##' @param theta Parameter vector
##' @param log whether to take return the logarithm of the prior probability density
##' @return Probability density
##' @author Sebastian Funk
dprior <- function(theta, log = FALSE)
{
    if ("N" %in% names(theta)) N <- round(theta[["N"]]) else N <- 1
    param_prior <- c()
    param_prior <- c(param_prior, dunif(theta[["pc"]], 0, 0.5, log = log))
    param_prior <- c(param_prior, dnorm(theta[["rc"]], 1/120, log = log))
    param_prior <- c(param_prior, dnorm(theta[["r1"]], 0.0019 * 30.42, log = log))
    param_prior <- c(param_prior, dnorm(theta[["r2"]], 0.0040 * 30.42, log = log))
    if ("alpha" %in% names(theta))
        param_prior <- c(param_prior, dunif(theta[["alpha"]], 0, 0.5, log = log))
    if ("delta" %in% names(theta))
        param_prior <- c(param_prior, dunif(theta[["delta"]], 0, 1, log = log))
    param_prior <- c(param_prior, dunif(theta[["screen1"]], 0.86, 0.98, log = log))
    param_prior <- c(param_prior, dnorm(theta[["screen2"]], 0.99, log = log))
    if (N == 1)
    {
        if ("lambda" %in% names(theta))
        {
            param_prior <- c(param_prior,
                             dunif(theta[["lambda"]], 0, 1e-2, log = log))
        }
        if ("beta" %in% names(theta))
        {
            param_prior <- c(param_prior,
                             dunif(theta[["beta"]], 0, 1e-2, log = log))
        }
        if ("p1" %in% names(theta))
        {
            param_prior <- c(param_prior, dunif(theta[["p1"]], 0, 100, log = log))
        }
        if ("p2" %in% names(theta))
        {
            param_prior <- c(param_prior, dunif(theta[["p2"]], 0, 100, log = log))
        }
    } else
    {
        lambda_names <-
            param_prior <-
                c(param_prior,
                  unlist(sapply(grep("^lambda\\.", names(theta), value = TRUE),
                                function(x)
                  {
                      dunif(theta[[x]], 10^(-5), 10^(0), log = log)
                  })))
        param_prior <-
            c(param_prior,
              unlist(sapply(grep("^p1\\.", names(theta), value = TRUE),
                            function(x)
              {
                  dunif(theta[[x]], 0, 100, log = log)
              })))
        param_prior <-
            c(param_prior,
              unlist(sapply(grep("^p2\\.", names(theta), value = TRUE),
                            function(x)
              {
                  dunif(theta[[x]], 0, 100, log = log)
              })))
    }

    if (log)
        return(sum(param_prior))
    else
        return(prod(param_prior))

}

##' Draw an initial condition sample from the prior density distribution of the trypanosomiasis model
##'
##' @param theta Parameter vector
##' @param village_number village number (if parameters for multiple villages are given
##' @return initial conditions vector
##' @author Sebastian Funk
rinit_consistent <- function(theta, village_number = 1)
{

    data(village_data)

    N <- village_screening[village.number == village_number, N]

    village_lambda <- paste("lambda", village_number, sep = ".")
    if (village_lambda %in% names(theta))
    {
        lambda <- theta[[village_lambda]]
    } else if ("lambda" %in% names(theta))
    {
        lambda <- theta[["lambda"]]
    } else
    {
        lambda <- 0
    }

    village_beta <- paste("beta", village_number, sep = ".")
    if (village_beta %in% names(theta))
    {
        beta <- theta[[village_beta]]
    } else if ("beta" %in% names(theta))
    {
        beta <- theta[["beta"]]
    } else
    {
        beta <- 0
    }

    if ("pc" %in% names(theta))
    {
        pc <- theta[["pc"]]
        rc <- theta[["rc"]]
    } else
    {
        pc <- 0
        rc <- 1
    }

    r1 <- theta[["r1"]]
    r2 <- theta[["r2"]]

    if ("p1" %in% names(theta))
    {
        p1 <- theta[["p1"]]
    } else
    {
        p1 <- 0
    }
    if ("p2" %in% names(theta))
    {
        p2 <- theta[["p2"]]
    } else
    {
        p2 <- 0
    }
    if ("delta" %in% names(theta))
    {
        delta <- theta[["delta"]]
    } else
    {
        delta <- 0
    }

    x1 <- pc / (1 - pc) * (r1 + p1) / rc
    x2 <- r1 / (r2 + p2)
    x3 <- x1 + x2 + 1
    x4 <- beta * (1 + delta * x1)
    x5 <- N * pc * lambda
    x6 <- pc * (x4 * N - x3 * lambda) - r1 - p1

    if (beta > 0)
    {
        ## beta > 0
        x7 <- - pc * x3 * x4
        p <- x6 / x7
        q <- x5 / x7

        I1_eq <- - p/2 + sqrt((p**2)/4 - q)
    } else
    {
        ## beta == 0
        I1_eq <- - x5 / x6
    }

    Ic_eq <- x1 * I1_eq
    I2_eq <- x2 * I1_eq

    stage1_detected <- village_screening[village.number == village_number,
                                         detected1_1]
    stage2_detected <- village_screening[village.number == village_number,
                                         detected2_1]


    ## generate random initial conditions that are consistent with the
    ## initial number of detected

    done <- FALSE

    prob_poisson_I1c <- ifelse(stage1_detected > 0, sum(sapply(seq_len(stage1_detected) - 1, function(x) dpois(x, lambda = Ic_eq + I1_eq))), 0)
    prob_poisson_I2 <- ifelse(stage2_detected > 0, sum(sapply(seq_len(stage2_detected) - 1, function(x) dpois(x, lambda = I2_eq))), 0)

    eq <- c(I1 = I1_eq, Ic = Ic_eq, I2 = I2_eq)
    
    if (suppressWarnings(is.finite(qpois(prob_poisson_I1c, I1_eq + Ic_eq) + qpois(prob_poisson_I2, I2_eq))))
    {
        initI1 <- 0
        initI2 <- 0

        while(!done || initI1c + initI2 > N)
    {
        initI1c <- qpois(runif(1, min = prob_poisson_I1c, max = 1), Ic_eq + I1_eq)
        initI2 <- qpois(runif(1, min = prob_poisson_I2, max = 1), I2_eq)
        done <- TRUE
    }

        ## rounding errors??
        ## if (stage1_detected > initI)
        ## {
        ##     initI <- stage1_detected
        ## }
        
        ## if (stage2_detected > initI2)
        ## {
        ##     initI2 <- stage2_detected
        ## }
        
        initS <- N - initI1c - initI2 + stage1_detected + stage2_detected

        initI1 <- rbinom(1, initI1c, 1 - pc)
        initIc <- initI1c - initI1

        initVec <- c(S = initS, Ic = initIc, I1 = initI1, I2 = initI2)
        ## prior density of initial state is density of initial
        ## observation times likelihood of that state
        ## logprior <- dinit(initVec, village_number, theta, TRUE) +
        logprior <-
            dpois(initIc, Ic_eq, log = TRUE) +
            dpois(initI1, I1_eq, log = TRUE) +
            dpois(initI2, I2_eq, log = TRUE)

        Ic_ind <- c(rep(1, initIc), rep(0, initI1))
        Ic_prob <- c(rep(theta[["alpha"]], initIc), rep(1, initI1))
        if (length(Ic_ind) == 0)
        {
            Ic_det <-  0
        } else
        {
            Ic_det <- sum(Ic_ind[sample(Ic_ind, stage1_detected, prob = Ic_prob)])
        }
        I1_det <- stage1_detected - Ic_det
        initIc <- initIc - Ic_det
        initI1 <- initI1 - I1_det
        initI2 <- initI2 - stage2_detected

        modInitVec <- c(S = initS, Ic = initIc, I1 = initI1, I2 = initI2)
    } else
    {
        logprior <- Inf
        modInitVec <- c(S = Inf, Ic = Inf, I1 = Inf, I2 = Inf)
    }

    return(list(state = modInitVec, logprior = logprior, eq = eq))
}

##' Draw an initial condition sample from the prior density distribution of the trypanosomiasis model
##'
##' @param theta Parameter vector
##' @param village_number village number (if parameters for multiple villages are given
##' @return initial conditions vector
##' @author Sebastian Funk
rinit <- function(theta, village_number = 1)
{
    sim_village_screening <- village_screening[village.number %in% village]

    N <- sim_village_screening[["N"]]

    lambda <- theta[["lambda"]]
    beta <- theta[["beta"]]

    if ("pc" %in% names(theta))
    {
        pc <- theta[["pc"]]
        rc <- theta[["rc"]]
    } else
    {
        pc <- 0
        rc <- 1
    }

    r1 <- theta[["r1"]]
    r2 <- theta[["r2"]]

    if ("p1" %in% names(theta))
    {
        p1 <- theta[["p1"]]
    } else
    {
        p1 <- 0
    }
    if ("p2" %in% names(theta))
    {
        p2 <- theta[["p2"]]
    } else
    {
        p2 <- 0
    }
    if ("delta" %in% names(theta))
    {
        delta <- theta[["delta"]]
    } else
    {
        delta <- 0
    }

    x1 <- pc / (1 - pc) * (r1 + p1) / rc
    x2 <- r1 / (r2 + p2)
    x3 <- x1 + x2 + 1
    x4 <- beta * (1 + delta * x1)
    x5 <- N * pc * lambda
    x6 <- pc * (x4 * N - x3 * lambda) - r1 - p1

    if (beta > 0)
    {
        ## beta > 0
        x7 <- - pc * x3 * x4
        p <- x6 / x7
        q <- x5 / x7

        I1_eq <- - p/2 + sqrt((p**2)/4 - q)
    } else
    {
        ## beta == 0
        I1_eq <- - x5 / x6
    }

    Ic_eq <- x1 * I1_eq
    I2_eq <- x2 * I1_eq

    testInit <- c(S = N - round(I1_eq) - round(Ic_eq) - round(I2_eq),
                  I1 = round(I1_eq), Ic = round(Ic_eq), I2 = round(I2_eq))

    initTraj <- sim_trajectory(theta, testInit, village = village_number,
                               stoptime = 12)
    init <- unlist(tail(initTraj, 1)[names(testInit)])

    return(init)
}

##' Evaluate the likelihood of a model state of the trypanosomiasis model
##'
##' @param state Model state
##' @param stage1_detected Number of stage 1 cases detected
##' @param stage2_detected Number of stage 2 cases detected
##' @param theta parameter vector
##' @param attendance screening attendance
##' @param log whether to return the logarithm
##' @return point likelihood
##' @author Sebastian Funk
likelihood <- function(state, stage1_detected, stage2_detected,
                       theta, attendance = 1, active = TRUE,
                       log = FALSE)
{
    if (log)
    {
        res <- 0
    } else
    {
        res <- 1
    }

    if (!any(is.na(state)) && !any(state < 0))
    {
        if (!missing(stage1_detected))
        {
            if (active) {
                detection_probability <- theta[["screen1"]] * attendance
            } else
            {
                detection_probability <- 1 ## sensitivity of passive screening
            }
            if ("Ic" %in% names(state))
            {
                stage1_ll <- sum(sapply(seq(0, stage1_detected), function(x)
                {
                    prob <-
                        dbinom(x, state[["Ic"]] + x,
                               theta[["alpha"]] * detection_probability)
                    prob <- prob * dbinom(stage1_detected - x, state[["I1"]] + x,
                                          detection_probability)
                }))
                if (log)
                {
                    res <- res + log(stage1_ll)
                } else
                {
                    res <- res * stage1_ll
                }
            } else
            {
                stage1_ll <- dbinom(stage1_detected, state[["I1"]] + stage1_detected,
                                    detection_probability, log = log)
                if (log)
                {
                    res <- res + stage1_ll
                } else
                {
                    res <- res * stage1_ll
                }
            }
        }
        if (!missing(stage2_detected))
        {
            if (active) {
                detection_probability <- theta[["screen1"]] * attendance
            } else
            {
                detection_probability <- 1 ## sensitivity of passive screening
            }
            stage2_ll <- dbinom(stage2_detected, state[["I2"]] + stage2_detected,
                                detection_probability, log = log)
            if (log)
            {
                res <- res + stage2_ll
            } else
            {
                res <- res * stage2_ll
            }
        }
    } else {
        if (log)
        {
            res <- -Inf
        } else
        {
            res <- 0
        }
    }

    return(res)
}

##' Evaluate the probability density of a given set of initial conditions of the trypanosomiasis model
##'
##' @param init initial conditions
##' @param village_data Village data
##' @param theta parameter vector
##' @param log whether to return the logarithm
##' @return Probability density
##' @author Sebastian Funk
dinit <- function(init, village_number, theta, log = FALSE)
{
    data(village_data)

    stage1_detected <- village_screening[village.number == village_number,
                                         detected1_1]
    stage2_detected <- village_screening[village.number == village_number,
                                         detected1_2]

    attendance <- min(village_screening[village.number == village_number,
                                        sigma.start], 1)

    return(likelihood(init, stage1_detected, stage2_detected,
                      theta, attendance, log = log))
}

##' Calculate the posterior density for a given set of parameters for the trypanosomiasis model
##'
##' @param theta parameter vector
##' @param village_data village data
##' @param cum_data cumulative data
##' @param nruns number of runs (to estimate the likelihood), 1 by default
##' @param log whether to return the logarithm of the posterior density
##' @param ...
##' @return posterior density
##' @author Sebastian Funk
traj_likelihood <- function(theta, village_number, nruns = 1, log = FALSE, ...)
{
    data(village_data)

    params <-
        grep(paste0("\\.", village_number, "$"), names(theta),
             value = TRUE)
    naked.params <- c()
    for (param in params)
    {
        naked.param <-
            sub(paste0("\\.", village_number, "$"), "", param)
        theta[[naked.param]] <- theta[[param]]
        naked.params <- c(naked.params, naked.param)
    }

    passive <- lapply(1:2, function(x)
    {
        village_cases[village.number == village.number & stage == x, cases]
    })

    res <- 0
    likelihoods <- c()
    final.attendance <- min(village_screening[village.number == village_number,
                                              sigma.end], 1)

    stoptime <- village_screening[village.number == village_number, stoptime]
    session.2.detected1 <-
        village_screening[village.number == village_number, detected1_2]

    for (j in seq_len(nruns))
    {
        init <- rinit(theta, village_number)
        log.init <- dinit(init, village_number, theta, TRUE)
        if (all(is.finite(init)))
        {
            ll <- 0
            run <- sim_trajectory(theta, init, village_number)

            if (any(c(run[, "I1"], run[, "I2"]) < 0)) {
                ll <- -Inf
            } else
            {
                ll <- ll + likelihood(run[nrow(run), ], session.2.detected1,
                                      theta = theta,
                                      attendance = final.attendance, log = TRUE)
            }

            for (stage in 1:2)
            {
                if (paste0("p", stage) %in% names(theta))
                {
                    passive_state <- run[nrow(run), paste0("Z", stage, "pass")]
                    names(passive_state) <- paste0("I", stage)
                    likelihood_options <- list(state = passive_state,
                                               detected = sum(passive[[stage]]),
                                               theta = theta, log = TRUE,
                                               active = FALSE)
                    detected_option <- which(names(likelihood_options) == "detected")
                    names(likelihood_options)[detected_option] <-
                        paste0("stage", stage, "_detected")
                    ll <- ll + do.call(likelihood, likelihood_options)
                }
            }

        } else {
            ll <- -Inf
        }
        likelihoods[j] <- log.init + ll
    }
    res <- sum(exp(likelihoods)) / nruns

    if (length(naked.params) > 0)
    {
        theta <- theta[-(which(names(theta) %in% naked.params))]
    }

    if (log)
    {
        return(log(res))
    } else
    {
        return(res)
    }
}

##' Calculate posterior probability density for a given set of parameters for the trypanosomiasis model, using a model that encompasses all villages
##'
##' @param theta parameter vector
##' @param village_number villages to simuate
##' @param ... 
##' @return posterior density
##' @author Sebastian Funk
param_posterior <- function(theta, init, village_number, ...)
{
    data(village_data)

    sim_village_screening <- village_screening[village.number %in% village]
    
    N <- sim_village_screening[["N"]]

    res <- list()

    attendance <-
        list(init = min(sim_village_screening[["sigma.start"]], 1),
             final = min(sim_village_screening[["sigma.end"]], 1))

    screening_prob <- list()

    passive_data <-
        village_cases[village.number == village_number,
                      list(cases = sum(cases)), by = stage]

    ## sample
    log.prior <- dprior(theta, log = TRUE)
    res <- list(theta = theta, log.prior = log.prior)
    if (missing(init)) init <- rinit(theta, village_number)
    log.init <- dinit(init, 1, theta, log = TRUE)
    res[["log.init"]] <- log.init
    res[["init"]] <- init

    ## work out initial detected
    if (is.finite(log.prior) && is.finite(log.init))
    {
        for (moment in c("init", "final"))
        {
            screening_prob[[moment]] <- 
                c(I1 = theta[["screen1"]] *
                      min(attendance[[moment]], 1),
                  Ic = theta[["screen1"]] *
                      min(attendance[[moment]] * theta[["alpha"]], 1), 
                  I2 = theta[["screen2"]] *
                      min(attendance[[moment]], 1))
        }

        data_detected <-
            list(init = c(I1c = sim_village_screening[["detected1_1"]],
                          I2 = sim_village_screening[["detected2_1"]]),
                 final = c(I1c = sim_village_screening[["detected1_2"]],
                           I2 = sim_village_screening[["detected.2_2"]]))
                
        init_detected <- c()
        for (compartment in names(screening_prob[["init"]]))
        {
            init_detected[compartment] <-
                rbinom(1, init[[compartment]],
                       screening_prob[["init"]][[compartment]])
        }
            
        modelled <- list(init = init)
        init <- init + c(sum(init_detected), -init_detected)
        traj <- sim_trajectory(theta, init, village = village_number)
        modelled[["final"]] <- unlist(tail(traj, 1))[names(init)]
        res[["final"]] <- modelled[["final"]]
        
        ## active screening posterior
        res[["log.posterior"]] <- log.prior + log.init
        for (moment in names(modelled))
        {
            I1c_seq <- seq_len(data_detected[[moment]][["I1c"]] + 1) - 1
            res[["log.posterior"]] <- res[["log.posterior"]] +
                sum(sapply(I1c_seq, function(x) {
                    dbinom(x, modelled[[moment]]["Ic"], 
                           screening_prob[[moment]][["Ic"]], log = TRUE) +
                    dbinom(x, modelled[[moment]]["I1"], 
                           screening_prob[[moment]][["I1"]], log = TRUE)
                })) +
                dbinom(data_detected[[moment]][["I2"]],
                       modelled[[moment]]["I2"],
                       screening_prob[[moment]][["I2"]], log = TRUE)
        }
        
        ## passive screening posterior
        for (passive.stage in 1:2)
        {
            if (grep(paste0("^p", passive.stage), names(theta)))
            {
                passive <-
                    unlist(tail(traj[[paste0("Z", passive.stage, "pass")]],
                                1))
                res[["log.posterior"]] <- res[["log.posterior"]] +
                    dpois(passive_data[stage == passive.stage, cases],
                          passive, log = TRUE)
            }
        }
    } else
    {
        traj <- data.frame(t(init))
        res[["log.posterior"]] <- log.prior + log.init
    }
           
    return(res)
}

##' Calculate summary statistics for a given set of parameters for the trypanosomiasis model, using a model that encompasses all villages
##'
##' @param theta parameter vector
##' @param nruns number of runs
##' @param villages villages to simuate
##' @param ... 
##' @return summary statistics
##' @author Sebastian Funk
param_sumstat_villages <- function(theta, init, nruns = 1, village_number, ...)
{
    data(village_data)

    sim_village_screening <- village_screening[village.number %in% village]
    
    N <- sim_village_screening[["N"]]

    res <- list()

    init.attendance <-  sim_village_screening[["sigma.start"]]
    final.attendance <- sim_village_screening[["sigma.end"]]
    
    init.stage1.prob <- theta[["screen1"]] * min(init.attendance, 1)
    init.stagec.prob <-
        theta[["screen1"]] * min(init.attendance * theta[["alpha"]], 1)
    init.stage2.prob <- theta[["screen2"]] * min(init.attendance, 1)
    
    final.stage1.prob <- theta[["screen1"]] * min(final.attendance, 1)
    final.stagec.prob <-
        theta[["screen1"]] * min(final.attendance * theta[["alpha"]], 1)
    final.stage2.prob <- theta[["screen2"]] * min(final.attendance, 1)

    passive_data <- village_cases[village.number == village_number]

    missing_init <- missing(init)

    sims <- lapply(seq_len(nruns), function(x)
    {
      log.prior <- dprior(theta, log = TRUE)
      if (missing_init) init <- rinit(theta, village_number)
      log.init <- dinit(init, 1, theta, log = TRUE)
      ret <- list(logprior = log.init + log.prior, init = init)
      if (all(is.finite(init)) && is.finite(log.init + log.prior))
      {
        ret[["traj"]] <-
          sim_trajectory(theta, init, village = village_number)
      } else
      {
        ret[["traj"]] <- data.frame(t(init))
      }
      ret
    })

    res[["logprior"]] <- sapply(sims, function(x) x[["logprior"]])

    if (is.finite(res[["logprior"]]))
    {
      res[["nneg"]] <- sum(sapply(sims, function(x)
      {
        any(x[["traj"]][["I1"]] < 0 | x[["traj"]][["I2"]] < 0)
      }))
      
      res[["init_active_stage1"]] <- sapply(sims, function(x)
      {
          rbinom(1, init[["Ic"]], init.stagec.prob) +
              rbinom(1, init[["I1"]], init.stage1.prob)
      })
            
      res[["init_active_stage2"]] <- sapply(sims, function(x)
      {
          rbinom(1, init[["I2"]], init.stage2.prob)
      })
      
      res[["final_active_stage1"]] <- sapply(sims, function(x)
      {
        rbinom(1, tail(x[["traj"]][["Ic"]], 1), final.stagec.prob) +
          rbinom(1, tail(x[["traj"]][["I1"]], 1), final.stage1.prob)
      })
            
      res[["final_active_stage2"]] <- sapply(sims, function(x)
      {
        rbinom(1, tail(x[["traj"]][["I2"]], 1), final.stage2.prob)
      })
      
      for (stage in 1:2)
      {
        if (grep(paste0("^p", stage), names(theta)))
        {
          passive_stage <- paste0("passive_stage", stage)
          res[[passive_stage]] <- sapply(sims, function(x)
          {
            tail(x[["traj"]][[paste0("Z", stage, "pass")]], 1)
          })
        }
      }
    } else
    {
      res[["nneg"]] <- NA_real_
      for (stage in 1:2)
      {
        passive_stage <- paste0("passive_stage", stage)
        res[[passive_stage]] <- NA_real_
        active_stage <- paste0("active_stage", stage)
        res[[active_stage]] <- NA_real_
      }
    }
    return(list(stat = res, traj = lapply(sims, function(x) {x[["traj"]]})))
}


##' Draw samples for the trypanosomiasis model with chronic carriers
##'
##' @param nsamples number of samples
##' @param seed random number seed
##' @param verbose whether to print out posteriors as it goes along
##' @param passive wheter to fit to passive detection data
##' @param calc.posterior whether to calculate the posterior. Alernatively, will calculate summary statistic
##' @param villages villages to simuate
##' @param progress.bar whether to show a progress bar (TRUE by default)
##' @param sample how to sample (lhs or prior)
##' @param ... further options for rprior
##' @return samples, posteriors, random numbers
##' @importFrom lhs randomLHS
##' @export
##' @author Sebastian Funk
chronic_carriers_sample <- function(nsamples = 1, seed,
                                    verbose = FALSE, passive = TRUE,
                                    calc.posterior = TRUE, demand.finite = FALSE, 
                                    village_number, progress.bar = TRUE,
                                    sample = c("lhs", "prior"),
                                    ...)
{
    data(village_data)

    sample <- match.arg(sample)

    samples <- list()
    likelihoods <- c()
    posteriors <- c()

    if (!missing(seed))
    {
        set.seed(seed)
    }

    if (sample == "lhs")
    {
      upper <- c(0.5, 0.5, 1, 1e-2, 1e-2)
      lower <- c(0, 0, 0, 0, 0)
      theta_names <- c("pc", "alpha", "delta", "lambda", "beta")
      if (passive)
      {
        upper <- c(upper, c(30, 30))
        lower <- c(lower, c(0, 0))
        theta_names <- c(theta_names, c("p1", "p2"))
      }
      
    }

    if (progress.bar && nsamples > 1) pb <- txtProgressBar(min = 0, max = nsamples - 1, style = 3)

    i <- 0
    accepted <- 0
    tosample <- nsamples
    while (tosample > 0)
    {
      if (sample == "lhs")
      {
        r <- lhs::randomLHS(tosample, length(theta_names))
        r <- t(apply(r, 1, function(x) { lower + (upper - lower) * x}))
        colnames(r) <- theta_names
      }
      i <- 0
      while(i < tosample)
      {
        i <- i + 1
        theta <- rprior(villages = village_number, ...)
        if (sample == "lhs") {
          theta[colnames(r)] <- r[i, ]
        }
        init <- rinit(theta, village_number)
        if (calc.posterior)
        {
          posterior <- param_posterior(theta, init = init,
                                       log = TRUE,
                                       village_number = village_number, 
                                       ...)
          
          if (!demand.finite || is.finite(posterior[["log.posterior"]]))
          {
            accepted <- accepted + 1
            samples[[accepted]] <- posterior
            if (verbose)
            {
              message(i, samples[[i]][["log.posterior"]], "\n")
            }
          }
        } else
        {
          samples[[i]] <-
            c(list(parameters = theta, init = init), 
              param_sumstat_villages(theta = theta, 
                                     init = init, 
                                     village_number = village_number,
                                     ...)[["stat"]])
          accepted <- accepted + 1
        }
        if (progress.bar && nsamples > 1) setTxtProgressBar(pb, accepted)
      }
      tosample <- nsamples - accepted
    }
    if (progress.bar && nsamples > 1) close(pb)

    return(samples)
}

##' run ABC-MCMC on the model for chronic carriers of trypanosomiasis
##'
##' @param start initial state of the chain
##' @param n_iterations number of iterations
##' @param sd standard deviation of each parameter
##' @param upper upper limit on parameters (if given)
##' @param lower lower limit on parameters (if given)
##' @param epsilon value of epsilon
##' @param data_summary summary statistic of the data
##' @param thin thinning (1 for no thinning)
##' @param return.traj whether to return trajectories
##' @param verbose whether to print verbose output
##' @param ... parameters to be passed to \code{param_sumstat_villages}
##' @return chain
##' @importFrom truncnorm rtruncnorm dtruncnorm
##' @author Sebastian Funk
chronic_carriers_abc_mcmc <- function(start, n_iterations, sd,
                                  upper, lower,
                                  epsilon, data_summary,
                                  thin = 1, return.traj = FALSE, 
                                  verbose = FALSE, ...)
{
    theta <- start
    prior_theta <- dprior(theta, log = TRUE)
    prior_init <- rinit(theta)$logprior
    prior_init_propose <- prior_init

    accepted <- 0

    chain <- list()

    if (missing(upper))
    {
        upper <- rep(Inf, length(theta))
    }
    if (missing(lower))
    {
        lower <- rep(-Inf, length(theta))
    }

    trajectories <- list()
    summary_statistics <- list()
    traj <- NULL
    stat <- NULL

    chain_counter <- 0
    for (i in seq_len(n_iterations))
    {
        is.accepted <- FALSE
        theta_propose <- rtruncnorm(length(theta[sd > 0]), lower[sd > 0], upper[sd > 0],
                                    theta[sd > 0], sd[sd > 0])
        names(theta_propose) <- names(theta[sd > 0])
        theta_propose <- c(theta_propose, theta[setdiff(names(theta),
                                                        names(theta_propose))])
        theta_propose <- theta_propose[names(theta)]

        prior_propose <- dprior(theta_propose, log = TRUE)

        if (is.finite(prior_propose)) {
            hastings_ratio <-
                ifelse(any(sd > 0),
                       log(dtruncnorm(theta[sd > 0], lower, upper,
                                      theta_propose[sd > 0], sd[sd > 0])) -
                       log(dtruncnorm(theta_propose[sd > 0], lower, upper,
                                      theta[sd > 0], sd[sd > 0])),
                       0)

            sumstats <-
                param_sumstat_villages(theta_propose, ...)
            prior_init_propose <- sumstats[["stat"]][["loginit"]]

            diff <- unlist(stat)[names(data_summary)] - data_summary
            pass <- all(abs(diff) <= epsilon)
            if (pass)
            {
                log.acceptance <- prior_propose - prior_theta + hastings_ratio + prior_init_propose - prior_init
                is.accepted <- (log(runif (1)) < log.acceptance)
                traj <- sumstats[["traj"]][[1]]
                stat <- sumstats[["stat"]]
            }
        }
        if (is.accepted)
        {
            accepted <- accepted + 1
            theta <- theta_propose
            prior_theta <- prior_propose
            prior_init <- prior_init_propose
            if (missing(epsilon)) posterior <- posterior_propose
        }
        if (i %% thin == 0)
        {
            chain_counter <- chain_counter + 1
            chain[[chain_counter]] <- theta
            if (return.traj)
            {
                trajectories[[chain_counter]] <- traj
            }
            summary_statistics[[chain_counter]] <- stat
        }
        if (verbose) cat(i, "acc:", is.accepted, accepted / i, "\n")
    }

    if (length(chain) > 0)
    {
        chain_names <- names(chain[[1]])
        df_chain <- data.frame(matrix(unlist(chain), ncol = length(chain_names), byrow = TRUE))
        colnames(df_chain) <- chain_names
    }

    ret <- list(acceptance.rate = accepted / n_iterations,
                trace = df_chain,
                summary.statistics = summary_statistics)
    if (return.traj)
    {
        ret[["trajectories"]] <- trajectories
    }
    return(ret)
}

sample_equilibrium <- function(theta, village)
{
  sim_village_screening <- village_screening[village.number %in% village]

  N <- sim_village_screening[["N"]]
  detectedI1c <- sim_village_screening[["detected1_1"]] ## number of people detected with stage 1 or chronic
  detectedI2<- sim_village_screening[["detected2_1"]] ## number of people detected with stage 2

  attendance <-  sim_village_screening[["sigma.start"]]

  stage1_prob <- theta[["screen1"]] * min(attendance, 1)
  stagec_prob <- theta[["screen1"]] * min(attendance * theta[["alpha"]], 1)
  stage2_prob <- theta[["screen2"]] * min(attendance, 1)

  ## number of detected in I1 that are chronic
  detectedIc <- sample(detectedI1c, 1) - 1

  ## sample stage2
  initI2_random <- runif(1, min = 0, max = 1 / stage2_prob)

  found <- FALSE
  prob_sum <- 0
  initI2 <- detectedI2

  found <- FALSE
  while (!found)
  {
    prob <- dbinom(detectedI2, initI2,  stage2_prob)
    prob_sum <- prob_sum + prob
    if (prob_sum > initI2_random) {
      found <- TRUE
    } else {
      initI2 <- initI2 + 1
    }
  }

  ## sample stage1
  found <- FALSE
  prob_sum <- 0
  initI1 <- detectedI1c - detectedIc

  ## normalisation
  while (!found)
  {
    if (initI1 < initI2 * theta[["r2"]] / theta[["r1"]])
    {
      prob <- dbinom(detectedI1c - detectedIc, initI1, stage1_prob)
      prob_sum <- prob_sum + prob
      initI1 <- initI1 + 1
    } else {
      found <- TRUE
    }
  }

  if (prob_sum < (1 / stage1_prob - 1e-15))
  {
    initI1_random <- runif(1, min = prob_sum, max = 1 / stage1_prob)
    
    found <- FALSE
    while (!found)
    {
      prob <- dbinom(detectedI1c - detectedIc, initI1,  stage1_prob)
      prob_sum <- prob_sum + prob
      if (prob_sum > initI1_random) {
        found <- TRUE
      } else {
        initI1 <- initI1 + 1
      }
    }
  }

  ## sample stagec
  initIc_random <- runif(1, min = 0, max = 1 / stagec_prob)

  found <- FALSE
  prob_sum <- 0
  initIc <- detectedIc

  while (!found)
  {
    prob <- dbinom(detectedIc, initIc, stagec_prob)
    prob_sum <- prob_sum + prob
    if (prob_sum > initIc_random || (initIc + initI1 + initI2 == N)) { 
      found <- TRUE
    } else {
      initIc <- initIc + 1
    }
  }

  loginit <- dbinom(detectedI1c - detectedIc, initI1, stage1_prob, log = TRUE) +
    dbinom(detectedIc, initIc, stagec_prob, log = TRUE) +
    dbinom(detectedI2, initI2, stage2_prob, log = TRUE)

  return(list(eq = c(S = N - initI1 - initIc - initI2, I1 = initI1, Ic = initIc,
                     I2 = initI2),
              state = c(S = N - initI1 - initIc - initI2 + detectedI1c +
                         detectedI2,
                       I1 = initI1 - detectedI1c + detectedIc,
                       Ic = initIc - detectedIc,
                       I2 = initI2 - detectedI2),
              logprior = loginit))
}

rprior_conditional <- function(village)
{
  theta <- c(rc = 1/120,
             r1 = 0.0019 * 30.42,
             r2 = 0.0040 * 30.42,
             screen1 = 0.95,
             screen2 = 0.99)
  theta["alpha"] <- runif(1, 0, 0.5)
  theta["delta"] <- runif(1, 0, 1)
  init <- sample_equilibrium(theta, village)

  found <- FALSE
  counter <- 0
  while (!found)
  {
    eq <- abs(runif(length(init[["eq"]]),
                    init[["eq"]] - 0.5, init[["eq"]] + 0.5))
    names(eq) <- names(init[["eq"]])

    x2 <- eq[["I2"]] / eq[["I1"]]
    theta["p2"] <- (theta[["r1"]] - x2 * theta[["r2"]]) / x2

    found <- (theta[["p2"]] > 0)
  }

  theta["p1"] <- runif(1, 0, theta[["p2"]])

  x1 <- eq[["Ic"]] / eq[["I1"]]
  if (is.finite(x1))
  {
    theta["pc"] <- x1 / (theta[["r1"]] + theta[["p1"]] + x1)
  } else
  {
    theta["pc"] <- 1
  }
  prop_background <- runif(1, 0, 1)

  lambda <- (theta[["rc"]] * eq[["Ic"]] + theta[["p1"]] * eq[["I1"]] +
             (theta[["p2"]] + theta[["r2"]]) * eq[["I2"]]) / eq[["S"]]

  theta["lambda"] <- lambda * prop_background
  theta["beta"] <- (lambda - theta[["lambda"]]) /
    (eq[["I1"]] + theta[["delta"]] * eq[["Ic"]])

  return(list(theta = theta, init = init))
}
