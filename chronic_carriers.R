library('truncnorm')

##' Produe a single simuation run
##'
##' @param theta parameter vector
##' @return data frame with trajctor
##' @author seb
##' @export
sim_trajectory <- function(theta, init, village, ...)
{
    stoptime <- village_screening[village_screening$village.number == village]$stoptime

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
        c(pc = runif(1, 0, 0.5), alpha = runif(1, 0, 1))
    param_vector["delta"] <- ifelse(chronic, runif(1, 0, 1), 0)

    if (length(villages) == 1)
    {
        param_vector <- c(param_vector,
                          lambda = ifelse(background,
                                          10^(runif(1, -5, 0)),
                                          0),
                          beta = ifelse(transmitted,
                                        10^(runif(1, -5, 0)),
                                        0),
                          p1 = ifelse(passive, runif(1, 0, 30), 0),
                          p2 = ifelse(passive, runif(1, 0, 30), 0))
    } else
    {
        village_vector <- c(ifelse(background,
                                   10^(runif(length(villages), -5,  0)),
                                   rep(0, length(villages))),
                            ifelse(transmitted,
                                   10^(runif(length(villages), -5, 0)),
                                   rep(0, length(villages))),
                            ifelse(passive,
                                   runif(length(villages), 0, 30),
                                   rep(0, length(villages))),
                            ifelse(passive,
                                   runif(length(villages), 0, 30),
                                   rep(0, length(villages))))
        names(village_vector) <-
            paste(rep(c("lambda", "beta", "p1", "p2"), each = length(villages)),
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
    N <- round(theta[["N"]])
    param_prior <- c()
    param_prior <- c(param_prior, dunif(theta[["pc"]], 0, 0.5, log = TRUE))
    param_prior <- c(param_prior, dnorm(theta[["rc"]], 1/120, log = TRUE))
    param_prior <- c(param_prior, dnorm(theta[["r1"]], 0.0019 * 30.42, log = TRUE))
    param_prior <- c(param_prior, dnorm(theta[["r2"]], 0.0040 * 30.42, log = TRUE))
    if ("alpha" %in% names(theta))
        param_prior <- c(param_prior, dunif(theta[["alpha"]], 0, 1, log = TRUE))
    if ("delta" %in% names(theta))
        param_prior <- c(param_prior, dunif(theta[["delta"]], 0, 1, log = TRUE))
    param_prior <- c(param_prior, dunif(theta[["screen1"]], 0.86, 0.98, log = TRUE))
    param_prior <- c(param_prior, dnorm(theta[["screen2"]], 0.99, log = TRUE))
    if (N == 1)
    {
        if ("lambda" %in% names(theta))
        {
            param_prior <- c(param_prior,
                             dunif(theta[["lambda"]], 10^(-5), 10^(0), log = TRUE))
        }
        if ("beta" %in% names(theta))
        {
            param_prior <- c(param_prior,
                             dunif(theta[["beta"]], 10^(-5), 10^(0), log = TRUE))
        }
        if ("p1" %in% names(theta))
        {
            param_prior <- c(param_prior, dunif(theta[["p1"]], 0, 30, log = TRUE))
        }
        if ("p2" %in% names(theta))
        {
            param_prior <- c(param_prior, dunif(theta[["p2"]], 0, 30, log = TRUE))
        }
    } else
    {
        lambda_names <-
            param_prior <-
                c(param_prior,
                  unlist(sapply(grep("^lambda\\.", names(theta), value = TRUE),
                                function(x)
                  {
                      dunif(theta[[x]], 10^(-5), 10^(0), log = TRUE)
                  })))
        param_prior <-
            c(param_prior,
              unlist(sapply(grep("^p1\\.", names(theta), value = TRUE),
                            function(x)
              {
                  dunif(theta[[x]], 0, 30, log = TRUE)
              })))
        param_prior <-
            c(param_prior,
              unlist(sapply(grep("^p2\\.", names(theta), value = TRUE),
                            function(x)
              {
                  dunif(theta[[x]], 0, 30, log = TRUE)
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
##' @param passive whether to accumulate infections for passive detection (Z1pass, Z2pass)
##' @return initial conditions vector
##' @author Sebastian Funk
rinit <- function(theta, village_number = 1)
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
    x4 <- beta * (1 + delta * x2)
    x5 <- N * pc * lambda
    x6 <- pc * (x4 * N - x3 * lambda) - r1 - p1

    if (beta > 0)
    {
        ## beta > 0
        x7 <- -pc * x3 * x4
        p <- - x6 / (2 * x7)
        q <- x5 / x7

        I1_eq <- p + sqrt(p**2 - q)
    } else
    {
        ## beta == 0
        I1_eq <- x5 / x6
    }

    Ic_eq <- x1 * I1_eq
    I2_eq <- x2 * I1_eq

    stage1_detected <- village_screening[village.number == village_number,
                                         detected1_1]
    stage2_detected <- village_screening[village.number == village_number,
                                         detected2_1]

    ## generate random initial conditions that are consistent with the
    ## initial number of detected

    prob_poisson_I <- ifelse(stage1_detected > 0, sum(sapply(seq_len(stage1_detected) - 1, function(x) dpois(x, lambda = Ic_eq + I1_eq))), 0)
    initI <- qpois(runif(1, min = prob_poisson_I, max = 1), Ic_eq + I1_eq)
    prob_poisson_I2 <- ifelse(stage2_detected > 0, sum(sapply(seq_len(stage2_detected) - 1, function(x) dpois(x, lambda = I2_eq))), 0)
    initI2 <- qpois(runif(1, min = prob_poisson_I2, max = 1), I2_eq)

    initS <- N - initI - initI2 + stage1_detected + stage2_detected

    initI1 <- rbinom(1, initI, 1 - pc)
    initIc <- initI - initI1

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
    Ic_det <- sum(Ic_ind[sample(Ic_ind, stage1_detected, prob = Ic_prob)])
    I1_det <- stage1_detected - Ic_det
    initIc <- initIc - Ic_det
    initI1 <- initI1 - I1_det
    initI2 <- initI2 - stage2_detected

    modInitVec <- c(S = initS, Ic = initIc, I1 = initI1, I2 = initI2)

    return(list(state = modInitVec, logprior = logprior))
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

##' Calculate summary statistic for a given set of parameters for the trypanosomiasis model, using a model that encompasses all villages
##'
##' @param theta parameter vector
##' @param nruns number of runs
##' @param villages villages to simuate
##' @param ... 
##' @return posterior density
##' @author Sebastian Funk
param_sumstat_villages <- function(theta, nruns = 1, villages, ...)
{
    data(village_data)

    if (missing(villages))
    {
        sim_village_screening <- village_screening
    } else
    {
        sim_village_screening <- village_screening[village.number %in% villages]
    }

    summary.statistics <- apply(sim_village_screening, 1, function(village)
    {
        village_number <- village[["village.number"]]
        res <- list(village.number = village_number)
        stoptime <- village[["stoptime"]]

        final.attendance <- min(village[["sigma.end"]], 1)
        passive_data <- village_cases[village.number == village_number]

        sims <- lapply(seq_len(nruns), function(x)
        {
            init <- rinit(theta)
            log.init <- init$logprior
            ret <- list(loginit = log.init)
            if (all(is.finite(init$state)))
            {
                ret[["traj"]] <-
                    sim_trajectory(theta, init$state, village = village_number)
            } else
            {
                ret[["traj"]] <- data.frame(t(init$state))
            }
            ret
        })

        res[["loginit"]] <- sapply(sims, function(x) x[["loginit"]])

        res[["nneg"]] <- sum(sapply(sims, function(x)
        {
            any(x[["traj"]][["I1"]] < 0 | x[["traj"]][["I2"]] < 0)
        }))

        res[["active_stage1"]] <- sapply(sims, function(x)
        {
            rbinom(1, tail(x[["traj"]][["Ic"]], 1),
                   theta[["alpha"]] * theta[["screen1"]] * final.attendance) +
                rbinom(1, tail(x[["traj"]][["I1"]], 1),
                       theta[["screen1"]] * final.attendance)
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
        return(list(stat = res, traj = lapply(sims, function(x) {x[["traj"]]})))
    })
    if (length(summary.statistics) == 1)
    {
        summary.statistics <- summary.statistics[[1]]
    }
    return(summary.statistics)
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
                                    calc.posterior = TRUE,
                                    villages, progress.bar = TRUE,
                                    sample = c("lhs", "prior"),
                                    ...)
{
    data(village_data)

    sample <- match.arg(sample)

    if (missing(villages))
    {
        villages <- seq_len(nrow(village_screening))
    }
    nb_villages <- length(villages)

    samples <- list()
    likelihoods <- c()
    posteriors <- c()

    if (!missing(seed))
    {
        set.seed(seed)
    }

    if (sample == "lhs")
    {
        if (passive)
        {
            repeat_vec_upper <- c(-2, 1, 2)
            repeat_vec_lower <- c(-4, 0, 0)
            repeat_names <- c("lambda", "p1", "p2")
        } else
        {
            repeat_vec_upper <- -2
            repeat_vec_lower <- -4
            repeat_names <- c("lambda")
        }

        r <- lhs::randomLHS(nsamples, nb_villages * length(repeat_names) + 3)
        upper <- c(1, 1, 1, rep(repeat_vec_upper, each = nb_villages))
        lower <- c(0, 0, 0, rep(repeat_vec_lower, each = nb_villages))
        r <- t(apply(r, 1, function(x) { lower + (upper - lower) * x}))
        theta_names <- c("pc", "alpha", "delta")
        if (nb_villages > 1)
        {
            theta_names <- c(theta_names,
                             paste(rep(repeat_names, each = nb_villages),
                                   villages, sep = "."))
        } else {
            theta_names <- c(theta_names, repeat_names)
        }
        colnames(r) <- theta_names
        r[, grep("^lambda", colnames(r), value = TRUE)] <-
            10^r[, grep("^lambda", colnames(r), value = TRUE)]
    }

    if (progress.bar && nsamples > 1) pb <- txtProgressBar(min = 0, max = nsamples - 1, style = 3)

    i <- 0
    while(i < nsamples)
    {
        i <- i + 1
        theta <- rprior(villages = villages, passive = passive, ...)
        if (sample == "lhs") {
            theta[colnames(r)] <- r[i, ]
        }
        if (calc.posterior)
        {
            posterior <- param_posterior_villages(theta, log = TRUE,
                                                  villages = villages,
                                                  ...)
            samples[[i]] <- list(parameters = theta)
            samples[[i]][["village_posteriors"]] <- posterior
            samples[[i]][["posterior"]] <- sum(posterior)

            if (verbose)
            {
                message(i, samples[[i]][["posterior"]], "\n")
            }
        } else
        {
            samples[[i]] <- list(parameters = theta,
                                 summary_statistics =
                                     param_sumstat_villages(theta = theta,
                                                            villages = villages,
                                                            ...)[["stat"]])
        }
        if (progress.bar && nsamples > 1) setTxtProgressBar(pb, i)
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
    traj <- NULL

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
            
            diff <- unlist(sumstats[["stat"]])[names(data_summary)] - data_summary
            pass <- all(diff <= epsilon)
            if (pass)
            {
                log.acceptance <- prior_propose - prior_theta + hastings_ratio + prior_init_propose - prior_init
                is.accepted <- (log(runif (1)) < log.acceptance)
                traj <- sumstats[["traj"]][[1]]
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
        }
        if (verbose) cat(i, "acc:", is.accepted, accepted / i, "\n")
    }

    if (length(chain) > 0)
    {
        chain_names <- names(chain[[1]])
        df_chain <- data.frame(matrix(unlist(chain), ncol = length(chain_names), byrow = TRUE))
        colnames(df_chain) <- chain_names
    }

    if (return.traj)
    {
        return(list(acceptance.rate = accepted / n_iterations, trace = df_chain, trajectories = trajectories))
    } else
    {
        return(list(acceptance.rate = accepted / n_iterations, trace = df_chain))
    }
}
