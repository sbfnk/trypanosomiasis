##' Produe a single simuation run
##'
##' @param theta parameter vector
##' @return data frame with trajctor
##' @author seb
##' @export
sim_trajectory <- function(theta, init, village)
{
    stoptime <- village_screening[village_screening$village.number == village]$stoptime
    if (missing(init)) init <- rinit(theta)

    chronic_options <- list(params = theta, init = init, times = seq(0, stoptime))
    for (stage in 1:2)
    {
        if (!paste0("p", stage) %in% names(theta))
        {
            chronic_options[[paste0("stage", stage, "_passive")]] <-
                passive[[stage]]
        }
    }

    return(do.call(chronic, chronic_options))
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
        c(pc = runif(1, 0, 1), alpha = runif(1, 0, 1))
    param_vector["delta"] <- ifelse(chronic, runif(1, 0, 1), 0)

    if (length(villages) == 1)
    {
        param_vector <- c(param_vector,
                          lambda = ifelse(background,
                                          10^(runif(1, -5, -2)),
                                          0),
                          beta = ifelse(transmitted,
                                        10^(runif(1, -5, -2)),
                                        0),
                          p1 = ifelse(passive, runif(1, 0, 52), 0),
                          p2 = ifelse(passive, runif(1, 0, 52), 0))
    } else
    {
        village_vector <- c(ifelse(background,
                                   10^(runif(length(villages), -5,  -2)),
                                   rep(0, length(villages))),
                            ifelse(transmitted,
                                   10^(runif(length(villages), -5, -2)),
                                   rep(0, length(villages))),
                            ifelse(passive,
                                   runif(length(villages), 0, 52),
                                   rep(0, length(villages))),
                            ifelse(passive,
                                   runif(length(villages), 0, 52),
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
    param_prior <- c(param_prior, dunif(theta[["pc"]], 0, 1, log = TRUE))
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
                             dunif(theta[["lambda"]], 10^(-5), 10^(-2), log = TRUE))
        }
        if ("beta" %in% names(theta))
        {
            param_prior <- c(param_prior,
                             dunif(theta[["lambda"]], 10^(-5), 10^(-2), log = TRUE))
        }
        if ("p1" %in% names(theta))
        {
            param_prior <- c(param_prior, dunif(theta[["p1"]], 0, 52, log = TRUE))
        }
        if ("p2" %in% names(theta))
        {
            param_prior <- c(param_prior, dunif(theta[["p2"]], 0, 52, log = TRUE))
        }
    } else
    {
        lambda_names <-
            param_prior <-
                c(param_prior,
                  unlist(sapply(grep("^lambda\\.", names(theta), value = TRUE),
                                function(x)
                  {
                      dunif(theta[[x]], 0, 0.05, log = TRUE)
                  })))
        param_prior <-
            c(param_prior,
              unlist(sapply(grep("^p1\\.", names(theta), value = TRUE),
                            function(x)
              {
                  dunif(theta[[x]], 0, 1, log = TRUE)
              })))
        param_prior <-
            c(param_prior,
              unlist(sapply(grep("^p2\\.", names(theta), value = TRUE),
                            function(x)
              {
                  dunif(theta[[x]], 0, 2, log = TRUE)
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
##' @param rand A list of combinations for how infections are distributed between chronic and symptomatic infections, with corresponding probabilities
##' @param equilibrium whether the initial conditions should be derived from equilibrium states
##' @param village_number village number (if parameters for multiple villages are given
##' @param passive whether to accumulate infections for passive detection (Z1pass, Z2pass)
##' @return initial conditions vector
##' @author Sebastian Funk
rinit <- function(theta, rand = NULL, equilibrium = TRUE, village_number = 1)
{

    data(village_data)

    N <- village_screening[village.number == village_number, N]

    village_lambda <- paste("lambda", village_number, sep = ".")
    if (!(village_lambda %in% names(theta)))
    {
        village_lambda <- "lambda"
    }

    if (equilibrium)
    {
        Ic_weight <- theta[["pc"]] / theta[["rc"]]
        I1_weight <- (1 - theta[["pc"]]) / theta[["r1"]]
        I2_weight <- (1 - theta[["pc"]]) / theta[["r2"]]

        sum_weights <- Ic_weight + I1_weight + I2_weight

        eq_I <- N * (1 - 1 / (1 + theta[[village_lambda]] * sum_weights))

        stage1_detected <- village_screening[village.number == village_number,
                                             detected1_1]
        stage2_detected <- village_screening[village.number == village_number,
                                             detected1_2]

        found.OK <- FALSE

        while(!found.OK) {
            initIc <- rpois(1, Ic_weight * eq_I / sum_weights)
            initI1 <- rpois(1, I1_weight * eq_I / sum_weights)
            initI2 <- rpois(1, I2_weight * eq_I / sum_weights)

            found.OK <- (initIc + initI1 >= stage1_detected) &&
                (initI2 >= stage2_detected)
        }
        Ic_ind <- c(rep(1, initIc), rep(0, initI1))
        Ic_prob <- c(rep(theta[["alpha"]], initIc), rep(1, initI1))
        Ic_det <- sum(Ic_ind[sample(Ic_ind, stage1_detected, prob = Ic_prob)])
        I1_det <- stage1_detected - Ic_det
        initIc <- initIc - Ic_det
        initI1 <- initI1 - I1_det
        initI2 <- initI2 - stage2_detected
    } else
    {
        Ic_weight <- theta[["pc"]] / theta[["rc"]]
        I1_weight <- (1 - theta[["pc"]]) / theta[["r1"]]

        proportion_chronic <- Ic_weight / (Ic_weight + I1_weight)

        attendance <- min(village_screening[village.number == village_number,
                                            sigma.start], 1)

        initI2 <-
            stage2_detected + rnbinom(1, stage2_detected,
                                      attendance * theta[["screen2"]])
        ## sample Ic and I1

        max.I <- N - stage2_detected
        if (is.null(rand))
        {
            rand <- list(dist = list(), prob = c())
            k <- 0
            for (i in seq_len(stage1_detected))
            {
                for (j in seq(stage1_detected, max.I))
                {
                    k <- k + 1
                    rand$prob[k] <-
                        dbinom(i, round(proportion_chronic * j),
                               attendance * theta[["alpha"]] * theta[["screen1"]])
                    rand$prob[k] <- rand$prob[k] *
                        dbinom(stage1_detected - i,
                               round((1 - proportion_chronic) * j),
                               attendance * theta[["screen1"]])
                    rand$dist[[k]] <- c(i = i, j = j)
                }
            }
        }
        pick <- sample(seq_len(k), 1, prob = rand$prob)

        initIc <-
            round(proportion_chronic * rand$dist[[pick]][["j"]]) -
            rand$dist[[pick]][["i"]]
        initI1 <-
            round((1 - proportion_chronic) * rand$dist[[pick]][["j"]]) -
            (stage1_detected - rand$dist[[pick]][["i"]])
    }

    initS <- N - initIc - initI1 - initI2

    initVec <- c(S = initS, Ic = initIc, I1 = initI1, I2 = initI2)

    return(initVec)
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

    if (!any(state < 0))
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
        if (is.finite(log.init))
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

##' Calculate the posterior density for a given set of parameters for the trypanosomiasis model, using a model that encompasses all villages
##'
##' @param theta parameter vector
##' @param villages villages to simuate
##' @param log whether to return the logarithm
##' @param ... further options for traj_likelihood
##' @return posterior density
##' @author Sebastian Funk
param_posterior_villages <- function(theta, villages, log = TRUE, ...)
{
    data(village_data)

    if (missing(villages))
    {
        sim_village_screening <- village_screening
    } else
    {
        sim_village_screening <- village_screening[village.number %in% villages]
    }

    res <- 0
    log.prior <- dprior(theta, log = TRUE)
    posteriors <- c()
    if (is.finite(log.prior))
    {
        village.posteriors <- apply(sim_village_screening, 1, function(village)
        {
            village_number <- village[["village.number"]]
            stoptime <- village[["stoptime"]]

            final.attendance <- min(village[["sigma.end"]], 1)

            posterior <- traj_likelihood(theta, village_number, log = TRUE, ...)
        })
        res <- village.posteriors
    } else
    {
        res <- rep(0, nrow(sim_village_screening))
    }

    if (log)
    {
        return(res)
    } else
    {
        return(exp(res))
    }
}

##' Calculate summary statistic for a given set of parameters for the trypanosomiasis model, using a model that encompasses all villages
##'
##' @param theta parameter vector
##' @param nruns number of runs
##' @param villages villages to simuate
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
            sim_trajectory(theta, village = village_number)
        })


        res[["nneg"]] <- sum(sapply(sims, function(x)
        {
            any(x[["I1"]] < 0 | x[["I2"]] < 0)
        }))

        res[["active_stage1"]] <- sapply(sims, function(x)
        {
            rbinom(1, x[["Ic"]],
                   theta[["alpha"]] * theta[["screen1"]] * final.attendance) +
                rbinom(1, x[["I1"]],
                       theta[["screen1"]] * final.attendance)
        })

        for (stage in 1:2)
        {
            if (grep(paste0("^p", stage), names(theta)))
            {
                passive_stage <- paste0("passive_stage", stage)
                res[[passive_stage]] <- sapply(sims, function(x)
                {
                    tail(x[[paste0("Z", stage, "pass")]], 1)
                })
            }
        }
        return(res)
    })
    if (length(summary.statistics) == 1)
    {
        summary.statistics <- unlist(summary.statistics)
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

    if (progress.bar) pb <- txtProgressBar(min = 0, max = nsamples - 1, style = 3)

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
                                                            ...))
        }
        if (progress.bar) setTxtProgressBar(pb, i)
    }
    if (progress.bar) close(pb)

    return(samples)
}

##' run MCMC on the model for chronic carriers of trypanosomiasis
##'
##' @param init initial state of the chain
##' @param n_iterations number of iterations
##' @param sd standard deviation of each parameter
##' @param upper upper limit on parameters (if given)
##' @param lower lower limit on parameters (if given)
##' @param epsilon value of epsilon (if ABC is run)
##' @param data_summary summary statistic of the data (if ABC is run)
##' @param verbose whether to print verbose output
##' @param ... parameters to be passed to \code{param_posterior_villages}
##' @return chain
##' @importFrom truncnorm rtruncnorm dtruncnorm
##' @author Sebastian Funk
chronic_carriers_mcmc <- function(init, n_iterations, sd,
                                  upper, lower,
                                  epsilon, data_summary,
                                  verbose = FALSE, ...)
{
    theta <- init
    prior_theta <- dprior(theta, log = TRUE)

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

    for (i in seq_len(n_iterations))
    {
        is.accepted <- FALSE
        theta_propose <- rtruncnorm(length(theta), lower, upper, theta, sd)
        names(theta_propose) <- names(theta)

        prior_propose <- dprior(theta_propose, log = TRUE)

        if (is.finite(prior_propose)) {
            hastings_ratio <-
                ifelse(any(sd > 0),
                       log(dtruncnorm(theta[sd > 0], lower, upper,
                                      theta_propose[sd > 0], sd[sd > 0])) -
                       log(dtruncnorm(theta_propose[sd > 0], lower, upper,
                                      theta[sd > 0], sd[sd > 0])),
                       0)

            if (missing(epsilon))
            {
                posterior_propose <-
                    param_posterior_villages(theta_propose, ...)
                if (is.finite(posterior_propose))
                {
                    log.acceptance <- posterior_propose - posterior + hastings_ratio
                    is.accepted <- (log(runif (1)) < log.acceptance)
                }
            } else
            {
                sumstats <-
                    param_sumstat_villages(theta_propose, ...)
                diff <- sumstats[names(data_summary)] - data_summary
                pass <- all(diff <= epsilon)
                if (pass)
                {
                    log.acceptance <- prior_propose - prior_theta + hastings_ratio
                    is.accepted <- (log(runif (1)) < log.acceptance)
                }
            }
        }
        if (is.accepted)
        {
            accepted <- accepted + 1
            theta <- theta_propose
            prior_theta <- prior_propose
            if (missing(epsilon)) posterior <- posterior_propose
        }
        chain[[i]] <- theta
        if (verbose) cat(i, "acc:", is.accepted, accepted / i, "\n")
    }

    return(list(acceptance.rate = accepted / n_iterations, trace = chain))
}
