library('data.table')
library('adaptivetau')
library('lhs')
library('reshape2')

sim_chronic_transitions <- list(
    c(S = -1, Ic = +1), # chronic infection
    c(Ic = -1, S = +1), # chronic recovery
    c(S = -1, I1 = +1), # pathogenic infection
    c(I1 = -1, I2 = +1), # progression into stage 2
    c(I1 = -1, S = +1, Z1pass = +1), # passive stage 1 detection
    c(I2 = -1, S = +1, Z2pass = +1), # passive stage 2 detection
    c(I2 = -1, S = +1) # stage 2 death
    )

sim_chronic_rates <- function(state, theta, t)
{
    # parameters
    pc <- theta[["pc"]]
    lambda <- theta[["lambda"]]
    rc <- theta[["rc"]]
    r1 <- theta[["r1"]]
    p1 <- theta[["p1"]]
    p2 <- theta[["p2"]]
    d <- theta[["d"]]

    S <- state[["S"]]
    Ic <- state[["Ic"]]
    I1 <- state[["I1"]]
    I2 <- state[["I2"]]
    Z1pass <- state[["Z1pass"]]
    Z2pass <- state[["Z2pass"]]

    return(c(
        pc * lambda * S, # chronic infection
        rc * Ic, # chronic recovery
        (1 - pc) * lambda * S, # pathogenic infection
        r1 * I1, # progression into stage 2
        p1 * I1, # passive stage 1 detection
        p2 * I2, # passitve stage 2 detection
        d * I2 # stage 2 death
        ))
}

rprior <- function(villages = 1)
{
    param_vector <- c(pc = runif(1, 0, 1), alpha = runif(1, 0, 1))
    if (villages == 1)
    {
        param_vector <- c(param_vector,
                          lambda = runif(1, 0, 0.05),
                          p1 = runif(1, 0, 1),
                          p2 = runif(1, 0, 2))
    } else
    {
        village_vector <- c(runif(villages, 0, 0.05),
                            runif(villages, 0, 1),
                            runif(villages, 0, 2))
        names(village_vector) <-
            paste(rep(c("lambda", "p1", "p2"), each = villages),
                  seq_len(villages), sep = ".")
        param_vector <- c(param_vector, village_vector)

    }
    param_vector <- c(param_vector,
                      rc = 1/120,
                      r1 = 0.0019 * 30.42,
                      d = 0.0040 * 30.42,
                      ## screen1 = runif(1, 0.86, 0.98),
                      screen1 = 0.95,
                      screen2 = 0.99,
                      N = villages)
    return(param_vector)
}

dprior <- function(theta, log = FALSE)
{
    N <- round(theta[["N"]])
    if (log)
    {
        res <- 0
        res <- res + dunif(theta[["pc"]], 0, 1, TRUE)
        res <- res + dnorm(theta[["rc"]], 1/120, log = TRUE)
        res <- res + dnorm(theta[["r1"]],  0.0019 * 30.42, log = TRUE)
        res <- res + dnorm(theta[["d"]], 0.0040 * 30.42, log = TRUE)
        res <- res + dunif(theta[["alpha"]], 0, 1, TRUE)
        res <- res + dunif(theta[["screen1"]], 0.86, 0.98, TRUE)
        res <- res + dnorm(theta[["screen2"]], 0.99, log = TRUE)
        if (N == 1)
        {
            res <- res + dunif(theta[["lambda"]], 0, 0.05, TRUE)
            res <- res + dunif(theta[["p1"]], 0, 1, TRUE)
            res <- res + dunif(theta[["p2"]], 0, 2, TRUE)
        } else
        {
            for (i in seq_len(N))
            {
                res <-
                    res + dunif(theta[[paste("lambda", i, sep = ".")]], 0, 0.05, TRUE)
                res <- res + dunif(theta[[paste("p1", i, sep = ".")]], 0, 1, TRUE)
                res <- res + dunif(theta[[paste("p2", i, sep = ".")]], 0, 2, TRUE)
            }
        }
    } else
    {
        res <- 1
        res <- res * dunif(theta[["pc"]], 0, 1, FALSE)
        res <- res * dnorm(theta[["rc"]], 1/120, log = FALSE)
        res <- res * dnorm(theta[["r1"]], 0.0019 * 30.42, log = FALSE)
        res <- res * dnorm(theta[["d"]], 0.0040 * 30.42, log = FALSE)
        res <- res * dunif(theta[["alpha"]], 0, 1, FALSE)
        res <- res * dunif(theta[["screen1"]], 0.86, 0.98, FALSE)
        res <- res * dnorm(theta[["screen2"]], 0.99, log = FALSE)
        if (N == 1)
        {
            res <- res * dunif(theta[["lambda"]], 0, 0.05, FALSE)
            res <- res * dunif(theta[["p1"]], 0, 1, FALSE)
            res <- res * dunif(theta[["p2"]], 0, 2, FALSE)
        } else
        {
            for (i in seq_len(N))
            {
                res <-
                    res * dunif(theta[[paste("lambda", i, sep = ".")]], 0, 0.05, FALSE)
                res <- res * dunif(theta[[paste("p1", i, sep = ".")]], 0, 1, FALSE)
                res <- res * dunif(theta[[paste("p2", i, sep = ".")]], 0, 2, FALSE)
            }
        }
    }

    return(res)

}

rinit <- function(theta, village_data, rand = NULL, equilibrium = TRUE)
{

    N <- village_data[["N"]]

    if (equilibrium)
    {
        Ic_weight <- theta[["pc"]] / theta[["rc"]]
        I1_weight <- (1 - theta[["pc"]]) / theta[["r1"]]
        I2_weight <- (1 - theta[["pc"]]) / theta[["d"]]

        sum_weights <- Ic_weight + I1_weight + I2_weight

        eq_I <- N * (1 - 1 / (1 + theta[["lambda"]] * sum_weights))

        initIc <- rpois(1, Ic_weight * eq_I / sum_weights)
        initI1 <- rpois(1, I1_weight * eq_I / sum_weights)
        initI2 <- rpois(1, I2_weight * eq_I / sum_weights)

        stage1_detected <- village_data[["detected1_1"]]
        stage2_detected <- village_data[["detected1_2"]]

        if (initIc + initI1 >= stage1_detected)
        {
            Ic_ind <- c(rep(1, initIc), rep(0, initI1))
            Ic_prob <- c(rep(theta[["alpha"]], initIc), rep(1, initI1))
            Ic_det <- sum(Ic_ind[sample(Ic_ind, stage1_detected, prob = Ic_prob)])
            I1_det <- stage1_detected - Ic_det
            initIc <- initIc - Ic_det
            initI1 <- initI1 - I1_det
        } else
        {
            initIc <- -1
            initI1 <- -1
        }
        if (initI2 >= stage2_detected)
        {
            initI2 <- initI2 - stage2_detected
        } else
        {
            initI2 <- -1
        }
    } else
    {
        Ic_weight <- theta[["pc"]] / theta[["rc"]]
        I1_weight <- (1 - theta[["pc"]]) / theta[["r1"]]

        proportion_chronic <- Ic_weight / (Ic_weight + I1_weight)

        attendance <- min(village_data[["sigma.start"]], 1)

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

    return(c(S = initS, Ic = initIc, I1 = initI1, I2 = initI2,
             Z1pass = 0, Z2pass = 0))
}

likelihood <- function(state, stage1_detected = NULL,
                       stage2_detected = NULL, theta, attendance = 1,
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
        if (!is.null(stage1_detected))
        {
            if ("Ic" %in% names(state))
            {
                stage1_ll <- sum(sapply(seq(0, stage1_detected), function(x)
                {
                    prob <-
                        dbinom(x, state[["Ic"]] + x,
                               theta[["alpha"]] * theta[["screen1"]] * attendance)
                    prob <- prob * dbinom(stage1_detected - x, state[["I1"]] + x,
                                          theta[["screen1"]] * attendance)
                    prob
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
                                    theta[["screen1"]] * attendance, log = log)
                if (log)
                {
                    res <- res + stage1_ll
                } else
                {
                    res <- res * stage1_ll
                }
            }
        }
        if (!is.null(stage2_detected))
        {
            stage2_ll <- dbinom(stage2_detected, state[["I2"]] + stage1_detected,
                                theta[["screen2"]] * attendance, log = log)
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

dinit <- function(init, village_data, theta, log = FALSE)
{
    stage1_detected <- village_data[["detected1_1"]]
    stage2_detected <- village_data[["detected1_2"]]

    attendance <- min(village_data[["sigma.start"]], 1)

    return(likelihood(init, stage1_detected, stage2_detected,
                      theta, attendance, log = log))
}

param_posterior <- function(theta, village_data, cum_data, nruns, log = FALSE, ...)
{
    res <- 0
    log.prior <- dprior(theta, log = TRUE)
    posteriors <- c()
    final.attendance <- min(village_data[["sigma.end"]], 1)
    if (is.finite(log.prior))
    {
        for (j in seq_len(nruns))
        {
            init <- rinit(theta, village_data, ...)
            log.init <- dinit(init, village_data, theta, TRUE)
            if (is.finite(log.init))
            {
                run <-
                    data.table(ssa.adaptivetau(init, sim_chronic_transitions,
                                               sim_chronic_rates, theta,
                                               village_data[["stoptime"]]))

                passive_state <- c(I1 = run[nrow(run), Z1pass],
                                   I2 = run[nrow(run), Z2pass])
                ll <- likelihood(passive_state, cum_data[1], cum_data[2], theta,
                                 log = TRUE)
                ll <- ll +
                    log(likelihood(run[nrow(run)],
                                   village_data[["detected1_2"]], theta = theta,
                                   attendance = final.attendance))
            } else {
                ll <- -Inf
            }
            posteriors[j] <- log.prior + log.init + ll
        }
        res <- sum(exp(posteriors)) / nruns
    } else
    {
        res <- 0
    }


    if (log)
    {
        return(log(res))
    } else
    {
        return(res)
    }
}

param_posterior_all <- function(theta, village_screening, cum_data, nruns,
                                log = FALSE, ...)
{
    res <- 0
    log.prior <- dprior(theta, log = TRUE)
    posteriors <- c()
    final.attendance <- min(village_data[["sigma.end"]], 1)
    if (is.finite(log.prior))
    {
        for (j in seq_len(nruns))
        {
            init <- rinit(theta, village_data, ...)
            log.init <- dinit(init, village_data, theta, TRUE)
            if (is.finite(log.init))
            {
                run <-
                    data.table(ssa.adaptivetau(init, sim_chronic_transitions,
                                               sim_chronic_rates, theta,
                                               village_data[["stoptime"]]))

                passive_state <- c(I1 = run[nrow(run), Z1pass],
                                   I2 = run[nrow(run), Z2pass])
                ll <- likelihood(passive_state, cum_data[1], cum_data[2], theta,
                                 log = TRUE)
                ll <- ll +
                    log(likelihood(run[nrow(run)],
                                   village_data[["detected1_2"]], theta = theta,
                                   attendance = final.attendance))
            } else {
                ll <- -Inf
            }
            posteriors[j] <- log.prior + log.init + ll
        }
        res <- sum(exp(posteriors)) / nruns
    } else
    {
        res <- 0
    }


    if (log)
    {
        return(log(res))
    } else
    {
        return(res)
    }
}
