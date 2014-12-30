library('data.table')
library('XLConnect')
library('adaptivetau')
library('compiler')
library('lhs')

enableJIT(1)

wb <- loadWorkbook("~/Research/Analysis/Sleeping_Sickness/Chronic/Village input data.xlsx")
village_screening <- data.table(readWorksheet(wb, "Sheet1"))

village_dirs <- list.files(path = path.expand("~/Research/Analysis/Sleeping_Sickness/Chronic/passive_screening_datasets/"))
village_ids <- as.integer(gsub("^village ", "", village_dirs))
village_ids <- village_ids[order(village_ids)]

village_cases <- data.table(village.number = integer(0), month = integer(0), cases = integer(0), stage = integer(0))
for (id in village_ids)
{
    for (stage in c(1, 2))
    {
        temp <- data.table(read.csv(paste("~/Research/Analysis/Sleeping_Sickness/Chronic/passive_screening_datasets/village ",
                                          id, "/ps", stage, "_", id, ".csv", sep = ""), header = F))
        setnames(temp, 1:2, c("month", "cases"))
        village_cases <- rbind(village_cases, data.table(village.number = rep(id, nrow(temp)),
                                                         month = temp[, month],
                                                         cases = temp[, cases],
                                                         stage = rep(stage, nrow(temp))))
    }
}

sim_chronic_transitions <- list(
    c(S = -1, Ic = +1), # chronic infection
    c(Ic = -1, S = +1), # chronic recovery
    c(S = -1, I1 = +1), # pathogenic infection
    c(I1 = -1, I2 = +1), # progression into stage 2
    c(I1 = -1, S = +1, Z1pass = +1), # passive stage 1 detection
    c(I2 = -1, S = +1, Z2pass = +2), # passive stage 2 detection
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

# village number 2

rprior <- function(vector)
{
    return(c(
        pc = runif(1, 0, 0.5),
        lambda = runif(1, 0, 0.05),
        rc = 1/120,
        r1 = 0.0019 * 30.42,
        p1 = runif(1, 0, 1),
        p2 = runif(1, 0, 1),
        d = 0.0040 * 30.42,
        alpha = runif(1, 0, 1),
        screen1 = runif(1, 0.86, 0.98),
        screen2 = 0.99
        ))
}

dprior <- function(theta, log = FALSE)
{
    if (log)
    {
        res <- 0
        res <- res + dunif(theta[["pc"]], 0, 0.5, TRUE)
        res <- res + dunif(theta[["lambda"]], 0, 0.05, TRUE)
        res <- res + dnorm(theta[["rc"]], 1/120, log = TRUE)
        res <- res + dnorm(theta[["r1"]],  0.0019 * 30.42, log = TRUE)
        res <- res + dunif(theta[["p1"]], 0, 1, TRUE)
        res <- res + dunif(theta[["p2"]], 0, 1, TRUE)
        res <- res + dnorm(theta[["d"]], 0.0040 * 30.42, log = TRUE)
        res <- res + dunif(theta[["alpha"]], 0, 1, TRUE)
        res <- res + dunif(theta[["screen1"]], 0.86, 0.98, TRUE)
        res <- res + dnorm(theta[["screen2"]], 0.99, log = TRUE)
    } else
    {
        res <- 1
        res <- res * dunif(theta[["pc"]], 0, 0.5, FALSE)
        res <- res * dunif(theta[["lambda"]], 0, 0.05, FALSE)
        res <- res * dnorm(theta[["rc"]], 1/120, log = FALSE)
        res <- res * dnorm(theta[["r1"]], 0.0019 * 30.42, log = FALSE)
        res <- res * dunif(theta[["p1"]], 0, 1, FALSE)
        res <- res * dunif(theta[["p2"]], 0, 1, FALSE)
        res <- res * dnorm(theta[["d"]], 0.0040 * 30.42, log = FALSE)
        res <- res * dunif(theta[["alpha"]], 0, 1, FALSE)
        res <- res * dunif(theta[["screen1"]], 0.86, 0.98, FALSE)
        res <- res * dnorm(theta[["screen2"]], 0.99, log = FALSE)
    }

    return(res)

}

rinit <- function(theta, village_data)
{

    N <- village_data[["N"]]

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

    if (initIc + initI1 > stage1_detected)
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

param_posterior <- function(theta, village_data, cum_data, nruns, log = FALSE)
{
    res <- 0
    log.prior <- dprior(theta)
    posteriors <- c()
    final.attendance <- min(village_data[["sigma.end"]], 1)
    for (j in seq_len(nruns))
    {
        init <- rinit(theta, village_data)
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
            ll <- ll + log(likelihood(run[nrow(run)],
                                      village_data[["detected1_2"]], theta = theta,
                                      attendance = final.attendance))
        } else {
            ll <- -Inf
        }
        posteriors[j] <- log.prior + log.init + ll
    }

    res <- sum(exp(posteriors)) / nruns

    if (log)
    {
        return(log(res))
    } else
    {
        return(res)
    }
}

samples <- list()
inits <- list()
distances <- c()
likelihoods <- c()
posteriors <- c()

set.seed(42)
n_samples <- 10000
n_runs <- 100

r <- randomLHS(n_samples, 6)
upper <- c(0.5, 0.05, 1, 1, 1, 0.98)
lower <- c(0, 0, 0, 0, 0, 0.86)
r <- t(apply(r, 1, function(x) { lower + (upper - lower) * x}))
colnames(r) <- c("pc", "lambda", "p1", "p2", "alpha", "screen1")
theta <- rprior()

cum_data <-
    village_cases[village.number == 2, list(cases = sum(cases)), by = stage]$cases

i <- 0
while(i < n_samples)
{
    theta[colnames(r)] <- r[i + 1, ]
    i <- i + 1
    samples[[i]] <- theta
    posteriors[i] <- param_posterior(theta, village_screening[village.number == 2],
                                     cum_data, n_runs, TRUE)
    cat(i, posteriors[i], "\n")
}

saveRDS(list(samples = samples, posteriors = posteriors), "cc_lhs_samples.rds")
top.10 <- order(posteriors, decreasing = TRUE)[1:10]

start.theta <- samples[[top.10[1]]]
chain <- list()
n_iterations <- 100000

theta <- start.theta
posterior <- posteriors[top.10[1]]
sdvec <- c(0.01, 0.001, 0, 0, 0.01, 0.01, 0, 0.01, 0.001, 0)
accepted <- 0

for (i in seq_len(n_iterations))
{
    theta_propose <- rnorm(length(theta), theta, sdvec)
    names(theta_propose) <- names(theta)
    log_prior <- dprior(theta_propose, log = TRUE)
    if (is.finite(log_prior))
    {
        init <- rinit(theta_propose, village_screening[village.number == 2])
        log_init <- dinit(init, village_screening[village.number == 2],
                          theta, TRUE)
        if (is.finite(log_init))
        {
            res <-
                data.table(ssa.adaptivetau(init, sim_chronic_transitions,
                                           sim_chronic_rates, theta,
                                           village_screening[village.number == 2,
                                                             stoptime]))

            sum_vec <-
                unname(unlist(res[nrow(res),
                                  list(sum(Z1pass), sum(Z2pass))]))

            ll <- sum(ifelse(any(sum_vec[1:2] - data_vec[1:2] < 0), -Inf, 0)) +
                log(stage1_likelihood(res[nrow(res)], data_vec[3], theta,
                                      final_attendance))
            if (is.finite(ll))
            {
                posterior_propose <- ll + log_prior + log_init
                log.acceptance <- posterior_propose - posterior
                is.accepted <- (log(runif (1)) < log.acceptance)
                if (is.accepted)
                {
                    accepted <- accepted + 1
                    theta <- theta_propose
                    posterior <- posterior_propose
                }
            }
        }
    }
    chain[[i]] <- list(theta = theta, posterior = posterior)
    cat(i, "acc:",  accepted / i, log_prior, log_init, ll, posterior_propose, "\n")
}
