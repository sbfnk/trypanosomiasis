library('compiler')
library('data.table')

args <- commandArgs(trailingOnly = TRUE)
village_no <- as.integer(args[1])

enableJIT(1)

source(path.expand("~/code/trypanosomiasis/chronic_carriers.R"))
village_data <- readRDS("village_data.rds")

village_cases <- village_data[["cases"]]
village_screening <- village_data[["screening"]]

theta <- rprior()

cum_data <-
    village_cases[village.number == village_no,
                  list(cases = sum(cases)), by = stage]$cases

res <- readRDS("cc_lhs_samples_eq.rds")
## samples <- readRDS(path.expand(paste("cc_mcmc_", village_no, ".rds", sep = "")))

samples <- res[["samples"]]
posteriors <- res[["posteriors"]]
nruns <- 100

dt_samples <- data.table(matrix(unlist(samples), ncol = 10, byrow = TRUE))
n_top10p <- round(sum(is.finite(posteriors))/10)
top10p <- order(posteriors, decreasing = TRUE)[1:n_top10p]

dt_top10p <- dt_samples[top10p]
setnames(dt_top10p, 1:ncol(dt_top10p), names(theta))
dt_top10p[, index := seq_len(nrow(dt_top10p))]
mdtt <- melt(dt_top10p, id.vars = "index")

start.theta <- samples[[top10p[1]]]
## max.index <- max(samples[, index])
 ##start.theta <- samples[index == max.index, value]
## names(start.theta) <- samples[index == max.index, variable]
chain <- list()
n_iterations <- 100000

theta <- start.theta
posterior <- param_posterior(theta, village_screening[village.number == village_no],
                             cum_data, nruns, TRUE)
#sdvec <- c(0.05, 0.001, 0, 0, 0.05, 0.5, 0, 0.1, 0, 0)
sdvec <- c(0.1, 0.0001, 0, 0, 0.04, 0.25, 0, 0.25, 0, 0)
accepted <- 0

for (i in seq_len(n_iterations))
{
    is.accepted <- FALSE
    theta_propose <- rnorm(length(theta), theta, sdvec / 2)
    names(theta_propose) <- names(theta)

    posterior_propose <-
        param_posterior(theta_propose, village_screening[village.number == village_no],
                        cum_data, nruns, TRUE)
    if (is.finite(posterior_propose))
    {
        log.acceptance <- posterior_propose - posterior
        is.accepted <- (log(runif (1)) < log.acceptance)
        if (is.accepted)
        {
            accepted <- accepted + 1
            theta <- theta_propose
            posterior <- posterior_propose
        }
    }
    chain[[i]] <- list(theta = theta, posterior = posterior)
    cat(i, "acc:", is.accepted, accepted / i, posterior, "\n")
}

dt_mcmc <- data.table(t(sapply(chain, function(x) {x[["theta"]]})))
dt_mcmc[, index := 1:nrow(dt_mcmc)]
mdtm <- melt(dt_mcmc, id.vars = "index")

saveRDS(mdtm, paste("cc_mcmc2_", village_no, ".rds", sep = ""))

## ggplot(mdtm, aes(x = value)) + geom_histogram() + facet_wrap(~variable)
## ggplot(mdtm, aes(x = index, y = value)) + geom_line() + facet_wrap(~variable, scales = "free")

