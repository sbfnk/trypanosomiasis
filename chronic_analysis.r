############################################################################
## Script for analysing chronic carriers of trypanosomiasis               ##
############################################################################

library('docopt')

"Script for analysing chronic carriers of trypanosomiasis
Usage: chronic_analysis.r [options]
Options:
-v --village=<village>                       number of the village
-e --epsilon=<epislon>                       epsilon
-n --nsamples=<nsamples>                     number of samples
-h --help                                    show this message" -> doc

opts <- docopt(doc)

if (opts[["help"]])
{
    print(opts)
    exit()
}

epsilon <- as.numeric(opts[["epsilon"]])
num_samples <- as.integer(opts[["nsamples"]])

if (length(epsilon) == 0)
{
    epsilon <- 0
}

if (length(num_samples) == 0)
{
    num_samples <- 1
}

village <- as.integer(opts[["village"]])
if (length(village) == 0)
{
    village <- 1
}

library('data.table')
library('trypR')
library('coda')

## lhs <- chronic_carriers_sample(nruns = 1, nsamples = 100, calc.posterior = FALSE,
##                                    villages = village, passive = TRUE, sample = "lhs")
prior <- chronic_carriers_sample(nsamples = 10000, calc.posterior = FALSE,
                                 villages = village, passive = TRUE, sample = "prior")

## lhs_new <- chronic_carriers_sample(nruns = 1, nsamples = 10000, calc.posterior = FALSE,
##                                    villages = village, passive = TRUE, sample = "lhs",
##                                    transmitted = TRUE)
## prior_new <- chronic_carriers_sample(nruns = 1, nsamples = 10000,
##                                      calc.posterior = FALSE,
##                                      villages = village, passive = TRUE,
##                                      sample = "prior", transmitted = TRUE)

## dt_lhs <- data.table(t(sapply(lhs, function(x)
## {
##     c(x[["parameters"]], unlist(x[["summary_statistics"]]))
## })))
dt <- data.table(t(sapply(prior, function(x)
{
    c(x[["parameters"]], unlist(x[["summary_statistics"]]))
})))
## dt_lhs_new <- data.table(t(sapply(lhs_new, function(x)
## {
##     c(x[["parameters"]], unlist(x[["summary_statistics"]]))
## })))
## dt_prior_new <- data.table(t(sapply(prior_new, function(x)
## {
##     c(x[["parameters"]], unlist(x[["summary_statistics"]]))
## })))

##dt_lhs[, sample := "lhs"]
## dt_prior[, sample := "prior"]

## dt <- rbind(dt_lhs, dt_prior)

## dt_lhs_new[, sample := "lhs"]
## dt_prior_new[, sample := "prior"]
## dt_new <- rbind(dt_lhs_new, dt_prior_new)

data_stage1_active <- village_screening[village.number == village, detected1_1]
data_stage1_passive <-
    village_cases[village.number == village & stage == 1, sum(cases)]
data_stage2_passive <-
    village_cases[village.number == village & stage == 2, sum(cases)]

dt[, active_stage1_diff := active_stage1 - data_stage1_active]
dt[, passive_stage1_diff := passive_stage1 - data_stage1_passive]
dt[, passive_stage2_diff := passive_stage2 - data_stage2_passive]
## dt_new[, active_stage1_diff := active_stage1 - data_stage1_active]
## dt_new[, passive_stage1_diff := passive_stage1 - data_stage1_passive]
## dt_new[, passive_stage2_diff := passive_stage2 - data_stage2_passive]

## epsilon <- 0

dt[, accept := (abs(active_stage1_diff) <= epsilon) & (abs(passive_stage1_diff) <= epsilon) & (abs(passive_stage2_diff) <= epsilon)]
## dt_new[, accept := (abs(active_stage1_diff) <= epsilon) & (abs(passive_stage1_diff) <= epsilon) & (abs(passive_stage2_diff) <= epsilon)]

prior_samples[, accept := (abs(active_stage1_diff) <= epsilon) & (abs(passive_stage1_diff) <= epsilon) & (abs(passive_stage2_diff) <= epsilon)]
prior_samples[, id := 1:nrow(prior_samples)]
accept_ids <- prior_samples[accept == TRUE, id]
current_accepted <- data.table(id = 1:nrow(prior_samples))
current_accepted[id >= min(accept_ids), accept_id := max(accept_ids[accept_ids <= id]), by = 1:nrow(current_accepted[id >= min(accept_ids)])]
prior_accepted <- prior_samples[current_accepted[!is.na(accept_id), accept_id]]

prior_parameters <- prior_accepted[, list(pc = pc, alpha = alpha, delta = delta, lambda = lambda, p1 = p1, p2 = p2, rc = rc, r1 = r1, r2 = r2, screen1 = screen1, screen2 = screen2, N = N)]
prior_sd <- apply(prior_parameters, 2, sd)
prior_zero <- prior_sd
prior_zero[] <- 0
prior_upper <- prior_sd


## only needed once, test for epsilon

## chronic_carriers_mcmc(init = unlist(prior_parameters[floor(runif(1, 1, 1:nrow(prior_parameters)))]), n_iterations = 100, sd = prior_zero, epsilon = 1, data_summary = c(active_stage1 = data_stage1_active, passive_stage1 = data_stage1_passive, passive_stage2 = data_stage2_passive), villages = village, verbose = TRUE)

x <- chronic_carriers_mcmc(init = unlist(prior_parameters[floor(runif(1, 1, 1:nrow(prior_parameters)))]), n_iterations = 100000, sd = prior_sd / 2, epsilon = 1, data_summary = c(active_stage1 = data_stage1_active, passive_stage1 = data_stage1_passive, passive_stage2 = data_stage2_passive), villages = village, verbose = TRUE, lower = prior_zero)

df <- data.frame(matrix(unlist(x), ncol = ncol(prior_parameters), byrow = TRUE))
colnames(df) <- colnames(prior_parameters)

saveRDS(df, paste0("village_", village, ".rds"))

## df_fixed <- data.frame(matrix(unlist(x), ncol = ncol(prior_parameters), byrow = TRUE))
## colnames(df_fixed) <- colnames(prior_parameters)
## mcmc_fixed <- mcmc(df_fixed)
## accRate <- 1 - min(rejectionRate(mcmc_fixed))
