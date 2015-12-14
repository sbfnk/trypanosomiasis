############################################################################
## Script for analysing chronic carriers of trypanosomiasis               ##
############################################################################

library('docopt')

"Script for analysing chronic carriers of trypanosomiasis
Usage: chronic_analysis.r [options]
Options:
-v --village=<village>                       number of the village
-n --nsamples=<nsamples>                     number of samples
-t --transmitted                             have the infection trasnmitted from cases
-b --background                              have a background force of infection
-c --chronic                                 let chronically infected cases infect
-l --lhs                                     lhs (as opposed to prior) sampling
-a --acceptance=<acceptance>                 require this many acceptances (for setting epsilon, default: 2)
-h --help                                    show this message" -> doc

opts <- docopt(doc)

if (opts[["help"]])
{
    print(opts)
    exit()
}

code_dir <- path.expand("~/code/trypanosomiasis/")
source(paste(code_dir, "chronic_carriers.R", sep = "/"))

num_samples <- as.integer(opts[["nsamples"]])

if (length(num_samples) == 0)
{
    num_samples <- 1
}

village <- as.integer(opts[["village"]])
if (length(village) == 0)
{
    village <- 1
}

require_acceptances <- as.integer(opts[["acceptance"]])
if (length(require_acceptances) == 0)
{
    require_acceptances <- 1
}

model_options <- list(transmitted = opts[["transmitted"]],
                      background = opts[["background"]],
                      chronic = opts[["chronic"]])

library('data.table')
library('trypR')

sample_options <- c(list(nsamples = num_samples,
                         calc.posterior = FALSE,
                         villages = village,
                         sample = ifelse(opts[["lhs"]], "lhs", "prior")
                         ),
                    model_options)

samples <- do.call(chronic_carriers_sample, sample_options)

dt <- data.table(t(sapply(samples, function(x)
{
    c(x[["parameters"]], unlist(x[["summary_statistics"]]))
})))

data_stage1_active <- village_screening[village.number == village, detected1_1]
data_stage1_passive <-
    village_cases[village.number == village & stage == 1, sum(cases)]
data_stage2_passive <-
    village_cases[village.number == village & stage == 2, sum(cases)]

dt[, active_stage1_diff := active_stage1 - data_stage1_active]
dt[, passive_stage1_diff := passive_stage1 - data_stage1_passive]
dt[, passive_stage2_diff := passive_stage2 - data_stage2_passive]

success <- FALSE

epsilon <- 0

while (!success)
{
    dt[, accept := (abs(active_stage1_diff) <= epsilon) & (abs(passive_stage1_diff) <= epsilon) & (abs(passive_stage2_diff) <= epsilon)]
    success <- (nrow(dt[accept == TRUE]) > 1)
    if (nrow(dt[accept == TRUE]) > (require_acceptances - 1)) 
    {
        success <- TRUE
    } else
    {
        epsilon <- epsilon + 1
    }
}

cat("First adaptation, epsilon =", epsilon, "\n")
## dt_new[, accept := (abs(active_stage1_diff) <= epsilon) & (abs(passive_stage1_diff) <= epsilon) & (abs(passive_stage2_diff) <= epsilon)]
dt[, id := 1:nrow(dt)]

accept_ids <- dt[accept == TRUE, id]

current_accepted <- data.table(id = 1:nrow(dt))
current_accepted[id >= min(accept_ids), accept_id := max(accept_ids[accept_ids <= id]), by = 1:nrow(current_accepted[id >= min(accept_ids)])]
prior_accepted <- dt[current_accepted[!is.na(accept_id), accept_id]]

prior_parameters <- copy(prior_accepted)
village_number <- which(names(prior_parameters) == "village.number")
prior_parameters <-
    prior_parameters[, -(village_number:ncol(prior_parameters)), with = FALSE]

## prior_sd <- apply(prior_parameters, 2, sd)
prior_sd <- c(0.05,0.1,0.1,0.001,0.001,10,5,0,0,0,0,0,0)
prior_zero <- prior_sd
prior_zero[] <- 0
prior_upper <- prior_sd


## success <- FALSE

## while (!success)
## {
##   mcmc_fixed <- chronic_carriers_mcmc(init = unlist(prior_parameters[floor(runif(1, 1, 1:nrow(prior_parameters)))]), n_iterations = 1000, sd = prior_zero, epsilon = epsilon, data_summary = c(active_stage1 = data_stage1_active, passive_stage1 = data_stage1_passive, passive_stage2 = data_stage2_passive), villages = village, verbose = TRUE)
##   if (mcmc_fixed$acceptance.rate > .2)
##   {
##       success <- TRUE
##   } else
##   {
##       epsilon <- epsilon + 1
##       cat("  increasing epsilon to", epsilon, "\n")
##   }
## }

## cat("Second adaption, epsilon =", epsilon, "\n")

cat(prior_sd, "\n")

mcmc_options <-
    c(list(init = unlist(prior_parameters[floor(runif(1, 1, 1:nrow(prior_parameters)))]),
           n_iterations = num_samples,
           sd = prior_sd,
           epsilon = epsilon,
           data_summary = c(active_stage1 = data_stage1_active,
                            passive_stage1 = data_stage1_passive,
                            passive_stage2 = data_stage2_passive),
           villages = village,
           verbose = TRUE,
           lower = prior_zero),
      model_options)
mcmc <- do.call(chronic_carriers_mcmc, mcmc_options)

df <- data.frame(matrix(unlist(mcmc$trace), ncol = ncol(prior_parameters), byrow = TRUE))
colnames(df) <- colnames(prior_parameters)

saveRDS(df, paste0(ifelse(opts[["transmitted"]], "tran_", ""),
                   ifelse(opts[["background"]], "back_", ""),
                   ifelse(opts[["chronic"]], "chro_", ""),
                   "village_", village, ".rds"))

## df_fixed <- data.frame(matrix(unlist(x), ncol = ncol(prior_parameters), byrow = TRUE))
## colnames(df_fixed) <- colnames(prior_parameters)
## mcmc_fixed <- mcmc(df_fixed)
## accRate <- 1 - min(rejectionRate(mcmc_fixed))
