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
-i --thin=<thin>                             thinning (1 for no thinning)
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

thin <- as.integer(opts[["thin"]])
if (length(thin) == 0)
{
    thin <- 1
}

require_acceptances <- as.integer(opts[["acceptance"]])
if (length(require_acceptances) == 0)
{
    require_acceptances <- 2
}

model_options <- list(transmitted = opts[["transmitted"]],
                      background = opts[["background"]],
                      chronic = opts[["chronic"]])

library('data.table')
library('trypR')

sample_options <- c(list(nsamples = num_samples,
                         calc.posterior = TRUE,
                         villages = village,
                         sample = ifelse(opts[["lhs"]], "lhs", "prior")
                         ),
                    model_options)

samples <- do.call(chronic_carriers_sample, sample_options)

dt <- data.table(t(sapply(samples, function(x) {unlist(x)})))

success <- FALSE

epsilon <- 0

while (!success)
{
  dt[, accept := (abs(active_stage1_diff) <= epsilon) & (abs(active_stage2_diff) <= epsilon) & (abs(passive_stage1_diff) <= epsilon) & (abs(passive_stage2_diff) <= epsilon)]
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
logprior_col <- which(names(prior_parameters) == "logprior")
prior_parameters <-
    prior_parameters[, -(logprior_col:ncol(prior_parameters)), with = FALSE]

## prior_sd <- apply(prior_parameters, 2, sd)
prior_sd <- c(pc = 0.1, alpha = 0.3, delta = ifelse(opts[["chronic"]], 0.3, 0), p1 = 0.01, p2 = 0.04, rc = 0, r1 = 0, r2 = 0, screen1 = 0, screen2 = 0, N = 0)



if (opts[["background"]])
{
    prior_sd["lambda"] = 1e-5
}
if (opts[["background"]])
{
    prior_sd["transmitted"] = 1e-5
}
prior_zero <- prior_sd
prior_zero[] <- 0
##prior_upper <- prior_sd


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

print(prior_sd)

mcmc_options <-
    c(list(start = unlist(prior_parameters[floor(runif(1, 1, 1:nrow(prior_parameters)))]),
           n_iterations = num_samples,
           sd = prior_sd / 5,
           epsilon = epsilon,
           data_summary = c(active_stage1 = data_stage1_active,
                            passive_stage1 = data_stage1_passive,
                            passive_stage2 = data_stage2_passive),
           village = village,
           verbose = TRUE,
           lower = prior_zero,
           thin = thin,
           return.traj = TRUE),
      model_options)
mcmc <- do.call(chronic_carriers_abc_mcmc, mcmc_options)

## df <- data.frame(matrix(unlist(mcmc$trace), ncol = ncol(prior_parameters), byrow = TRUE))
## colnames(df) <- names(mcmc$trace[[1]])

saveRDS(mcmc, paste0(ifelse(opts[["transmitted"]], "tran_", ""),
                     ifelse(opts[["background"]], "back_", ""),
                     ifelse(opts[["chronic"]], "chro_", ""),
                     "village_", village, ".rds"))

## df_fixed <- data.frame(matrix(unlist(x), ncol = ncol(prior_parameters), byrow = TRUE))
## colnames(df_fixed) <- colnames(prior_parameters)
## mcmc_fixed <- mcmc(df_fixed)
## accRate <- 1 - min(rejectionRate(mcmc_fixed))
