############################################################################
## Script for analysing chronic carriers of trypanosomiasis               ##
############################################################################

library('docopt')

"Script for analysing chronic carriers of trypanosomiasis
Usage: chronic_analysis.r [options]
Options:
-v --village=<village>                       number of the village
-h --help                                    show this message" -> doc

opts <- docopt(doc)

if (opts[["help"]])
{
    print(opts)
    exit()
}

village_string <- opts[["nsamples"]]
villages <- as.integer(unlist(strsplit(village_string, ",")))

library('data.table')
library('trypR')

village <- villages[1]

lhs <- chronic_carriers_sample(nruns = 1, nsamples = 100, calc.posterior = FALSE,
                                   villages = village, passive = TRUE, sample = "lhs")
prior <- chronic_carriers_sample(nruns = 1, nsamples = 10000, calc.posterior = FALSE,
                                 villages = village, passive = TRUE, sample = "prior")

lhs_new <- chronic_carriers_sample(nruns = 1, nsamples = 10000, calc.posterior = FALSE,
                                   villages = village, passive = TRUE, sample = "lhs",
                                   transmitted = TRUE)
prior_new <- chronic_carriers_sample(nruns = 1, nsamples = 10000,
                                     calc.posterior = FALSE,
                                     villages = village, passive = TRUE,
                                     sample = "prior", transmitted = TRUE)

dt_lhs <- data.table(t(sapply(lhs, function(x)
{
    c(x[["parameters"]], unlist(x[["summary_statistics"]]))
})))
dt_prior <- data.table(t(sapply(prior, function(x)
{
    c(x[["parameters"]], unlist(x[["summary_statistics"]]))
})))
dt_lhs_new <- data.table(t(sapply(lhs_new, function(x)
{
    c(x[["parameters"]], unlist(x[["summary_statistics"]]))
})))
dt_prior_new <- data.table(t(sapply(prior_new, function(x)
{
    c(x[["parameters"]], unlist(x[["summary_statistics"]]))
})))

dt_lhs[, sample := "lhs"]
dt_prior[, sample := "prior"]

dt <- rbind(dt_lhs, dt_prior)

dt_lhs_new[, sample := "lhs"]
dt_prior_new[, sample := "prior"]
dt_new <- rbind(dt_lhs_new, dt_prior_new)

data_stage1_active <- village_screening[village.number == village, detected1_1]
data_stage1_passive <-
    village_cases[village.number == village & stage == 1, sum(cases)]
data_stage2_passive <-
    village_cases[village.number == village & stage == 2, sum(cases)]

dt[, active_stage1_diff := active_stage1 - data_stage1_active]
dt[, passive_stage1_diff := passive_stage1 - data_stage1_passive]
dt[, passive_stage2_diff := passive_stage2 - data_stage2_passive]
dt_new[, active_stage1_diff := active_stage1 - data_stage1_active]
dt_new[, passive_stage1_diff := passive_stage1 - data_stage1_passive]
dt_new[, passive_stage2_diff := passive_stage2 - data_stage2_passive]

threshold <- 0

dt[, accept := (abs(active_stage1_diff) <= threshold) & (abs(passive_stage1_diff) <= threshold) & (abs(passive_stage2_diff) <= threshold)]
dt_new[, accept := (abs(active_stage1_diff) <= threshold) & (abs(passive_stage1_diff) <= threshold) & (abs(passive_stage2_diff) <= threshold)]

good <- dt[accept == TRUE]
good_new <- dt_new[accept == TRUE]

prior_samples <- dt[sample == "prior"]
prior_samples[, sample := NULL]
saveRDS(prior_samples, "prior_samples.rds")
