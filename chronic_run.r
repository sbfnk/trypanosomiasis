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

sample_options <- c(list(nsamples = 1,
                         calc.posterior = TRUE,
                         village_number = village,
                         sample = "prior", 
                         demand.finite = TRUE, 
                         n = 1 
                         ),
                    model_options)

first_sample <- do.call(chronic_carriers_sample, sample_options)

sample_options[["nsamples"]] <- num_samples - 1
sample_options[["demand.finite"]] <- FALSE
sample_options[["sample"]] <- ifelse(opts[["lhs"]], "lhs", "prior")

samples <- do.call(chronic_carriers_sample, sample_options)
dt <- rbindlist(lapply(c(first_sample, samples), function(x) {data.frame(t(unlist(x)))}), fill = TRUE)
saveRDS(dt, paste0("prior_",
                   ifelse(opts[["transmitted"]], "tran_", ""),
                   ifelse(opts[["background"]], "back_", ""),
                   ifelse(opts[["chronic"]], "chro_", ""),
                   "village_", village, "_prior.rds"))

prior <- dt[is.finite(log.posterior)]
prior_sd <- apply(dt, 2, sd)
prior_sd <- prior_sd[is.finite(prior_sd)]
prior_sd <- prior_sd[grep("^theta", names(prior_sd))]
prior_sd <- prior_sd[prior_sd > 0]
names(prior_sd) <- sub("^theta\\.", "", names(prior_sd))
prior <- prior[, names(prior_sd), with = FALSE]

setnames(prior, colnames(prior), sub("^theta\\.", "", colnames(prior)))
cov_prior <- cov(prior)

start <- unlist(dt[log.posterior == max(log.posterior)])
start <- start[grep("^theta", names(start))]
names(start) <- sub("^theta.", "", names(start))

prior_zero <- start
prior_zero[] <- 0

print(prior_sd)

mcmc_posterior <- function(theta)
{
  theta["delta"] = theta["alpha"]
  post <- do.call(param_posterior,
                  c(list(theta, village_number = village, n = 10),
                    model_options))
  return(list(log.density = post[["log.posterior"]],
              trace = c(theta, log.prior = post[["log.prior"]])))
}

mcmc_options <-
  list(target = mcmc_posterior,
       init.theta = start, 
       proposal.sd = prior_sd,
       n.iterations = num_samples,
       limits = list(lower = prior_zero),
       adapt.size.start = 100,
       adapt.shape.start = 500,
       verbose = TRUE) 
mcmc <- do.call(mcmcMH, mcmc_options)



## ## df <- data.frame(matrix(unlist(mcmc$trace), ncol = ncol(prior_parameters), byrow = TRUE))
## ## colnames(df) <- names(mcmc$trace[[1]])

## saveRDS(mcmc, paste0(ifelse(opts[["transmitted"]], "tran_", ""),
##                      ifelse(opts[["background"]], "back_", ""),
##                      ifelse(opts[["chronic"]], "chro_", ""),
##                      "village_", village, ".rds"))

## ## df_fixed <- data.frame(matrix(unlist(x), ncol = ncol(prior_parameters), byrow = TRUE))
## ## colnames(df_fixed) <- colnames(prior_parameters)
## ## mcmc_fixed <- mcmc(df_fixed)
## ## accRate <- 1 - min(rejectionRate(mcmc_fixed))

system.time(init_test <- lapply(1:10000,  function(x) {
theta <- rprior(villages = village, background = TRUE, transmitted = TRUE, chronic = TRUE)
init <- rinit(theta)
return(list(theta = theta, init = init, log.init = dinit(init, village, theta, TRUE)))}))

param_bounds <- list(
  pc = c(0, 0.5), 
  alpha = c(0, 0.5), 
  delta = c(0, 0.5), 
  p1 = c(0, 30), 
  p2 = c(0, 30), 
  lambda = c(0, 0.01), 
  beta = c(0, 0.01))

bounds_to_grid <- function(bounds, n)
{
  grid <- expand.grid(lapply(param_bounds, function(x)
  {
    if (x[1] == 0) {
      x[1] <- 10**(log10(x[2]) - n + 1)
    }
    return(10**seq(log10(x[1]), log10(x[2]), length.out = n))
  }))
  return(grid)
}

system.time(scan <- apply(grid, 1, function(x)
{
  theta <- rprior(village, transmitted = TRUE)
  theta[names(x)] <- x
  log.posterior <- param_posterior(theta, village_number = village)
  return(log.posterior)
}))

dt <- rbindlist(lapply(scan, function(x) {data.frame(t(unlist(x)))}), fill = TRUE)
id <- melt(dt, id.vars = "log.init")
