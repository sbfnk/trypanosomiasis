library('lhs')

source(path.expand("~/code/trypanosomiasis/chronic_carriers.R"))
village_data <- readRDS("village_data.rds")

village_cases <- village_data[["cases"]]
village_screening <- village_data[["screening"]]

args <- commandArgs(trailingOnly = TRUE)
village_no <- as.integer(args[1])

samples <- list()
inits <- list()
distances <- c()
likelihoods <- c()
posteriors <- c()

## set.seed(42)
n_samples <- 10000
n_runs <- 100

r <- randomLHS(n_samples, 5)
upper <- c(1, 0.05, 1, 2, 1)
lower <- c(0, 0, 0, 0, 0)
r <- t(apply(r, 1, function(x) { lower + (upper - lower) * x}))
colnames(r) <- c("pc", "lambda", "p1", "p2", "alpha")
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

saveRDS(list(samples = samples, posteriors = posteriors, seed = .Random.seed),
        paste("cc_lhs_eq_", village_no, ".rds", sep = ""))

