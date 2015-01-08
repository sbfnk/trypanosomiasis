library('lhs')

source(path.expand("~/code/trypanosomiasis/chronic_carriers.R"))
village_data <- readRDS("village_data.rds")

village_cases <- village_data[["cases"]]
village_screening <- village_data[["screening"]]
nb_villages <- nrow(village_screening)

samples <- list()
inits <- list()
distances <- c()
likelihoods <- c()
posteriors <- c()

## set.seed(42)
n_samples <- 100000
n_runs <- 100

r <- randomLHS(n_samples, nb_villages * 3 + 2)
## upper <- c(1, 0.05, 1, 2, 1)
## lower <- c(0, 0, 0, 0, 0)
upper <- c(1, 1, rep(c(0.05, 1, 2), each = nb_villages))
lower <- c(0, 0, rep(c(0, 0, 0), each = nb_villages))
r <- t(apply(r, 1, function(x) { lower + (upper - lower) * x}))
## colnames(r) <- c("pc", "lambda", "p1", "p2", "alpha")
colnames(r) <- c("pc", "alpha",
                 paste(rep(c("lambda", "p1", "p2"), each = nb_villages),
                       seq_len(nb_villages), sep = "."))

theta <- rprior(villages = nb_villages)

cum_data <-
    village_cases[, list(cases = sum(cases)), by = list(village.number, stage)]

i <- 0
while(i < n_samples)
{
    theta[colnames(r)] <- r[i + 1, ]
    i <- i + 1
    samples[[i]] <- theta
    posteriors[i] <- param_posterior_all(theta, village_screening,
                                         cum_data, n_runs, TRUE)
    cat(i, posteriors[i], "\n")
}

saveRDS(list(samples = samples, posteriors = posteriors, seed = .Random.seed),
        paste("cc_lhs_eq_", village_no, ".rds", sep = ""))

