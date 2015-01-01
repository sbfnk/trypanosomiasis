samples <- list()
inits <- list()
distances <- c()
likelihoods <- c()
posteriors <- c()

set.seed(42)
n_samples <- 10000
n_runs <- 100

r <- randomLHS(n_samples, 5)
upper <- c(0.5, 0.05, 1, 1, 1)
lower <- c(0, 0, 0, 0, 0)
r <- t(apply(r, 1, function(x) { lower + (upper - lower) * x}))
colnames(r) <- c("pc", "lambda", "p1", "p2", "alpha")
theta <- rprior()

cum_data <-
    village_cases[village.number == 2, list(cases = sum(cases)), by = stage]$cases

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

i <- 0
while(i < n_samples)
{
    theta[colnames(r)] <- r[i + 1, ]
    i <- i + 1
    samples[[i]] <- theta
    posteriors[i] <- param_posterior(theta, village_screening[village.number == 2],
                                     cum_data, n_runs, TRUE, rand = rand)
    cat(i, posteriors[i], "\n")
}

saveRDS(list(samples = samples, posteriors = posteriors), "cc_lhs_samples.rds")

