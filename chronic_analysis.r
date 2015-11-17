library('data.table')
library('trypR')

village <- 1

system.time(lhs_samples <- chronic_carriers_lhs(passive = FALSE, nsamples = 1000, nruns = 5, verbose = FALSE, village = village))

theta <- rprior(villages = 1, passive = FALSE)

current.posterior <- traj_likelihood(theta, 1, 5, log = TRUE, passive = FALSE)

acceptance <- 0
for (i in 1:10000)
{
    proposal.posterior <- traj_likelihood(theta, 1, 5, log = TRUE, passive = FALSE)
    if (runif(1) < exp(proposal.posterior - current.posterior))
    {
        current.posterior <- proposal.posterior
        acceptance <- acceptance + 1
    }
    cat(acceptance / i, "\n")
}


lhs_posteriors <- sapply(seq_along(lhs_samples), function(x)
{
    lhs_samples[[x]][["posterior"]]
})

## param_posterior_villages(theta = samples[[92]]$parameters, nruns = 100, log = TRUE)

top.10 <- order(lhs_posteriors, decreasing = TRUE)[1:10]

top10_lhs_posteriors <- lhs_posteriors[top.10]
top10_lhs_samples<- sapply(top.10, function(x)
{
    lhs_samples[[x]][["parameters"]]
})
