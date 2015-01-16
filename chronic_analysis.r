library('data.table')
library('trypR')
library('compiler')

enableJIT(1)

{
    timing_lhs <- system.time(lhs_samples <- chronic_carriers_lhs(passive = FALSE, nsamples = 100))
    timing_prior <- system.time(prior_samples <- chronic_carriers_prior(passive = FALSE, nsamples = 10, finite = TRUE, verbose = TRUE))
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
