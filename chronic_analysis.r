library('data.table')
library('trypR')

samples <- chronic_carriers_lhs(passive = FALSE, nsamples = 100, finite = FALSE)

posteriors <- sapply(seq_along(samples), function(x)
{
    samples[[x]][["posterior"]]
})

## param_posterior_villages(theta = samples[[92]]$parameters, nruns = 100, log = TRUE)

n.inf <- apply(posteriors, 2, function(x) {sum(!is.finite(x))})
finite.runs <- which(n.inf == 0)
top.10 <- finite.runs[order(posteriors)][1:10]
finite.samples <- sapply(top.10, function(x)
{
    samples[[x]][["parameters"]]
})
