library('data.table')

chains <- NULL
for (i in seq_len(26)) {
    chain <- readRDS(paste("cc_mcmc_short_", i, ".rds", sep = ""))[["chain"]]
    chain[, village := i]
    if (is.null(chains))
    {
        chains <- chain
    } else
    {
        chains <- rbind(chains, chain)
    }
    p <- ggplot(chain, aes(x = value))
    p <- p + geom_histogram()
    p <- p + facet_wrap(~variable, scales = "free")
    p <- p + theme(legend.position = "bottom",
                   axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste("chain_village_", i, ".pdf", sep = ""))
}

