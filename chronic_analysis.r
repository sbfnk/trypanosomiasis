library('data.table')
library('coda')
library('modeest')
library('cowplot')
library('bde')

## m <- matrix(c(1,0,
##               10,0,
##               11,0,
##               12,8,
##               13,1,
##               14,2,
##               15,6,
##               16,9,
##               17,11,
##               18,5,
##               19,1,
##               2,1,
##               20,1,
##               21,7,
##               22,3,
##               23,2,
##               24,0,
##               25,4,
##               26,1,
##               3,0,
##               4,0,
##               5,5,
##               6,1,
##               7,1,
##               8,1,
##               9,0), ncol = 2, byrow = TRUE)
## runs <- data.table(m)
## setnames(runs, 1:2, c("run", "epsilon"))

limits <- list(pc = c(0, 0.5), alpha = c(0, 1), delta = c(0, 1), p1 = c(0, 1/7), p2 = c(0, 1/7), beta = c(0, 0.05), lambda = c(10^(-5), 10^(-2)))

## for (model in c("tran_back_chro", "back"))
for (model in c("tran_back_chro"))
{
    file_path <- path.expand("~/Data/Trypanosomiasis/")
    file_pattern <- paste(model, ".*", sep = "_")
    file_nos <- as.integer(gsub("[^0-9]", "",
                                list.files(path = file_path, pattern = file_pattern)))

    mcmc <- list()
    i <- 0
    ## file_nos <- runs[epsilon <= 1, run]

    for (file_no in file_nos)
    {
        i <- i + 1
        file_name <-
            path.expand(paste0(file_path, "/", model, "_village_", file_no, ".rds"))
        mcmc[[i]] <- readRDS(file_name)
        mcmc[[i]][["trace"]] <- data.table(mcmc[[i]][["trace"]])
        mcmc[[i]][["trace"]][, village := file_no]
    }

    traces <- rbindlist(lapply(mcmc, function(x) {x[["trace"]]}))
    acceptance.rates <- data.table(file_no = file_nos, acceptance.rate = sapply(mcmc, function(x) {x[["acceptance.rate"]]}))
    traces[, id := rep(seq(0, nrow(traces) / length(file_nos) - 1), length(file_nos))]

    ## burn-in
    traces <- traces[id > 99]
    ## ggplot(thin.traces, aes(x = pc))+geom_density()+facet_wrap(~village)

    ## Kernel density estimation

    densities <- list()
    hists <- list()

    densm <- list()
    histsm <- list()

    village_density <- list()
    village_hist <- list()

    dt <- list()
    dtm <- list()

    for (param in c("pc", "alpha", "delta", "p1", "p2", "beta", "lambda"))
    {

        densities[[param]] <- list()
        hists[[param]] <- list()
        
        min <- floor(traces[, min(get(param))])
        max <- ceiling(traces[, max(get(param))])
        
        i <- 0
        for (village.number in file_nos)
        {
            i <- i + 1
            densities[[param]][[i]] <- bde(traces[village == village.number, get(param)],
                                           estimator = "vitale", lower.limit = limits[[param]][1],
                                           upper.limit = limits[[param]][2])
            hists[[param]][[i]] <- hist(traces[village == village.number, get(param)],
                                        breaks = seq(limits[[param]][1], limits[[param]][2],
                                                     length.out = 101),
                                        plot = FALSE)
        }

        densm[[param]] <- sapply(densities[[param]], function(x) {x@densityCache})
        histsm[[param]] <- sapply(hists[[param]], function(x) log(x[["density"]]))

        dt[[param]] <- data.table(x = seq(limits[[param]][1], limits[[param]][2],
                                          length.out = 101), y = apply(densm[[param]], 1, prod))
        dtm[[param]] <- data.table(x = hists[[param]][[1]]$breaks[-length(hists[[param]][[1]]$breaks)] + 0.0005,
                                   y = exp(apply(histsm[[param]], 1, sum)))
        dt[[param]][, y.norm := y * 100 / (sum(y) * (limits[[param]][2] - limits[[param]][1]))]
        dtm[[param]][, y.norm := y * 100 / (sum(y) * (limits[[param]][2] - limits[[param]][1]))]

        village_density[[param]] <- data.table(melt(densm[[param]]))
        setnames(village_density[[param]], c("Var1", "Var2"), c("param.value", "village"))
        village_density[[param]][, param.value := seq(limits[[param]][1], limits[[param]][2],
                                                      length.out = 101)[param.value]]
        village_density[[param]][, param := param]

        village_hist[[param]] <- data.table(melt(densm[[param]]))
        setnames(village_hist[[param]], c("Var1", "Var2"), c("param.value", "village"))
        village_hist[[param]][, param.value := seq(limits[[param]][1], limits[[param]][2],
                                            length.out = 101)[param.value]]
        village_hist[[param]][, param := param]
    }
}

village_density <- rbindlist(village_density)
village_hist <- rbindlist(village_hist)

for (param in c("pc", "alpha", "delta", "p1", "p2", "beta", "lambda"))
{
    p <- ggplot(dt[[param]], aes(x = x, y = y.norm))+geom_line()
    ggsave(paste0("density_", param, ".pdf"))
##    p <- ggplot(dtm[[param]], aes(x = x, y = y.norm))+geom_bar(stat = "identity")
##    ggsave(paste0("histogram_", param, ".pdf"))
}


## check caterpillars
p <- list()

for (this.param in c("pc", "alpha", "delta", "p1", "p2", "lambda", "beta"))
{
    p[[this.param]] <- list()
    p[[this.param]][["trace"]] <- ggplot(traces, aes_string(x = "id", y = this.param))
    p[[this.param]][["trace"]] <- p[[this.param]][["trace"]] + geom_line()
    p[[this.param]][["trace"]] <- p[[this.param]][["trace"]] + facet_wrap(~village,
                                                                scales = "free")
    ggsave(paste0("caterpillars_", this.param, ".pdf"), p[[this.param]][["trace"]],
           width = 14, height = 14)

    p[[this.param]][["density"]] <- ggplot(village_density[param == this.param], aes(x = param.value, y = value))
    p[[this.param]][["density"]] <- p[[this.param]][["density"]] + geom_line()
    p[[this.param]][["density"]] <- p[[this.param]][["density"]] + facet_wrap(~village,
                                                                    scales = "free")
    ggsave(paste0("densities_", this.param, ".pdf"), p[[this.param]][["density"]],
           width = 14, height = 14)
    p[[this.param]][["histogram"]] <- ggplot(village_hist[param == this.param], aes(x = param.value, y = value))
    p[[this.param]][["histogram"]] <- p[[this.param]][["histogram"]] + geom_line()
    p[[this.param]][["histogram"]] <- p[[this.param]][["histogram"]] + facet_wrap(~village,
                                                                        scales = "free")
    ggsave(paste0("histograms_", this.param, ".pdf"), p[[this.param]][["histogram"]],
           width = 14, height = 14)
}
