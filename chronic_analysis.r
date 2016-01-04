library('data.table')
library('coda')
library('modeest')
library('cowplot')

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

for (model in c("tran_back_chro", "back"))
{
    file_path <- path.expand("~/Data/Trypanosomiasis/")
    file_pattern <- paste(model, ".*", sep = "_")
    file_nos <- as.integer(gsub("[^0-9]", "",
                                list.files(path = file_path, pattern = file_pattern)))

    traces <- list()
    i <- 0
    ## file_nos <- runs[epsilon <= 1, run]

    for (file_no in file_nos)
    {
        i <- i + 1
        file_name <-
            path.expand(paste0(file_path, "/", model, "_village_", file_no, ".rds"))
        traces[[i]] <- data.table(readRDS(file_name))
        traces[[i]][, village := file_no]
    }

    traces <- rbindlist(traces)
    traces[, id := rep(seq(0, nrow(traces) / length(file_nos) - 1), length(file_nos))]

    thin.traces <- traces[id %% 10 == 0]
    ## ggplot(thin.traces, aes(x = pc))+geom_density()+facet_wrap(~village)

    ## Kernel density estimation

    densities <- list()
    hists <- list()

    densm <- list()
    histsm <- list()

    dt <- list()
    dtm <- list()

    for (param in c("pc", "alpha", "delta", "p1", "p2", "lambda"))
    {

        densities[[param]] <- list()
        hists[[param]] <- list()
        
        min <- floor(traces[, min(get(param))])
        max <- ceiling(traces[, max(get(param))])
        
        i <- 0
        for (village.number in file_nos)
        {
            i <- i + 1
            densities[[param]][[i]] <- density(traces[village == village.number, get(param)],
                                      from = -0.2, to=1.2, kernel = "gaussian")
            hists[[param]][[i]] <- hist(traces[village == village.number, get(param)],
                               breaks = seq(min, max, length.out = 1000), plot = FALSE)
        }

        densm[[param]] <- sapply(densities[[param]], function(x) x[["y"]])
        histsm[[param]] <- sapply(hists[[param]], function(x) log(x[["density"]]))

        dt[[param]] <- data.table(x = densities[[param]][[1]]$x, y = apply(densm[[param]], 1, prod))
        dtm[[param]] <- data.table(x = hists[[param]][[1]]$breaks[-length(hists[[param]][[1]]$breaks)] + 0.0005,
                          y = exp(apply(histsm[[param]], 1, sum)))
        dt[[param]][, y.norm := y / 427]
        dtm[[param]][, y.norm := y / (sum(y) * 1000)]
        
        ## p <- ggplot(dt, aes(x = x, y = y.norm))+geom_line()
        ## p <- ggplot(dtm, aes(x = x, y = y.norm))+geom_bar(stat = "identity")
        ## ggsave(paste0("density_", param, ".pdf"))
    }
}

## check caterpillars
thin.traces <- traces[id %% 10 == 0]
thin.traces <- thin.traces[id > 9999]
p <- list()

for (param in c("pc", "alpha", "delta", "p1", "p2", "lambda", "beta"))
{
    p[[param]] <- list()
    p[[param]][["trace"]] <- ggplot(thin.traces, aes_string(x = "id", y = param))
    p[[param]][["trace"]] <- p[[param]][["trace"]] + geom_line()
    p[[param]][["trace"]] <- p[[param]][["trace"]] + facet_wrap(~village,
                                                                scales = "free")
    ggsave(paste0("caterpillars_", param, ".pdf"), p[[param]][["trace"]],
           width = 14, height = 14)
    p[[param]][["density"]] <- ggplot(thin.traces, aes_string(x = param))
    p[[param]][["density"]] <- p[[param]][["density"]] + geom_density(adjust = 2)
    p[[param]][["density"]] <- p[[param]][["density"]] + facet_wrap(~village,
                                                                    scales = "free")
    ggsave(paste0("densities_", param, ".pdf"), p[[param]][["density"]],
           width = 14, height = 14)
    p[[param]][["histogram"]] <- ggplot(thin.traces, aes_string(x = param))
    p[[param]][["histogram"]] <- p[[param]][["histogram"]] + geom_histogram()
    p[[param]][["histogram"]] <- p[[param]][["histogram"]] + facet_wrap(~village,
                                                                        scales = "free")
    ggsave(paste0("histograms_", param, ".pdf"), p[[param]][["histogram"]],
           width = 14, height = 14)
}
