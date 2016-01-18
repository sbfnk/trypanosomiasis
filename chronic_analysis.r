library('dplyr')
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

limits <- list(pc = c(0, 0.5), alpha = c(0, 1), delta = c(0, 1), p1 = c(0, 30), p2 = c(0, 30), beta = c(10^(-5), 10^(0)), lambda = c(10^(-5), 10^(0)))

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
        mcmc[[i]][["trace"]] <- mcmc[[i]][["trace"]] %>% mutate(village = file_no)
    }

    traces <- rbind_all(lapply(mcmc, function(x) {x[["trace"]]}))
    acceptance.rates <- data.frame(file_no = file_nos, acceptance.rate = sapply(mcmc, function(x) {x[["acceptance.rate"]]}))
    traces <- traces %>% mutate(id = rep(seq(0, nrow(traces) / length(file_nos) - 1), length(file_nos)))

    ## burn-in
    traces <- traces %>% filter(id > 99)
    ## ggplot(thin.traces, aes(x = pc))+geom_density()+facet_wrap(~village)

    ## Kernel density estimation

    densities <- list()
    hists <- list()

    village_density <- list()
    village_hist <- list()

    dt <- list()
    dtm <- list()

    for (param in c("pc", "alpha", "delta", "p1", "p2", "beta", "lambda"))
    {

        densities[[param]] <- list()
        hists[[param]] <- list()
        
        min <- traces %>% select(get(param)) %>% min %>% floor
        max <- traces %>% select(get(param)) %>% max %>% ceiling
        
        i <- 0
        for (village.number in file_nos)
        {
            village_param <- traces %>% filter(village == village.number) %>% .[[param]]
            i <- i + 1
            densities[[param]][[i]] <- bde(village_param, 
                                           estimator = "vitale", lower.limit = limits[[param]][1],
                                           upper.limit = limits[[param]][2])
            hists[[param]][[i]] <- hist(village_param, 
                                        breaks = seq(limits[[param]][1], limits[[param]][2],
                                                     length.out = 101),
                                        plot = FALSE)
        }

        densm_x <- sapply(densities[[param]], function(x) {x@lower.limit + x@dataPointsCache * x@upper.limit})
        densm_y <- sapply(densities[[param]], function(x) {x@densityCache})
        histsm_x <- sapply(hists[[param]], function(x) log(x[["mids"]]))
        histsm_y <- sapply(hists[[param]], function(x) log(x[["density"]]))

        dt[[param]] <- data.frame(x = densm_x[, 1], y = apply(densm_y, 1, prod))
        dtm[[param]] <- data.frame(x = histsm_x[, 1],y = exp(apply(histsm_y, 1, sum)))

        ## normalise
        dt[[param]] <- dt[[param]] %>%
            mutate(y.norm = y * (length(x) - 1) / (sum(y) * (max(x) - min(x))))
        dtm[[param]] <- dtm[[param]] %>%
            mutate(y.norm = y * (length(x) - 1) / (sum(y) * (max(x) - min(x))))

        village_density[[param]] <-
            densm_y %>% as.data.frame %>% gather %>%
            mutate(village = as.integer(key), param = param, param.value = rep(densm_x[, 1], dim(densm_x)[2]), density = value) %>%
            select(param, village, param.value, density)

        village_hist[[param]] <-
            histsm_y %>% as.data.frame %>% gather %>%
            mutate(village = as.integer(key), param = param, param.value = rep(histsm_x[, 1], dim(densm_x)[2]), density = value) %>%
            select(param, village, param.value, density)
    }
}

village_density <- rbind_all(village_density)
village_hist <- rbind_all(village_hist)

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

    p[[this.param]][["density"]] <- ggplot(village_density %>% filter(param == this.param),
                                           aes(x = param.value, y = density))
    p[[this.param]][["density"]] <- p[[this.param]][["density"]] + geom_line()
    p[[this.param]][["density"]] <- p[[this.param]][["density"]] + facet_wrap(~village,
                                                                    scales = "free")
    ggsave(paste0("densities_", this.param, ".pdf"), p[[this.param]][["density"]],
           width = 14, height = 14)
    p[[this.param]][["histogram"]] <- ggplot(village_hist %>% filter(param == this.param),
                                             aes(x = param.value, y = density))
    p[[this.param]][["histogram"]] <- p[[this.param]][["histogram"]] + geom_line()
    p[[this.param]][["histogram"]] <- p[[this.param]][["histogram"]] + facet_wrap(~village,
                                                                        scales = "free")
    ggsave(paste0("histograms_", this.param, ".pdf"), p[[this.param]][["histogram"]],
           width = 14, height = 14)
}
