library('dplyr')
library('tidyr')
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

limits <- list(pc = c(0, 0.5), alpha = c(0, 1), delta = c(0, 1), p1 = c(0, 30), p2 = c(0, 30), beta = c(0, 10^(0)), lambda = c(0, 10^(0)))

mcmc <- list()
village_density <- list()
densities <- list()
dt <- list()
traces <- list()

## for (model in c("tran_back_chro", "tran_back", "back"))
for (model in c("tran_back_chro", "tran_back"))
{
    file_path <- path.expand("~/Data/Trypanosomiasis/")
    file_pattern <- paste(model, ".*", sep = "_")
    file_nos <- as.integer(gsub("[^0-9]", "",
                                list.files(path = file_path, pattern = file_pattern)))

    mcmc[[model]] <- list()
    i <- 0
    ## file_nos <- runs[epsilon <= 1, run]

    for (file_no in file_nos)
    {
        i <- i + 1
        file_name <-
            path.expand(paste0(file_path, "/", model, "_village_", file_no, ".rds"))
        mcmc[[model]][[i]] <- readRDS(file_name)
        mcmc[[model]][[i]][["trace"]] <- mcmc[[model]][[i]][["trace"]] %>%
            mutate(village = file_no)
    }

    traces[[model]] <- rbind_all(lapply(mcmc[[model]], function(x) {x[["trace"]]})) %>%
        mutate(model = model)
    acceptance.rates <- data.frame(file_no = file_nos, acceptance.rate = sapply(mcmc[[model]], function(x) {x[["acceptance.rate"]]}))
    traces[[model]] <- traces[[model]] %>% mutate(id = rep(seq(0, n() / length(file_nos) - 1), length(file_nos)))

    ## burn-in
    traces[[model]] <- traces[[model]] %>% filter(id > 99)
    ## ggplot(thin.traces[[model]], aes(x = pc))+geom_density()+facet_wrap(~village)

    ## Kernel density estimation

    densities[[model]] <- list()
    village_density[[model]] <- list()

    dt[[model]] <- list()
    dtm[[model]] <- list()

    for (param in c("pc", "alpha", "delta", "p1", "p2", "beta", "lambda"))
    {
        if (traces[[model]] %>% select(get(param)) %>% unique %>% nrow > 1)
        {

            densities[[model]][[param]] <- list()

            min <- traces[[model]] %>% select(get(param)) %>% min %>% floor
            max <- traces[[model]] %>% select(get(param)) %>% max %>% ceiling

            i <- 0
            for (village.number in file_nos)
            {
                village_param <- traces[[model]] %>%
                    filter(village == village.number) %>% .[[param]]
                i <- i + 1
                densities[[model]][[param]][[i]] <-
                    bde(village_param,
                        estimator = "vitale", lower.limit = limits[[param]][1],
                        upper.limit = limits[[param]][2])
            }

            densm_x <- sapply(densities[[model]][[param]],
                              function(x) {x@lower.limit +
                                               x@dataPointsCache * x@upper.limit})
            densm_y <- sapply(densities[[model]][[param]], function(x) {x@densityCache})

            dt[[model]][[param]] <-
                data.frame(x = densm_x[, 1], y = apply(densm_y, 1, prod))
            ## normalise
            dt[[model]][[param]] <- dt[[model]][[param]] %>%
                mutate(y.norm = y * (length(x) - 1) / (sum(y) * (max(x) - min(x))))

            village_density[[model]][[param]] <-
                densm_y %>% as.data.frame %>% gather %>%
                mutate(file = as.integer(key), param = param,
                       param.value = rep(densm_x[, 1], dim(densm_x)[2]),
                       density = value, model = model) %>%
                mutate(village = file_nos[file]) %>%
                select(model, param, village, param.value, density)
        }

    }
    village_density[[model]] <- rbind_all(village_density[[model]])
}

village_density <- rbind_all(village_density)

p <- list()
for (model in c("tran_back_chro", "tran_back", "back"))
{
    for (param in names(dt[[model]]))
    {
        p <- ggplot(dt[[model]][[param]], aes(x = x, y = y.norm))+geom_line()
        ggsave(paste0("density_", model, "_", param, ".pdf"))
    }


    ## check caterpillars
    p[[model]] <- list()

    for (this.param in names(dt[[model]]))
    {
        p[[this.param]] <- list()
        p[[model]][[this.param]][["trace"]] <-
            ggplot(traces[[model]], aes_string(x = "id", y = this.param))
        p[[model]][[this.param]][["trace"]] <-
            p[[model]][[this.param]][["trace"]] + geom_line()
        p[[model]][[this.param]][["trace"]] <-
            p[[model]][[this.param]][["trace"]] + facet_wrap(~village, scales = "free")
        ggsave(paste0("caterpillars_", model, "_", this.param, ".pdf"),
               p[[model]][[this.param]][["trace"]], width = 14, height = 14)

        p[[model]][[this.param]][["density"]] <-
            ggplot(village_density %>% filter(param == this.param),
                   aes(x = param.value, y = density))
        p[[model]][[this.param]][["density"]] <-
            p[[model]][[this.param]][["density"]] + geom_line()
        p[[model]][[this.param]][["density"]] <-
            p[[model]][[this.param]][["density"]] +
            facet_wrap(~village, scales = "free")
        ggsave(paste0("densities_", model, "_", this.param, ".pdf"),
               p[[model]][[this.param]][["density"]],
               width = 14, height = 14)
    }
}
