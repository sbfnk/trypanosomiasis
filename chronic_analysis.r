library('data.table')
library('coda')
library('modeest')
library('ggthemr')
ggthemr('fresh')

m <- matrix(c(1,0,
              10,0,
              11,0,
              12,8,
              13,1,
              14,2,
              15,6,
              16,9,
              17,11,
              18,5,
              19,1,
              2,1,
              20,1,
              21,7,
              22,3,
              23,2,
              24,0,
              25,4,
              26,1,
              3,0,
              4,0,
              5,5,
              6,1,
              7,1,
              8,1,
              9,0), ncol = 2, byrow = TRUE)
runs <- data.table(m)
setnames(runs, 1:2, c("run", "epsilon"))

traces <- list()

## tran_back_chro
file_nos <- as.integer(gsub("[^0-9]", "", list.files(path = (path.expand("~/Data/Trypanosomiasis/")), pattern = "tran_back_chro_.*")))

for (file_no in file_nos)
{
    traces[[file_no]] <- data.table(readRDS(path.expand(paste0("~/Data/Trypanosomiasis/tran_back_chro_village_", file_no, ".rds"))))
}

traces <- list()
i <- 0
file_nos <- runs[epsilon <= 1, run]

for (file_no in file_nos)
{
    i <- i + 1
    traces[[i]] <- data.table(readRDS(path.expand(paste0("~/Data/Trypanosomiasis/village_", file_no, ".rds"))))
    traces[[i]][, village := file_no]
}

traces <- rbindlist(traces)

## Kernel density estimation

for (param in c("pc", "alpha", "delta"))
{
    densities <- list()

    i <- 0
    for (village.number in file_nos)
    {
        i <- i + 1
        densities[[i]] <- density(traces[village == village.number, get(param)],
                                  from = -0.1, to=1.1)
    }

    densm <- sapply(densities, function(x) x[["y"]])
    dt <- data.table(x = densities[[1]]$x, y = apply(densm, 1, prod))
    dt[, y.norm := y / 427]

    p <- ggplot(dt, aes(x = x, y = y.norm))+geom_line()
    ggsave(paste0("density_", param, ".pdf"))
}
