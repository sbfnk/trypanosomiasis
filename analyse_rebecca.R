library('XLConnect')
library('data.table')
library('ggplot2')
library('ggthemr')

path <- path.expand("~/Research/Francesco/MSc\ project\ Rebecca\ Oettle/Results/")

files.list <- list.files(path, pattern = "^Village.+.xlsx")

village_data <- list()

for (file in files.list)
{
    village_number <- sub("Village([0-9]+).+$", "\\1", file)
    
    wb <- loadWorkbook(paste0(path, "/", file), create = FALSE)
    dt <- data.table(readWorksheet(wb, sheet = "Sheet3"))
    dt[, village := village_number]

    setnames(dt, colnames(dt), tolower(colnames(dt)))

    if ("bathrun.number" %in% colnames(dt))
    {
        setnames(dt, "bathrun.number", "batchrun.number")
    }
    dt <- dt[, colnames(dt)[!grepl("^col", colnames(dt))], with = FALSE]

    village_data[[village_number]] <- copy(dt)
}

dt <- rbindlist(village_data)
dt[, village := as.integer(as.character(village))]
setkey(dt, village)
saveRDS(dt, "village_results.rds")

plot_marginal <- function(df, var)
{
    ggthemr('fresh')
    
    dt <- data.table(df)
    func_dt <- dt[, list(fit = mean(fit)), by = var]
    p <- ggplot(func_dt, aes_string(x = var, y = "fit"))
    p <- p + geom_line()
    return(p)
}


