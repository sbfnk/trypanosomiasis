library('XLConnect')
library('data.table')
library('ggplot2')
library('ggthemr')
library('reshape2')

ggthemr('fresh')

path <- path.expand("~/Research/Francesco/MSc\ project\ Rebecca\ Oettle/Results/")

files.list <- list.files(path, pattern = "Village.+.xlsx")

village_data <- list()

for (file in files.list)
{
    village_number <- sub("Village([0-9]+).+$", "\\1", file)
    
    wb <- loadWorkbook(paste0(path, "/", file), create = FALSE)
    df <- readWorksheet(wb, sheet = "Sheet3")
    dt <- data.table(df)
    setnames(dt, names(dt), tolower(names(dt)))
    dt <- dt[, 1:8, with = FALSE]
    if ("bathrun.number" %in% colnames(dt))
    {
        setnames(dt, "bathrun.number", "batchrun.number")
    }
    dt[, village.number := village_number]

    village_data[[village_number]] <- copy(dt)
    
}

dt <- rbindlist(village_data)

fit_dt <- data.table(value = numeric(0), fit = numeric(0), param = character(0))

for (param in setdiff(names(dt), c("batchrun.number", "fit", "village.number")))
{
    param_dt <- dt[, list(fit = mean(fit)), by = c(param, "village.number")]
    param_dt <- param_dt[, list(fit = sum(fit)), by = param]
    setnames(param_dt, param, "value")
    param_dt[, param := param]
    fit_dt <- rbind(fit_dt, param_dt)
}

p <- ggplot(fit_dt, aes(x = value, y = fit))
p <- p + geom_line()
p <- p + facet_wrap(~ param, scales = "free")

final_p <- ggplot(fit_dt[param %in% c("alpha", "pc")], aes(x = value, y = fit))
final_p <- final_p + geom_line()
final_p <- final_p + facet_wrap(~ param, scales = "free")

lambda_dt <- dt[, list(fit = mean(fit)), by = c("lambda", "village.number")]
lambda_p <- ggplot(lambda_dt, aes(x = lambda, y = fit))+ geom_line()+ facet_wrap(~village.number, scales = "free")


pc_dt <- dt[, list(fit = mean(fit)), by = c("pc", "village.number")]
pc_p <- ggplot(pc_dt, aes(x = pc, y = fit))+ geom_line()+ facet_wrap(~village.number, scales = "free")

lambda1_dt <- dt[, list(fit = mean(fit)), by = c("lambda", "village.number")]
lambda_p <- ggplot(alpha_dt, aes(x = lambda, y = fit))+ geom_line()+ facet_wrap(~village.number, scales = "free")

