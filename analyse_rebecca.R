library('XLConnect')
library('data.table')
library('ggplot')
library('ggthemr')

path <- path.expand("~/Research/Francesco/MSc\ project\ Rebecca\ Oettle/Results/")

files.list <- list.files(path, pattern = "Village.+.xlsx")

village_data <- list()

for (file in files.list)
{
    village_number <- sub("Village([0-9]+).+$", "\\1", file)
    
    wb <- loadWorkbook(paste0(path, "/", file), create = FALSE)
    df <- readWorksheet(wb, sheet = "Sheet3")

    village_data[[village_number]] <- copy(df)
}


plot_marginal <- function(df, var)
{
    ggthemr('fresh')
    
    dt <- data.table(df)
    func_dt <- dt[, list(Fit = mean(Fit)), by = var]
    p <- ggplot(func_dt, aes_string(x = var, y = "Fit"))
    p <- p + geom_line()
    return(p)
}

plot_marginal(dt, "Stage1")
