library('XLConnect')
library('data.table')

wb <- loadWorkbook("~/Research/Francesco/MSc\ project\ Rebecca\ Oettle/Results/Village10_Result1.xlsx", create = FALSE)
df <- readWorksheet(wb, sheet = "Sheet3")
dt <- data.table(df)
