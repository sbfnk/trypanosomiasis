library('XLConnect')

wb <- loadWorkbook(path.expand("~/Research/Analysis/Sleeping_Sickness/Chronic/Village input data.xlsx"))
village_screening <- data.table(readWorksheet(wb, "Sheet1"))

village_dirs <- list.files(path = path.expand("~/Research/Analysis/Sleeping_Sickness/Chronic/passive_screening_datasets/"))
village_ids <- as.integer(gsub("^village ", "", village_dirs))
village_ids <- village_ids[order(village_ids)]

village_cases <- data.table(village.number = integer(0), month = integer(0), cases = integer(0), stage = integer(0))
for (id in village_ids)
{
    for (stage in c(1, 2))
    {
        temp <- data.table(read.csv(paste("~/Research/Analysis/Sleeping_Sickness/Chronic/passive_screening_datasets/village ",
                                          id, "/ps", stage, "_", id, ".csv", sep = ""), header = F))
        setnames(temp, 1:2, c("month", "cases"))
        village_cases <- rbind(village_cases, data.table(village.number = rep(id, nrow(temp)),
                                                         month = temp[, month],
                                                         cases = temp[, cases],
                                                         stage = rep(stage, nrow(temp))))
    }
}

saveRDS(list(cases = village_cases, screening = village_screening), "village_data.rds")

sim_chronic_transitions <- list(
    c(S = -1, Ic = +1), # chronic infection
    c(Ic = -1, S = +1), # chronic recovery
    c(S = -1, I1 = +1), # pathogenic infection
    c(I1 = -1, I2 = +1), # progression into stage 2
    c(I1 = -1, S = +1, Z1pass = +1), # passive stage 1 detection
    c(I2 = -1, S = +1, Z2pass = +1), # passive stage 2 detection
    c(I2 = -1, S = +1) # stage 2 death
    )

