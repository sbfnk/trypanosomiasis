library(filehash)
library(data.table)

data <- read.csv("results_20121025.csv", header=T, sep=',', nrow=1500000)
db <- dumpDF(data, dbName="reservoirs")
rm(data)

renv <- db2env(db)
(hosts <- ls(renv)[c(73,43,86,97,1,11,21,33,62,107,117,127)])
(domains <- ls(renv)[c(73,75,83,31,32,137,96)])

species_data <- data.frame(db[hosts], db[domains[-1]], N=renv$Human_N,
  xi=renv$G._palpalis_gambiense_xi, groups=renv$groups, habitat=renv$habitat)
rm(renv)

data <- read.csv("results_human_N.csv", header=T, sep=',')
db <- dumpDF(data, dbName="human_N")
rm(data)

renv <- db2env(db)
human_n_data <- data.frame(db[hosts], db[domains[-1]], N=renv$Human_N,
  xi=renv$G._palpalis_gambiense_xi, groups=renv$groups, habitat=renv$habitat)

rm(renv)

data <- read.csv("results_habitats.csv", header=T, sep=',')
db <- dumpDF(data, dbName="habitats")
rm(data)

renv <- db2env(db)
(hosts <- ls(renv)[c(91,59,107,119,1,15,29,45,79,130,144,156)])
(domains <- ls(renv)[c(91,93,101,43,44,167,118)])

habitat_data <- data.frame(db[hosts], db[domains[-1]], N=renv$Human_N,
                           xi=renv$G._palpalis_gambiense_xi,
                           tau=renv$G._palpalis_gambiense_tau,
                           groups=renv$groups, habitat=renv$habitat)

rm(renv)


data <- read.csv("/data/Tryps/results_habitats_new.csv", header=T, sep=',')
db <- dumpDF(data, dbName="habitats_new")
rm(data)

renv <- db2env(db)
(hosts <- ls(renv)[c(91,59,107,119,1,15,29,45,79,130,144,156)])
(domains <- ls(renv)[c(91,93,101,43,44,167,118)])

habitat_data_new<- data.frame(db[hosts], db[domains[-1]],
                                N=renv$Human_N,
                                xi=renv$G._palpalis_gambiense_xi,
                                groups=renv$groups, habitat=renv$habitat)

rm(renv)
