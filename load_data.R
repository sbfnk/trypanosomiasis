library(filehash)
library(data.table)

data <- read.csv("results_20121025.csv", header=T, sep=',', nrow=1500000)
db <- dumpDF(data, dbName="reservoirs")
rm(data)

renv <- db2env(db)
## (hosts <- ls(renv)[c(67,39,79,89,1,10,19,30,57,98,107,116)])
## (domains <- ls(renv)[c(67,69,77,28,29,125,88)])
(hosts <- ls(renv)[c(73,43,86,97,1,11,21,33,62,107,117,127)])
(domains <- ls(renv)[c(73,75,83,31,32,137,96)])
species_data <- data.frame(db[hosts], db[domains[-1]], N=renv$Human_N,
  xi=renv$G._palpalis_gambiense_xi, groups=renv$groups, habitat=renv$habitat)

data <- read.csv("results_single.csv", header=T, sep=',')
db <- dumpDF(data, dbName="single")
rm(data)

renv <- db2env(db)
species_data_single <- data.frame(db[hosts], db[domains[-1]], N=renv$Human_N,
  xi=renv$G._palpalis_gambiense_xi, groups=renv$groups, habitat=renv$habitat)

rm(renv)

