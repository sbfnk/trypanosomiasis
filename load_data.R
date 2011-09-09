library(filehash)

data <- read.csv("results.csv", header=T, sep=',', nrow=1499375)
db <- dumpDF(data, dbName="reservoirs")
rm(data)

renv <- db2env(db)
(hosts <- ls(renv)[c(67,39,79,89,1,10,19,30,57,98,107,116)])
(domains <- ls(renv)[c(67,69,77,28,29,125,88)])
species_data <- data.frame(db[hosts], db[domains[-1]], N=renv$Human_N,
  xi=renv$G._palpalis_gambiense_xi, groups=renv$groups, habitat=renv$habitat)
