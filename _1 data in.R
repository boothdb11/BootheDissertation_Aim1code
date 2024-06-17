
#### Note: read in the csv data file, make necessary cleaning and correction, and then save the data
####  as a R data object.

rm(list=ls())

library(tidyr)

setwd("C:/Users/dbboothe/Desktop/TBdatafor dissertation/Aim1code")

dat.ori <- read.csv('./aim1datasetf.csv', stringsAsFactors = FALSE)

unique(dat.ori$subpop)
dat.ori$subpop <- trimws(dat.ori$subpop)
unique(dat.ori$subpop)

dat.ori$subpop[dat.ori$subpop=='totalcasecount'] <- 'total'

dat.inc <- dat.ori[, c('Year','subpop', 'case')] %>% spread(subpop, case)
dat.pop <- dat.ori[, c('Year','subpop', 'popsize')] %>% spread(subpop, popsize)

#visually checking if the summation equal to the total given in the dataset
data.frame(rowSums(dat.inc[,2:9]), dat.inc$total); rowSums(dat.inc[,2:9])==dat.inc$total
data.frame(rowSums(dat.pop[,2:9]), dat.pop$total); rowSums(dat.pop[,2:9])==dat.pop$total #rounding error shows difference of 1

save(dat.ori, dat.inc, dat.pop, file='./inc pop_12-2-2023.RData')

#graph showing scatter plot of population by year 
sapply(dat.pop, class)

for (i in 1:ncol(dat.pop[,-1]))
  plot(dat.pop[,1], dat.pop[,-1][,i], xlab = "Year", ylab = "Population", main=colnames(dat.pop[,-1])[i])
