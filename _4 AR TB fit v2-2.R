
#### 3. use the functions from StMoMo to fit the data and to generate the predictions for future years

rm(list=ls())

library(tidyr)
library(ggplot2)

setwd("C:/Users/dbboothe/Desktop/TBdatafor dissertation/Aim1code")

load('./inc pop_12-2-2023.RData')

####
inc <- as.matrix(dat.inc[,c(-1, -ncol(dat.inc))]); inc[which(inc==0, arr.ind=TRUE)] <- 0.001
colnames(inc) <- colnames(inc)
rownames(inc) <- dat.inc[,1]

ages <- 1:ncol(inc)
years <- as.integer(rownames(inc))

pop <- as.matrix(dat.pop[,c(-1, -ncol(dat.pop))])
colnames(pop) <- ages #make the names as integers for functions from the R packages to run with
rownames(pop) <- years 

library(demography)
library(StMoMo)

dat <- demogdata(t(inc/pop), t(pop), ages=ages, years=years, 
                 type="mortality", label='US', name='total') 

dat <- StMoMoData(dat) 

model <- fit(lc(), dat, verbose=F) #fit LC model on the data of 2009-2021
model #names(model) will provide details including Deviance, loglik things and etc.

plot(model)

#### forecast_rates func
forecast_rates <- function(model, until=2035, n_sim=1000) {
  h <- until - max(model$years)
  m <- model
  
  dat <- t(m$data$Dxt/m$data$Ext) #observed TB rate for the bynow years
  est <- exp(t(fitted(m))) #fitted TB rate for the bynow years
  
  fore <- forecast(m, h=h, kt.method="iarima") #need to study kt.order=order_f 
  #fore <- forecast(m, h=h, kt.method="mrwd")
  ktf <- c(fore$kt.f$mean)                     #, kt.order=c(1,1,1)
  fore <- t(fore$rates) #this shows the TB reates of risk groups for all future year
  colnames(fore) <- colnames(dat) <- colnames(est) <- m$ages
  
  sim <- simulate(m, nsim=n_sim, h=h, kt.method="iarima") #, kt.order=order_f, here the func is simulate.fitStMoMo
  #sim <- simulate(m, nsim=n_sim, h=h, kt.method="mrwd")
  kts <- sim$kt.s$sim[1,,]                              #, kt.order=c(1,1,1)    
  sim <- lapply(1:n_sim, function(x) t(sim$rates[,,x])) #simulated TB rates for future years
  
  res <- list(
    Data=dat, Forecast=fore, Fitted=est, 
    Boot=sim, kt.fore=ktf, kt.sim=kts)
  res
  
  #Data: observed TB rates; Forecast: predicted TB rates; Fitted: fitted TB rates;
  #Boot: n_sim bootstrap samples; kt.fore: kt paras estimated for future yrs; 
  #kt.sim: kt paras estimated for bt samples
}

# Forecast future incidence rates
fore <- forecast_rates(model, 2035, n_sim=1000) 


#### 
load('./dat.fupop.RData')

simulate_incidence <- function(fore, pop.proj) {
  nms <- dimnames(fore$Forecast)
  
  mc <- list()
  
  for (i in 1:length(fore$Boot)) {
    p <- pop.proj #future pop 
    r <- fore$Boot[[i]] #TB rate
    n <- rpois(nrow(r)*ncol(r), as.matrix(p*r)) #TB numbers wi randomness added from Pois
    n <- matrix(n, nrow(r), ncol(r))
    dimnames(p) <- dimnames(n) <- nms
    mc[[i]] <- list(p=p, n=n)
  }
  
  dimnames(pop.proj) <- nms
  
  list(
    Data=list(p=pop, n=fore$Data*pop),
    Fitted=list(p=pop, n=fore$Fitted*pop),
    Forecast=list(p=pop.proj, n=fore$Forecast*pop.proj),
    Boot=mc
  )
  
  #Data: the dataset having both observed population and TB case numbers 
  #Fitted: the dataset having both esimated population and TB case numbers
  #Forecast: the dataset having the projected population and predicated TB case numbers
  #Boot: each element has a dataset containing both bootstrp population and TB case numbers
}

# Sample future incidence cases based on rates and population sizes
sims <- simulate_incidence(fore, dat.fupop[,c(-1, -ncol(dat.fupop))])

summary_with_age <- function(inc, agp, agl) {
  n.total <- inc$n
  p.total <- inc$p
  
  n.sums <- rowSums(n.total)
  p.sums <- rowSums(p.total)
  
  res <- data.frame(Year=as.numeric(rownames(n.total)), 
                    TotalN=n.sums, TotalR=n.sums/p.sums*1e5)
  
  
  #agp <- list(Chd=1, Young=2:3, Mid=4:5, Old=6:7, nonusborn=8)
  #agl <- c(Chd="0-17", Young="18-44", Mid="45-64", Old="65+", nonusborn='nonusborn')
  
  year.end <- 2035
  
  for (a in names(agp)) {
    ind <- agp[[a]]
    
    if (length(ind)>1)
    {
      n <- rowSums(n.total[,ind])
      pop <- rowSums(p.total[,ind])
      res[paste0(a, "N")] <- n
      res[paste0(a, "R")] <- (n/pop) * 1e5
      res[paste0(a, "Pr")] <- n/n.sums  #proportion of TB cases of this group out of total TB cases in a yr
    } else
    {
      n <- n.total[,ind]
      pop <- p.total[,ind]
      res[paste0(a, "N")] <- n
      res[paste0(a, "R")] <- (n/pop) * 1e5
      res[paste0(a, "Pr")] <- n/n.sums  #proportion of TB cases of this group out of total TB cases in a yr
    }}
    
  rownames(res) <- NULL
  res
}

agp <- list(Agegrp1=1, Agegrp2=2, Agegrp3=3, Agegrp4=4,
            Agegrp5=5, Agegrp6=6, Agegrp7=7, Nonusborn=8)

agl <- c(Agegrp1="0-17", Agegrp2="18-24", Agegrp3="25-44", 
         Agegrp4="45-54", Agegrp5="55-64", Agegrp6="65-74",
         Agegrp7="75+", nonusborn='nonusborn')

est.agegrp <- summary_with_age(inc=sims$Forecast, agp, agl)

summarise_all_age <- function(sims, agp, agl, q=0.95) {
  fcst <- summary_with_age(sims$Forecast, agp)
  
  mc <- array(0, c(dim(fcst), length(sims$Boot)))
  
  for (i in 1:length(sims$Boot)) {
    summ <- summary_with_age(sims$Boot[[i]], agp)
    mc[,, i] <- as.matrix(summ)
  }
  
  dimnames(mc)[[1]] <- rownames(fcst)
  dimnames(mc)[[2]] <- colnames(fcst)
  
  mc.summary <- array(0, c(dim(fcst), 3))
  mc.summary[,, 1] <- apply(mc, c(1, 2), mean)
  mc.summary[,, 2] <- apply(mc, c(1, 2), function(x) quantile(x, .5-q/2))
  mc.summary[,, 3] <- apply(mc, c(1, 2), function(x) quantile(x, .5+q/2))
  
  dimnames(mc.summary)[[1]] <- rownames(fcst)
  dimnames(mc.summary)[[2]] <- colnames(fcst)
  dimnames(mc.summary)[[3]] <- c("Mean", "Lower", "Upper")
  
  list(
    Data=summary_with_age(sims$Data, agp),
    Fitted=summary_with_age(sims$Fitted, agp),
    Forecast=fcst,
    Boot=mc.summary
  )
  
  #Data: caclulated rate by group from the bynow observed data
  #Fitted: fitted rate by group form the bynow data
  #Forecast: 
}

summ.agegrp <- summarise_all_age(sims, agp, agl) 
#this summary is by risk groups (originally by age groups, that's how the functions are named)

year <- 2035

summ <- summ.agegrp

dat <- summ$Data[c("Year", "TotalR")]
dat <- subset(dat, Year <= year)
ftt <- summ$Fitted[c("Year", "TotalR")]
ftt <- subset(ftt, Year <= year)
fcst <- data.frame(Year=summ$Forecast$Year, summ$Boot[,"TotalR",])
fcst <- subset(fcst, Year <= year)


endTB <- data.frame(Year=c(2020, 2025, 2030, 2035),
                    Roles=c("Milestone", "Milestone", "Goal(SDG)", "Goal"),
                    Lab=c("20%", "50%", "80%", "90%"),
                    Noti=(1 - c(0.2, 0.5, 0.8, 0.9)) * 45.65)


gnr <- ggplot(data=dat, aes(x=Year)) +
  geom_point(aes(y=TotalR, shape="Data")) +
  geom_line(aes(y=Mean, linetype="Forecast"), data=fcst) +
  geom_ribbon(data=fcst, aes(ymin=Lower, ymax=Upper), alpha=0.3) +
  geom_line(aes(y=TotalR, linetype="Fitted"), data=ftt) +
  #geom_point(data=endTB, aes(y=Noti, shape="End TB Strategy"), size=2) +
  #geom_point(aes(y=1, shape="Reduction"), size=0) +
  #geom_text(data=endTB, aes(y=Noti, label=Lab, vjust=1.5, hjust=1)) +
  scale_linetype_manual(element_blank(), values=c("Fitted"=1, "Forecast"=2), drop=F) +
  #scale_shape_manual(element_blank(), values=c(Data=16, "End TB Strategy"=2, "Reduction"=17), drop=F) +
  guides(linetype=guide_legend(order=1), shape=guide_legend(order=2), size="none") +
  scale_y_continuous("Incidence rate per 100,000", limits=c(1.5, 3.75)) +
  scale_x_continuous("Year", breaks=seq(2005, year, by=5)) +
  theme_minimal() #+ 
#theme(legend.position="bottom")

gnr

###
dat <- summ$Data[c("Year", "Agegrp1R")]
dat <- subset(dat, Year <= year)
ftt <- summ$Fitted[c("Year", "Agegrp1R")]
ftt <- subset(ftt, Year <= year)
fcst <- data.frame(Year=summ$Forecast$Year, summ$Boot[,"Agegrp1R",])
fcst <- subset(fcst, Year <= year)

gnr <- ggplot(data=dat, aes(x=Year)) +
  geom_point(aes(y=Agegrp1R, shape="Data")) +
  geom_line(aes(y=Mean, linetype="Forecast"), data=fcst) +
  geom_ribbon(data=fcst, aes(ymin=Lower, ymax=Upper), alpha=0.3) +
  geom_line(aes(y=Agegrp1R, linetype="Fitted"), data=ftt) +
  scale_linetype_manual(element_blank(), values=c("Fitted"=1, "Forecast"=2), drop=F) +
  guides(linetype=guide_legend(order=1), shape=guide_legend(order=2), size="none") +
  scale_y_continuous("Incidence rate per 100,000", limits=c(0, 7)) +
  scale_x_continuous("Year", breaks=seq(2005, year, by=5)) +
  theme_minimal()

gnr

### age group 2
dat <- summ$Data[c("Year", "Agegrp2R")]
dat <- subset(dat, Year <= year)
ftt <- summ$Fitted[c("Year", "Agegrp2R")]
ftt <- subset(ftt, Year <= year)
fcst <- data.frame(Year=summ$Forecast$Year, summ$Boot[,"Agegrp2R",])
fcst <- subset(fcst, Year <= year)

gnr <- ggplot(data=dat, aes(x=Year)) +
  geom_point(aes(y=Agegrp2R, shape="Data")) +
  geom_line(aes(y=Mean, linetype="Forecast"), data=fcst) +
  geom_ribbon(data=fcst, aes(ymin=Lower, ymax=Upper), alpha=0.3) +
  geom_line(aes(y=Agegrp2R, linetype="Fitted"), data=ftt) +
  scale_linetype_manual(element_blank(), values=c("Fitted"=1, "Forecast"=2), drop=F) +
  guides(linetype=guide_legend(order=1), shape=guide_legend(order=2), size="none") +
  scale_y_continuous("Incidence rate per 100,000", limits=c(-1, 3)) +
  scale_x_continuous("Year", breaks=seq(2005, year, by=5)) +
  theme_minimal()

gnr

### age group 3
dat <- summ$Data[c("Year", "Agegrp3R")]
dat <- subset(dat, Year <= year)
ftt <- summ$Fitted[c("Year", "Agegrp3R")]
ftt <- subset(ftt, Year <= year)
fcst <- data.frame(Year=summ$Forecast$Year, summ$Boot[,"Agegrp3R",])
fcst <- subset(fcst, Year <= year)

gnr <- ggplot(data=dat, aes(x=Year)) +
  geom_point(aes(y=Agegrp3R, shape="Data")) +
  geom_line(aes(y=Mean, linetype="Forecast"), data=fcst) +
  geom_ribbon(data=fcst, aes(ymin=Lower, ymax=Upper), alpha=0.3) +
  geom_line(aes(y=Agegrp3R, linetype="Fitted"), data=ftt) +
  scale_linetype_manual(element_blank(), values=c("Fitted"=1, "Forecast"=2), drop=F) +
  guides(linetype=guide_legend(order=1), shape=guide_legend(order=2), size="none") +
  scale_y_continuous("Incidence rate per 100,000", limits=c(0, 3)) +
  scale_x_continuous("Year", breaks=seq(2005, year, by=5)) +
  theme_minimal()

gnr

### age group 4
dat <- summ$Data[c("Year", "Agegrp4R")]
dat <- subset(dat, Year <= year)
ftt <- summ$Fitted[c("Year", "Agegrp4R")]
ftt <- subset(ftt, Year <= year)
fcst <- data.frame(Year=summ$Forecast$Year, summ$Boot[,"Agegrp4R",])
fcst <- subset(fcst, Year <= year)

gnr <- ggplot(data=dat, aes(x=Year)) +
  geom_point(aes(y=Agegrp4R, shape="Data")) +
  geom_line(aes(y=Mean, linetype="Forecast"), data=fcst) +
  geom_ribbon(data=fcst, aes(ymin=Lower, ymax=Upper), alpha=0.3) +
  geom_line(aes(y=Agegrp4R, linetype="Fitted"), data=ftt) +
  scale_linetype_manual(element_blank(), values=c("Fitted"=1, "Forecast"=2), drop=F) +
  guides(linetype=guide_legend(order=1), shape=guide_legend(order=2), size="none") +
  scale_y_continuous("Incidence rate per 100,000", limits=c(-1, 4)) +
  scale_x_continuous("Year", breaks=seq(2005, year, by=5)) +
  theme_minimal()

gnr

### age group 5
dat <- summ$Data[c("Year", "Agegrp5R")]
dat <- subset(dat, Year <= year)
ftt <- summ$Fitted[c("Year", "Agegrp5R")]
ftt <- subset(ftt, Year <= year)
fcst <- data.frame(Year=summ$Forecast$Year, summ$Boot[,"Agegrp5R",])
fcst <- subset(fcst, Year <= year)

gnr <- ggplot(data=dat, aes(x=Year)) +
  geom_point(aes(y=Agegrp5R, shape="Data")) +
  geom_line(aes(y=Mean, linetype="Forecast"), data=fcst) +
  geom_ribbon(data=fcst, aes(ymin=Lower, ymax=Upper), alpha=0.3) +
  geom_line(aes(y=Agegrp5R, linetype="Fitted"), data=ftt) +
  scale_linetype_manual(element_blank(), values=c("Fitted"=1, "Forecast"=2), drop=F) +
  guides(linetype=guide_legend(order=1), shape=guide_legend(order=2), size="none") +
  scale_y_continuous("Incidence rate per 100,000", limits=c(0, 5)) +
  scale_x_continuous("Year", breaks=seq(2005, year, by=5)) +
  theme_minimal()

gnr

### age group 6
dat <- summ$Data[c("Year", "Agegrp6R")]
dat <- subset(dat, Year <= year)
ftt <- summ$Fitted[c("Year", "Agegrp6R")]
ftt <- subset(ftt, Year <= year)
fcst <- data.frame(Year=summ$Forecast$Year, summ$Boot[,"Agegrp6R",])
fcst <- subset(fcst, Year <= year)

gnr <- ggplot(data=dat, aes(x=Year)) +
  geom_point(aes(y=Agegrp6R, shape="Data")) +
  geom_line(aes(y=Mean, linetype="Forecast"), data=fcst) +
  geom_ribbon(data=fcst, aes(ymin=Lower, ymax=Upper), alpha=0.3) +
  geom_line(aes(y=Agegrp6R, linetype="Fitted"), data=ftt) +
  scale_linetype_manual(element_blank(), values=c("Fitted"=1, "Forecast"=2), drop=F) +
  guides(linetype=guide_legend(order=1), shape=guide_legend(order=2), size="none") +
  scale_y_continuous("Incidence rate per 100,000", limits=c(-1, 5)) +
  scale_x_continuous("Year", breaks=seq(2005, year, by=5)) +
  theme_minimal()

gnr

### age group 7
dat <- summ$Data[c("Year", "Agegrp7R")]
dat <- subset(dat, Year <= year)
ftt <- summ$Fitted[c("Year", "Agegrp7R")]
ftt <- subset(ftt, Year <= year)
fcst <- data.frame(Year=summ$Forecast$Year, summ$Boot[,"Agegrp7R",])
fcst <- subset(fcst, Year <= year)

gnr <- ggplot(data=dat, aes(x=Year)) +
  geom_point(aes(y=Agegrp7R, shape="Data")) +
  geom_line(aes(y=Mean, linetype="Forecast"), data=fcst) +
  geom_ribbon(data=fcst, aes(ymin=Lower, ymax=Upper), alpha=0.3) +
  geom_line(aes(y=Agegrp7R, linetype="Fitted"), data=ftt) +
  scale_linetype_manual(element_blank(), values=c("Fitted"=1, "Forecast"=2), drop=F) +
  guides(linetype=guide_legend(order=1), shape=guide_legend(order=2), size="none") +
  scale_y_continuous("Incidence rate per 100,000", limits=c(0, 11)) +
  scale_x_continuous("Year", breaks=seq(2005, year, by=5)) +
  theme_minimal()

gnr

### nonusborn
dat <- summ$Data[c("Year", "NonusbornR")]
dat <- subset(dat, Year <= year)
ftt <- summ$Fitted[c("Year", "NonusbornR")]
ftt <- subset(ftt, Year <= year)
fcst <- data.frame(Year=summ$Forecast$Year, summ$Boot[,"NonusbornR",])
fcst <- subset(fcst, Year <= year)

gnr <- ggplot(data=dat, aes(x=Year)) +
  geom_point(aes(y=NonusbornR, shape="Data")) +
  geom_line(aes(y=Mean, linetype="Forecast"), data=fcst) +
  geom_ribbon(data=fcst, aes(ymin=Lower, ymax=Upper), alpha=0.3) +
  geom_line(aes(y=NonusbornR, linetype="Fitted"), data=ftt) +
  scale_linetype_manual(element_blank(), values=c("Fitted"=1, "Forecast"=2), drop=F) +
  guides(linetype=guide_legend(order=1), shape=guide_legend(order=2), size="none") +
  scale_y_continuous("Incidence rate per 100,000", limits=c(10, 35)) +
  scale_x_continuous("Year", breaks=seq(2005, year, by=5)) +
  theme_minimal()

gnr
