# Lisbon workshop May 2022
# example: Joint longitudinal - Survival

library(INLA)
set.seed(1)
# data generation - one longitudinal marker
nsujet=500 # number of individuals
# L1 (longitudinal continuous)
b_0=0.1 # intercept
b_1=-0.1 # slope
b_e=0.1 # residual error
# S1 (survival)
phi_b0=1 # random intercept association with survival
gap=1 # gap between measurements
followup=5 # follow-up time
mestime=seq(0,followup,gap) # measurement times
time=rep(mestime, nsujet) # time column
nmesindiv=followup/gap+1 # number of individual measurements
nmesy= nmesindiv*nsujet # total number of measurements
id<-rep(1:nsujet, each=nmesindiv) # individual id
# random effect
b_i <- rnorm(nsujet, 0, 0.1)
b_i0 <- rep(b_i, each=nmesindiv) # random intercept
# linear predictor
linPredL1 <- b_i0 + b_0 + b_1 * time
# continuous outcome L1
L1 <- rnorm(nmesy, linPredL1, b_e)
lon <- data.frame(L1, id, time)
## generation of exponential death times
u <- runif(nsujet) # uniform distribution for survival times generation
baseScale=0.1
deathTimes <- -(log(u) / (baseScale * exp(b_i*phi_b0)))
d <- as.numeric(deathTimes<followup) # deathtimes indicator
## censoring individuals at end of follow-up (not at random)
deathTimes[deathTimes>=followup]=followup
surv <- data.frame(id=1:nsujet,deathTimes, d) # survival times dataset
## removing longi measurements after death
ind <- rep(NA, nsujet*length(mestime))
for (i in 1:nsujet){
  for(j in 1:length(mestime)){
    if(lon[(i-1)*length(mestime)+j, "time"]<=surv[i,"deathTimes"]) ind[(i-1)*length(mestime)+j]=1
  }
}
lon <- lon[!is.na(ind),]
summary(surv)
summary(lon)










