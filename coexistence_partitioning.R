# Create figure 2b -- low density growth rates for observed environment

# load the data
load("model.dat.output.RData")

# First determine how common each environmental type is

# what about for what we actually see in terms of the number of years in each env. condition
# first read in the data
rain <- read_csv("Data/PRISM_brownsvalley_long.csv", skip = 10) %>%
  mutate(ppt = `ppt (inches)`*2.54*10) %>%
  separate(Date, c("year", "month")) %>%
  mutate(year = as.numeric(year),
         month = as.numeric(month)) %>%
  mutate(year = ifelse(month == 12 | month == 11 | month == 10 | month == 9, year + 1, year)) %>%
  mutate(season = "Early",
         season = ifelse(month == 2 | month == 3 | month == 4, "Late", season)) %>%
  filter(month != 5, month != 6, month!= 7, month != 8)

## Summarize by year 
## Using 50% as the cutoff 
rainsummary <-  rain %>%
  group_by(year, season) %>%
  summarize(ppt = sum(ppt)) %>%
  spread(season, ppt) %>%
  mutate(Total = Early + Late) 

rainsummary <- rainsummary %>%
  mutate(raintype = "controlRain",
         raintype = ifelse(Early < quantile(rainsummary$Early, .5), "fallDry", raintype),
         raintype = ifelse(Late < quantile(rainsummary$Late, .5), "springDry", raintype),
         raintype = ifelse(Total < quantile(rainsummary$Total, .5), "consistentDry", raintype)) 

fall.dry <- length(which(rainsummary$raintype == "fallDry")) / nrow(rainsummary)
spring.dry <- length(which(rainsummary$raintype == "springDry")) / nrow(rainsummary)
consistent.dry <- length(which(rainsummary$raintype == "consistentDry")) / nrow(rainsummary)
control.rain <- length(which(rainsummary$raintype == "controlRain")) / nrow(rainsummary)

# Determine equilibrium conditions for each species in isolation 
pop.equilibrium <- function (N0, s, g, a_intra, lambda) {
  # to run for only a single timestep
    N <- s*(1-g)*N0 + N0*(lambda*g)/(1+a_intra*N0)
  return(N)
}

avena <- subset(model.dat, species=="Avena")
erodium <- subset(model.dat, species=="Erodium")

## Set germination and survival fractions from the literature
as <- .4
ag <- .9
es <- .82
eg <- .6

# use the timeseries of environmental conditions for environmental variability
# for avena
N0 <- 550
time <- length(rainsummary$raintype)
N_avena <- rep(NA, time)
N_avena[1] <- N0

for (t in 1:time) {
  params <- subset(avena, treatment==rainsummary$raintype[t])
  N_avena[t+1] <- pop.equilibrium(N0=N_avena[t], s=as, g=ag, a_intra=params$aiA, lambda=params$lambda)
}

# check output
plot(seq(1:(time+1)), N_avena, type="l")


# for erodium
N0 <- 70
time <- length(rainsummary$raintype)
N_erodium <- rep(NA, time)
N_erodium[1] <- N0

for (t in 1:time) {
  params <- subset(erodium, treatment==rainsummary$raintype[t])
  N_erodium[t+1] <- pop.equilibrium(N0=N_erodium[t], s=es, g=eg, a_intra=params$aiE, lambda=params$lambda)
}

# check output
plot(seq(1:(time+1)), N_erodium, type="l")

# say that resident community has reached equilibrium after a 50 year "burn in" period
# and invade each species into the community one at a time

pop.invade <- function (N0, resident, s, g, a_inter, lambda) {
  # to run for only a single timestep
  N <- s*(1-g)*N0 + N0*(lambda*g)/(1+a_inter*resident)
  return(N)
}

pop.resident <- function (N0, resident, s, g, a_intra, a_inter, lambda) {
  # to run for only a single timestep
  N <- s*(1-g)*resident + resident*(lambda*g)/(1+a_intra*resident+a_inter*N0)
  return(N)
}

# invade avena first
avena_invade <- rep (NA, 72)
erodium_resident <- rep (NA, 72)
temp <- 1
for (t in 50:time) {
  params <- subset(avena, treatment==rainsummary$raintype[t])
  params_resident <- subset(erodium, treatment==rainsummary$raintype[t])
  avena_invade[temp] <- pop.invade(N0=1, resident=N_erodium[t], s=as, g=ag, a_inter=params$aiE, lambda=params$lambda)
  
  # sanity check that the resident isn't affected
  erodium_new <- pop.resident(N0=1, resident=N_erodium[t], s=es, g=eg, 
                                         a_intra=params_resident$aiE, a_inter=params_resident$aiA, 
                                         lambda=params_resident$lambda)
  erodium_resident[temp] <- erodium_new/N_erodium[t]
  
  temp  <- temp + 1 
}

# then have erodium invade into avena
erodium_invade <- rep (NA, 72)
avena_resident <- rep (NA, 72)
temp <- 1
for (t in 50:time) {
  params <- subset(erodium, treatment==rainsummary$raintype[t])
  params_resident <- subset(avena, treatment==rainsummary$raintype[t])
  
  erodium_invade[temp] <- pop.invade(N0=1, resident=N_avena[t], s=es, g=eg, a_inter=params$aiA, lambda=params$lambda)
 
   # sanity check that the resident isn't affected
  avena_new <- pop.resident(N0=1, resident=N_avena[t], s=as, g=ag, 
                              a_intra=params_resident$aiA, a_inter=params_resident$aiE, 
                              lambda=params_resident$lambda)
  avena_resident[temp] <- avena_new/N_avena[t]
  temp  <- temp + 1 
}

# ------------------------------------------------------------------------------------
# These are the values for Lauren Hallett to use for figure 2b
avena_invader <- log(avena_invade)
erodium_invader <- log(erodium_invade)
