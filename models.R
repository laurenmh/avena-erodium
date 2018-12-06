library(tidyverse)
library(minpack.lm)
library(nlstools)
library(grid)
library(gridExtra)

source("seeds_datasummary.R")

## FORMAT DATA FOR POPULATION MODELS ##

av1 <- togdat %>%
  filter(species == "Avena") %>%
  mutate(AVseedin = seed,
         AVseedout = seedsout) %>%
  dplyr::select(plot, subplot, density, treatment, shelterBlock, AVseedin, AVseedout) 

er1 <- togdat %>%
  filter(species == "Erodium") %>%
  mutate(ERseedin = seed,
         ERseedout = seedsout) %>%
  dplyr::select(plot, subplot, density, treatment, shelterBlock, ERseedin, ERseedout)

dat <- left_join(av1, er1) %>%
  mutate(rm = ifelse(AVseedout/AVseedin > 60 & treatment == "fallDry", 1, 0)) %>%
  # remove D3 for being unrealistically high density
  filter(density != "D3")

treatments <- unique(dat$treatment)

## Set germination and survival fractions from the literature
as <- .4
ag <- .9
es <- .82
eg <- .6

##### ERODIUM MODEL ####
## With seed bank 

# no seed bank in the experiment, so we don't want to include the seedbank term when estimating parameters
# measured the seeds produced
#m1 <- as.formula(log(ERseedout +1) ~  log(eg*(ERseedin+1)*exp(log(lambda/(1+aiE*(ERseedin+1)*eg+aiA*(AVseedin+1)*ag)))))

# same as above, just different formalization
m1E <- as.formula(log(ERseedout +1) ~  log(eg*(ERseedin+1)*exp(log(lambda)-log((1+aiE*(ERseedin+1)*eg+aiA*(AVseedin+1)*ag)))))

# old m1 version
#m1 <- as.formula(log(ERseedout +1) ~  log((ERseedin*eg +1))*((lambda)/(1+aiE*log((ERseedin*eg + 1)) + aiA*log(AVseedin*ag + 1))))

treatments <- unique(togdat$treatment)

ERoutput <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(ERoutput) = c("estimate", "se", "t", "p", "params", "treatment", "species")
for (i in 1:length(treatments)){
  m1out <- nlsLM(m1E, start=list(lambda=1, aiE = .01, aiA=.01),
                 lower = c(0, 0, 0), upper = c(200, 1, 1),
                 control=nls.lm.control(maxiter=500), trace=T,
                 data = subset(dat, !is.na(ERseedout) & treatment == treatments[i]))
  
  outreport <- as.data.frame(summary(m1out)$coef[1:3, 1:4])
  names(outreport) = c("estimate", "se", "t", "p")
  outreport$params <- row.names(outreport)
  outreport$treatment <- treatments[i]
  outreport$species <- "Erodium"
  ERoutput <- rbind(ERoutput, outreport)
}


#### AVENA MODEL ###
## With seed bank 

# new model, log transformed
m1A <- as.formula(log(AVseedout +1) ~  log(ag*(AVseedin+1)*exp(log(lambda)-log((1+aiE*(ERseedin+1)*eg+aiA*(AVseedin+1)*ag)))))


# old model
#m1 <- as.formula(log(AVseedout +1) ~  log(AVseedin*ag +1)*((lambda)/(1+aiE*log(ERseedin*eg + 1) + aiA*log(AVseedin*ag + 1))))

treatments <- unique(togdat$treatment)

AVoutput <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(AVoutput) = c("estimate", "se", "t", "p", "params", "treatment", "species")
for (i in 1:length(treatments)){
  
  m1out <- nlsLM(m1A, start=list(lambda=1, aiE = .01, aiA=.01), 
                 control=nls.lm.control(maxiter=500), 
                 lower = c(0, 0, 0), upper = c(200, 1, 1),
                 trace=T,
                 data = subset(dat, !is.na(AVseedout) & treatment == treatments[i] & rm == 0))
  outreport <- as.data.frame(summary(m1out)$coef[1:3, 1:4])
  names(outreport) = c("estimate", "se", "t", "p")
  outreport$params <- row.names(outreport)
  outreport$treatment <- treatments[i]
  outreport$species <- "Avena"
  AVoutput <- rbind(AVoutput, outreport)
}

## PUT THE TWO OUTPUTS TOGETHER
model.dat <- rbind(ERoutput, AVoutput) %>%
  tbl_df() %>%
  dplyr::select(estimate, params, treatment, species) %>%
  spread(params, estimate)

model.dat
# ggplot(model.dat, aes(x=treatment, y=aiA)) + geom_bar(stat="identity") + facet_wrap(~species)
# ggplot(model.dat, aes(x=treatment, y=aiE)) + geom_bar(stat="identity") + facet_wrap(~species)


parameter_table <- rbind(ERoutput, AVoutput) %>%
  tbl_df() %>%
  mutate_if(is.numeric, round, 2) %>%
  mutate(estimateout = paste(estimate, "±", se)) %>%
  dplyr::select(estimateout, params, treatment, species) %>%
  spread(params, estimateout) %>%
  select(treatment, species, lambda, aiA, aiE) %>%
  arrange(species, treatment)

