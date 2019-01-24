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

# model
m1 <- as.formula(log(ERseedout +1) ~  log(eg*(ERseedin+1)*exp(log(lambda)-log((1+aiE*(ERseedin+1)*eg+aiA*(AVseedin+1)*ag)))))


m1out <- nlsLM(m1, start=list(lambda=1, aiE = .01, aiA=.01),
                 lower = c(0, 0, 0), upper = c(30, 1, 1),
                 control=nls.lm.control(maxiter=500), trace=T,
                 data = subset(dat, !is.na(ERseedout) & treatment == "fallDry"))



m1out_big <- nlsLM(m1, start=list(lambda=158, aiE = .01, aiA=.01),
                lower = c(0, 0, 0), upper = c(158, 1, 1),
                control=nls.lm.control(maxiter=500), trace=T,
                data = subset(dat, !is.na(ERseedout) & treatment == "fallDry"))
summary(m1out_big)
# logLik(m1out_big)
# 
# lrtest( m1out_big, m1out)

# ## attempts at nls_multstart 
# library(nls.multstart)
# 
# growthfunc <- function(lambda, aiE, aiA, ERseedin, AVseedin){
#   as <- .4
#   ag <- .9
#   es <- .82
#   eg <- .6
#   newcol <- log(eg*(ERseedin+1)*exp(log(lambda)-log((1+aiE*(ERseedin+1)*eg+aiA*(AVseedin+1)*ag))))
#   return(newcol)
# }
# 
# 
# dat$newcol <- log(dat$ERseedout + 1)
# subdat <- subset(dat, !is.na(ERseedout) & treatment == "fallDry")
# 
# a <- nls_multstart(newcol ~ growthfunc(lambda, aiE, aiA, ERseedin, AVseedin), 
#                    data = subdat,
#                    iter = 500,
#                    start_lower=c(lambda=0, aiE = 0, aiA=0),
#                    start_upper = c(lambda=200, aiE = 1, aiA=1),
#                    lower = c(lambda=0, aiE = 0, aiA=0),
#                    upper = c(lambda=200, aiE =2, aiA=1))
# summary(a)
# 
