library(tidyverse)
library(minpack.lm)
library(nlstools)
library(grid)
library(gridExtra)


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
m1 <- as.formula(log(ERseedout +1) ~  log(eg*(ERseedin+1)*exp(log(lambda)-log((1+aiE*(ERseedin+1)*eg+aiA*(AVseedin+1)*ag)))))

# old m1 version
#m1 <- as.formula(log(ERseedout +1) ~  log((ERseedin*eg +1))*((lambda)/(1+aiE*log((ERseedin*eg + 1)) + aiA*log(AVseedin*ag + 1))))

treatments <- unique(togdat$treatment)

ERoutput <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(ERoutput) = c("estimate", "se", "t", "p", "params", "treatment", "species")
for (i in 1:length(treatments)){
  m1out <- nlsLM(m1, start=list(lambda=1, aiE = .01, aiA=.01),
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
m1 <- as.formula(log(AVseedout +1) ~  log(ag*(AVseedin+1)*exp(log(lambda)-log((1+aiE*(ERseedin+1)*eg+aiA*(AVseedin+1)*ag)))))


# old model
#m1 <- as.formula(log(AVseedout +1) ~  log(AVseedin*ag +1)*((lambda)/(1+aiE*log(ERseedin*eg + 1) + aiA*log(AVseedin*ag + 1))))

treatments <- unique(togdat$treatment)

AVoutput <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(AVoutput) = c("estimate", "se", "t", "p", "params", "treatment", "species")
for (i in 1:length(treatments)){
  
  m1out <- nlsLM(m1, start=list(lambda=1, aiE = .01, aiA=.01), 
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




### RUN THE MODELS FOR EACH SPECIES AT EQUILIBRIUM OR INVADING FOR A CONSISTENT CLIMATE CONDITION ###

### CREATE A FUNCTION THAT RUNS THE MODEL
# new growth function

growth = function(N, par.dat, t.num){
  for (i in 1:(t.num-1)){
    N$Na[i+1] = as*(1-ag)*N$Na[i]  + ag*N$Na[i]*(par.dat$Alambda[i]/(1 + par.dat$AaiE[i]*eg*N$Ne[i] + par.dat$AaiA*ag*N$Na[i]))
    
    N$Ne[i+1] = es*(1-eg)*N$Ne[i] + eg*N$Ne[i]*(par.dat$Elambda[i]/(1 + par.dat$EaiE[i]*eg*N$Ne[i] + par.dat$EaiA[i]*ag*N$Na[i]))
    
  }
  return(N)
}

### CREATE A FUNCTION THAT CALCULATES POPULATION CHANGE OVER TIME AND GRWR ####

grwr = function(par.dat, t.num) {
  N1 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N1) = c("Na", "Ne")
  N1[1,] = c(1,0)
  Aequil <- growth(N1, par.dat, t.num)
  
  N2 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N2) = c("Na", "Ne")
  N2[1,] = c(0,1)
  Eequil <- growth(N2, par.dat, t.num)
  
  N3 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N3) = c("Na", "Ne")
  N3[1,] = c(Aequil[t.num, 1],1)
  Einvade <- growth(N3, par.dat, t.num)
  Egrwr <- Einvade[2,2]/Einvade[1,2]
  Agrwc <- Einvade[2,1]/Einvade[1,1]
  Einvade$invader <- "Erodium"
  Einvade$grwr <- Egrwr
  Einvade$time <- as.numeric(row.names(Einvade))
  Einvade$Cgrwc <- Agrwc
  
  N4 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N4) = c("Na", "Ne")
  N4[1,] = c(1, Eequil[t.num, 2])
  Ainvade <- growth(N4, par.dat, t.num)
  
  Agrwr <- Ainvade[2,1]/Ainvade[1,1]
  Egrwc <- Ainvade[2,2]/Ainvade[1,2]
  
  Ainvade$invader <- "Avena"
  Ainvade$grwr <- Agrwr
  Ainvade$time <- as.numeric(row.names(Ainvade))
  Ainvade$Cgrwc <- Egrwc
  
  
  out <- rbind(Ainvade, Einvade) 
  
  return(out)
}  


### CREATE A FUNCTION THAT SETS CONSISTENT PARAMETERS ACROSS TIMESTEPS ####
consistent.par <- function(model.dat, j.num, t.num){
  par = as.data.frame(matrix(NA, nrow=t, ncol=6))
  colnames(par) = c("Alambda", "AaiE", "AaiA", "Elambda", "EaiE", "EaiA")
  trtselect = rep(j.num, t.num)
  
  for (i in 1:(t.num)){
    myparams <- subset(model.dat, treatment == treatments[trtselect[i]])
    par$Alambda[i] <- as.numeric(subset(myparams, species == "Avena")[5])
    par$AaiE[i] <- as.numeric(subset(myparams, species == "Avena")[4])
    par$AaiA[i] <- as.numeric(subset(myparams, species == "Avena")[3])
    par$Elambda[i] <- as.numeric(subset(myparams, species == "Erodium")[5])
    par$EaiE[i] <- as.numeric(subset(myparams, species == "Erodium")[4])
    par$EaiA[i] <- as.numeric(subset(myparams, species == "Erodium")[3])
  }
  return(par)
}


### Run the model for the consistent environment ###
t = 100

conDryPar <- consistent.par(model.dat, 1, t)
fallDryPar <- consistent.par(model.dat, 2, t)
springDryPar <- consistent.par(model.dat, 4, t)
conRainPar <- consistent.par(model.dat, 3, t)


conDry <- grwr(conDryPar, t) %>%
  mutate(treatment = "consistentDry")
fallDry <- grwr(fallDryPar, t) %>%
  mutate(treatment = "fallDry")
springDry <- grwr(springDryPar, t) %>%
  mutate(treatment = "springDry")
conRain <- grwr(conRainPar, t) %>%
  mutate(treatment = "controlRain")

consistent.out <- rbind(conDry, fallDry, springDry, conRain) %>%
  mutate(treatment = as.factor(treatment)) %>%
  mutate(treatment=ordered(treatment, levels = c( consistentDry="consistentDry", fallDry="fallDry",springDry="springDry", controlRain="controlRain"))) %>%
  mutate(treatment = recode(treatment, consistentDry = "Consistent dry", fallDry = "Fall dry",springDry = "Spring dry", controlRain = "Consistent wet"))

consistent.grwr.out <- consistent.out %>%
  dplyr::select(invader, grwr, treatment, Cgrwc) %>%
  unique()  %>%
  mutate(grwrChesson = log(grwr)-log(Cgrwc))

# calculate average GRWR for consistant conditions (e.g. 25% each)
GRWR.avg.avena <- mean(subset(consistent.grwr.out$grwrChesson, consistent.grwr.out$invader =="Avena"))
GRWR.avg.erodium <- mean(subset(consistent.grwr.out$grwrChesson, consistent.grwr.out$invader =="Erodium"))
# shoot, erodium doesn't coexist

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

obs.avena <- subset(consistent.grwr.out$grwrChesson, consistent.grwr.out$invader =="Avena")
GRWR.obs.avena <- consistent.dry*obs.avena[1] + fall.dry*obs.avena[2] + spring.dry*obs.avena[3] + control.rain*obs.avena[4]
obs.erodium <- subset(consistent.grwr.out$grwrChesson, consistent.grwr.out$invader =="Erodium")
GRWR.obs.erodium  <- consistent.dry*obs.erodium [1] + fall.dry*obs.erodium [2] + spring.dry*obs.erodium [3] + control.rain*obs.erodium [4]

overall.grwr <- c(GRWR.avg.avena, GRWR.avg.erodium, GRWR.obs.avena, GRWR.obs.erodium)
treatments <- c("Consistent", "Consistent", "Observed", "Observed")
species <- c("Avena", "Erodium", "Avena", "Erodium")
time.avg <- data.frame(cbind(species, treatments, overall.grwr)) %>%
  mutate(overall.grwr = as.numeric(as.character(overall.grwr)))

# ------------------------------------------------------------------------------------------------
# Lauren H, can we make this graph pretty?
# Absolutely! Although the Avena values so dwarf the Erodium ones, we might want to put it in the sup? Can't log because Erodium is less than 1 in consistent 
ggplot(time.avg, aes(x=treatments, y=overall.grwr, fill=species)) + geom_bar(stat = "identity", position = "dodge") + theme_classic() + 
  labs(x = "Long-Term Conditions", y = "Average growth rate when rare", fill = "Species") + theme(text = element_text(size = 24)) +
  geom_hline(yintercept  =0) + scale_fill_manual(values = c("grey80", "grey30")) + theme(legend.position = "none")  #+ scale_y_log10()
ggsave(here("Figs", "fig2b_avg-grwr-consistent-observed.pdf"), width = 8, height = 6)
ggsave(here("Figs", "fig2b_avg-grwr-consistent-observed.jpg"), width = 8, height = 6)

# Break it up by panels to avoid this issue
ggplot(time.avg, aes(x=treatments, y=overall.grwr, fill=species)) + geom_bar(stat = "identity", position = "dodge") + theme_classic() + 
  labs(x = "Long-term conditions", y = "Average growth rate when rare") + theme(text = element_text(size = 24)) +
  geom_hline(yintercept  =0) + scale_fill_manual(values = c("grey80", "grey30")) + theme(legend.position = "none")  + facet_wrap(~species ,scales = "free") +
  theme(strip.background = element_blank(), text = element_text(size = 16), 
        strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(here("Figs", "fig2b_avg-grwr-consistent-observed-faceted.pdf"), width = 8, height = 4)
ggsave(here("Figs", "fig2b_avg-grwr-consistent-observed-faceted.jpg"), width = 8, height = 4)


# ------------------------------------------------------------------------------------------------

consistent.out2 <- consistent.out %>%
  gather(species, count, Na:Ne) %>%
  tbl_df() %>%
  mutate(count = ifelse(count <= 0.1, 0.1, count))

## Simulation graph ##
#pdf("Consistentsimulation.pdf", width = 12, height = 6)
ggplot(consistent.out2, aes(x=time, y=(count), color = species)) + geom_line(size = 4) + facet_grid(invader ~ treatment) + 
  theme_bw() +  theme(text = element_text(size = 24), legend.position = "none", strip.background = element_blank(),
                      #strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16),
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_log10(limits=c(.1, 1200), breaks = c( 1, 10, 100, 1000)) +
  labs(y=expression(paste("Count (individuals/m"^"2",")")), x = "Time step") + scale_color_manual(values = c("grey80", "grey30"))
ggsave(here("Figs", "fig3_consistent-simulation.pdf"), width = 12, height = 8)
ggsave(here("Figs", "fig3_consistent-simulation.jpg"), width = 12, height = 8)


#pdf("ddGRWRfromModel.pdf", width = 8, height = 6)
ggplot(consistent.grwr.out, aes(x=treatment, y=grwrChesson, fill=invader)) + geom_bar(stat = "identity", position = "dodge") + theme_classic() + 
  labs(x = "Rainfall treatment", y = expression("Growth rate when rare (r"[i]* "-r"[r]*")"), fill = "Species") + theme(text = element_text(size = 24)) +
  geom_hline(yintercept  =0) + scale_fill_manual(values = c("grey80", "grey30")) + theme(legend.position = "none")
ggsave(here("Figs", "fig2_GRWR-from-model.pdf"), width = 8, height = 6)
ggsave(here("Figs", "fig2_GRWR-from-model.jpg"), width = 8, height = 6)


# ## CHECK MODEL FITS
# model.resid <- nlsResiduals(m1out)
# plot(model.resid)
# plot(model.resid, which=5)
# 
# predicted <- m1out$m$predict()
# actual <- log(subset(dat, !is.na(ERseedout) & treatment == treatments[i])$ERseedout)
# 
# pa <- cbind(predicted, actual) %>%
#   data.frame() %>%
#   tbl_df() 
# 
# plot(log(pa$actual), (pa$predicted))
# 
# l <- lm( pa$predicted ~ log(pa$actual +1))
# abline(l)
# 
# m1.resid <- nlsResiduals(m1out)
# plot(m1.resid)
# plot(m1.resid, which=5)
# plot(m1.resid, which=6)
# 
# preview(m1, dat, start=list(lambda=1, aii=.01, p1=.1))
# preview(m1, dat, start=list(lambda=1, aii=.01, p1=.1), variable = 13)
# preview(m1, dat, start=list(lambda=1, aii=.01, p1=.1), variable = 8)


### CREATE A FUNCTION THAT CALCULATES POPULATION CHANGE OVER TIME WHEN A SPECIES IS ALONE ####

grwr_nocomp = function(par.dat, t.num) {
  N1 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N1) = c("Na", "Ne")
  N1[1,] = c(1,0)
  Ainvade <- growth(N1, par.dat, t.num)
  Agrwr <- Ainvade[2,1]/Ainvade[1,1]
  Ainvade$invader <- "Avena"
  Ainvade$grwr <- Agrwr
  Ainvade$time <- as.numeric(row.names(Ainvade))
  
  
  N2 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N2) = c("Na", "Ne")
  N2[1,] = c(0,1)
  Einvade <- growth(N2, par.dat, t.num)
  Egrwr <- Einvade[2,2]/Einvade[1,2]
  Einvade$invader <- "Erodium"
  Einvade$grwr <- Egrwr
  Einvade$time <- as.numeric(row.names(Einvade))
  
  
  out <- rbind(Ainvade, Einvade) 
  
  return(out)
}  

### Run the model for each species grown alone in the consistent environment ###
t = 100
 
# Set parameters
conDryPar <- consistent.par(model.dat, 1, t)
fallDryPar <- consistent.par(model.dat, 2, t)
springDryPar <- consistent.par(model.dat, 4, t)
conRainPar <- consistent.par(model.dat, 3, t)

# Run for each scenario
conDry_GRWR_nocomp <- grwr_nocomp(conDryPar, t) %>%
  mutate(treatment = "consistentDry")
fallDry_GRWR_nocomp <- grwr_nocomp(fallDryPar, t)  %>%
  mutate(treatment = "fallDry")
springDry_GRWR_nocomp <- grwr_nocomp(springDryPar, t)  %>%
  mutate(treatment = "springDry")
conRain_GRWR_nocomp <- grwr_nocomp(conRainPar, t)  %>%
  mutate(treatment = "controlRain")

GRWR_nocomp <- rbind(conDry_GRWR_nocomp, fallDry_GRWR_nocomp, springDry_GRWR_nocomp, conRain_GRWR_nocomp) %>%
  mutate(grwalone = grwr) %>%
  dplyr::select(invader, grwalone, treatment) %>%
  unique()

# ggplot(GRWR_nocomp, aes(x=treatment, y=grwalone, color = invader)) + geom_point()

GRWR_withcomp <- rbind(conDry, fallDry, springDry, conRain) %>%
  dplyr::select(invader, grwr, treatment) %>%
  unique() 

# Combine with GRWR with a competitor; compare growth rates
GRWRtog <- left_join(GRWR_nocomp, GRWR_withcomp) %>%
  mutate(species = invader) %>%
  mutate(compeffect = log(grwalone/grwr))

# Visualize simulated lambda (ie, growth rate when alone and rare) with competitive effect
ggplot(GRWRtog, aes(y=log(grwalone), x=compeffect, color = invader)) + geom_point(size = 4) + geom_smooth(method = "lm", se =F) + 
  labs(y= "Environment", x = "Competition", color = "Species") +
  theme_classic() + theme(text = element_text(size = 20))


# Put it together with fited lambda for covariances
outdat <- left_join(GRWRtog, model.dat) %>%
  group_by(species) %>%
  mutate(Ravg = mean(grwalone),
         Ejx = grwalone/Ravg,
         Cjx = compeffect/Ravg,
         Ejx2 = log(grwalone))

GRWRav <- outdat %>%
  filter(invader == "Avena")

cov(GRWRav$grwalone, GRWRav$grwr)
cov(GRWRav$Ejx2, GRWRav$compeffect)


GRWRero <- outdat %>%
  filter(invader == "Erodium")

cov(GRWRero$grwalone, GRWRero$grwr)
cov(GRWRero$grwalone, GRWRero$compeffect)
cov(GRWRero$Ejx2, GRWRero$compeffect)
cov(GRWRero$Ejx, GRWRero$compeffect)



ggplot(outdat, aes(y=Ejx, x=compeffect, color = species)) + geom_point(pch =0, size = 4 ) + geom_smooth(method = "lm", se = F) + theme_classic()

ggplot(outdat, aes(y=Ejx, x=compeffect, color = species)) + geom_point() 


#ggplot(outdat, aes(x=lambda, y=grwalone)) + geom_point() + xlim(0,10)
