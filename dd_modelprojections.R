library(tidyverse)
library(minpack.lm)
library(nlstools)
library(grid)
library(gridExtra)


## FORMAT DATA FOR POPULATION MODELS
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
  mutate(rm = ifelse(AVseedout/AVseedin > 60 & treatment == "fallDry", 1, 0))

treatments <- unique(dat$treatment)

## Set germination and survival fractions from the literature
as <- .4
ag <- .9
es <- .82
eg <- .6

##### ERODIUM MODELS ####

## With seed bank ##
m1 <- as.formula(log(ERseedout +1) ~  log((ERseedin*eg +1))*((lambda)/(1+aiE*.5*log((ERseedin*eg + 1)) + aiA*log(AVseedin*ag + 1))))

treatments <- unique(togdat$treatment)

ERoutput <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(ERoutput) = c("estimate", "se", "t", "p", "params", "treatment", "species")
for (i in 1:length(treatments)){
  m1out <- nlsLM(m1, start=list(lambda=1, aiE = .01, aiA=.01),
                 control=nls.lm.control(maxiter=500), trace=T,
                 data = subset(dat, !is.na(ERseedout) & treatment == treatments[i]))
  
  outreport <- as.data.frame(summary(m1out)$coef[1:3, 1:4])
  names(outreport) = c("estimate", "se", "t", "p")
  outreport$params <- row.names(outreport)
  outreport$treatment <- treatments[i]
  outreport$species <- "Erodium"
  ERoutput <- rbind(ERoutput, outreport)
}


#### AVENA MODELS ###

#### AVENA MODELS with seed bank ###

m1 <- as.formula(log(AVseedout +1) ~  log(AVseedin*ag +1)*((lambda)/(1+aiE*log(ERseedin*eg + 1) + aiA*log(AVseedin*ag + 1))))

treatments <- unique(togdat$treatment)

AVoutput <- as.data.frame(matrix(nrow = 0, ncol = 7))
names(AVoutput) = c("estimate", "se", "t", "p", "params", "treatment", "species")
for (i in 1:length(treatments)){
  
  m1out <- nlsLM(m1, start=list(lambda=1, aiE = .01, aiA=.01), 
                 control=nls.lm.control(maxiter=500), trace=T,
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

# ggplot(model.dat, aes(x=treatment, y=aiA)) + geom_bar(stat="identity") + facet_wrap(~species)
# ggplot(model.dat, aes(x=treatment, y=aiE)) + geom_bar(stat="identity") + facet_wrap(~species)








### RUN THE MODELS FOR EACH SPECIES AT EQUILIBRIUM OR INVADING FOR A CONSISTENT CLIMATE CONDITION ###

### CREATE A FUNCTION THAT RUNS THE MODEL
growth = function(N, par.dat, t.num){
  for (i in 1:(t.num-1)){
    N$Na[i+1] = as*(1-ag)*N$Na[i] + exp(log(N$Na[i]*ag+1)*par.dat$Alambda[i]/(1 + par.dat$AaiE[i]*log((N$Ne[i] + 1)*eg) +  par.dat$AaiA[i]*log(N$Na[i]*ag + 1))) -1
    #  N$Na[i+1] =ifelse(N$Na[i +1] == 1, 0, N$Na[i + 1])
    N$Ne[i+1] = es*(1-eg)*N$Ne[i] + exp(log((N$Ne[i]*eg + 1))*par.dat$Elambda[i]/(1 + par.dat$EaiE[i]*log((N$Ne[i]*eg + 1)) + par.dat$EaiA[i]*log(N$Na[i]*ag + 1))) -1
    # N$Ne[i+1] =ifelse(N$Ne[i +1] == 1, 0, N$Ne[i + 1])
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
  Einvade$invader <- "Erodium"
  Einvade$grwr <- Egrwr
  Einvade$time <- as.numeric(row.names(Einvade))
  
  N4 = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
  colnames(N4) = c("Na", "Ne")
  N4[1,] = c(1, Eequil[t.num, 2])
  Ainvade <- growth(N4, par.dat, t.num)
  
  Agrwr <- Ainvade[2,1]/Ainvade[1,1]
  
  Ainvade$invader <- "Avena"
  Ainvade$grwr <- Agrwr
  Ainvade$time <- as.numeric(row.names(Ainvade))
  
  
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
  mutate(treatment = revalue(treatment, c( consistentDry = "Consistent dry", fallDry = "Fall dry",springDry = "Spring dry", controlRain = "Consistent wet")))

consistent.grwr.out <- consistent.out %>%
  dplyr::select(invader, grwr, treatment) %>%
  unique() 


consistent.out2 <- consistent.out %>%
  gather(species, count, Na:Ne) %>%
  tbl_df() %>%
  mutate(count = ifelse(count <= 0.1, 0.1, count))

## Simulation graph ##
pdf("Consistentsimulation.pdf", width = 12, height = 6)
ggplot(consistent.out2, aes(x=time, y=(count), color = species)) + geom_line(size = 4) + facet_grid(invader ~ treatment) + 
  theme_bw() +  theme(text = element_text(size = 24), legend.position = "none", strip.background = element_blank(),
                      #strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16),
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_log10(limits=c(.1, 1200), breaks = c( 1, 10, 100, 1000)) +
  labs(y=expression(paste("Count (individuals/m"^"2",")")), x = "Time step") + scale_color_manual(values = c("tan3", "darkgreen"))
dev.off()


pdf("ddGRWRfromModel.pdf", width = 8, height = 6)
ggplot(consistent.grwr.out, aes(x=treatment, y=log(grwr), fill=invader)) + geom_bar(stat = "identity", position = "dodge") + theme_classic() + 
  labs(x = "Rainfall treatment", y = "Growth rate when rare (logged)", fill = "Species") + theme(text = element_text(size = 24)) +
  geom_hline(yintercept  =0) + scale_fill_manual(values = c("tan3", "darkgreen")) + theme(legend.position = "none")
dev.off()


# 
# 
# ## Check model fit
# model.resid <- nlsResiduals(m1out)
# plot(model.resid)
# plot(model.resid, which=5)
# 
# predicted <- m1out$m$predict()
# actual <- log(subset(dat, !is.na(ERseedout) & treatment == treatments[i])$ERseedout)
# 
# 
# pa <- cbind(predicted, actual) %>%
#   data.frame() %>%
#   tbl_df() 
# 
# plot(log(pa$actual), (pa$predicted))
# 
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
# 
# 
# 
# 

## 

### CREATE A FUNCTION THAT CALCULATES POPULATION CHANGE OVER TIME WHEN ALONE! ####

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

### Run the model for the consistent environment ###
t = 100
 
conDryPar <- consistent.par(model.dat, 1, t)
fallDryPar <- consistent.par(model.dat, 2, t)
springDryPar <- consistent.par(model.dat, 4, t)
conRainPar <- consistent.par(model.dat, 3, t)

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

ggplot(GRWR_nocomp, aes(x=treatment, y=grwalone, color = invader)) + geom_point()

GRWR_withcomp <- rbind(conDry, fallDry, springDry, conRain) %>%
  dplyr::select(invader, grwr, treatment) %>%
  unique() 

GRWRtog <- left_join(GRWR_nocomp, GRWR_withcomp) %>%
  mutate(species = invader) %>%
  mutate(compeffect = log(grwalone/grwr))

ggplot(GRWRtog, aes(y=log(grwalone), x=compeffect, color = invader)) + geom_point(size = 4) + geom_smooth(method = "lm", se =F) + 
  labs(y= "Environment", x = "Competition", color = "Species") +
  theme_classic() + theme(text = element_text(size = 20))


pdf("EnvironmentCompetitionCovariance.pdf", width = 10, height = 8)
ggplot(GRWRtog, aes(x=grwalone, y=grwr, color = invader)) + geom_point(size = 4) + geom_smooth(method = "lm", se =F) + 
  labs(x= "Growth rate when rare and alone", y = "Growth rate when rare in competition", color = "Species") +
  theme_classic() + theme(text = element_text(size = 20)) + geom_abline(slope = 1)
dev.off()

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

outdat <- left_join(GRWRtog, model.dat) %>%
  group_by(species) %>%
  mutate(Ravg = mean(grwalone),
         Ejx = grwalone/Ravg,
         Cjx = compeffect/Ravg,
         Ejx2 = log(grwalone))

resinvader 

ggplot(outdat, aes(y=Ejx, x=compeffect, color = species)) + geom_point(pch =0, size = 4 ) + geom_smooth(method = "lm", se = F) + theme_classic()

ggplot(outdat, aes(y=Ejx, x=compeffect, color = species)) + geom_point() 


ggplot(outdat, aes(x=lambda, y=grwalone)) + geom_point() + xlim(0,10)
