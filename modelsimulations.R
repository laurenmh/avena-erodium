source("models.R")

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
