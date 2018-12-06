library(tidyverse)
library(here)
source("modelprojections.R")
source("rainfallsimulations.R")

### CREATE A FUNCTION THAT SETS VARIABLE PARAMETERS ACROSS TIMESTEPS ####
variable.par <- function(model.dat, trtselect, t.num){
  par = as.data.frame(matrix(NA, nrow=t.num, ncol=7))
  colnames(par) = c("Alambda", "AaiE", "AaiA", "Elambda", "EaiE", "EaiA", "treatment")
  for (i in 1:(t.num)){
    myparams <- subset(model.dat, treatment == trtselect[i])
    par$Alambda[i] <- as.numeric(subset(myparams, species == "Avena")[5])
    par$AaiE[i] <- as.numeric(subset(myparams, species == "Avena")[4])
    par$AaiA[i] <- as.numeric(subset(myparams, species == "Avena")[3])
    par$Elambda[i] <- as.numeric(subset(myparams, species == "Erodium")[5])
    par$EaiE[i] <- as.numeric(subset(myparams, species == "Erodium")[4])
    par$EaiA[i] <- as.numeric(subset(myparams, species == "Erodium")[3])
    par$treatment[i] <- as.character(trtselect[i])
    
  }
  return(par)
}

## set terms for the projection
varPar = variable.par(model.dat, rainsummary$raintype, t.num = nrow(rainsummary))
t.num = nrow(rainsummary)

N = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
colnames(N) = c("Na", "Ne")
# calculate average equilibrium to figure out where to set this at!
N[1,] =c(400, 100)

vartreatment <- varPar$treatment
vartreatmentprev <- c(NA, vartreatment[1:t.num-1])

## run the projection
variable.out <- cbind(growth(N, varPar, t.num), vartreatment, vartreatmentprev) %>%
  tbl_df() 

variable.out <- variable.out %>%
  mutate(time = as.numeric(row.names(variable.out)))

## Clean up the labels
variable.out2 <- variable.out %>%
  gather(species, count, Na:Ne) %>%
  tbl_df() %>%
  mutate( species2 = ifelse(species == "Ne", "Erodium", "Avena")) %>%
  mutate(vartreatment=ordered(vartreatment, levels = c(fallDry="fallDry", consistentDry="consistentDry", springDry="springDry", controlRain="controlRain"))) %>%
  mutate(vartreatmentprev=ordered(vartreatmentprev, levels = c(fallDry="fallDry", consistentDry="consistentDry", springDry="springDry", controlRain="controlRain")))

# Plot it
ggplot(variable.out2, aes(x=time, y=(count), color = species)) + geom_line(size = 2) + 
  theme_bw() +  theme(text = element_text(size = 24), legend.position = "none",
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  scale_y_log10() + #scale_y_log10(limits=c(.1, 1200), breaks = c(1, 10, 100, 1000)) +
  labs(y=expression(paste("Count (individuals/m"^"2",")")), x = "Time step")  + scale_color_manual(values = c("grey80", "grey30")) # +
ggsave(here("Figs", "fig3.pdf"), width = 8, height = 6)

