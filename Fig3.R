library(tidyverse)
library(here)
library(cowplot)
source("modelsimulations.R")
source("rainfallsimulation.R")

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

rainsummary$time <- seq(1:121)
variable.out3 <- left_join(variable.out2, rainsummary)

# # Plot it
# ggplot(variable.out3, aes(x=time, y=(count), color = species)) + geom_line(size = 2) + 
#   theme_bw() +  theme(text = element_text(size = 24), legend.position = "none",
#                       panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  scale_y_log10() + #scale_y_log10(limits=c(.1, 1200), breaks = c(1, 10, 100, 1000)) +
#   labs(y=expression(paste("Count (individuals/m"^"2",")")), x = "Time step")  + scale_color_manual(values = c("grey80", "grey30")) # +
# #ggsave(here("Figs", "fig3.pdf"), width = 8, height = 6)


ggplot(subset(variable.out3), aes(x=year, y=(count), color = species)) + geom_line(size = 1.2) +
  theme_bw() +  theme(text = element_text(size = 20), legend.position = "none",
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  scale_y_log10() + #scale_y_log10(limits=c(.1, 1200), breaks = c(1, 10, 100, 1000)) +
  labs(y=expression(paste("Count (individuals/m"^"2",")")), x = "Year")  +
  scale_color_manual(values = c("tan3", "darkgreen"))
ggsave("fig3-talks.pdf", width = 8, height = 6)



# Plot it
a <-  ggplot(subset(variable.out3), aes(x=year, y=(count), color = species)) + geom_line(size = 1.2) +
  theme_bw() +  theme(text = element_text(size = 20), legend.position = "none",
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
   scale_y_log10() + #scale_y_log10(limits=c(.1, 1200), breaks = c(1, 10, 100, 1000)) +
  labs(y=expression(paste("Count (individuals/m"^"2",")")), x = "Year")  +
  scale_color_manual(values = c("grey80", "grey30")) + 
  annotate("text", x= 1895, y = 1100, label= "(a)", size = 7) # +

b <- ggplot(subset(variable.out3), aes(x=year, y=Total)) +
  geom_line(size = 1.2) +
  theme_bw() +  theme(text = element_text(size = 20), 
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(y="Rainfall (mm)", x = "Year")   + 
  annotate("text", x= 1895, y = 1200, label= "(b)", size = 7)

pdf(here::here("Figs", "fig3_v2.pdf"), width = 8, height = 8)
plot_grid(a + theme(axis.text.x = element_blank(), 
                    axis.title.x = element_blank(), axis.ticks.x = element_blank()),
          b, rel_heights = c(3/5, 2/5), nrow = 2, align = "v")
dev.off()


#ggsave(here("Figs", "fig3.pdf"), width = 8, height = 6)
