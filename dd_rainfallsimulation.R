### Historical projects for Browns Valley ###
library(tidyverse)
library(lubridate)

## Pull in the prism data and clean
rain <- read_csv("PRISM_brownsvalley_long.csv", skip = 10) %>%
  mutate(ppt = `ppt..inches.`*2.54*10) %>%
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
  mutate(Total = Early + Late) %>%
  mutate(raintype = "controlRain",
         raintype = ifelse(Early < quantile(rainsummary$Early, .5), "fallDry", raintype),
         raintype = ifelse(Late < quantile(rainsummary$Late, .5), "springDry", raintype),
         raintype = ifelse(Total < quantile(rainsummary$Total, .5), "consistentDry", raintype)) 


## Visualize the rainfall scenarios
#pdf("Rainfal history.pdf", width = 10, height = 8)
ggplot(rainsummary, aes(x=year, y=Total))  + geom_line()+
  geom_point(aes(color = raintype), size = 4) + theme_bw() + labs(x="Year", y="Total rainfall (mm)")
#dev.off()

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
N[1,] =c(1002.06748, 10)

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

pdf("Variable_Projection.pdf", width = 8, height = 6)
ggplot(variable.out2, aes(x=time, y=(count), color = species)) + geom_line(size = 4) + 
  theme_bw() +  theme(text = element_text(size = 24), legend.position = "none",
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  scale_y_log10() + #scale_y_log10(limits=c(.1, 1200), breaks = c(1, 10, 100, 1000)) +
  labs(y=expression(paste("Count (individuals/m"^"2",")")), x = "Time step")  + scale_color_manual(values = c("tan3", "darkgreen")) # +
dev.off()  
#geom_hline(yintercept = 1) + geom_hline(yintercept = 10)#+ 
  #ggtitle("Variable rainfall")


simco <- variable.out2 %>%
  dplyr::filter(time > 20) %>%
  mutate(species = species2) %>%
  dplyr::select(-species2)





###  ERODIUM MONOCULTURE
## set terms for the projection
varPar = variable.par(model.dat, rainsummary$raintype, t.num = nrow(rainsummary))
t.num = nrow(rainsummary)

N = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
colnames(N) = c("Na", "Ne")
N[1,] =c(0, 30)

### CREATE A FUNCTION THAT RUNS THE MODEL for just ERO
growthEROmono = function(N, par.dat, t.num){
  for (i in 1:(t.num-1)){
    N$Na[i+1] = 0
    #  N$Na[i+1] =ifelse(N$Na[i +1] == 1, 0, N$Na[i + 1])
    N$Ne[i+1] = es*(1-eg)*N$Ne[i] + exp(log((N$Ne[i]*eg + 1))*par.dat$Elambda[i]/(1 + par.dat$EaiE[i]*log((N$Ne[i]*eg + 1)) + par.dat$EaiA[i]*log(N$Na[i]*ag + 1))) -1
    # N$Ne[i+1] =ifelse(N$Ne[i +1] == 1, 0, N$Ne[i + 1])
  }
  return(N)
}

eroMONO <- growthEROmono(N, varPar, t.num) %>%
  tbl_df()
eroMONO$timestep <- as.numeric(row.names(eroMONO))
  
ggplot(eroMONO, aes(x=timestep, y=Ne)) + geom_line()  
ggplot(subset(variable.out2, species == "Ne"), aes(x=time, y=count)) + geom_line()  



###  AVENA MONOCULTURE
## set terms for the projection
varPar = variable.par(model.dat, rainsummary$raintype, t.num = nrow(rainsummary))
t.num = nrow(rainsummary)

N = as.data.frame(matrix(NA, nrow=t.num, ncol=2))
colnames(N) = c("Na", "Ne")
N[1,] =c(1000, 0)


### CREATE A FUNCTION THAT RUNS THE MODEL for just ERO
### CREATE A FUNCTION THAT RUNS THE MODEL for just ERO

### CREATE A FUNCTION THAT RUNS THE MODEL
growthAVmono = function(N, par.dat, t.num){
  for (i in 1:(t.num-1)){
    N$Na[i+1] = as*(1-ag)*N$Na[i] + exp(log(N$Na[i]*ag+1)*par.dat$Alambda[i]/(1 + par.dat$AaiE[i]*log((N$Ne[i] + 1)*eg) +  par.dat$AaiA[i]*log(N$Na[i]*ag + 1))) -1
    #  N$Na[i+1] =ifelse(N$Na[i +1] == 1, 0, N$Na[i + 1])
    N$Ne[i+1] = 0
    # N$Ne[i+1] =ifelse(N$Ne[i +1] == 1, 0, N$Ne[i + 1])
  }
  return(N)
}


avMONO <- growthAVmono(N, varPar, t.num) %>%
  tbl_df()
avMONO$timestep <- as.numeric(row.names(avMONO))


ggplot(avMONO, aes(x=timestep, y=Na)) + geom_line()  
ggplot(subset(variable.out2, species == "Na"), aes(x=time, y=count)) + geom_line()  

