library(nlme)
library(multcomp)
library(xtable)
source("recruitment_datasummary.R")
source("seeds_datasummary.R")

### CODE FOR MIXED-EFFECT MODELS OF RECRUITMENT AND FECUNDITY ###

# Put recruitment and fecundity data together
# Need to have run dd_recruitment_datasummary and dd_seeds_datasummary previously

alltog <- left_join(recruittog[c("plot", "subplot", "propgerm", "species")], togdat) %>%
  mutate(growth = seedsout/seed) %>%
  mutate(treatment1 = as.character(treatment))  %>%
  mutate(fallrain = ifelse(treatment1 == "controlRain" | treatment1 == "springDry", "yes", "no")) %>%
  mutate(denstrt = as.factor(paste(density, treatment1, sep="_")),
         densprop = as.factor(paste(density, prop, sep = "_")),
         denstrtprop = as.factor(paste(density, prop, treatment1, sep ="_")))

### ERODIUM propgerm as a function of density, prop and treatment ###
ggplot(subset(alltog, species == "Erodium"), aes(x=prop, y=propgerm)) + geom_point() + 
  facet_grid(density ~ treatment, scales = "free") + geom_smooth(method = "lm", se = F)


ggplot(subset(alltog, species == "Erodium"), aes(x=prop, y=propgerm)) + geom_point() + 
  facet_grid(density ~ fallrain, scales = "free") + geom_smooth(method = "lm", se = F)

## Effect of fall rain on germination rates 
l <- lm(propgerm~prop*density + prop*treatment1, data = subset(alltog, species == "Erodium", !is.na(propgerm)))
summary(l)    

l <- lme(propgerm ~ prop*treatment1 + prop*density, random = ~1|shelterBlock, data = subset(alltog, species == "Erodium" & !is.na(propgerm)))
l1 <- summary(l)[5]

l2 <- tidy(l, effects = "fixed")
xtable(l2)

## Now just do yes versus no for fall rain
l <- lm(propgerm~prop*density + prop*fallrain, data = subset(alltog, species == "Erodium", !is.na(propgerm)))
summary(l)    

l <- lme(propgerm ~ prop*fallrain + prop*density, random = ~1|shelterBlock, data = subset(alltog, species == "Erodium" & !is.na(propgerm)))
summary(l)


### AVENA propgerm as a function of density, prop and treatment ###
ggplot(subset(alltog, species == "Avena"), aes(x=prop, y=propgerm)) + geom_point() + 
  facet_grid(density ~ treatment, scales = "free") + geom_smooth(method = "lm", se = F)


ggplot(subset(alltog, species == "Avena"), aes(x=prop, y=propgerm)) + geom_point() + 
  facet_wrap(~density, scales = "free") + geom_smooth(method = "lm", se = F)

## Effect of fall rain on germination rates 
l <- lm(propgerm~prop*density + prop*treatment1, data = subset(alltog, species == "Avena", !is.na(propgerm)))
summary(l)    

l <- lme(propgerm ~ prop*treatment1 + prop*density, random = ~1|shelterBlock, data = subset(alltog, species == "Avena" & !is.na(propgerm)))
summary(l)

l2 <- tidy(l, effects = "fixed")
xtable(l2)
## Now just do yes versus no for fall rain
l <- lm(propgerm~prop*density + prop*fallrain, data = subset(alltog, species == "Avena", !is.na(propgerm)))
summary(l)    

l <- lme(propgerm ~ prop*fallrain + prop*density, random = ~1|shelterBlock, data = subset(alltog, species == "Avena" & !is.na(propgerm)))
summary(l)

############################
#### NOW DO SEEDS! #########
############################

### ERODIUM R as a function of density, prop and treatment ###
ggplot(subset(alltog, species == "Erodium"), aes(x=prop, y=R)) + geom_point() + 
  facet_grid(density ~ treatment, scales = "free") + geom_smooth(method = "lm", se = F)


ggplot(subset(alltog, species == "Erodium"), aes(x=prop, y=R)) + geom_point() + 
  facet_grid(density ~ fallrain, scales = "free") + geom_smooth(method = "lm", se = F)


## Effect of fall rain on germination rates 
l <- lm(R~prop*density + prop*treatment1, data = subset(alltog, species == "Erodium", !is.na(propgerm)))
summary(l)    


l <- lme(R ~ denstrtprop, random = ~1|shelterBlock, data = subset(alltog, species == "Erodium" & !is.na(propgerm)))
summary(l)   


summary(glht(l, linfct=mcp(denstrtprop="Tukey")))



l <- lme(R ~ prop*density*treatment1, random = ~1|shelterBlock, data = subset(alltog, species == "Erodium" & !is.na(propgerm)))
summary(l)
l2 <- tidy(l, effects = "fixed")
xtable(l2)

## Now just do yes versus no for fall rain
l <- lm(R~prop*density + prop*fallrain, data = subset(alltog, species == "Erodium", !is.na(propgerm)))
summary(l)    

l <- lme(R ~ prop*fallrain + prop*density, random = ~1|shelterBlock, data = subset(alltog, species == "Erodium" & !is.na(propgerm)))
summary(l)



### AVENA propgerm as a function of density, prop and treatment ###
ggplot(subset(alltog, species == "Avena"), aes(x=prop, y=R)) + geom_point() + 
  facet_grid(density ~ treatment, scales = "free") + geom_smooth(method = "lm", se = F)


## Effect of fall rain on germination rates 
l <- lm(R~prop*density + prop*treatment1, data = subset(alltog, species == "Avena", !is.na(propgerm)))
summary(l)    

l <- lme(R ~ prop*treatment1 + prop*density, random = ~1|shelterBlock, data = subset(alltog, species == "Avena" & !is.na(R)))
summary(l)

## Now just do yes versus no for fall rain
l <- lm(R~prop*density + prop*fallrain, data = subset(alltog, species == "Avena", !is.na(R)))
summary(l)    

l <- lme(R ~ prop*fallrain + prop*density, random = ~1|shelterBlock, data = subset(alltog, species == "Avena" & !is.na(R)))
summary(l)



############################################
#### NOW DO GERM EFFECTS ON SEEDS! #########
##############################################

### ERODIUM R as a function of density, prop and treatment ###
ggplot(subset(alltog, species == "Erodium"), aes(x=prop, y=R)) + geom_point() + 
  facet_grid(density ~ treatment, scales = "free") + geom_smooth(method = "lm", se = F)


ggplot(subset(alltog, species == "Erodium"), aes(x=prop, y=R)) + geom_point() + 
  facet_grid(density ~ fallrain, scales = "free") + geom_smooth(method = "lm", se = F)

## Effect of fall rain on germination rates 
l <- lm(R~propgerm*density + propgerm*treatment1, data = subset(alltog, species == "Erodium", !is.na(propgerm)))
summary(l)    

l <- lme(R ~ prop*treatment1 +prop*density, random = ~1|shelterBlock, data = subset(alltog, species == "Erodium" & !is.na(propgerm)))
summary(l)
l2 <- tidy(l, effects = "fixed")
xtable(l2)


## Now just do yes versus no for fall rain
l <- lm(R~prop*density + prop*fallrain, data = subset(alltog, species == "Erodium", !is.na(propgerm)))
summary(l)    

l <- lme(R ~ prop*fallrain + prop*density, random = ~1|shelterBlock, data = subset(alltog, species == "Erodium" & !is.na(propgerm)))
summary(l)



### AVENA propgerm as a function of density, prop and treatment ###
ggplot(subset(alltog, species == "Avena"), aes(x=prop, y=R)) + geom_point() + 
  facet_grid(density ~ treatment, scales = "free") + geom_smooth(method = "lm", se = F)


## Effect of fall rain on germination rates 
l <- lm(R~prop*density + prop*treatment1, data = subset(alltog, species == "Avena", !is.na(propgerm)))
summary(l)    

l <- lme(R ~ prop*treatment1 + prop*density, random = ~1|shelterBlock, data = subset(alltog, species == "Avena" & !is.na(R)))
summary(l)
l2 <- tidy(l, effects = "fixed")
xtable(l2)

## Now just do yes versus no for fall rain
l <- lm(R~prop*density + prop*fallrain, data = subset(alltog, species == "Avena", !is.na(R)))
summary(l)    

l <- lme(R ~ prop*fallrain + prop*density, random = ~1|shelterBlock, data = subset(alltog, species == "Avena" & !is.na(R)))
summary(l)




############ NOW GERMINATION EFFECTS ON RECRUITMENT #############
ggplot(subset(alltog, species == "Erodium"), aes(x=propgerm, y=R, color = treatment)) + geom_point() + 
  geom_smooth(method = "lm", se = F)

ggplot(subset(alltog, species == "Avena"), aes(x=propgerm, y=R, color = treatment)) + geom_point() + 
   geom_smooth(method = "lm", se = F)





ggplot(subset(alltog, species == "Erodium"), aes(x=propgerm, y=R)) + geom_point() + 
  facet_grid(density ~ treatment, scales = "free") + geom_smooth(method = "lm", se = F)


ggplot(subset(alltog, species == "Erodium"), aes(x=propgerm, y=R,  color = treatment)) + geom_point() + 
  geom_smooth(method = "lm", se = F)

ggplot(subset(alltog, species == "Erodium"), aes(x=prop, y=propgerm)) + geom_point() + geom_smooth(method = "lm", se = F) +
  facet_grid(density~treatment)




ggplot(subset(alltog, species == "Erodium"), aes(x=propgerm, y=R,  color = treatment)) + geom_point() + 
  geom_smooth(method = "lm", se = F) + facet_wrap(~treatment)

# prop germination led to the most seeds by a lot, then a partially sig effect of density and a density/trt interaction for one
l <- lm(R~propgerm*density*treatment, data = subset(alltog, species == "Erodium"))
summary(l)        


ggplot(subset(alltog, species == "Avena"), aes(x=propgerm, y=R)) + geom_point() + 
  facet_grid(density ~ treatment, scales = "free") + geom_smooth(method = "lm", se = F)



ggplot(subset(alltog, species == "Avena"), aes(x=propgerm, y=R, color =treatment)) + geom_point() + 
 geom_smooth(method = "lm", se = F)


ggplot(subset(alltog, species == "Avena"), aes(x=propgerm, y=R)) + geom_point() + 
  facet_wrap(~treatment , scales = "free") + geom_smooth(method = "lm", se = F)


ggplot(subset(alltog, species == "Avena"), aes(x=propgerm, y=R)) + geom_point() + 
  facet_grid(density ~ treatment, scales = "free") + geom_smooth(method = "lm", se = F)


## a higher density led to more seeds for avena, but otherwise nothing!
l <- lm(R~propgerm*density*treatment, data = subset(alltog, species == "Avena"))
summary(l)  



ggplot(subset(alltog, species == "Avena"), aes(x=prop, y=R)) + geom_point() + 
  facet_grid(density ~ treatment, scales = "free") + geom_smooth(method = "lm", se = F)

l <- lm(R~prop*density*treatment, data = subset(alltog, species == "Avena"))
summary(l)    
