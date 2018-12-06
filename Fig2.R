source("dd_modelprojections.R")
source("coexistence_partitioning.R")

library(tidyverse)
library(gridExtra)

calcSE<-function(x){
  x <- x[is.na(x)==F]
  sd(x)/sqrt(length(x))
}

## data for consistent from modelprojections
dat2a <- consistent.grwr.out %>%
  select(invader, treatment, grwrChesson) %>%
  mutate(sdgrwr = NA, segrwr = NA)
  
## data for historic from coexistence partitioning
dat2b <- as.data.frame(cbind(Avena = avena_invader, Erodium = erodium_invader)) %>%
  gather(invader, grwr, Avena:Erodium) %>%
  mutate(treatment = "Historic") %>%
  group_by(treatment, invader) %>%
  summarize(grwrChesson = mean(grwr), sdgrwr = sd(grwr), segrwr = calcSE(grwr)) %>%
  as.data.frame()

dat2 <- as.data.frame(rbind(dat2a, dat2b))

# Figure 2

  ggplot(dat2, aes(x=treatment, y=grwrChesson, fill=invader)) +
  geom_bar(stat = "identity", position = "dodge") + theme_classic() + 
    geom_errorbar(aes(ymin = grwrChesson - segrwr, ymax = grwrChesson + segrwr),position = position_dodge(width=0.9), width=0.2) + 
  labs(x = "Rainfall condition", y = expression("Growth rate when rare"), fill = "Species") + 
  theme(text = element_text(size = 24)) +
  geom_hline(yintercept  =0) + scale_fill_manual(values = c("grey80", "grey30")) + 
  theme(legend.position = "none") + 
  geom_vline(xintercept = 4.5, linetype = "dashed", lwd = 1.2) 


  ### old #

# Figure 2a
fig2a <- ggplot(consistent.grwr.out, aes(x=treatment, y=grwrChesson, fill=invader)) + geom_bar(stat = "identity", position = "dodge") + theme_classic() + 
  labs(x = "Rainfall treatment", y = expression("Growth rate when rare (r"[i]* "-r"[r]*")"), fill = "Species") + theme(text = element_text(size = 24)) +
  geom_hline(yintercept  =0) + scale_fill_manual(values = c("grey80", "grey30")) + theme(legend.position = "none")


# Figure 2b
# These are the values for Lauren Hallett to use for figure 2b
dat2b <- as.data.frame(cbind(Avena = avena_invader, Erodium = erodium_invader)) %>%
  gather(species, grwr, Avena:Erodium) %>%
  group_by(species) %>%
  summarize(grwrChesson = mean(grwr), sdgrwr = sd(grwr)) 

fig2b <- ggplot(dat2b, aes(x=species, y=meangrwr, fill=species)) + geom_bar(stat = "identity", position = "dodge") + 
  theme_classic() + 
  labs(x = "Species", y = "Average growth rate when rare") + 
  theme(text = element_text(size = 24)) +
  geom_hline(yintercept  =0) + scale_fill_manual(values = c("grey80", "grey30")) + theme(legend.position = "none")  +
  #facet_wrap(~species ,scales = "free") +
  theme(strip.background = element_blank(), text = element_text(size = 16), 
        strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

grid.arrange(fig2a, fig2b, nrow = 1)

# Break it up by panels to avoid this issue
ggplot(time.avg, aes(x=treatments, y=overall.grwr, fill=species)) + geom_bar(stat = "identity", position = "dodge") + theme_classic() + 
  labs(x = "Long-term conditions", y = "Average growth rate when rare") + theme(text = element_text(size = 24)) +
  geom_hline(yintercept  =0) + scale_fill_manual(values = c("grey80", "grey30")) + theme(legend.position = "none")  + facet_wrap(~species ,scales = "free") +
  theme(strip.background = element_blank(), text = element_text(size = 16), 
        strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
