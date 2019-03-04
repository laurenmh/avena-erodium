source("modelsimulations.R")
source("invader_resident_comparison.R")

library(tidyverse)
library(gridExtra)

calcSE<-function(x){
  x <- x[is.na(x)==F]
  sd(x)/sqrt(length(x))
}

### CONSISTENT 

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
  ggsave(here("Figs", "fig2.pdf"), width = 10, height = 6)
  ggsave(here("Figs", "fig2.jpg"), width = 10, height = 6)
  
