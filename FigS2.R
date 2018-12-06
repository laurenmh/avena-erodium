library(tidyverse)
source("modelsimulations.R")

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
