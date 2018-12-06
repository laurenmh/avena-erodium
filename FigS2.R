library(tidyverse)
source("modelsimulations.R")

consistent.out2 <- consistent.out %>%
  gather(species, count, Na:Ne) %>%
  tbl_df() %>%
  mutate(count = ifelse(count <= 0.1, 0.1, count))

## Simulation graph ##
#pdf("Consistentsimulation.pdf", width = 12, height = 6)
p <- ggplot(consistent.out2, aes(x=time, y=(count), color = species)) + geom_line(size = 2) + 
  facet_grid(invader ~ treatment) + 
  theme_bw() +  
  theme(strip.background = element_blank(), legend.position = "none",text = element_text(size = 16), 
        strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_y_log10(limits=c(.1, 3000), breaks = c( 1, 10, 100, 1000)) +
  labs(y=expression(paste("Count (individuals/m"^"2",")")), x = "Time step") + scale_color_manual(values = c("grey80", "grey30"))
#ggsave(here("Figs", "fig3_consistent-simulation.pdf"), width = 12, height = 8)
#ggsave(here("Figs", "fig3_consistent-simulation.jpg"), width = 12, height = 8)



#Use grid.text

pdf(here("Figs", "fig2s.pdf"), width = 8, height = 5)
p
grid.text(c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"),x = c(0.29,0.5,.715,.925 ),
          y = c(rep(.89, 4), rep(.48, 4)),
          gp=gpar(fontsize=16))
dev.off()
