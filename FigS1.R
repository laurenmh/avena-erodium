source("recruitment_datasummary.R")

## Graph of raw recruitment values 

## BW version
r <- ggplot(subset(recruittog, species == "Avena"), aes(x=(prop/10), y=(propgerm))) + 
  geom_point(size = 3, color = "grey80") + facet_grid(density~treatment,  scale="free") +
  geom_smooth(method="lm", color ="grey80", lwd = 1.5, se = F) + 
  theme_bw() + ylab("Percent recruitment") + 
  geom_point(dat = subset(recruittog, species == "Erodium"), size = 3, color = "grey30") +
  geom_smooth(dat = subset(recruittog, species == "Erodium"), 
              method="lm", color = "grey30", lwd = 1.5, se = F) +
  xlab("Seeding ratio")  +
  scale_x_continuous(limits = c(0, 1), breaks = c(.1, .5, .9, 1), 
                     labels = c(".1", ".5", ".9", "1")) + 
  scale_y_continuous(limits = c(0,1.15), breaks = c(0, .25, .5, .75, 1)) + 
  theme(strip.background = element_blank(), text = element_text(size = 16), 
        strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ggsave(here("Figs", "figS1_avena-erodium-germination.pdf"), width = 8, height = 8)
# ggsave(here("Figs", "figS1_avena-erodium-germination.jpg"), width = 8, height = 8)

pdf(here("Figs", "figS1.pdf"), width = 8, height =5)
r
grid.text(c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"),x = c(0.27,0.49,.71,.925 ),
          y = c(rep(.895, 4), rep(.48, 4)),
          gp=gpar(fontsize=16))
dev.off()
