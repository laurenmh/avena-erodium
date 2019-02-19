library(tidyverse)

source("invader_resident_comparison.R")

params <- c("a", "b", "c", "d", "e")

dat <- as.data.frame(t(rbind(params, ir_avena_results_weighted, ir_erodium_results_weighted)))
names(dat)[2:3] = c("Avena", "Erodium")

dat2 <- dat %>%
  gather(species, value, Avena:Erodium) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(mycol = ifelse(params == "a", "a", "b"))


quartz(width=5, height=5)
p <- ggplot(dat2, aes(x=params, y=value, fill = mycol)) + geom_bar(stat = "identity") + 
  facet_wrap(~species) + 
  scale_x_discrete("Parameters", labels = c("a" = expression(bar("r")[i]) ,
                                            "b" = expression(bar(Delta)[i]^0),
                                            "c" = expression(bar(Delta)[i]^alpha),
                                            "d" = expression(bar(Delta)[i]^lambda),
                                            "e" = expression(bar(Delta)[i]^{alpha*lambda}))) +
  theme_bw() + 
  theme(legend.position = "none", strip.background = element_blank(), 
        text = element_text(size = 16), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 20)) +
  scale_fill_manual(values = c( "grey40", "grey70")) + 
  ylab("Partitioning of growth rate when rare") + geom_hline(yintercept = 0)


# add panel labels
g <- ggplotGrob(p)
#Use grid.text

pdf(here("Figs", "fig4_ir.pdf"), width = 8, height = 5)
p
grid.text(c("(a)", "(b)"), x = c(0.495,.95),
          y = c(.88,.88),
          gp=gpar(fontsize=16))
dev.off()



# adjust equation to be deltas for the invader resident 
