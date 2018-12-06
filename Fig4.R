library(tidyverse)

source("coexistence_partitioning_weighted.R")

# double check it all works
avena_LDGR <- mean(avena_invader)
erodium_LDGR <- mean(erodium_invader)

avena_results <- c(avena_LDGR, avena_epilson_0, avena_epsilon_alpha, 
                   avena_epsilon_lambda, avena_epsilon_interaction)
erodium_results <- c(erodium_LDGR, erodium_epilson_0, erodium_epsilon_alpha, 
                     erodium_epsilon_lambda, erodium_epsilon_interaction)


params <- c("a", "b", "c", "d", "e")

dat <- as.data.frame(t(rbind(params, avena_results, erodium_results)))
names(dat)[2:3] = c("Avena", "Erodium")

dat2 <- dat %>%
  gather(species, value, Avena:Erodium) %>%
  mutate(value = as.numeric(value)) %>%
  mutate(mycol = ifelse(params == "a", "a", "b"))


ggplot(dat2, aes(x=params, y=value, fill = mycol)) + geom_bar(stat = "identity") + 
  facet_wrap(~species) + 
  scale_x_discrete("Parameters", labels = c("a" = expression(bar("r")[i]) ,
                                            "b" = expression(bar(epsilon)[i]^0),
                                            "c" = expression(bar(epsilon)[i]^alpha),
                                            "d" = expression(bar(epsilon)[i]^lambda),
                                            "e" = expression(bar(epsilon)[i]^{alpha*lambda}))) +
  theme_bw() + 
  theme(legend.position = "none", strip.background = element_blank(), 
        text = element_text(size = 16), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 20)) +
  scale_fill_manual(values = c( "grey40", "grey70")) + 
  ylab("Low density growth rate partitioning") 
ggsave(here("Figs", "Fig4.pdf"), width = 8, height = 5)
