alltog <- merge(recruittog, togdat2) %>%
  tbl_df() %>%
  mutate(relgerm = AV2/(ERO2+AV2),
         relgerm = ifelse(species == "Erodium", ERO2/(AV2+ERO2), relgerm),
         relgerm = ifelse(!is.finite(relgerm), 1, relgerm)) %>%
  mutate(compgerm = AV2,
         compgerm = ifelse(species == "Avena", ERO2, compgerm)) %>%
  tbl_df()

ggplot(subset(alltog, species == "Avena" & R != 66), aes(x=relgerm, y=(R))) + 
  geom_point(size = 2, color = "grey80")+ facet_grid(density~treatment,  scale="free_y") +
  geom_smooth(method="lm", color ="grey80", lwd = 1, se = F) + theme_bw() +
  ylab("Per capita population growth rate") + 
  geom_point(dat = subset(alltog, species == "Erodium"), size = 2, color = "grey30") +
  geom_smooth(dat = subset(alltog, species == "Erodium"), method="lm", color = "grey30", lwd = 1, se = F) +
  xlab("Stem frequency")  + geom_hline(yintercept=1) +
  scale_x_continuous(limits = c(0, 1), breaks = c(.1, .5, .9, 1), labels = c(".1", ".5", ".9", "1")) + 
  theme(strip.background = element_blank(), text = element_text(size = 16), 
        strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(subset(alltog, species == "Avena" & R != 66), aes(x=compgerm, y=(R))) + 
  geom_point(size = 2, color = "grey80")+ facet_grid(density~treatment,  scale="free_y") +
  geom_smooth(method="lm", color ="grey80", lwd = 1, se = F) + theme_bw() +
  ylab("Per capita population growth rate") + 
  geom_point(dat = subset(alltog, species == "Erodium"), size = 2, color = "grey30") +
  geom_smooth(dat = subset(alltog, species == "Erodium"), method="lm", color = "grey30", lwd = 1, se = F) +
  xlab("Seeding ratio")  + geom_hline(yintercept=1) +
 # scale_x_continuous(limits = c(0, 1), breaks = c(.1, .5, .9, 1), labels = c(".1", ".5", ".9", "1")) + 
  theme(strip.background = element_blank(), text = element_text(size = 16), 
        strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

