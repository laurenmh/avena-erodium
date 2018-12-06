library(tidyverse)
library(nlme)
library(multcomp)

### FORMATTING THE RECRUITMENT DATA ###

# read data in
mydat<-read_csv("Data/dd_stem_count_20141216.csv")

trtconvert <- mydat %>%
  dplyr::select(plot, subplot, number, AVtrt, EROtrt, density)

#set up conversion key for different density measurements
subplots <- as.data.frame(sort(unique(mydat$subplot)))
names(subplots)[1]="subplot"
AVcor <- as.data.frame(c(1,1,1,1,0,1,2.5,2.5,2.5,0,1,25,25,25,0))
names(AVcor)[1]="AVcor"
EROcor <- as.data.frame(c(1,0,1,1,1,2.5,0,2.5,1,2.5, 25,0,25,1,25))
names(EROcor)[1]="EROcor"
subkey<-cbind(subplots, AVcor, EROcor)

mydat <- tbl_df(merge(mydat, subkey)) %>%
  mutate(ERO2=ERO*EROcor, AV2=AV*AVcor) %>%
  mutate(plot=extract_numeric(plot))

#add the shelter conversion key
shelterkey<-read.csv("Data/Shelter_key.csv")
mydat2<-merge(mydat, shelterkey)%>%
  tbl_df() %>%
gather(spptrt, prop, AVtrt:EROtrt)%>%
  mutate(prop=extract_numeric(prop))%>%
  mutate(treatment=ordered(treatment, levels = c(fallDry="fallDry", consistentDry="consistentDry", springDry="springDry", controlRain="controlRain")))
  
##look at each species by density, proportion and treatment
erodat <-mydat2 %>%
  filter(spptrt=="EROtrt")%>%
  mutate(seed=20*(prop/10), seed=ifelse(density=="D2", 200*(prop/10), seed), seed=ifelse(density=="D3", 2000*(prop/10), seed))%>%
  mutate(ERO2=ifelse(ERO2>seed, seed, ERO2))%>%
  mutate(propgerm=ERO2/seed) %>%
  mutate(ratio = prop/10) %>%
  mutate(species = "Erodium") %>%
  mutate(trts = as.factor(treatment)) %>%
  mutate(bigfact = as.factor(paste(ratio, treatment, density, sep="_"))) %>%
  mutate(littlefact = as.factor(paste(ratio, treatment, sep = "_")))
# 
# ero <- erodat %>%
#   ggplot(aes(x=(prop/10), y=propgerm))+ geom_point()+ facet_grid(density~treatment, scale="free")
# ero+geom_smooth(method="lm") +theme_bw() +ylab("Proportion recruited") + xlab("Seeding ratio") + ggtitle("Erodium") #+ ylim(0,.5)
# 
# l <- lme(propgerm ~ ratio*treatment + density, random = ~1|shelterBlock, data = subset(erodat, !is.na(propgerm)))
# summary(l)
# 
# l2 <- lme(propgerm ~ littlefact + density, random = ~1|shelterBlock, data = subset(erodat, !is.na(propgerm)))
# summary(l2)
# glht(l2, linfct = mcp(density = "Tukey")) 

avdat <-mydat2 %>%
  filter(spptrt=="AVtrt")%>%
  mutate(seed=20*(prop/10), seed=ifelse(density=="D2", 200*(prop/10), seed), seed=ifelse(density=="D3", 2000*(prop/10), seed))%>%
  mutate(AV2=ifelse(AV2>seed, seed, AV2))%>%
  mutate(propgerm=AV2/seed) %>%
  mutate(species = "Avena")

# 
# av <- avdat %>%
#   ggplot(aes(x=(prop/10), y=propgerm))+ geom_point()+ facet_grid(density~treatment, scale="free")
# av+geom_smooth(method="lm") + theme_bw()+ylab("Proportion recruited") + xlab("Seeding ratio") + ggtitle("Avena")
#   
# l <- lme(propgerm ~ ratio*treatment*density, random = ~1|shelterBlock, data = subset(erodat, !is.na(propgerm)))
# summary(l)

# Put it all together
recruittog <- rbind(erodat[c("plot", "subplot", "prop", "propgerm", "density", "treatment", "species")], 
                    avdat[c("plot", "subplot", "prop", "propgerm", "density", "treatment", "species")])  %>%
  filter(density != "D3") %>%
  mutate(treatment=ordered(treatment, levels = c( consistentDry="consistentDry", fallDry="fallDry",springDry="springDry", controlRain="controlRain"))) %>%
  mutate(treatment = recode(treatment, consistentDry = "Consistent dry", fallDry = "Fall dry",  springDry = "Spring dry", controlRain = "Consistent wet")) %>%
  mutate(density = ordered(density, levels = c(D1 = "D1", D2 = "D2"))) %>%
  mutate(density = recode(density, D1 = "Low density", D2 = "High density"))


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
