library(tidyverse)
library(here)

## DECEMBER DATA
mydat <- read_csv("Data/dd_stem_count_20141216.csv")


## CREATE A KEY TO ADD TREATMENT INFO TO SUBSEQUENT DATA
trtconvert <- mydat %>%
  dplyr::select(plot, subplot, number, AVtrt, EROtrt, density)

# ## FEBRUARY DATA
# mydat_feb <- read.csv("dd_stem_count_20150221.csv") %>%
#   tbl_df() %>%
#   mutate(plot = paste("Plot", Plot, sep="_"),
#          subplot = Treatment,
#          number = Subplot) %>%
#   select(-Plot, -Subplot, -Treatment, -AVquadrat, -EROquadrat, -Notes)%>%
#   select(date, plot, subplot, number, EROcount, AVcount)
# 
# mydat_feb2 <- tbl_df(merge(mydat_feb, trtconvert))
# names(mydat_feb2)[5:6] = c("ERO", "AV")
# 
# 
# ## MARCH DATA
# mydat_mar <- read.csv("dd_stem_count_20150320.csv") %>%
#   tbl_df() %>%
#   mutate(plot=paste("Plot", plot, sep="_")) %>%
#   select(date, plot, subplot, number, EROcount, AVcount)
# 
# mydat_mar2 <- tbl_df(merge(mydat_mar, trtconvert))
# names(mydat_mar2)[5:6] = c("ERO", "AV")


## APRIL DATA
mydat_apr <- read_csv("Data/dd_stem_count_20150423.csv") %>%
mutate(AVcor=1, AVcor=ifelse(AVquadrat2 == "10x25cm", 2.5, AVcor),
       AVcor=ifelse(AVquadrat2 == "5x5cm", 25, AVcor)) %>%
mutate(EROcor=1, EROcor=ifelse(EROquadrat2 == "10x25cm", 2.5, EROcor),
       EROcor=ifelse(EROquadrat2 == "5x5cm", 25, EROcor)) %>%
  mutate(AV_seeds=(AVgoodseed + AVmoldseed)*2*AVcor,
         AV_seedsgood=AVgoodseed*2*AVcor) %>%
    mutate(AV_propgood = AV_seedsgood/AV_seeds) %>%
  mutate(ERO_seeds=(EROstork + EROseed)*5*EROcor) %>%
  mutate(AV_stemtog=NA, AV_stemtog = ifelse(is.na(AVstemsTotal) == F, AVstemsTotal, AV_stemtog),
                AV_stemtog = ifelse(is.na(AVstemswseed) == F, AVstemswseed, AV_stemtog), 
    AV_stems = (AV_stemtog)*AVcor,
         ERO_stems = EROcount*EROcor) %>%
  mutate(AV_r = AV_seeds/AV_stems,
         ERO_r = ERO_seeds/ERO_stems) %>%
  dplyr::select(plot, subplot, AV_seeds, ERO_seeds, AV_stems, ERO_stems, AV_r, ERO_r, AV_seedsgood, AV_propgood)

trtconvert2 <- trtconvert %>%
  mutate(plot = extract_numeric(plot))

mydat_apr2 <- tbl_df(merge(mydat_apr, trtconvert2))



# ADD SHELTER CONVERSION KEY
shelterkey <- read_csv("Data/Shelter_key.csv")
mydat2 <- merge(mydat_apr2, shelterkey) %>%
  tbl_df() %>%
gather(spptrt, prop, AVtrt:EROtrt)%>%
  mutate(prop=extract_numeric(prop))%>%
  mutate(treatment=ordered(treatment, levels = c(consistentDry="consistentDry", fallDry="fallDry", springDry="springDry", controlRain="controlRain")))

# #just look at the monocultures
# mono <- mydat2 %>%
#   filter(prop==10) %>%
#   select(plot, subplot, density, ERO_seeds, AV_seedsgood, treatment, shelterBlock, spptrt, prop)%>%
#   gather(spp, count, ERO_seeds:AV_seedsgood) %>%
#   filter((spptrt=="AVtrt"&spp=="AV_seedsgood")|(spptrt=="EROtrt"&spp=="ERO_seeds")) %>%
#   mutate(seed=20, seed=ifelse(density=="D2", 200, seed), seed=ifelse(density=="D3", 2000, seed))%>%
#   mutate(propgerm=count/seed)%>%
#   ggplot(aes(x=density, y=propgerm, group=1))+geom_point() + facet_grid(spp~treatment) +geom_smooth(stat="smooth")
# 
# mono+theme_bw()


## PRELIMINARY VISUALIZATIONS

# Erodium
ero <- mydat2 %>%
  filter(spptrt=="EROtrt") %>%
  mutate(seed=20*(prop/10), seed=ifelse(density=="D2", 200*(prop/10), seed), seed=ifelse(density=="D3", 2000*(prop/10), seed)) %>%
  mutate(ERO_stemcorrect=1, ERO_stemcorrect=ifelse(ERO_stems > seed, seed/ERO_stems, ERO_stemcorrect)) %>%
  mutate(ERO_seeds2=ERO_seeds*ERO_stemcorrect) %>%
  mutate(R=ERO_seeds/seed) %>%
  ggplot(aes(x=(prop/10), y=R))+ geom_point()+ facet_grid(density~treatment, scale="free") 

ero + geom_smooth(method="lm") +theme_bw() +ylab("Per capita growth rate") + xlab("Seeding ratio") + ggtitle("Erodium") #+ ylim(0,.5)

# Avena
av <- mydat2 %>%
  filter(spptrt=="AVtrt") %>%
 # filter(density=="D1") %>%
  mutate(seed=20*(prop/10), seed=ifelse(density=="D2", 200*(prop/10), seed), seed=ifelse(density=="D3", 2000*(prop/10), seed)) %>%
  mutate(AV_stemcorrect=1, AV_stemcorrect=ifelse(AV_stems > seed, seed/AV_stems, AV_stemcorrect))%>%
  mutate(AV_seeds2=AV_seeds*AV_stemcorrect) %>%
  mutate(AV_seedsgood2=AV_seedsgood*AV_stemcorrect) %>%
  mutate(R=AV_seedsgood2/seed)%>%
  ggplot(aes(x=(prop/10), y=R))+ geom_point()+ facet_grid(density~treatment,  scale="free")
av + geom_smooth(method="lm") + theme_bw()+ylab("Per capita growth rate") + xlab("Seeding ratio") + ggtitle("Avena")
  

## Prop viable
av_propgood <- mydat2 %>%
  filter(spptrt=="AVtrt")%>%
  # filter(density=="D1") %>%
  mutate(seed=20*(prop/10), seed=ifelse(density=="D2", 200*(prop/10), seed), seed=ifelse(density=="D3", 2000*(prop/10), seed))%>%
  mutate(AV_stemcorrect=1, AV_stemcorrect=ifelse(AV_stems > seed, seed/AV_stems, AV_stemcorrect))%>%
  mutate(AV_seeds2=AV_seeds*AV_stemcorrect) %>%
  mutate(AV_seedsgood2=AV_seedsgood*AV_stemcorrect) %>%
  mutate(R=AV_seedsgood2/seed)%>%
  ggplot(aes(x=(prop/10), y=AV_propgood))+ geom_point()+ facet_grid(density~treatment,  scale="free")
av_propgood + geom_smooth(method="lm") + theme_bw()+ylab("Proportion of viable Avena seeds") + xlab("Seeding ratio")

# 
# 
# pdf("ddExperiment_growthrate_20150423.pdf")
#  ero+geom_smooth(method = "lm", formula = y ~ poly(x,2)) +theme_bw() + xlab("Seeding ratio") +ylab("Per capita growth rate") + 
#   ggtitle("Erodium") + geom_hline(y=1)
# av+geom_smooth(method = "lm") + theme_bw()+ xlab("Seeding ratio")+ylab("Per capita growth rate")  + ggtitle("Avena") + geom_hline(y=1)
# dev.off()

## APPEND DATA
# Avena
avdat <- mydat2 %>%
  mutate(density = as.character(density)) %>%
  filter(spptrt=="AVtrt")%>%
  mutate(seed=20*(prop/10), 
         seed=ifelse(density=="D2", 200*(prop/10), seed), 
         seed=ifelse(density=="D3", 2000*(prop/10), seed)) %>%
  mutate(AV_stemcorrect=1, AV_stemcorrect=ifelse(AV_stems > seed, seed/AV_stems, AV_stemcorrect))%>%
  mutate(AV_seeds2=AV_seeds*AV_stemcorrect) %>%
  mutate(seedsout=AV_seedsgood*AV_stemcorrect) %>%
  mutate(R=seedsout/seed) %>%
  mutate(species="Avena") %>%
  dplyr::select(plot, R, density, treatment, prop, shelterBlock, seed, species, seedsout, subplot)

# Erodium
erodat <- mydat2 %>%
  filter(spptrt=="EROtrt")%>%
  mutate(seed=20*(prop/10), seed=ifelse(density=="D2", 200*(prop/10), seed), seed=ifelse(density=="D3", 2000*(prop/10), seed))%>%
  mutate(ERO_stemcorrect=1, ERO_stemcorrect=ifelse(ERO_stems > seed, seed/ERO_stems, ERO_stemcorrect))%>%
  mutate(seedsout=ERO_seeds*ERO_stemcorrect) %>%
  mutate(R=seedsout/seed) %>%
  mutate(species="Erodium") %>%
  dplyr::select(plot, R, density, treatment, prop, shelterBlock, seed, species, seedsout, subplot)

# Together
togdat <- rbind(avdat, erodat) %>%
  # remove D3 for being unrealistically high density
  filter(density != "D3")


## RAW VISUAL
togdat2 <- togdat %>%
  mutate(treatment = as.character(treatment)) %>%
  #   mutate(treatment = revalue(treatment, c(fallDry = "Fall dry", consistentDry = "Consistent dry", springDry = "Spring dry", controlRain = "Consistent wet"))) %>%
  mutate(treatment=ordered(treatment, levels = c( consistentDry="consistentDry", fallDry="fallDry",springDry="springDry", controlRain="controlRain"))) %>%
  mutate(treatment = recode(treatment, consistentDry = "Consistent dry", fallDry = "Fall dry",  springDry = "Spring dry", controlRain = "Consistent wet")) %>%
  mutate(density = ordered(density, levels = c(D1 = "D1", D2 = "D2"))) %>%
  mutate(density = recode(density, D1 = "Low density", D2 = "High density"))

ggplot(subset(togdat2, species == "Avena" & R != 66), aes(x=(prop/10), y=(R)))+ geom_point(size = 4, color = "tan3")+ facet_grid(density~treatment,  scale="free") +
  geom_smooth(method="lm", color ="tan3", lwd = 2, se = F) + theme_bw() + ylab("Per capita population growth rate") + 
  geom_point(dat = subset(togdat2, species == "Erodium"), size = 4, color = "darkgreen") +
  geom_smooth(dat = subset(togdat2, species == "Erodium"), method="lm", color = "darkgreen", lwd = 2, se = F) +
  xlab("Seeding ratio")  + geom_hline(yintercept=1) + scale_x_continuous(limits = c(0, 1), breaks = c(.1, .5, .9, 1)) + 
  theme(strip.background = element_blank(), text = element_text(size = 16), 
        strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## BW version
ggplot(subset(togdat2, species == "Avena" & R != 66), aes(x=(prop/10), y=(R)))+ geom_point(size = 3, color = "grey80")+ facet_grid(density~treatment,  scale="free") +
  geom_smooth(method="lm", color ="grey80", lwd = 1.5, se = F) + theme_bw() + ylab("Per capita population growth rate") + 
  geom_point(dat = subset(togdat2, species == "Erodium"), size = 3, color = "grey30") +
  geom_smooth(dat = subset(togdat2, species == "Erodium"), method="lm", color = "grey30", lwd = 1.5, se = F) +
  xlab("Seeding ratio")  + geom_hline(yintercept=1) + scale_x_continuous(limits = c(0, 1), breaks = c(.1, .5, .9, 1), labels = c(".1", ".5", ".9", "1")) + 
  theme(strip.background = element_blank(), text = element_text(size = 16), 
        strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(here("Figs", "fig1_avena-erodium-percapita.pdf"), width = 8, height = 8)
ggsave(here("Figs", "fig1_avena-erodium-percapita.jpg"), width = 8, height = 8)

