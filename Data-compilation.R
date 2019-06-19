## all data compiled for repository

source("models.R")
source("recruitment_datasummary.R")


dat2 <- dat %>%
  dplyr::select(-density, -rm)

recruittog2 <- recruittog %>%
  filter(species == "Erodium") %>%
  mutate(ERprop = prop, AVprop = 10-prop,
         ERstems = ERO2, AVstems = AV2) %>%
  dplyr::select(plot, subplot, density,ERprop, AVprop, ERstems, AVstems)

tog <- left_join(dat2, recruittog2) %>%
  mutate(ERprop = ERprop/10, AVprop = AVprop/10) %>%
  mutate(ERstem = round(ERstems, 0),
         AVstem = round(AVstems, 0),
         AVseedout = round(AVseedout, 0),
         ERseedout = round(ERseedout, 0)) %>%
  dplyr::select(shelterBlock, plot, subplot, treatment, density, ERprop, AVprop, 
                ERseedin, ERstem, ERseedout,
          AVseedin, AVstem, AVseedout) 

write_csv(tog, here("Data", "Hallett_avena-erodium.csv"))