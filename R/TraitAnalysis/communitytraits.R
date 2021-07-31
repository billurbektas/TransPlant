#### COMMUNITY TRAITS ANALYSIS (CWM) ####
#Analysing CWM trends in functional traits

# Source traits drakeplan from runsource_traitplan.R first before running this script. Have removed originSiteID from trait drake plant 'merge_traits', should I add that back in?

# Source ggplot theme
source('./R/theme_ggplot.R')

traitdat <- alltraits %>% select(Region, destSiteID, SpeciesName, Plant_Veg_Height_cm, Leaf_Area_cm2, SLA_cm2_g, C_percent, N_percent, P_percent, Seed_Mass, Stem_density) %>%
  rename(plant_height = Plant_Veg_Height_cm, leaf_area = Leaf_Area_cm2, SLA = SLA_cm2_g , leaf_C = C_percent , leaf_N = N_percent, leaf_P = P_percent, seed_mass = Seed_Mass, SD = Stem_density)

# Create ODT treatment and nest data
dd <- dat %>% select(Region, originSiteID, destSiteID, Treatment) %>% 
  distinct() %>% 
  filter(Treatment == "Warm") %>% 
  select(-Treatment) %>% 
  mutate(comm = pmap(.l = list(R = Region, O = originSiteID, D = destSiteID), .f = function(R, O, D){
    bind_rows(
      originControls = dat %>% filter(Region == R, destSiteID == O, Treatment == "LocalControl"),
      destControls = dat %>% filter(Region == R, destSiteID == D, Treatment == "LocalControl"),
      warmed =  dat %>% filter(Region == R, destSiteID == D, Treatment == "Warm"),
      .id = "ODT") 
  })) 

# Unnest species and lengthen dataframe
traitdat <- dd %>% select(comm) %>% unnest(comm) 

# Calculate coverage of traits per plot
coverage <- traitdat %>% group_by(Region, Year, destSiteID, destPlotID) %>% summarise_at(.vars=c('leaf_C', 'leaf_P', 'leaf_N', 'plant_height', 'SLA', 'leaf_area', 'seed_mass', 'SD'), ~mean(is.na(.))) #using all data
coverage %>% group_by(Region) %>% summarise_at(.vars=c('leaf_C', 'leaf_P', 'leaf_N', 'plant_height', 'SLA', 'leaf_area', 'seed_mass', 'SD'), ~sum(.==1))
coverage %>% group_by(Region) %>% summarise_at(.vars=c('leaf_C', 'leaf_P', 'leaf_N', 'plant_height', 'SLA', 'leaf_area', 'seed_mass', 'SD'), ~sum(.>0.5))
coverage_f <- traitdat %>% filter(Rel_Cover > 0.05) %>% group_by(Region, Year, destSiteID, destPlotID) %>% summarise_at(.vars=c('leaf_C', 'leaf_P', 'leaf_N', 'plant_height', 'SLA','leaf_area', 'seed_mass', 'SD'), ~mean(is.na(.))) #first subsetting plants <10% cover
coverage %>% group_by(Region) %>% summarise_at(.vars=c('leaf_C', 'leaf_P', 'leaf_N', 'plant_height', 'SLA', 'leaf_area', 'seed_mass', 'SD'), ~sum(.==1))
# this removed ~50 plots (500 if >0.1) 

# Scaling traits and removing species with zero trait values
traitdat <- traitdat %>% select(Region, Year, originSiteID, destSiteID, ODT, destPlotID, species, Rel_Cover, plant_height:SD) %>% tibble() %>%
                 filter(Rel_Cover > 0.1) 
traitdat_f <- traitdat %>% filter(!rowSums(is.na(.[,9:16]))==8)

# Calculate CWMs
CWM <- traitdat_f %>% group_by(Region, Year, destSiteID, ODT, destPlotID) %>%
  summarise(across(.cols= plant_height:SD, .fns = list(CWM = ~weighted.mean({.}, Rel_Cover, na.rm=T)))) %>%
  pivot_longer(plant_height_CWM:SD_CWM, names_to="CWM") %>%
  filter(!is.na(value), value>0, value != 'Inf')

# Plot change in CWMs over time
colour_odt <- c("#A92420", "#016367", "#FBC00E")

CWM %>%
  group_by(Region) %>% mutate(Year_diff=max(Year - min(Year))) %>%
  ungroup() %>%
  filter(Year_diff > 0, !is.na(value)) %>%
  group_by(Region, destSiteID, ODT, destPlotID, CWM) %>%
  nest() %>% 
  mutate(npoints = map_int(data, nrow)) %>%
  filter(npoints > 2) %>%
  mutate(model = map(data, ~lm(log(value) ~ Year, data = .x))) %>% 
  mutate(tidied = map(model, broom::tidy)) %>% 
  unnest(tidied) %>%
  filter(term == 'Year') %>% #, CWM != "leaf_C_CWM") %>% #ODT == "warmed"
  ggplot(aes(x=estimate, y=ODT, col = ODT)) + 
  geom_pointrange(aes(xmin=min(estimate), xmax=max(estimate)), adjust = 3, alpha=0.2) + #used to be geom_density
  geom_errorbar(aes(xmin = min(estimate), xmax = max(estimate)), width = 0.2) +
  scale_color_manual(values = colour_odt) + 
  xlim(-0.5, 0.5) +
  TP_theme() + 
  labs(x=expression(Delta*Slope)) +
  labs(color="Treatment", y="Number of communities") +
  facet_wrap(~CWM, scales = 'free') +
  geom_vline(xintercept=0, linetype="dashed", 
             color = "grey", size=1) 

# Plot CWM of last year per treatment
colour_odt <- c("#016367", "#FBC00E", "#A92420")

CWM %>%
  group_by(Region) %>% mutate(Year=max(Year)) %>%
  ungroup() %>%
  mutate(ODT = factor(ODT, levels=c("originControls","warmed", "destControls"))) %>%
  filter(CWM=="SLA_CWM", !is.na(value)) %>%
  ggplot(aes(y=log(value), x=ODT, fill=ODT)) + 
  geom_boxplot(alpha=0.6) + 
  scale_fill_manual(values = colour_odt) + 
  TP_theme() + 
  labs(x="Treatment") +
  labs(color="Treatment", y="log(SLA)") +
  facet_wrap(~Region, scales='free') 

CWM %>%
  group_by(Region) %>% mutate(Year=max(Year)) %>%
  ungroup() %>%
  mutate(ODT = factor(ODT, levels=c("originControls","warmed", "destControls"))) %>%
  filter(CWM=="plant_height_CWM", !is.na(value)) %>%
  ggplot(aes(y=log(value), x=ODT, fill=ODT)) + 
  geom_boxplot(alpha=0.6) + 
  scale_fill_manual(values = colour_odt) + 
  TP_theme() + 
  labs(x="Treatment") +
  labs(color="Treatment", y="log(Plant Height)") +
  facet_wrap(~Region, scales='free') 


