library(codyn)
library(broom)
library(brms)
library(modelr)
library(tidybayes)
source('R/theme_ggplot.R')

#### Make plant community cover relative abundance (at the moment includes categories like other, rock, moss, etc.)
dat2 <- dat %>%  
  select(Region, Year, originSiteID, destSiteID, destBlockID, Elevation, Treatment, destPlotID, SpeciesName, Rel_Cover) %>%
  filter(!is.na(Rel_Cover)) %>% #not creating other cover as biomass, doesn't exist
  mutate(Cover = Rel_Cover) %>%
  select(-Rel_Cover) %>%
  group_by(Region, Year, originSiteID, destSiteID, destBlockID, Elevation, Treatment, destPlotID) %>%
  mutate(Total_Cover = sum(Cover, na.rm=T), Rel_Cover = Cover / Total_Cover) %>%
  ungroup()
dd %>% filter(Region == "IT_MatschMazia2") %>% .$comm %>% .[[1]] %>% group_by(originSiteID, destSiteID, Year, Treatment) %>% tally

# Create ODT and nest communities
dd <- dat2 %>% select(Region, originSiteID, destSiteID, Treatment) %>% 
  distinct() %>% 
  filter(Treatment == "Warm") %>% 
  select(-Treatment) %>% 
  mutate(comm = pmap(.l = list(R = Region, O = originSiteID, D = destSiteID), .f = function(R, O, D){
    bind_rows(
      originControls = dat2 %>% filter(Region == R, destSiteID == O, Treatment == "LocalControl"),
      destControls = dat2 %>% filter(Region == R, destSiteID == D, Treatment == "LocalControl"),
      warmed =  dat2 %>% filter(Region == R, destSiteID == D, Treatment == "Warm"),
      .id = "ODT") 
  })) %>%
  mutate(specrich = map(comm, ~ {.} %>% group_by(Year, ODT, destPlotID) %>% summarize(SR=n_distinct(SpeciesName))), 
         colonisation = map(comm, ~turnover(.x, time.var= "Year", species.var= "SpeciesName", abundance.var= "Rel_Cover", replicate.var="destPlotID", metric = "appearance")),
         extinction = map(comm, ~turnover(.x, time.var= "Year", species.var= "SpeciesName", abundance.var= "Rel_Cover", replicate.var="destPlotID", metric = "disappearance")),
         turnover = map(comm, ~turnover(.x, time.var= "Year", species.var= "SpeciesName", abundance.var= "Rel_Cover", replicate.var="destPlotID", metric = "total"))) 

# calculate coldyn metrics (can add this later, puts it all in one for easier manipulation)
# test <- dd %>% 
#   filter(!Region %in% c("FR_Lautaret", "IN_Kashmir", "US_Colorado")) %>%
#   mutate(comm_sim = map(comm, ~.x %>% select(ODT, destPlotID) %>% distinct())) %>%
#   mutate(dat = pmap(.l = list(C=colonisation, E=extinction, To=turnover), .f = function(C, E, To){
#     C <- C %>% rename(value = appearance)
#     E <- E %>% rename(value = disappearance)
#     To <- To %>% rename(value = total)
#     bind_rows('C' = C,'E'= E, 'T' = To, .id='comdyn')})) %>% 
#   mutate(dat = map2(dat, comm_sim, ~left_join(.x, .y, by = "destPlotID"))) %>%
#   mutate(dat = map(dat, ~.x %>% mutate(Year_0 = Year-min(Year)))) %>%
#   unnest(dat) %>% 
#  filter(ODT=='warmed') %>%
#  filter(comdyn%in% c('C', 'E', 'T')) 

#### Plot SR patterns ####
#dest control, origin control, warmed
colour_odt <- c("#A92420", "#016367", "#FBC00E")

SR <- dd %>% select(Region, originSiteID, destSiteID, specrich) %>%
  mutate(dat = map(specrich, ~.x %>% mutate(Year_0 = Year-min(Year)))) %>%
  unnest(specrich) 

SR %>%
  ggplot(aes(x = Year, y = SR)) + 
  scale_colour_manual(values = colour_odt) + 
  geom_line(aes(group=destPlotID, color=ODT), alpha=0.2) +
  geom_smooth(aes(group = ODT, color=ODT), method = "lm", se = F) +
  facet_wrap(~Region, nrow=2) +
  TP_theme() + 
  scale_x_continuous(breaks=c(2010,2013,2016)) +
  labs(color = "Treatment", y = 'Species Richness') +
  labs(title = 'Species Richness over time', color = "Treatment") 

SR_m <- SR %>% 
  mutate(Gradient = paste(originSiteID, destSiteID, sep="_")) %>%
  mutate(Years = scale(Year))

mS_1 = brm(
  SR ~ 0 + ODT + ODT:Years + (Years|Region/Gradient), 
  data = SR_m,
  family=poisson(link = "log") # data is counts
)
summary(mS_1)
posterior_samples(mS_1)
save(mS_1, file="./R/ModelOutputs/fit_mS_spec.Rdata")


#### Plot colonisation patterns ####
#Need to remove IN_Kashmir, FR_Lautaret and US_Colorado as only 2 years of data
CO <- dd %>% filter(!Region %in% c("FR_Lautaret", "IN_Kashmir", "US_Colorado")) %>%
  mutate(comm_sim = map(comm, ~.x %>% select(ODT, destPlotID) %>% distinct())) %>%
  mutate(dat = map2(colonisation, comm_sim, ~left_join(.x, .y, by = "destPlotID"))) %>%
  mutate(dat2 = map(dat, ~.x %>% mutate(Year_0 = Year-min(Year)))) %>%
  unnest(dat2) 

CO %>%
  ggplot(aes(x = Year, y = appearance)) + 
  scale_colour_manual(values = colour_odt) + 
  geom_point(aes(group=destPlotID, color=ODT), alpha=0.1) +
  geom_smooth(aes(group = ODT, color=ODT), method = "lm", se = F) +
  facet_wrap(~Region, nrow=2) +
  TP_theme() + 
  scale_x_continuous(breaks=c(2010,2013,2016)) + 
  labs(title = 'Colonisation over time', color = "Treatment") 

CO %>%
  ggplot(aes(x = Year, y = appearance)) + 
  scale_colour_manual(values = colour_odt) + 
  #geom_point(aes(group=destPlotID, color=ODT), alpha=0.1) +
  geom_smooth(aes(group = ODT, color=ODT), method = "lm", se = F) +
  facet_wrap(~Region, nrow=2) +
  TP_theme() + 
  scale_x_continuous(breaks=c(2010,2013,2016)) + 
  labs(title = 'Colonisation over time', color = "Treatment") 

CO_m <- CO %>% 
  mutate(Gradient = paste(originSiteID, destSiteID, sep="_")) %>%
  mutate(col = 0.00001+(1-2*0.00001)*(appearance-min(appearance))/(max(appearance)-min(appearance)))  %>%
  mutate(Years = scale(Year))

mC_1 = brm(
  col ~0 + ODT:Years:Region + (1|Gradient), 
  data = CO_m,
  family=Beta(link = "logit", link_phi = "log")  # data is counts
)

mC_2 = brm(
  col ~ 0 + Years:ODT:Region + (1|Gradient), 
  data = CO_m,
  family=Beta(link = "logit", link_phi = "log")  # data is counts
)

save(mC_1, file="./R/ModelOutputs/fit_mC_col_all.Rdata")
save(mC_2, file="./R/ModelOutputs/fit_mC_col_site.Rdata")


#### Plot extinction patterns ####
EX <- dd %>% filter(!Region %in% c("FR_Lautaret", "IN_Kashmir", "US_Colorado")) %>%
  mutate(comm_sim = map(comm, ~.x %>% select(ODT, destPlotID) %>% distinct())) %>%
  mutate(dat = map2(extinction, comm_sim, ~left_join(.x, .y, by = "destPlotID"))) %>%
  mutate(dat2 = map(dat, ~.x %>% mutate(Year_0 = Year-min(Year)))) %>%
  unnest(dat2) 

EX %>%
  ggplot(aes(x = Year, y = disappearance)) + 
  scale_colour_manual(values = colour_odt) + 
  geom_point(aes(group=destPlotID, color=ODT), alpha=0.1) +
  geom_smooth(aes(group = ODT, color=ODT), method = "lm", se = F) +
  facet_wrap(~Region, nrow=2) +
  TP_theme() + 
  scale_x_continuous(breaks=c(2010,2013,2016)) + 
  labs(title = 'Extinction over time', color = "Treatment") 

EX_m <- EX %>% 
  mutate(Gradient = paste(originSiteID, destSiteID, sep="_")) %>%
  mutate(ext = 0.00001+(1-2*0.00001)*(disappearance-min(disappearance))/(max(disappearance)-min(disappearance)))  %>%
  mutate(Years = scale(Year))

mE_1 = brm(
  ext ~ 0 + Years*ODT + (1|Region/Gradient), 
  data = EX_m,
  family=Beta(link = "logit", link_phi = "log")  # data is counts
)

mE_2 = brm(
  ext ~ 0 + Years:ODT:Region + (1|Gradient), 
  data = EX_m,
  family=Beta(link = "logit", link_phi = "log")  # data is counts
)

save(mE_1, file="./R/ModelOutputs/fit_mE_ext_all.Rdata")
save(mE_2, file="./R/ModelOutputs/fit_mE_ext_site.Rdata")

#### Plot turnover patterns ####
TU <- dd %>% filter(!Region %in% c("FR_Lautaret", "IN_Kashmir", "US_Colorado")) %>%
  mutate(comm_sim = map(comm, ~.x %>% select(ODT, destPlotID) %>% distinct())) %>%
  mutate(dat = map2(turnover, comm_sim, ~left_join(.x, .y, by = "destPlotID"))) %>%
  mutate(dat2 = map(dat, ~.x %>% mutate(Year_0 = Year-min(Year)))) %>%
  unnest(dat2) 

TU %>%
  ggplot(aes(x = Year, y = total)) + 
  scale_colour_manual(values = colour_odt) + 
  geom_point(aes(group=destPlotID, color=ODT), alpha=0.1) +
  geom_smooth(aes(group = ODT, color=ODT), method = "lm", se = F) +
  facet_wrap(~Region, nrow=2) +
  TP_theme() + 
  scale_x_continuous(breaks=c(2010,2013,2016)) + 
  labs(title = 'Turnover over time', color = "Treatment") 

TU %>%
  ggplot(aes(x = Year, y = total)) + 
  scale_colour_manual(values = colour_odt) + 
  geom_line(aes(group=destPlotID, color=ODT), alpha=0.2) +
  geom_smooth(aes(group = ODT, color=ODT), method = "lm", se = F) +
  #geom_smooth(aes(color=darkgrey), method = "lm", se = F) +
  #facet_wrap(~Region, nrow=2) +
  TP_theme() + 
  scale_x_continuous(breaks=c(2010,2013,2016)) + 
  labs(title = 'Turnover over time', color = "Treatment") 

TU_m <- TU %>% 
  mutate(Gradient = paste(originSiteID, destSiteID, sep="_")) %>%
  mutate(tur = 0.00001+(1-2*0.00001)*(total-min(total))/(max(total)-min(total)))  %>%
  mutate(Years = scale(Year))

mT_1 = brm(
  tur ~ 0 + Years:ODT:Region + (1|Gradient), 
  data = TU_m,
  family=Beta(link = "logit", link_phi = "log")  # data is counts
)

mT_2 = brm(
  tur ~ 0 + Years:ODT:Region + (1|Gradient), 
  data = TU_m,
  family=Beta(link = "logit", link_phi = "log")  # data is counts
)

save(mT_1, file="./R/ModelOutputs/fit_mT_tur_all.Rdata")
save(mT_2, file="./R/ModelOutputs/fit_mT_tur_site.Rdata")

#### Plot C, E and T for only transplanted turfs ####
colour_comdyn <- c("#49BEB7", "black", "black")
#"#7496D2", "#CECD7B"
#c("#fe875d", "#49BEB7", "black")

dd %>%
  filter(!Region %in% c("FR_Lautaret", "IN_Kashmir", "US_Colorado")) %>%
  mutate(comm_sim = map(comm, ~.x %>% select(ODT, destPlotID) %>% distinct())) %>%
  mutate(dat = pmap(.l = list(C=colonisation, E=extinction, To=turnover), .f = function(C, E, To){
    C <- C %>% rename(value = appearance)
    E <- E %>% rename(value = disappearance)
    To <- To %>% rename(value = total)
    bind_rows('C' = C,'E'= E, 'T' = To, .id='comdyn')})) %>% 
  mutate(dat = map2(dat, comm_sim, ~left_join(.x, .y, by = "destPlotID"))) %>%
  unnest(dat) %>% 
  filter(ODT == "warmed", comdyn %in% c('C', 'E')) %>% 
  ggplot(aes(x=Year, y=value, color=comdyn, group = interaction(destPlotID, comdyn))) + 
  scale_colour_manual(values = colour_comdyn) + 
  geom_line(alpha=0.2) +
  geom_smooth(aes(group=comdyn, color=comdyn), method='lm', se=F) +
  facet_wrap(~Region, nrow=2) +
  TP_theme() + 
  scale_x_continuous(breaks=c(2011,2013,2015,2017)) + 
  labs(color="Process", y="Proportional change", x='Year') 

#On one graph, with years scaled to 0 +

dd_T <- dd %>%
  filter(!Region %in% c("FR_Lautaret", "IN_Kashmir", "US_Colorado")) %>%
  mutate(comm_sim = map(comm, ~.x %>% select(ODT, destPlotID) %>% distinct())) %>%
  mutate(dat = pmap(.l = list(C=colonisation, E=extinction, To=turnover), .f = function(C, E, To){
    C <- C %>% rename(value = appearance)
    E <- E %>% rename(value = disappearance)
    To <- To %>% rename(value = total)
    bind_rows('C' = C,'E'= E, 'T' = To, .id='comdyn')})) %>% 
  mutate(dat = map2(dat, comm_sim, ~left_join(.x, .y, by = "destPlotID"))) %>%
  mutate(dat = map(dat, ~.x %>% mutate(Year_0 = Year-min(Year)))) %>%
  unnest(dat) %>% 
  filter(ODT == "warmed") %>%
  filter(comdyn%in% c('C', 'E', 'T'))  

dd_T %>% 
  ggplot(aes(x=Year_0, y=value, color=comdyn, group = interaction(Region, destPlotID, comdyn))) + 
  scale_colour_manual(values = colour_comdyn) + 
  geom_line(alpha=0.2) +
  geom_smooth(aes(group=interaction(Region, comdyn), color=comdyn), method='lm', se=F) +
  #geom_smooth(aes(group=comdyn, color=comdyn), method='lm', se=F) +
  TP_theme() + 
  facet_grid(~comdyn) +
  labs(color="Process", y="Proportional change in communities", x='Duration (years)') 

library(nlme)
library(emmeans)
#Colonisation
dd_C <- dd_T %>% filter(comdyn=='C')
m1<-lme(value ~ Year_0*Region, random = ~1|originSiteID, method = "ML", data=dd_C) 
mnull<-lme(value ~ Year_0+ Region, random = ~1|originSiteID, method = "ML", data=dd_C) 
anova(m1, mnull)  #***, 99.582
mnull1<-lme(value ~ Year_0, random = ~1|originSiteID, method = "ML", data=dd_C) 
anova(mnull, mnull1) #***, 42.692
mnull2<-lme(value ~ Region, random = ~1|originSiteID, method = "ML", data=dd_C) 
anova(mnull, mnull2) #**, 8.674

##Extinction
dd_C <- dd_T %>% filter(comdyn=='E')
m1<-lme(value ~ Year_0*Region, random = ~1|originSiteID, method = "ML", data=dd_C) 
mnull<-lme(value ~ Year_0+ Region, random = ~1|originSiteID, method = "ML", data=dd_C) 
anova(m1, mnull)  #***, 100.610
mnull1<-lme(value ~ Year_0, random = ~1|originSiteID, method = "ML", data=dd_C) 
anova(mnull, mnull1) #***, 59.958
mnull2<-lme(value ~ Region, random = ~1|originSiteID, method = "ML", data=dd_C) 
anova(mnull, mnull2) 

##Turnover
dd_C <- dd_T %>% filter(comdyn=='T')
m1<-lme(value ~ Year_0*Region, random = ~1|originSiteID, method = "ML", data=dd_C) 
mnull<-lme(value ~ Year_0+ Region, random = ~1|originSiteID, method = "ML", data=dd_C) 
anova(m1, mnull)  #***, 50.267
mnull1<-lme(value ~ Year_0, random = ~1|originSiteID, method = "ML", data=dd_C) 
anova(mnull, mnull1) #***, 68.218
mnull2<-lme(value ~ Region, random = ~1|originSiteID, method = "ML", data=dd_C) 
anova(mnull, mnull2) 

#With region as fixed effect
#Colonisation
dd_C <- dd_T %>% filter(comdyn=='C')
m1<-lme(value ~ Year_0*Region, random = ~1|originSiteID, method = "ML", data=dd_C) 
summary(m1)
anova(m1)
mnull<-lme(value ~ Region, random = ~1|originSiteID, method = "ML", data=dd_C) 
anova(m1, mnull)

##Extinction
dd_C <- dd_T %>% filter(comdyn=='E')
m1<-lme(value ~ Year_0*Region, random = ~1|originSiteID, method = "ML", data=dd_C) 
summary(m1)
anova(m1)

##Turnover
dd_C <- dd_T %>% filter(comdyn=='T')
m1<-lme(value ~ Year_0*Region, random = ~1|originSiteID, method = "ML", data=dd_C) 
summary(m1) 
anova(m1)


#### Extract elevation and compare ####
library(broom)
colour_comdyn <- c("#fe875d", "#356288", "black")
dd_E <- dd %>%
  mutate(comm_sim = map(comm, ~.x %>% select(Elevation, ODT, destPlotID) %>% 
                          distinct() %>%
                          mutate(Elev_diff=max(Elevation - min(Elevation))))) %>% 
  mutate(dat = pmap(.l = list(C=colonisation, E=extinction, To=turnover), .f = function(C, E, To){
    C <- C %>% rename(value = appearance)
    E <- E %>% rename(value = disappearance)
    To <- To %>% rename(value = total)
    bind_rows('C' = C,'E'= E, 'T' = To, .id='comdyn')})) %>% 
  mutate(dat = map2(dat, comm_sim, ~left_join(.x, .y, by = "destPlotID"))) %>%
  mutate(dat = map(dat, ~.x %>% filter(ODT == "warmed"))) %>%
  unnest(dat) %>% 
  group_by(Region, destSiteID, Elevation, Elev_diff, comdyn) %>%
  nest() %>% 
  mutate(model = map(data, ~lm(value ~ Year, data = .x))) %>% 
  mutate(tidied = map(model, tidy)) %>% 
  unnest(tidied) %>%
  filter(term == 'Year') 

dd_E %>%
  ggplot(aes(x=Elev_diff, y=abs(estimate), color=comdyn)) + 
  geom_point() + 
  scale_colour_manual(values = colour_comdyn) + 
  geom_smooth(aes(group=comdyn, color=comdyn), method='lm', se=F) +
  TP_theme() + 
  labs(color="Process", y="Absolute slope over time", x=expression(Delta*Elevation)) 

#Colonisation
dd_C <- dd_E %>% filter(comdyn=='C')
m1<-lme(abs(estimate) ~ Elev_diff, random = ~1|Region, method = "ML", data=dd_C) 
summary(m1)
anova(m1)

##Extinction
dd_C <- dd_E %>% filter(comdyn=='E')
m1<-lme(abs(estimate) ~ Elev_diff, random = ~1|Region, method = "ML", data=dd_C) 
summary(m1)
anova(m1)

##Turnover
dd_C <- dd_E %>% filter(comdyn=='T')
m1<-lme(abs(estimate) ~ Elev_diff, random = ~1|Region, method = "ML", data=dd_C) 
summary(m1) 
anova(m1)

## Elevation but not slope estimates, just to raw values
dd_E <- dd %>%
  mutate(comm_sim = map(comm, ~.x %>% select(Elevation, ODT, destPlotID) %>% 
                          distinct() %>%
                          mutate(Elev_diff=max(Elevation - min(Elevation))))) %>% 
  mutate(dat = pmap(.l = list(C=colonisation, E=extinction, To=turnover), .f = function(C, E, To){
    C <- C %>% rename(value = appearance)
    E <- E %>% rename(value = disappearance)
    To <- To %>% rename(value = total)
    bind_rows('C' = C,'E'= E, 'T' = To, .id='comdyn')})) %>% 
  mutate(dat = map2(dat, comm_sim, ~left_join(.x, .y, by = "destPlotID"))) %>%
  mutate(dat = map(dat, ~.x %>% filter(ODT == "warmed"))) %>%
  unnest(dat) 

dd_E %>%
  ggplot(aes(x=Elev_diff, y=value, color=comdyn)) + 
  geom_point() + 
  scale_colour_manual(values = colour_comdyn) + 
  geom_smooth(aes(group=comdyn, color=comdyn), method='lm', se=F) +
  TP_theme() + 
  labs(color="Process", y="Proportional change", x=expression(Delta*Elevation)) 


#Colonisation
dd_C <- dd_E %>% filter(comdyn=='C')
m1<-lme(value ~ Elev_diff, random = ~1|Region, method = "ML", data=dd_C) 
summary(m1)
anova(m1)

dd_C <- dd_E %>% filter(comdyn=='C')
m1<-lm(value ~ Elev_diff, data=dd_C) 
summary(m1) 

##Extinction
dd_C <- dd_E %>% filter(comdyn=='E')
m1<-lme(value ~ Elev_diff, random = ~1|Region, method = "ML", data=dd_C) 
summary(m1)
anova(m1)

dd_C <- dd_E %>% filter(comdyn=='E')
m1<-lm(value ~ Elev_diff, data=dd_C) 
summary(m1) 

##Turnover
dd_C <- dd_E %>% filter(comdyn=='T')
m1<-lme(value ~ Elev_diff, random = ~1|Region, method = "ML", data=dd_C) 
summary(m1) 
anova(m1)

dd_C <- dd_E %>% filter(comdyn=='T')
m1<-lm(value ~ Elev_diff, data=dd_C) 
summary(m1) 



#### Extract year since establishment and compare ####
dd %>%
  mutate(comm_sim = map(comm, ~.x %>% mutate(Year_diff=max(Year - min(Year))) %>%
                                select(Elevation, Year_diff, ODT, destPlotID) %>% 
                                distinct())) %>% 
  mutate(dat = pmap(.l = list(C=colonisation, E=extinction, To=turnover), .f = function(C, E, To){
    C <- C %>% rename(value = appearance)
    E <- E %>% rename(value = disappearance)
    To <- To %>% rename(value = total)
    bind_rows('C' = C,'E'= E, 'T' = To, .id='comdyn')})) %>% 
  mutate(dat = map2(dat, comm_sim, ~left_join(.x, .y, by = "destPlotID"))) %>%
  mutate(dat = map(dat, ~.x %>% filter(ODT == "warmed"))) %>%
  unnest(dat) %>% 
  group_by(Region, Year_diff, comdyn) %>%
  nest() %>% 
  mutate(model = map(data, ~lm(value ~ Year, data = .x))) %>% 
  mutate(tidied = map(model, tidy)) %>% 
  unnest(tidied) %>%
  filter(term == 'Year') %>%
  ggplot(aes(x=Year_diff, y=estimate, color=comdyn)) + 
  geom_point() + 
  scale_colour_manual(values = colour_comdyn) + 
  geom_smooth(method = "nls", method.args = list(formula = "y~(a/(b+x))", algorithm = 'plinear', start=list(a=1,b=0.5)), aes(group=comdyn, color=comdyn), se=F) +
  TP_theme() + 
  labs(y="Absolute slope over time", x=expression(Delta*Time)) +
  labs(color="Process") 

## Compare distributions of change in these three

Fig <- dd %>%
  mutate(comm_sim = map(comm, ~.x %>% mutate(Year_diff=max(Year - min(Year))) %>%
                          select(Elevation, Year_diff, ODT, destPlotID) %>% 
                          distinct())) %>% 
  mutate(dat = pmap(.l = list(C=colonisation, E=extinction, To=turnover), .f = function(C, E, To){
    C <- C %>% rename(value = appearance)
    E <- E %>% rename(value = disappearance)
    To <- To %>% rename(value = total)
    bind_rows('C' = C,'E'= E, 'T' = To, .id='comdyn')})) %>% 
  mutate(dat = map2(dat, comm_sim, ~left_join(.x, .y, by = "destPlotID"))) %>%
  mutate(dat = map(dat, ~.x %>% filter(ODT == "warmed"))) %>%
  unnest(dat) %>% 
  group_by(Region, Year_diff, comdyn) %>%
  nest() %>% 
  mutate(model = map(data, ~lm(value ~ Year, data = .x))) %>% 
  mutate(tidied = map(model, tidy)) %>% 
  unnest(tidied) %>%
  filter(term == 'Year') 

Fig %>%
  ggplot(aes(x=estimate, fill=comdyn)) + 
  geom_density(alpha=0.6) + 
  scale_colour_manual(values = colour_comdyn) +
  TP_theme() + 
  scale_fill_manual(values=colour_comdyn) +
  geom_vline(xintercept = 0) +
  facet_grid(~comdyn, nrows=3) +
  labs(y="Number of sites", x=expression(Delta*Slope)) +
  labs(color="Process") 



#### Plot colonisation of invaders vs. local species ####
dd %>% filter(!Region %in% c("FR_Lautaret", "IN_Kashmir", "US_Colorado")) %>%
  mutate(comm_sim = map(comm, ~.x %>% mutate(Year_diff=max(Year - min(Year))) %>%
                          select(Elevation, Year_diff, ODT, destPlotID) %>% 
                          distinct())) %>% 
  mutate(dat = pmap(.l = list(C=colonisation, E=extinction, To=turnover), .f = function(C, E, To){
    C <- C %>% rename(value = appearance)
    E <- E %>% rename(value = disappearance)
    To <- To %>% rename(value = total)
    bind_rows('C' = C,'E'= E, 'T' = To, .id='comdyn')})) %>% 
  mutate(dat = map2(dat, comm_sim, ~left_join(.x, .y, by = "destPlotID"))) #%>%
  mutate(dat = map(dat, ~.x %>% filter(ODT == "warmed"))) %>%
  unnest(dat) %>% 
  group_by(Region, Year_diff, comdyn) %>%
  nest() %>% 
  unnest(dat) %>%
  ggplot(aes(x = Year, y = appearance)) + 
  scale_colour_manual(values = colour_odt) + 
  geom_line(aes(group=destPlotID, color=ODT), alpha=0.2) +
  geom_smooth(aes(group = ODT, color=ODT), method = "lm", se = F) +
  facet_wrap(~Region, nrow=2) +
  TP_theme() + 
  scale_x_continuous(breaks=c(2010,2013,2016)) + 
  labs(title = 'Colonisation over time', color = "Treatment") 



