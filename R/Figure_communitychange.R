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

#### A. COMMUNITY DISSIMILARITY ####

dat_fil <- dat2 %>% filter(!is.na(Rel_Cover)) %>% #4 NAs in DE_Susalps, correct that in cleaning file. We are missing data from 2016 from some sites (EB-FE and EB-BT and EB-GW)
  filter(!(Region=="DE_Susalps" & Year==2016))

# RUN PCA AND EXTRACT DISTANCES 
dd <- dat_fil %>% 
  select(Region, originSiteID, destSiteID, Treatment) %>% 
  distinct() %>% 
  filter(Treatment == "Warm") %>% 
  select(-Treatment) %>% 
  mutate(comm = pmap(.l = list(R = Region, O = originSiteID, D = destSiteID), .f = function(R, O, D){
    bind_rows(
      originControls = dat_fil %>% filter(Region == R, destSiteID == O, Treatment == "LocalControl"),
      destControls = dat_fil %>% filter(Region == R, destSiteID == D, Treatment == "LocalControl"),
      warmed =  dat_fil %>% filter(Region == R, destSiteID == D, Treatment == "Warm"),
      .id = "ODT")
  })) %>% 
  mutate(comm_wide = map(comm, ~{
    .x %>% select(ODT, Year, SpeciesName, Rel_Cover, destPlotID) %>% 
      pivot_wider(names_from = SpeciesName, values_from = Rel_Cover, values_fill = list(Rel_Cover = 0), values_fn = list(Rel_Cover = sum))
  })) %>% 
  mutate(comm_meta  = map(comm_wide, select, ODT, destPlotID, Year), 
         comm_spp = map(comm_wide, select, -ODT, -destPlotID, -Year)) %>% 
  select(-comm_wide) %>% 
  mutate(PCA = map(comm_spp, ~rda(sqrt(.x)))) %>% 
  mutate(scores = map(.x = PCA, ~scores(.x, choices=c(1,2), display='sites'))) %>%
  mutate(scores = map(.x = scores, ~data.frame(PCA1=.x[,1], PCA2=.x[,2]))) %>%
  mutate(scores = map2(.x = comm_meta, .y = scores, bind_cols))


dd2 <- dd %>% 
  select(-(comm:PCA)) %>% 
  unnest(scores) %>% 
  group_by(Region, originSiteID, destSiteID, ODT, Year) %>% 
  summarise_at(vars(matches("PC")), .funs = mean) %>% 
  group_by(Region, originSiteID, destSiteID, Year) %>% 
  nest() %>% 
  mutate(distances = map(data, ~dist(select(.x, matches("PC"))))) %>% 
  mutate(distances = map(distances, ~tibble(what = c("Low_High", "Low_TP", "High_TP"), dist = as.vector(.x)))) %>% 
  unnest(distances)


## Individual region PCA
#dest control, origin control, warmed
colour_odt <- c("#D1362F", "#27223C", "#DEB18B")
#shape_odt <- c(16,16,25)

p1 <-dd %>% filter(Region == "CH_Calanda") %>% 
  pmap(function(scores, Region, originSiteID, ...){
    scores %>% arrange(Year) %>% 
      ggplot(aes(x = PCA1, y = PCA2, colour  = ODT, scale_fill_manual(values = colour_otd), group = destPlotID)) +
      geom_point(mapping = aes(size = Year == min(Year))) +
      scale_colour_manual(values = colour_odt) + 
      scale_fill_manual(values = colour_odt) +
      geom_path() +
      coord_equal() +
      TP_theme() +
      labs(title = "Calanda, CH", size = "First Year", color = "Treatment") 
  })

p2 <-dd %>% filter(Region == "US_Montana") %>% 
  pmap(function(scores, Region, originSiteID, ...){
    scores %>% arrange(Year) %>% 
      ggplot(aes(x = PCA1, y = PCA2, colour  = ODT, scale_fill_manual(values = colour_otd), group = destPlotID)) +
      geom_point(mapping = aes(size = Year == min(Year))) +
      scale_colour_manual(values = colour_odt) + 
      scale_fill_manual(values = colour_odt) +
      geom_path() +
      coord_equal() +
      TP_theme() +
      labs(title = "Montana, US", size = "First Year", color = "Treatment") 
  })

p3 <-dd %>% filter(Region == "CN_Gongga") %>% 
  pmap(function(scores, Region, originSiteID, ...){
    scores %>% arrange(Year) %>% 
      ggplot(aes(x = PCA1, y = PCA2, colour  = ODT, scale_fill_manual(values = colour_otd), group = destPlotID)) +
      geom_point(mapping = aes(size = Year == min(Year))) +
      scale_colour_manual(values = colour_odt) + 
      scale_fill_manual(values = colour_odt) +
      geom_path() +
      coord_equal() +
      TP_theme() +
      labs(title = "Gongga, CN", size = "First Year", color = "Treatment") 
  })

p1a <- p1[[1]]
p2a <- p2[[1]]
p3a <- p3[[1]]

# p1a + p2a + p3a
library(ggpubr)
ggarrange(p1a, p2a,p3a, ncol=3, nrow=1, common.legend = TRUE, legend="bottom", align = 'hv')

### Plot centroid distance over time ####
#colour_cd <- c("#A92420", "darkgrey", "#016367")
colour_cd <- c("#49BEB7", "black", "black")

dd2 %>%  filter(!Region %in% c("FR_Lautaret", "US_Colorado", "IN_Kashmir")) %>%
  filter(what != "Low_High") %>%
  ggplot(aes(x = Year, y = dist, color = what)) + 
  TP_theme() +
  geom_point() +
  scale_colour_manual(values = colour_cd) + 
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ Region) + 
  labs(color = "Treatment Comparisons", y="Distance between centroids", x='Year') 


## Plot by duration of time on one graph
colour_cd <- c("#49BEB7", "black", "black")
ddcent <- dd2 %>%  filter(!Region %in% c("FR_Lautaret", "US_Colorado", "IN_Kashmir")) %>%
  group_by(Region) %>%
  mutate(Year_0 = Year-min(Year)) 

ddcent %>% filter(what != "Low_High") %>%
  ggplot(aes(x = Year_0, y = dist, color = what, group=interaction(Region, what))) + 
  TP_theme() +
  geom_point() +
  scale_colour_manual(values = colour_cd) + 
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(~what) + 
  labs(color = "Treatment Comparisons", y="Distance between centroids", x='Duration (years)') 

ddcent %>% group_by(Region, originSiteID, destSiteID, what) %>%
  do(fitdd = tidy(lm(dist ~ Year_0, data = .))) %>% 
  unnest(fitdd) %>%
  filter(term=="Year_0", what!="Low_High") %>%
  group_by(Region, what) %>%
  summarize(mid = mean(estimate, na.rm=T), sd=sd(estimate, na.rm=T)) %>%
  ggplot() + 
  scale_colour_manual(values = colour_odt) + 
  scale_size_manual(values = c(4,2,2)) +
  geom_point(aes(y=what, x=mid, color=what, size=what)) +
  geom_linerange(aes(xmin=mid-sd, xmax=mid+sd, 
                     y=what, color=what)) +
  geom_vline(aes(xintercept=0)) +
  facet_wrap(~Region, nrow=3) +
  TP_theme() + 
  labs(color = "Treatment", y = 'Colonisation rate') +
  labs(title = 'Colonisation rate', color = "Treatment") 

ddcent %>% group_by(Region, originSiteID, destSiteID, what) %>%
  do(fitdd = tidy(lm(dist ~ Year_0, data = .))) %>% 
  unnest(fitdd) %>%
  filter(term=="Year_0", what!="Low_High") %>%
  group_by(Region, what) %>%
  summarize(mid = mean(estimate, na.rm=T), sd=sd(estimate, na.rm=T)) %>%
  ggplot() + 
  scale_colour_manual(values = colour_odt) + 
  scale_fill_manual(values = colour_odt) + 
  geom_density(aes(x=mid, color=what, fill=what)) +
  geom_vline(aes(xintercept=0)) +
  TP_theme() + 
  labs(color = "Treatment", y = 'Colonisation rate') +
  labs(title = 'Colonisation rate', color = "Treatment") 

#### B. POPULATIONDYNAMICS ####

# Use test from C_D.R
#Graph A: Turnover (all)
CO <- test %>% filter(!Region %in% c("FR_Lautaret", "IN_Kashmir", "US_Colorado", "CH_Calanda", "CH_Calanda2")) %>% 
  select(Region, originSiteID, destSiteID, comm_inv, colonisation, colonisation_inv, colonisation_res) %>%
  mutate(comm_sim = map(comm_inv, ~.x %>% select(destPlotID) %>% distinct())) %>%
  mutate(dat = map2(colonisation, comm_sim, ~left_join(.x, .y, by = "destPlotID"))) %>%
  mutate(dat2 = pmap(.l=list(a=dat, b=colonisation_inv, c=colonisation_res), function(a,b,c) bind_rows(list(a, b,c), .id="sp_pool"))) %>%
  unnest(dat2) %>%
  mutate(type = recode(sp_pool, "1"="all", "2"="invader", "3"="resident")) 

CO %>% 
  ggplot(aes(x = Year, y = appearance)) + 
  scale_colour_manual(values = colour_odt) + 
  geom_point(aes(group=destPlotID, color=type), alpha=0.05) +
  geom_smooth(aes(group = type, color=type), method = "lm", se = F) +
  facet_wrap(~Region, nrow=3) +
  TP_theme() + 
  scale_x_continuous(breaks=c(2010,2013,2016)) + 
  labs(title = 'Colonisation over time', color = "Treatment") 

CO %>% filter(type != "all") %>%
  ggplot(aes(x = Year, y = appearance, color = type)) + 
  TP_theme() +
  scale_x_continuous(breaks=c(2010,2013,2016, 2019)) + 
  geom_line(aes(alpha=0.8)) +
  scale_colour_manual(values = colour_odt) + 
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(~type) + 
  labs(color = "Treatment Comparisons", y="Turnover", x='Duration (years)') 

CO_m <- CO %>% mutate(Gradient = paste(originSiteID, destSiteID, sep="_")) %>%
  mutate(col = 0.00001+(1-2*0.00001)*(appearance-min(appearance))/(max(appearance)-min(appearance)))  %>%
  mutate(Years = scale(Year)) #%>%
  #filter(type != "all")

mC = brm(
  col ~ 0 + type + type:Years + (1|Region), 
  data = CO_m,
  family=Beta(link = "logit", link_phi = "log") # data scaled between 0 and 1
)

plot(mC)
mC

CO_m %>%
  data_grid(type=type, Years = seq_range(Years, n = 101)) %>%
  add_predicted_draws(mC, scale = "response", re_formula = NA) %>%
  ggplot(aes(x = Years, y = col, color = ordered(type), fill = ordered(type))) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50), alpha = 1/4) +
  geom_point(data = CO_m) +
  facet_grid(~type) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2")

antilogit = function(x) 1 / (1 + exp(-x))
antilogit(0.3)
antilogit(0.17)
antilogit(0.21)

pl <-posterior_linpred(mC)
str(pl)
plot(pl)

lme(col ~ type*Years, random ~Region, data=CO_m)

#### RATE INTERVAL ####
# Create ODT and nest communities
ddrar <- dat2 %>% select(Region, originSiteID, destSiteID, Treatment) %>% 
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
         turnover = map(comm, ~turnover(.x, time.var= "Year", species.var= "SpeciesName", abundance.var= "Rel_Cover", replicate.var="destPlotID", metric = "total")),
         ratechange = map(comm, ~ {.} %>% mutate(plot_odt = paste(ODT, destPlotID, sep="___")) %>% rate_change(., time.var= "Year", species.var= "SpeciesName", abundance.var= "Rel_Cover", replicate.var="plot_odt")),
         rateinterval = map(comm, ~ {.} %>% mutate(plot_odt = paste(ODT, destPlotID, sep="___")) %>% rate_change_interval(., time.var= "Year", species.var= "SpeciesName", abundance.var= "Rel_Cover", replicate.var="plot_odt")))

# Graph A

colour_odt <- c("#D1362F", "#27223C", "#FCD16B")

RC <- ddrar %>% select(Region, originSiteID, destSiteID, ratechange) %>%
  mutate(rates = map(ratechange, ~ {.} %>% separate(col=plot_odt, into=c("ODT", "destPlotID"), sep = "___"))) %>% #sep = "([.?:])"
  unnest(rates) 

RI <- ddrar %>% select(Region, originSiteID, destSiteID, rateinterval) %>%
  mutate(rates = map(rateinterval, ~ {.} %>% separate(col=plot_odt, into=c("ODT", "destPlotID"), sep = "___"))) %>% #sep = "([.?:])"
  unnest(rates) 
  

RI %>%
  ggplot(aes(x = interval, y = distance)) + 
  scale_colour_manual(values = colour_odt) + 
  #geom_line(aes(group=destPlotID, color=ODT), alpha=0.2) +
  geom_smooth(aes(group = ODT, color=ODT), method = "lm", se = F) +
  facet_wrap(~Region, nrow=2) +
  TP_theme() + 
  scale_x_continuous(breaks=c(2010,2013,2016)) +
  labs(color = "Treatment", y = 'Species Richness') +
  labs(title = 'Species Richness over time', color = "Treatment") 


RC %>%
  ggplot(aes(x=rate_change)) + 
  scale_fill_manual(values = colour_odt) + 
  scale_color_manual(values = colour_odt) + 
  geom_density(aes(fill=ODT, color=ODT), alpha=0.5) +
  #facet_wrap(~Region, nrow=2) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey", size=3) +
  TP_theme() 
  # labs(color = "Treatment", y = 'Species Richness') +
  # labs(title = 'Species Richness over time', color = "Treatment") 
