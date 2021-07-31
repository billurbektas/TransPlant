dat_fil <- dat %>% filter(!is.na(Rel_Cover)) %>% #4 NAs in DE_Susalps, correct that in cleaning file. We are missing data from 2016 from some sites (EB-FE and EB-BT and EB-GW)
  filter(!(Region=="DE_Susalps" & Year==2016))

source('R/theme_ggplot.R')
library(LearnGeom) #for trig functions

#### RUN PCA AND EXTRACT DISTANCES ####
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

#### RUN PCA AND EXTRACT STANDARDIZED DISTANCES ####
# Get single point per region

dd2 <- dd %>% 
  select(-(comm:PCA)) %>% 
  unnest(scores) %>% 
  group_by(Region, originSiteID, destSiteID, ODT, Year) %>%
  summarise_at(vars(matches("PC")), .funs = mean) %>% 
  group_by(Region, originSiteID, destSiteID) %>% 
  filter(Year == max(Year)) %>%
  select(-Year) %>%
  nest() %>% 
  #mutate(distances = map(data, ~dist(select(.x, matches("PC"))))) %>% 
  #mutate(distances = map(distances, ~tibble(what = c("Low_High", "Low_TP", "High_TP"), dist = as.vector(.x)))) #%>% 
  mutate(point = map(data, function(x) {
    a=as.numeric(x[2,2:3]) # a vertice (origin)
    b=as.numeric(x[1,2:3]) # b vertice (dest)
    c=as.numeric(x[3,2:3]) # c vertice (warmed)
    dis <- DistancePoints(a, b) # distance between origin and dest
    # Line <- CreateLinePoints(a, b)
    # pro <- ProjectPoint(c, Line)
    c <- abs(c/dis) #the length of intersecting line between c and the a-b segment
    # pro <- abs(pro/dis)
    ppro <- c[1] # the projected intersection point of c along that a-b segment
    return(tibble(Treatment = c("origin", "destination", "warmed", "origin", "destination", "proj_warmed"), value = c(1,1,1,2,2,2), Xcoord = c(0, 1, c[1], 0, c[1], c[1]), Ycoord = c(0, 0, c[2], 0, 0, c[2])))
  })) %>% 
  unnest(point)

# Plot triangles for each site
dd2 %>%
  ggplot(aes(x = Xcoord, y = Ycoord)) + 
  TP_theme() +
  geom_polygon(aes(fill = as.factor(value), group = as.factor(value))) +
  scale_fill_manual(values = c("light grey", "blue")) +
  facet_wrap(~Region)  +
  xlim(0,1) 

# Get point per plot

dd3 <- dd %>% 
  #Add species richness per warmed plots (average of first year)
  mutate(SR_ave = map(comm, function(x) x %>% 
                        filter(ODT=="warmed") %>%
                        group_by(destPlotID) %>%
                        summarise(SR = n_distinct(SpeciesName)))) %>% 
  #Calculate bray-curtis dissimilarity between low and high species pools
  mutate(BC = map(comm, function(x) x %>% 
                        filter(!ODT %in% c("warmed"), Year==min(Year)) %>%
                        select(ODT, SpeciesName, Rel_Cover) %>%
                        group_by(ODT, SpeciesName) %>%
                        dplyr::summarize(Rel_Cover = sum(Rel_Cover, na.rm=T)) %>%
                        ungroup() %>%
                        mutate(Rel_Cover = ifelse(Rel_Cover>0, 1, 0)) %>%
                        pivot_wider(names_from = SpeciesName, values_from = Rel_Cover, values_fill = 0) %>%
                        select(-ODT) %>%
                        vegdist(., method="bray", na.rm=T) %>%
                        as.numeric(.))) %>%
  select(-(comm:PCA)) %>% 
  mutate(odt_ave = map(scores, function(x) x %>% 
                         group_by(ODT, Year) %>%
                         summarise_at(vars(matches("PC")), .funs = mean) %>% 
                         filter(Year == max(Year)) %>%
                         select(-Year))) %>% 
  mutate(dis = map(odt_ave, function(x) {
    a=as.numeric(x[2,2:3]) #origin
    b=as.numeric(x[1,2:3]) #destination
    c=as.numeric(x[3,2:3]) #warmed plots
    dis = DistancePoints(a, b)
    return(dis)}),
  avepoints = map(odt_ave, function(x) {
      a=as.numeric(x[2,2:3]) #origin
      b=as.numeric(x[1,2:3]) #destination
      c=as.numeric(x[3,2:3]) #warmed plots
      dis = DistancePoints(a, b)
      c <- abs(c/dis) #the length of intersecting line between c and the a-b segment
      return(tibble(Treatment = c("origin", "destination", "warmed"), PCA1 = c(0, 1, c[1]), PCA2 = c(0, 0, c[2])))
    })) %>%
  unnest(dis) %>%
  mutate(odt_s = map2(scores, dis, function(.x,.y) {
    newpoint = .x %>% 
      mutate(PCA1 = abs(as.numeric(PCA1)/.y), PCA2=abs(as.numeric(PCA2)/.y)) %>%
      filter(ODT == "warmed", Year == max(Year)) %>%
      select(destPlotID, PCA1, PCA2)
    return(newpoint)
  })) %>%
  select(-scores, -odt_ave, -dis) %>%
  unnest(BC)  # unnest site-level data
  

# Get plot-level data
dd_ind_odt <- dd3 %>% select(-avepoints, -SR_ave) %>%
  unnest(odt_s) #unnest plot-scale obs, however we may have >1 as the scaling is between average values per site, not plot-level (because this is impossible)

dd_ind_SR <- dd3 %>% select(Region:SR_ave) %>% unnest(SR_ave) 

dd_ind <- left_join(dd_ind_odt, dd_ind_SR)

# Get average per Site to overplot
dd_ind_ave <- dd_ind %>% group_by(Region, originSiteID, destSiteID) %>% summarize(pca1 = mean(PCA1), pca2=mean(PCA2), se1 = sd(PCA1)/sqrt(n()), se2 = sd(PCA2)/sqrt(n()))

dd_group <- dd3 %>% unnest(avepoints) # group obs

ggplot(dd_ind, aes(x = PCA1, y = PCA2, color = Region)) +
  geom_point(alpha = .4) +
  geom_point(data = dd_group, size = 4) +
  theme_bw() +
  guides(color = guide_legend("Region")) +
  labs(
    x = "Distance from origin to destination",
    y = "Distance to point"
  )

#### ADD CLIMATE DATA ####

# Read in temperature and compute cumulative temp for ordering graph

temp <- read.csv('./climate/worlclim2_processedclimate.csv')

#Calculate annual cumulative warming
Temp <- temp %>% mutate(YearRange = YearMax-YearEstablished) %>%
  select(Gradient, destSiteID, T_ann_cor, T_sum_cor, V_ann, P_ann, YearRange, PlotSize_m2) 

# Get unique sets of gradients
sites <- dd_ind %>% select(Region, originSiteID, destSiteID) %>% distinct()

#Join temp data to originSiteID
df1 <- left_join(sites, Temp, by=c("Region"="Gradient", "originSiteID"="destSiteID")) %>% 
  mutate(T_ann_origin = T_ann_cor, T_sum_origin = T_sum_cor, V_ann_origin = V_ann, P_ann_origin = P_ann) %>% 
  select(-T_ann_cor, -T_sum_cor, -V_ann, -P_ann, -YearRange, -PlotSize_m2)

#Join temp data to destSiteID
df2 <- left_join(df1, Temp, by=c("Region"="Gradient", "destSiteID"="destSiteID")) %>% 
  mutate(T_ann_dest = T_ann_cor, T_sum_dest = T_sum_cor, V_ann_dest = V_ann, P_ann_dest = P_ann) %>%
  select(-T_ann_cor, -T_sum_cor, -V_ann, -P_ann)

#Calculate origin - destSiteID temps and *Year for cumulative 
df3 <- df2 %>% select(Region,  originSiteID, destSiteID, T_ann_dest, T_ann_origin, T_sum_dest, T_sum_origin, P_ann_dest, P_ann_origin, YearRange, PlotSize_m2) %>% 
  mutate(T_warm = T_ann_dest-T_ann_origin, P_warm = P_ann_dest-P_ann_origin, T_swarm = T_sum_dest-T_sum_origin) %>% #US_Arizona has cooler temps for annual in low site.
  mutate(T_warm_cum = T_warm*YearRange, P_warm_cum = P_warm*YearRange, T_swarm_cum = T_swarm*YearRange)

dd_clim_sum <- df3 %>% 
  select(Region, originSiteID, destSiteID, T_warm, T_warm_cum, T_swarm, T_swarm_cum, P_warm, P_warm_cum, YearRange, PlotSize_m2) %>%
  mutate(Gradient = paste(originSiteID, destSiteID, sep='_')) %>%
  mutate(T_warm = abs(T_warm), T_warm_cum=abs(T_warm_cum), T_swarm = abs(T_swarm), T_swarm_cum = abs(T_swarm_cum)) %>% # to fix issue of GW and FE giving - values to BT (issue with BT coords?)
  mutate(T_warms = scale(T_warm), T_warm_cums = scale(log(T_warm_cum)), T_swarms = scale(log(T_swarm)), T_swarm_cums = scale(log(T_swarm_cum)), P_warms = scale(P_warm), P_warm_cums = scale(P_warm_cum),
         YearRanges = scale(YearRange), PlotSize_m2s = scale(PlotSize_m2))

# Standardize response (PCA1) between 0 and 1 (0,1) for beta distribution
dd_ind_stand <- dd_ind %>% mutate(PCA1 = 0.00001+(1-2*0.00001)*(PCA1-min(PCA1))/(max(PCA1)-min(PCA1)),
                                  BCs = scale(BC), 
                                  SRs = scale(SR))

# Merge climate data to dd_ind (standardized) based on Region, originSiteID and destSiteID
dd_plot <- left_join(dd_ind_stand, dd_clim_sum)

# Get Region averages for dd_ind and dd_clim_sum
dd_ind_ave <- dd_ind_stand %>% 
  select(Region, originSiteID, destSiteID, destPlotID, BC, SR, BCs, SRs, PCA1, PCA2) %>%
  group_by(Region) %>%
  summarize_each(BC:PCA2, funs=c(mean, function(x) sd(x)/sqrt(n()))) 

dd_clim_ave <- dd_clim_sum %>% 
  select(Region, T_warm:PlotSize_m2, T_warms:PlotSize_m2s) %>%
  group_by(Region) %>%
  summarize_each(T_warm:PlotSize_m2s, funs=c(mean))

dd_site <- left_join(dd_ind_ave, dd_clim_ave)



#### PLOT RESULTS ####

## All transplanted turfs
dd_plot %>%
  group_by(Gradient) %>%
  data_grid(T_warm_cums = seq_range(T_warm_cums, n = 51)) %>%
  add_fitted_draws(m4) %>%
  ggplot(aes(x = T_warm_cums, y = PCA1)) +
  stat_lineribbon(aes(y = .value)) +
  geom_point(data = dd_Temp_sum) +
  scale_fill_brewer(palette = "Greys") +
  #scale_color_brewer(palette = "Set2") +
  TP_theme() 

## Region- scale using climate data
#### PLOT AVERAGE SITE VALUES BY TEMPERATURE ####

#order by Temp (black to red) and P_warm (black to blue) <- VPD doesn't vary much
#dd_clim_ave$SiteID <- dd_clim_ave$Gradient

dd_site$YearRange[2] <- 4
dd_site$SiteID <- dd_site$Region
TS <- dd_site$T_swarm_cums
TS[2] <- TS[1]
TSA <- dd_site$T_swarms
TSA[2] <- TSA[1]
PR <- dd_site$P_warm_cums
PR[2] <- PR[1]
PRA <- dd_site$P_warms
PRA[2] <- PRA[1]

TS <- (TS-min(TS)) / (max(TS)-min(TS))
TSA <- (TSA-min(TSA)) / (max(TSA)-min(TSA))
PR <- (PR-min(PR)) / (max(PR)-min(PR))
PRA <- (PRA-min(PRA)) / (max(PRA)-min(PRA))

TScol <- rgb(TS*0.9, TS/2, 1-TS*0.9) #R, G, B colors, between 0 to 1
TSAcol <- rgb(TSA, 0.3, 1-TSA)
PRcol <- rgb(0, 0, PR) 
PRAcol <- rgb(0, 0, PRA) 

names(TScol) <- as.character(dd_site$SiteID)
names(TSAcol) <- as.character(dd_site$SiteID)
names(PRcol) <- as.character(dd_site$SiteID)
names(PRAcol) <- as.character(dd_site$SiteID)

# add in site names
#order by T_warm_cum
dd_site2 <- dd_site %>% arrange(T_warm_cum)
#dd_site2$SiteID <- c('COL', 'SWE', 'MON', 'ARI', 'FRA1', 'SWI1', 'NOR1', 'FRA2', 'GER1', 'CHI1', 'CHI2', 'IND', 'CHI3', 'NOR2', 'SWI2', 'NOR3', 'NOR4', 'GER2', 'ITA')
colfunc <- colorRampPalette(c("#DEB18B", "#AB2330"))
cols <- colfunc(22)
cols <- data.frame(SiteID=dd_site2$SiteID, cols=cols)


#### For change in summer temp (dest-origin) 
p1 <-  dd_site %>% arrange(T_swarm_cums) %>%    
  mutate(SiteID=factor(SiteID, levels=SiteID)) %>%
  ggplot(aes(x = PCA1_fn1, y = PCA2_fn1, group=SiteID, fill = SiteID)) +
  geom_errorbar(aes(xmin=PCA1_fn1-PCA1_fn2, xmax=PCA1_fn1+PCA1_fn2), width=.01, alpha=1, col=cols$cols) +
  geom_errorbar(aes(ymin=PCA2_fn1-PCA2_fn2, ymax=PCA2_fn1+PCA2_fn2), width=.01, alpha=1, col=cols$cols) +
  geom_point(alpha = 1, size=5, colour="black",pch=21) + 
  scale_fill_manual(values=cols$cols) +
  TP_theme() +
  xlim(0,1) +
  ylim(0,1) +
  guides(color = guide_legend("Cumulative warming (degC)")) +
  labs(
    x = "Distance from origin to destination",
    y = "Deviation from origin-destination"
  )

p2 <-  dd_site %>% arrange(T_swarm_cums) %>%
  mutate(SiteID=factor(SiteID, levels=SiteID)) %>%
  ggplot(aes(x = SiteID, y=PCA2_fn1, fill = SiteID)) +
  geom_errorbar(aes(ymin=PCA2_fn1-PCA2_fn2, ymax=PCA2_fn1+PCA2_fn2), width=.01, alpha=1, col=cols$cols) +
  geom_point(alpha = 1, size=5, colour="black",pch=21) + 
  scale_colour_manual(values=cols$cols) +
  TP_theme() +
  ylim(0,1) +
  labs(
    x = "Region",
    y = ""
  )

p3 <- dd_site %>% arrange(T_swarm_cums) %>%
  mutate(SiteID=factor(SiteID, levels=SiteID)) %>%
  ggplot(aes(x = PCA1_fn1, y=SiteID, color = SiteID)) +
  geom_point(alpha = 3, size=3) +
  geom_errorbar(aes(xmin=PCA1_fn1-PCA1_fn2, xmax=PCA1_fn1+PCA1_fn2), width=.01) +
  scale_colour_manual(values=TScol) +
  TP_theme() +
  xlim(0,1) +
  labs(
    x = "Distance from origin to destination",
    y = "Region"
  )


library(cowplot)

legend <- get_legend(
  # create some space to the left of the legend
  p1 + theme(legend.position = "right")
)


plot_grid(p1 + theme(legend.position="none"), 
          legend,
          p3 + theme(legend.position="none"), 
          p2 +  theme(legend.position="none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
          align = "vh", nrow = 2, 
          rel_heights = c(1/2, 1/2, 1/2),
          rel_widths = c(1/2, 1/2, 1/2),
          labels = c('A', 'Legend', 'B', 'C'))


## Add in temperature or SR or distance between groups correlations
p2 <-  dd_site %>% arrange(T_swarm_cums) %>%
  mutate(SiteID=factor(SiteID, levels=SiteID)) %>%
  ggplot(aes(x = BCs_fn1, y=PCA1_fn1, fill = SiteID)) +
  geom_errorbar(aes(ymin=PCA1_fn1-PCA1_fn2, ymax=PCA1_fn1+PCA1_fn2), width=.01, alpha=1, col=cols$cols) +
  geom_point(alpha = 1, size=5, colour="black",pch=21) +
  scale_fill_manual(values=cols$cols) +
  TP_theme() +
  #ylim(0,1) +
  labs(
    x = "BC Dissimilarity",
    y = "Distance from origin to destination"
  )

p3 <- dd_site %>% arrange(T_swarm_cums) %>%
  mutate(SiteID=factor(SiteID, levels=SiteID)) %>%
  ggplot(aes(x = T_swarms, y=PCA1_fn1, fill = SiteID)) +
  geom_errorbar(aes(ymin=PCA1_fn1-PCA1_fn2, ymax=PCA1_fn1+PCA1_fn2), width=.01, alpha=1, col=cols$cols) +
  geom_point(alpha = 1, size=5, colour="black",pch=21) + 
  scale_fill_manual(values=cols$cols) +
  TP_theme() +
  #ylim(0,1) +
  labs(
    x = "Average Temperature (degC)",
    y = "Distance from origin to destination"
  )

p4 <- dd_site %>% mutate(SiteID=factor(SiteID, levels=SiteID)) %>%
  ggplot(aes(x = YearRange, y=PCA1_fn1, fill = SiteID)) +
  geom_smooth(aes(x=YearRange, y=PCA1_fn1)) + 
  geom_errorbar(aes(ymin=PCA1_fn1-PCA1_fn2, ymax=PCA1_fn1+PCA1_fn2), width=.01, alpha=1, col=cols$cols) +
  geom_point(alpha = 1, size=5, colour="black",pch=21) + 
  scale_fill_manual(values=cols$cols) +
  TP_theme() +
  #ylim(0,1) +
  labs(
    x = "Time (years)",
    y = "Distance from origin to destination"
  )


plot_grid(p1 + theme(legend.position="none"), 
          legend,
          p2 + theme(legend.position="none"), 
          p3 +  theme(legend.position="none"),
          align = "vh", nrow = 2, 
          rel_heights = c(1/2, 1/2, 1/2),
          rel_widths = c(1/2, 1/2, 1/2),
          labels = c('A', 'Legend', 'B', 'C'))

#### RUNNING BRMS MODEL AGAINST ALL PLOTS (REGION RANDOM EFFECT) ####
library(brms)
library(modelr)
library(tidybayes)

m1 = brm(
  PCA1 ~ T_warm_cums + P_warm_cums + (1|Gradient), 
  data = dd_plot,
  iter = 4000,
  family=Beta(link = "logit", link_phi = "log") # data scaled between 0 and 1
)

m2 = brm(
  PCA1 ~ YearRange + T_warms + P_warms + (1|Gradient), 
  data = dd_plot,
  iter = 4000,
  family=Beta(link = "logit", link_phi = "log") # data scaled between 0 and 1
)

m3 = brm(
  PCA1 ~ YearRange + (1|Gradient), 
  data = dd_plot,
  iter = 4000,
  family=Beta(link = "logit", link_phi = "log") # data scaled between 0 and 1
)

m4 = brm(
  PCA1 ~ T_warm_cums + (1|Gradient), 
  data = dd_plot,
  iter = 4000,
  family=Beta(link = "logit", link_phi = "log") # data scaled between 0 and 1
)

m5 = brm(
  PCA1 ~ SR + (1|Gradient), 
  data = dd_plot,
  iter = 4000,
  family=Beta(link = "logit", link_phi = "log") # data scaled between 0 and 1
)

m6 = brm(
  PCA1 ~ BC + (1|Gradient), 
  data = dd_plot,
  iter = 4000,
  family=Beta(link = "logit", link_phi = "log") # data scaled between 0 and 1
)

m7 = brm(
  PCA1 ~ YearRange + SRs + BCs + PlotSize_m2s + T_swarms  + P_warms + T_swarm_cums + P_warm_cums + (1|Region:Gradient), 
  data = dd_plot,
  #iter = 4000,
  family=Beta(link = "logit", link_phi = "log") # data scaled between 0 and 1
)

summary(m7)


#### BRMS MODEL PCA1 ~ ALL CLIMATE DATA
# for gaussian use: PCA1_fn1|resp_se(PCA1_fn2, sigma = TRUE)

dd_site[2,14:19] <- dd_site[1,14:19]

dd_site[2,21:29] <- dd_site[1,21:29]

m1 = brm(
  PCA1_fn1|mi(PCA1_fn2) ~ PlotSize_m2s + YearRanges, # plot size and year range no effect
  data = dd_site,
  #iter = 10000,
  #cores=3, 
  #chains=3,
  #control = list(adapt_delta=0.99, max_treedepth = 20), 
  family=Beta(link = "logit", link_phi = "log") 
)

m2 = brm(
  PCA1_fn1|mi(PCA1_fn2) ~ YearRanges + T_warms + P_warms, # none significant
  data = dd_site,
  family=Beta(link = "logit", link_phi = "log") 
)

m3 = brm(
  PCA1_fn1|mi(PCA1_fn2) ~ T_warm_cums + P_warm_cums, # t-warm sig, but negative !?
  data = dd_site,
  family=Beta(link = "logit", link_phi = "log") 
)

m4 = brm(
  PCA1_fn1|mi(PCA1_fn2) ~  T_warm + P_warm + T_warm_cum + P_warm_cum, # use 0+ if unscaled
  data = dd_site,
  family=Beta(link = "logit", link_phi = "log") 
)


m5 = brm(
  PCA1_fn1|mi(PCA1_fn2) ~  PlotSize_m2s_fn1*T_warm_cums_fn1 + (1|Gradient), # no interaction of plotsize, but negative tcum
  data = dd_site,
  iter=8000,
  control = list(adapt_delta=0.99, max_treedepth = 15), 
  family=Beta(link = "logit", link_phi = "log") 
)

m6 = brm(
  PCA1_fn1|mi(PCA1_fn2) ~  SRs_fn1 + BCs_fn1 + PlotSize_m2s + T_swarms + P_warms + T_swarm_cums + P_warm_cums, # no interaction of plotsize, but negative tcum
  data = dd_site,
  family=Beta(link = "logit", link_phi = "log") 
)

m6 = brm(
  PCA1_fn1|mi(PCA1_fn2) ~  T_swarms, # BCs_fn1 positive, Tswarms negative
  data = dd_site,
  family=Beta(link = "logit", link_phi = "log") 
)

summary(m6) # BC +, Tswarms neg, Tswarm_cum

# Plotting for SNF report

mT = brm(
  PCA1_fn1|mi(PCA1_fn2) ~  T_swarms, # BCs_fn1 positive, Tswarms negative
  data = dd_site,
  family=Beta(link = "logit", link_phi = "log") 
)

newdf <- dd_site %>%
  distinct(PCA1_fn2, T_swarms) %>%
  mutate(Region = "fake")  %>%
  add_fitted_draws(mT, allow_new_levels = TRUE) %>%
  mutate(SiteID = "fake")

pT <- dd_site %>% arrange(T_swarm_cums) %>%
  mutate(SiteID=factor(SiteID, levels=SiteID)) %>%
  ggplot(aes(x = T_swarms, y=PCA1_fn1, col = SiteID)) +
  stat_lineribbon(data=newdf, aes(x=T_swarms, y=.value), fill= "#DEB18B", col = "#541F12", alpha=0.8, .width = c(.95)) +
  geom_errorbar(aes(ymin=PCA1_fn1-PCA1_fn2, ymax=PCA1_fn1+PCA1_fn2), width=.01, alpha=1, col=cols$cols) +
  geom_point(alpha = 1, size=5, col=cols$cols) + 
  TP_theme() +
  #ylim(0,1) +
  labs(
    x = expression("Average Regional Temperature " ( degree*C)),
    y = ""
  )


mB = brm(
  PCA1_fn1|mi(PCA1_fn2) ~ BCs_fn1, # BCs_fn1 positive, Tswarms negative
  data = dd_site,
  family=Beta(link = "logit", link_phi = "log") 
)

newdf <- dd_site %>%
  distinct(PCA1_fn2, BCs_fn1) %>%
  mutate(Region = "fake")  %>%
  add_fitted_draws(mB, allow_new_levels = TRUE) %>%
  mutate(SiteID = "fake")

pB <-  dd_site %>% arrange(T_swarm_cums) %>%
  mutate(SiteID=factor(SiteID, levels=SiteID)) %>%
  ggplot(aes(x = BCs_fn1, y=PCA1_fn1, fill = SiteID)) +
  stat_lineribbon(data=newdf, aes(x=BCs_fn1, y=.value), fill= "#DEB18B", col = "#541F12", alpha=0.8, .width = c(.95)) +
  geom_errorbar(aes(ymin=PCA1_fn1-PCA1_fn2, ymax=PCA1_fn1+PCA1_fn2), width=.01, alpha=1, col=cols$cols) +
  geom_point(alpha = 1, size=5, col=cols$cols) + 
  TP_theme() +
  #ylim(0,1) +
  labs(
    x = "Initial dissimilarity between high-low controls",
    y = ""
  )


mY = brm(
  PCA1_fn1|mi(PCA1_fn2) ~ YearRanges, # BCs_fn1 positive, Tswarms negative
  data = dd_site,
  family=Beta(link = "logit", link_phi = "log") 
)

newdf <- dd_site %>%
  distinct(PCA1_fn2, YearRanges) %>%
  mutate(Region = "fake")  %>%
  add_fitted_draws(mY, allow_new_levels = TRUE) %>%
  mutate(SiteID = "fake")

pY <- dd_site %>% mutate(SiteID=factor(SiteID, levels=SiteID)) %>%
  ggplot(aes(x = YearRanges, y=PCA1_fn1, fill = SiteID)) +
  #stat_lineribbon(data=newdf, aes(x=YearRanges, y=.value), fill= "#DEB18B", col = "#541F12", alpha=0.8, .width = c(.95)) +
  geom_errorbar(aes(ymin=PCA1_fn1-PCA1_fn2, ymax=PCA1_fn1+PCA1_fn2), width=.01, alpha=1, col=cols$cols) +
  geom_point(alpha = 1, size=5, col=cols$cols) + 
  TP_theme() +
  #ylim(0,1) +
  labs(
    x = "Time (years)",
    y = "Transplanted community composition convergence"
  )

prow <- plot_grid(pY + theme(legend.position="none"), 
          pT + theme(legend.position="none"), 
          pB +  theme(legend.position="none"),
          #legend,
          align = "vh", nrow = 1,
          hjust = -1,
          #rel_heights = c(1/2, 1/2, 1/2),
          #rel_widths = c(1/2, 1/2, 1/2),
          labels = c('A', 'B', 'C'))

legend_b <- get_legend(pB + theme(legend.position="bottom"))

# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
p
prow
  