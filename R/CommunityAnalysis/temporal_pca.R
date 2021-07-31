dat_fil <- dat %>% filter(!is.na(Rel_Cover)) %>% #4 NAs in DE_Susalps, correct that in cleaning file. We are missing data from 2016 from some sites (EB-FE and EB-BT and EB-GW)
  filter(!(Region=="DE_Susalps" & Year==2016))
source('R/theme_ggplot.R')
library("patchwork")

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

#### PLOT PCA ####
colour_otd <- c("orange", "blue","green3")
 
pmap(dd, function(scores, Region, originSiteID, ...){
  scores %>% arrange(Year) %>% 
    ggplot(aes(x = PC1, y = PC2, colour  = ODT, scale_fill_manual(values = colour_otd), group = destPlotID)) +
    geom_point(mapping = aes(shape = Year == min(Year))) +
    scale_colour_manual(values = colour_otd) + 
    geom_path() +
    coord_equal() +
    TP_theme()+
    labs(title = paste(Region, originSiteID), shape = "First Year")
  })%>%
    wrap_plots() +
  plot_layout(guides = 'collect')
  
## Individual region PCA
#dest control, origin control, warmed
colour_odt <- c("#A92420", "#016367", "#FBC00E")
#shape_odt <- c(16,16,25)

p1 <- dd %>% filter(Region == "DE_Susalps") %>% #Insert desired region name
  pmap(function(scores, Region, originSiteID, ...){
    scores %>% arrange(Year) %>% 
    ggplot(aes(x = PC1, y = PC2, colour  = ODT, scale_fill_manual(values = colour_otd), group = destPlotID)) +
      geom_point(mapping = aes(size = Year == min(Year))) +
      scale_colour_manual(values = colour_odt) + 
      scale_fill_manual(values = colour_odt) +
      geom_path() +
      coord_equal() +
      TP_theme() +
      labs(title = "IT_MatschMazia", size = "First Year", color = "Treatment") 
  })

# Plotting 3x plots for display
# p2a <- p2[[2]]
# p1a <- p1[[3]]
# p3a <- p3[[2]]
# p1a + p2a + p3a
# library(ggpubr)
# ggarrange(p1a, p2a, p3a, ncol=3, nrow=1, common.legend = TRUE, legend="bottom", align = 'hv')
#Next!: remove rare species by region, max cover e.g. 1% or n occur ==1


### Plot centroid distance over time ####
colour_cd <- c("#A92420", "darkgrey", "#016367")
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

library(nlme)
library(emmeans)
#Controls
dd_C <- ddcent %>% filter(what=='Low_High')
m1<-lme(dist ~ Year_0*Region, random = ~1|originSiteID, method = "ML", data=dd_C) 
summary(m1)
anova(m1) #*, -0.006

##To low
dd_C <- ddcent %>% filter(what=='Low_TP')
m1<-lme(dist ~ Year_0*Region, random = ~1|originSiteID, method = "ML", data=dd_C) 
summary(m1)
anova(m1) #***, -0.034
mnull<-lme(dist ~ Year_0+ Region, random = ~1|originSiteID, method = "ML", data=dd_C) 
anova(m1, mnull)
mnull1<-lme(dist ~ Year_0, random = ~1|originSiteID, method = "ML", data=dd_C) 
anova(mnull, mnull1)
mnull2<-lme(dist ~ Region, random = ~1|originSiteID, method = "ML", data=dd_C) 
anova(mnull, mnull2)


##To alpine
dd_C <- ddcent %>% filter(what=='High_TP')
m1<-lme(dist ~ Year_0*Region, random = ~1|originSiteID, method = "ML", data=dd_C) 
summary(m1) 
anova(m1) #***, 0.028
mnull<-lme(dist ~ Year_0+ Region, random = ~1|originSiteID, method = "ML", data=dd_C) 
anova(m1, mnull)
mnull1<-lme(dist ~ Year_0, random = ~1|originSiteID, method = "ML", data=dd_C) 
anova(mnull, mnull1)
mnull2<-lme(dist ~ Region, random = ~1|originSiteID, method = "ML", data=dd_C) 
anova(mnull, mnull2)



#TRY TO DO: 
#Order by duration
#change shape to up and down triangle

### Plot centroid distance over time lags...is this necessary? ####
#Gives an indication of directionality and rate of change

#### RUN PCA AND EXTRACT STANDARDIZED DISTANCES ####
# Should we do an NMDS for this? This is so we can rotate distances (instead of it being orientated to the maximum variation)
library(LearnGeom) #for trig functions

dd2 <- dd %>% 
  filter(!Region %in% c("IT_MatschMazia2")) %>%
  select(-(comm:PCA)) %>% 
  unnest(scores) %>% 
  #group_by(Region, originSiteID, destSiteID, ODT, Year) %>%
  group_by(Region, ODT, Year) %>%
  summarise_at(vars(matches("PC")), .funs = mean) %>% 
  group_by(Region) %>% 
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
  
# Plot sites on one graph with error bars

dd2 <- dd %>% 
  filter(!Region %in% c("IT_MatschMazia2")) %>%
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
      filter(ODT == "warmed") %>%
      select(destPlotID, PCA1, PCA2)
    return(newpoint)
  })) 

dd_ind <- dd2 %>% unnest(odt_s) #individual obs

dd_ind_ave <- dd_ind %>% group_by(Region, originSiteID, destSiteID) %>% summarize(pca1 = mean(PCA1), pca2=mean(PCA2), se1 = sd(PCA1)/sqrt(n()), se2 = sd(PCA2)/sqrt(n()))

dd_group <- dd2 %>% unnest(avepoints) # group obs

ggplot(dd_ind, aes(x = PCA1, y = PCA2, color = Region)) +
  geom_point(alpha = .4) +
  geom_point(data = dd_group, size = 4) +
  theme_bw() +
  guides(color = guide_legend("Region")) +
  labs(
    x = "Distance from origin to destination",
    y = "Distance to point"
  )

# This includes gradients (can remove those as well)
ggplot(dd_ind_ave, aes(x = pca1, y = pca2, color = Region)) +
  geom_point(alpha = 3, size=3) +    
  geom_errorbar(aes(ymin=pca2-se2, ymax=pca2+se2), width=.01) +
  geom_errorbar(aes(xmin=pca1-se2, xmax=pca1+se2), width=.01) +
  geom_line() +
  TP_theme() +
  xlim(0,1) +
  ylim(0,1) +
  guides(color = guide_legend("Region")) +
  labs(
    x = "Distance from origin to destination",
    y = "Distance to point"
  )

# If this needs to be changed to nmds, use monoMDS and vegdist(method="bray") to do this

# Read in temperature and compute cumulative temp for ordering graph

temp <- read.csv('./climate/worlclim2_processedclimate.csv')

#Calculate annual cumulative warming
Temp <- temp %>% select(Gradient, destSiteID, T_ann_cor, T_sum_cor, V_ann, P_ann, YearRange) 

dd_temp <- left_join(dd_ind_ave, Temp, by=c("Region"="Gradient", "originSiteID"="destSiteID")) %>% 
  mutate(T_ann_origin = T_ann_cor, T_sum_origin = T_sum_cor, V_ann_origin = V_ann, P_ann_origin = P_ann) %>% select(-T_ann_cor, -T_sum_cor, -V_ann, -P_ann, -YearRange)
dd_Temp <- left_join(dd_temp, Temp, by=c("Region"="Gradient", "destSiteID"="destSiteID")) %>% 
  mutate(T_ann_dest = T_ann_cor, T_sum_dest = T_sum_cor, V_ann_dest = V_ann, P_ann_dest = P_ann, T_warm = T_ann_dest-T_ann_origin, T_warm_cum = T_warm*YearRange, T_swarm = T_sum_dest-T_sum_origin, T_swarm_cum = T_swarm*YearRange, V_warm = V_ann_dest-V_ann_origin, V_warm_cum = V_warm*YearRange, P_warm = P_ann_dest-P_ann_origin, P_warm_cum = P_warm*YearRange) %>% 
  select(-T_ann_cor, -V_ann, -P_ann)

dd_clim_ave <- dd_Temp %>% group_by(Region) %>% summarize_each(pca1:P_warm_cum, funs=mean)

# Create colors by values
#order by Temp (black to red) and P_warm (black to blue) <- VPD doesn't vary much
dd_clim_ave$SiteID <- substr(dd_clim_ave$Region, 1,4)
TS <- log(dd_clim_ave$T_swarm_cum)
TS_ave <- log(dd_clim_ave$T_sum_cor)
PR <- log(dd_clim_ave$P_warm_cum+1)
PR_ave <- log(dd_clim_ave$P_ann_dest)
VPD <- log(dd_clim_ave$V_warm_cum+1)
VPD_ave <- log(dd_clim_ave$V_ann_dest)

TS <- (TS-min(TS)) / (max(TS)-min(TS))
TSA <- (TS_ave-min(TS_ave)) / (max(TS_ave)-min(TS_ave))
PR <- (PR-min(PR)) / (max(PR)-min(PR))
PRA <- (PR_ave-min(PR_ave)) / (max(PR_ave)-min(PR_ave))
VPD <- (VPD-min(VPD)) / (max(VPD)-min(VPD))
VPDave <- (VPD_ave-min(VPD_ave)) / (max(VPD_ave)-min(VPD_ave))
TScol <- rgb(TS, 0, 0) #R, G, B colors, between 0 to 1
TSAcol <- rgb(TSA, 0, 0)
PRcol <- rgb(0, 0, PR) 
PRAcol <- rgb(0, 0, PRA) 
Vcol <- rgb(0, 0, VPD) 
VAcol <- rgb(0, 0, VPDave) 
names(TScol) <- as.character(dd_clim_ave$SiteID)
names(TSAcol) <- as.character(dd_clim_ave$SiteID)
names(PRcol) <- as.character(dd_clim_ave$SiteID)
names(PRAcol) <- as.character(dd_clim_ave$SiteID)
names(Vcol) <- as.character(dd_clim_ave$SiteID)
names(VAcol) <- as.character(dd_clim_ave$SiteID)

#### For change in summer temp (dest-origin) 
p1 <-  dd_clim_ave %>% arrange(T_swarm_cum) %>%    
  mutate(SiteID=factor(SiteID, levels=SiteID)) %>%
  ggplot(aes(x = pca1, y = pca2, color = SiteID)) +
  geom_point(alpha = 3, size=3) +    
  geom_errorbar(aes(ymin=pca2-se2, ymax=pca2+se2), width=.01) +
  geom_errorbar(aes(xmin=pca1-se1, xmax=pca1+se1), width=.01) +
  geom_line() +
  scale_colour_manual(values=col1) +
  TP_theme() +
  xlim(0,1) +
  ylim(0,1) +
  guides(color = guide_legend("Region")) +
  labs(
    x = "",
    y = "Distance to point"
  )

p2 <-  dd_clim_ave %>% arrange(T_swarm_cum) %>%    
  mutate(SiteID=factor(SiteID, levels=SiteID)) %>%
  ggplot(aes(x = SiteID, y=pca2, color = SiteID)) +
  geom_point(alpha = 3, size=3) +   
  geom_errorbar(aes(ymin=pca2-se2, ymax=pca2+se2), width=.01) +
  geom_line() +
  scale_colour_manual(values=col1) +
  TP_theme() +
  ylim(0,1) +
  labs(
    x = "Region",
    y = ""
  )

p3 <- dd_clim_ave %>% arrange(T_swarm_cum) %>%    
  mutate(SiteID=factor(SiteID, levels=SiteID)) %>%
  ggplot(aes(x = pca1, y=SiteID, color = SiteID)) +
  geom_point(alpha = 3, size=3) +   
  geom_errorbar(aes(xmin=pca1-se1, xmax=pca1+se1), width=.01) +
  geom_line() +
  scale_colour_manual(values=col1) +
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
                    p2 + theme(legend.position="none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)), 
                    p3 +  theme(legend.position="none"), 
                    align = "vh", nrow = 2, 
                    rel_heights = c(1/2, 1/2, 1/2),
                    rel_widths = c(1/2, 1/2, 1/2),
                    labels = c('A', 'B', 'C'))

#### For annual precip (VPD and p_warm didn't show much)

p1 <-  dd_clim_ave %>% arrange(P_warm) %>%    
  mutate(SiteID=factor(SiteID, levels=SiteID)) %>%
  ggplot(aes(x = pca1, y = pca2, color = SiteID)) +
  geom_point(alpha = 3, size=3) +    
  geom_errorbar(aes(ymin=pca2-se2, ymax=pca2+se2), width=.01) +
  geom_errorbar(aes(xmin=pca1-se1, xmax=pca1+se1), width=.01) +
  geom_line() +
  scale_colour_manual(values=PRAcol) +
  TP_theme() +
  xlim(0,1) +
  ylim(0,1) +
  guides(color = guide_legend("Region")) +
  labs(
    x = "",
    y = "Distance to point"
  )

p2 <-  dd_clim_ave %>% arrange(P_ann_dest) %>%   
  mutate(SiteID=factor(SiteID, levels=SiteID)) %>%
  ggplot(aes(x = SiteID, y=pca2, color = SiteID)) +
  geom_point(alpha = 3, size=3) +   
  geom_errorbar(aes(ymin=pca2-se2, ymax=pca2+se2), width=.01) +
  geom_line() +
  scale_colour_manual(values=PRAcol) +
  TP_theme() +
  ylim(0,1) +
  labs(
    x = "Region",
    y = ""
  )

p3 <-  dd_clim_ave %>% arrange(P_ann_dest) %>%    
  mutate(SiteID=factor(SiteID, levels=SiteID)) %>%
  ggplot(aes(x = pca1, y=SiteID, color = SiteID)) +
  geom_point(alpha = 3, size=3) +   
  geom_errorbar(aes(xmin=pca1-se1, xmax=pca1+se1), width=.01) +
  geom_line() +
  scale_colour_manual(values=PRAcol) +
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
          p2 + theme(legend.position="none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)), 
          p3 +  theme(legend.position="none"), 
          align = "vh", nrow = 2, 
          rel_heights = c(1/2, 1/2, 1/2),
          rel_widths = c(1/2, 1/2, 1/2),
          labels = c('A', 'B', 'C'))


#### For VPD (dest-origin)

p1 <-  dd_clim_ave %>% arrange(V_warm_cum) %>%    
  mutate(SiteID=factor(SiteID, levels=SiteID)) %>%
  ggplot(aes(x = pca1, y = pca2, color = SiteID)) +
  geom_point(alpha = 3, size=3) +    
  geom_errorbar(aes(ymin=pca2-se2, ymax=pca2+se2), width=.01) +
  geom_errorbar(aes(xmin=pca1-se1, xmax=pca1+se1), width=.01) +
  geom_line() +
  scale_colour_manual(values=Vcol) +
  TP_theme() +
  xlim(0,1) +
  ylim(0,1) +
  guides(color = guide_legend("Region")) +
  labs(
    x = "",
    y = "Distance to point"
  )

p2 <-  dd_clim_ave %>% arrange(V_warm_cum) %>%   
  mutate(SiteID=factor(SiteID, levels=SiteID)) %>%
  ggplot(aes(x = SiteID, y=pca2, color = SiteID)) +
  geom_point(alpha = 3, size=3) +   
  geom_errorbar(aes(ymin=pca2-se2, ymax=pca2+se2), width=.01) +
  geom_line() +
  scale_colour_manual(values=Vcol) +
  TP_theme() +
  ylim(0,1) +
  labs(
    x = "Region",
    y = ""
  )

p3 <-  dd_clim_ave %>% arrange(V_warm_cum) %>%    
  mutate(SiteID=factor(SiteID, levels=SiteID)) %>%
  ggplot(aes(x = pca1, y=SiteID, color = SiteID)) +
  geom_point(alpha = 3, size=3) +   
  geom_errorbar(aes(xmin=pca1-se1, xmax=pca1+se1), width=.01) +
  geom_line() +
  scale_colour_manual(values=Vcol) +
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
          p2 + theme(legend.position="none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)), 
          p3 +  theme(legend.position="none"), 
          align = "vh", nrow = 2, 
          rel_heights = c(1/2, 1/2, 1/2),
          rel_widths = c(1/2, 1/2, 1/2),
          labels = c('A', 'B', 'C'))

#### For VPD average

p1 <-  dd_clim_ave %>% arrange(V_ann_dest) %>%    
  mutate(SiteID=factor(SiteID, levels=SiteID)) %>%
  ggplot(aes(x = pca1, y = pca2, color = SiteID)) +
  geom_point(alpha = 3, size=3) +    
  geom_errorbar(aes(ymin=pca2-se2, ymax=pca2+se2), width=.01) +
  geom_errorbar(aes(xmin=pca1-se1, xmax=pca1+se1), width=.01) +
  geom_line() +
  scale_colour_manual(values=Vcol) +
  TP_theme() +
  xlim(0,1) +
  ylim(0,1) +
  guides(color = guide_legend("Region")) +
  labs(
    x = "",
    y = "Distance to point"
  )

p2 <-  dd_clim_ave %>% arrange(V_ann_dest) %>%   
  mutate(SiteID=factor(SiteID, levels=SiteID)) %>%
  ggplot(aes(x = SiteID, y=pca2, color = SiteID)) +
  geom_point(alpha = 3, size=3) +   
  geom_errorbar(aes(ymin=pca2-se2, ymax=pca2+se2), width=.01) +
  geom_line() +
  scale_colour_manual(values=Vcol) +
  TP_theme() +
  ylim(0,1) +
  labs(
    x = "Region",
    y = ""
  )

p3 <-  dd_clim_ave %>% arrange(V_ann_dest) %>%    
  mutate(SiteID=factor(SiteID, levels=SiteID)) %>%
  ggplot(aes(x = pca1, y=SiteID, color = SiteID)) +
  geom_point(alpha = 3, size=3) +   
  geom_errorbar(aes(xmin=pca1-se1, xmax=pca1+se1), width=.01) +
  geom_line() +
  scale_colour_manual(values=Vcol) +
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
          p2 + theme(legend.position="none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)), 
          p3 +  theme(legend.position="none"), 
          align = "vh", nrow = 2, 
          rel_heights = c(1/2, 1/2, 1/2),
          rel_widths = c(1/2, 1/2, 1/2),
          labels = c('A', 'B', 'C'))

#### RUNNING BRMS MODEL WITH CLIMATE FOR POSITION OF X AXIS ####
library(brms)
library(modelr)
library(tidybayes)

m_ddclim = brm(
  pca1 ~ T_sum_cor*YearRange, 
  data = dd_clim_ave
)

dd_clim_ave %>%
  group_by(YearRange) %>%
  data_grid(T_sum_cor = seq_range(T_sum_cor, n = 51)) %>%
  add_fitted_draws(m_ddclim) %>%
  ggplot(aes(x = T_sum_cor, y = pca1, color = ordered(YearRange))) +
  stat_lineribbon(aes(y = .value)) +
  geom_point(data = dd_clim_ave) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2") +
  TP_theme() 

