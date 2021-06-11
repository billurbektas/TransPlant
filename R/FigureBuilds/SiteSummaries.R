#### Creating site-specific information for Tables 1 and 2 of results ####

#Get total number of sites, num species, etc. (Table 1)
Tb1 <- dat %>% group_by(Region) %>%
  summarise(nSites = n_distinct(destSiteID),
            nPlots = n_distinct(destPlotID),
            nSpec = n_distinct(SpeciesName),
            nObs = n(),
            YearRange = max(Year)-min(Year),
            ElevRange = max(Elevation)-min(Elevation))

View(Tb1)

# For lat and lon information (supp Site table with lat/lon)
alldat <- tibble::lst(NO_Ulvhaugen, NO_Lavisdalen, NO_Gudmedalen, NO_Skjellingahaugen, 
                      CH_Lavey, CH_Calanda, 
                      US_Colorado, US_Montana, US_Arizona,
                      CN_Damxung, IN_Kashmir, CN_Gongga, CN_Heibei, 
                      DE_Grainau, DE_Susalps, DE_TransAlps, FR_AlpeHuez, SE_Abisko, FR_Lautaret, IT_MatschMazia1, IT_MatschMazia2)
meta <- alldat %>% 
  map_df("meta", .id='Region') %>%
  ungroup() 


meta %>% group_by(Region) %>% summarize(mElev = mean(dist(Elevation)), mPlot = mean(PlotSize_m2)) %>% View()

View(meta)
write.csv(meta, '~/Desktop/Subsite_summaries.csv')
