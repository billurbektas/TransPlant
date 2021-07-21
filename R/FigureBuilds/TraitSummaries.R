# Table of trait coverage for sites and TRY

#using alltraits from trait_drakeplan.R
ARI <- alltraits %>% filter(Region == "US_Arizona")

try_sum <- ARI %>% 
  select(Region, destSiteID.y, SpeciesName, TLeaf_Area_cm2:TStem_density) %>%
  #pivot_longer(cols=c(TLeaf_Area_cm2:TStem_density), names_to="TRYTraits") %>%
  dplyr::group_by(Region, destSiteID.y) %>% 
  dplyr::summarize(nNum = n(),
            nSp=n_distinct(SpeciesName),
            nLA = length(!is.na(TLeaf_Area_cm2))/n(),
            nC = length(!is.na(TC_percent))/n(), 
            nN = length(!is.na(TN_percent))/n(), 
            nP = length(!is.na(TP_percent))/n(),
            nH = length(!is.na(TPlant_Veg_Height_cm))/n(), 
            nSM = length(!is.na(TSeed_Mass))/n(), 
            nSLA = length(!is.na(TSLA_cm2_g))/n(), 
            nSD = length(!is.na(TStem_density))/n())
            
try_sum <- alltraits %>% 
  select(Region, destSiteID, SpeciesName, TLeaf_Area_cm2:TStem_density) %>%
  pivot_longer(cols=c(TLeaf_Area_cm2:TStem_density), names_to="TRYTraits") %>%
  dplyr::group_by(Region, destSiteID, TRYTraits) %>% 
  dplyr::summarize(n = n(), 
                   nSp = n_distinct(SpeciesName))
percov = sum(!is.na(value))/(n/8))

site_traits_means %>%
  group_by(Region, destSiteID) %>%
  dplyr::summarize(n=n(), nSp=n_distinct(SpeciesName)) %>% View()
  