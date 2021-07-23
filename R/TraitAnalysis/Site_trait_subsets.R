#Site specific trait coverage

#US_ARIZONA
ARI <- alltraits %>% filter(Region == "US_Arizona") %>% ungroup() %>%
  select(Region, SpeciesName, submitted_name, TLeaf_Area_cm2:TStem_density) %>%
  distinct() %>% 
  filter(!grepl('NID', submitted_name), !is.na(SpeciesName)) 

ARI_missing <- ARI %>% filter(rowAny(
  across(
    .cols = TLeaf_Area_cm2:TSLA_cm2_g,
    .fns = ~ is.na(.x)
  ))
) 

splist_ARI <- dat %>% filter(Region == "US_Arizona") %>% select(Region, SpeciesName) %>% distinct() %>% View()


ARI <- ARI %>% select(submitted_name, TLeaf_Area_cm2:TStem_density) 

rowAny <- function(x) rowSums(x) > 0 
 
ARI %>% 
  filter(rowAny(
    across(
      .cols = TLeaf_Area_cm2:TSLA_cm2_g,
      .fns = ~ is.na(.x)
    ))
  ) %>% select(-destSiteID) %>%
  distinct() %>% 
  filter(!grepl('NID', submitted_name), !is.na(SpeciesName)) 
  

#DE_SUSALPS AND TRANSALPS


# 
#   mutate(Plant_Veg_Height_cm = ifelse(!is.na(Plant_Veg_Height_cm), Plant_Veg_Height_cm, TPlant_Veg_Height_cm),
#          Leaf_Area_cm2 = ifelse(!is.na(Leaf_Area_cm2), Leaf_Area_cm2, TLeaf_Area_cm2),
#          SLA_cm2_g = ifelse(!is.na(SLA_cm2_g), SLA_cm2_g, TSLA_cm2_g),
#          C_percent = ifelse(!is.na(C_percent), C_percent, TC_percent),
#          N_percent = ifelse(!is.na(N_percent), N_percent, TN_percent),
#          P_percent = TP_percent,
#          Seed_Mass = TSeed_Mass,
#          Stem_density = TStem_density)

SUS <- alltraits %>% filter(Region %in% c("DE_Susalps", "DE_TransAlps")) %>% ungroup() %>%
  select(SpeciesName, submitted_name, TLeaf_Area_cm2:TStem_density) %>%
  distinct() %>% 
  filter(!grepl('NID', submitted_name), !is.na(SpeciesName)) 

SUS_missing <- SUS %>% filter(rowAny(
  across(
    .cols = TLeaf_Area_cm2:TSLA_cm2_g,
    .fns = ~ is.na(.x)
  ))
) 


rowAny <- function(x) rowSums(x) > 0 

TRYSUS <- SUS %>% 
  filter(rowAny(
    across(
      .cols = TLeaf_Area_cm2:TSLA_cm2_g,
      .fns = ~ is.na(.x)
    ))
  ) %>% distinct()  

splist_SUS <- dat %>% filter(Region %in% c("DE_Susalps", "DE_TransAlps")) %>% select(Region, SpeciesName) %>% distinct()

allSUS <- left_join(splist_SUS, TRYSUS, by="SpeciesName")

write.csv(TRYSUS, '~/Desktop/try_sus_trans_alp.csv')

# Add in site level data
try_fil <- alltraits %>% ungroup() %>%
  select(Region, -destSiteID, SpeciesName, Plant_Veg_Height_cm:N_percent) %>%
  group_by(SpeciesName) %>%
  summarise(across(Plant_Veg_Height_cm:N_percent, ~mean(.x, na.rm = TRUE)))


try_fil <- left_join(SUS, try_fil)

try_fil_missing <- try_fil %>% 
  filter(rowAny(
    across(
      .cols = TLeaf_Area_cm2:TSLA_cm2_g,
      .fns = ~ is.na(.x)
    ))
  ) %>% distinct()  

write.csv(try_fil, '~/Desktop/try_site_data.csv')
