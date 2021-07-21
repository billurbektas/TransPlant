#Site specific trait coverage

#US_ARIZONA
ARI <- alltraits %>% filter(Region == "US_Arizona")

splist_ARI <- dat %>% filter(Region == "US_Arizona") %>% select(Region, SpeciesName) %>% distinct() %>% View()


ARI <- ARI %>% select(submitted_name, TLeaf_Area_cm2:TStem_density) 

rowAny <- function(x) rowSums(x) > 0 
 
ARI %>% 
  filter(rowAny(
    across(
      .cols = TLeaf_Area_cm2:TSLA_cm2_g,
      .fns = ~ is.na(.x)
    ))
  ) %>% View()

#DE_SUSALPS AND TRANSALPS
SUS <- alltraits %>% filter(Region %in% c("DE_Susalps", "DE_TransAlps"))

SUS <- SUS %>% select(submitted_name, TLeaf_Area_cm2:TStem_density) 

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
