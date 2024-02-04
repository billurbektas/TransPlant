#### Read in climate data for all sites ####
library(GSODR)
library(remotes)
library(sf)

install.packages("remotes")
remotes::install_github("adamhsparks/GSODRdata")
nearest_stations(46.20299633094538, 7.062330594198718, 5)
#1 within 5 km for Lavey, 1 within 75 km for Heibei...

#not many stations, start with CRU data as more recent grids

# Create lat lon metadata
temp <- read.csv('./climate/worlclim2_processedclimate.csv')
#missing Calanda2 data, adding here
temp[58:59,] <- temp[17:18,]
temp[58:59, 2] <- "CH_Calanda2"
temp[58:59, "YearEstablished"] <- 2017 # Check this!
temp[58:59, "YearMin"] <- 2017
temp[58:59, "YearMax"] <- 2020
temp[58:59, "YearRange"] <- 3
temp[58:59, "PlotSize_m2"] <- 1 # Check this!
head(temp)

metadata = read_excel(path = 'data/metadata/site_overview.xlsx', 
                       sheet='Site level', skip = c(2))
metadat = data.frame(gradient = metadata[,1], site = metadata[,2], lat =  metadata[,5], long = metadata[,6], elev=metadata[4], year1=metadata[7], yearn=metadata[8])
colnames(metadat) = c('gradient', 'site', 'lat', 'long', 'elev', 'year1', 'yearn') 
metadat = metadat %>% group_by(gradient,site) %>% summarize(lat=mean(as.numeric(lat)), long=mean(as.numeric(long)),
                                                             yearn=mean(yearn), year1=mean(year1), 
                                                             elev=elev, Yearrange=yearn-year1)
coords = data.frame(x=metadat$long,y=metadat$lat)
points = vect(coords, proj4string = pre@crs)

# Read in CRU TS data

library(raster)
library(ncdf4)

pre = brick("data/climate/cru_ts4.07.1901.2022.pre.dat.nc", varname="pre") #Precipitation
tmp <- brick("./cru_ts4.04.1901.2019.tmp.dat.nc", varname="tmp") # Mean monthly temperature
pre.sites <- extract(pre, points)
tmp.sites <- extract(tmp, points)

# Read in elevation data
#elv <- brick("./climate/CRU_TS/halfdesg.elv", varname="elv") #not working
#adiabatic lapse rate 9.8 deg C / km (dry)

# Change column names
years <- 1901:2019
month <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
colnames(pre.sites) <- c(paste(rep(years, each=12), rep(month, times=119), sep="_"))
colnames(tmp.sites) <- c(paste(rep(years, eachhe=12), rep(month, times=119), sep="_"))

# Join and create long-form data
temp <- temp %>% dplyr::select(Gradient, destSiteID, Latitude, Longitude, PlotSize_m2, Elevation, YearRange, YearEstablished, YearMin, YearMax)
prec <- cbind.data.frame(temp, pre.sites)
temp <- cbind.data.frame(temp, tmp.sites)
pre.sites <- prec %>% pivot_longer(cols='1901_Jan':'2019_Dec') %>% separate(name, sep='_', into = c("year", "month")) %>% mutate(year = as.numeric(year)) %>% filter(year>2000)
temp.sites <- temp %>% pivot_longer(cols='1901_Jan':'2019_Dec')%>% separate(name, sep='_', into = c("year", "month")) %>% mutate(year = as.numeric(year)) %>% filter(year>2000)

# Get annual mean summer temperature
avetemp <- temp.sites %>% filter(month %in% c("Jun", "Jul", "Aug")) %>% group_by(Gradient, destSiteID, Latitude, Longitude, PlotSize_m2, Elevation, YearRange, YearEstablished, YearMin, YearMax, year) %>%
  dplyr::summarize(avetemp = mean(value, na.rm=T)) 

# Get summed summer precipitation
sumprec <-  pre.sites %>% filter(month %in% c("Jun", "Jul", "Aug")) %>% group_by(Gradient, destSiteID, Latitude, Longitude, PlotSize_m2, Elevation, YearRange, YearEstablished, YearMin, YearMax, year) %>%
  dplyr::summarize(sumprec = sum(value, na.rm=T)) 

# join data
annclim <- left_join(avetemp, sumprec)

# Note this is uncorrected data! At the summer year level

write.csv(annclim, './climate/cru_ts_processedclimate.csv')


