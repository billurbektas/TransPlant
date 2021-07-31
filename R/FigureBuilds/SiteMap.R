#### MAKE A MAP ####
#10.04.2019
library(tidyverse)
library(readxl)
library(RColorBrewer)
library(maps)
library(mapproj)

metadata <- read_excel(path = '~/Dropbox/projects/SNF Experiment/Transplant_Incline/data/overview/site_overviews/MetadataOverview_Gradient_Site.xlsx', 
           sheet='Site level', skip = c(2))
metadat <- data.frame(gradient = metadata[,1], site = metadata[,2], lat =  metadata[,5], long = metadata[,6], elev=metadata[4], year1=metadata[7], yearn=metadata[8])
colnames(metadat) <- c('gradient', 'site', 'lat', 'long', 'elev', 'year1', 'yearn') 
metadat <- metadat %>% group_by(gradient) %>% dplyr::summarize(nSites=n(), lat=mean(as.numeric(lat)), long=mean(as.numeric(long)),
                                                        yearn=mean(yearn), year1=mean(year1), 
                                                        Elevrange=(max(elev)-min(elev)), Yearrange=yearn-year1)


world_map <- map_data("world")

#Creat a base plot with gpplot2
p <- ggplot() + coord_fixed() +
  xlab("") + ylab("")

#Add map to base plot
base_world_messy <- p + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
                                    colour="white", fill="#dddddd") +
                          scale_y_continuous(breaks = (-2:2) * 30) +
                          scale_x_continuous(breaks = (-4:4) * 45) +
                          coord_map("ortho",ylim=c(25,180), orientation=c(61, 0, 0))

base_world_messy

cleanup <- 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white', colour = 'white'), 
        axis.line = element_line(colour = "white"), legend.position="none",
        axis.ticks=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.line.y = ggplot2::element_blank())

base_world <- base_world_messy + cleanup


map_data <- base_world +
  geom_point(data=metadat, 
             aes(x=as.numeric(long), y=as.numeric(lat), fill=Yearrange, size=Elevrange), 
             pch=21, alpha=I(0.9)) + 
  theme(legend.position="bottom", 
        legend.text.align = 0) + # omit plot title saying 'color'
  scale_fill_distiller(palette ="Greys", direction = 1) +
  labs(size="Max Transplant Downslope (m)", fill="Year Since Established") +
  ggtitle('TransPlant Network Sites') 

map_data

# Zooming into map for Europe
library(sf)
#zoom_to <- c(13.38, 52.52)  # Berlin
zoom_to <- c(8.8, 46.3)  # Ticino
zoom_level <- 5
lon_span <- 360 / 2^zoom_level
lat_span <- 180 / 2^zoom_level

lon_bounds <- c(zoom_to[1] - lon_span / 2, zoom_to[1] + lon_span / 2)
lat_bounds <- c(zoom_to[2] - lat_span / 2, zoom_to[2] + lat_span / 2)

base_world +
  geom_point(data=metadat, 
             aes(x=as.numeric(long), y=as.numeric(lat), fill=Yearrange, size=Elevrange), 
             pch=21, alpha=I(0.9)) + 
  theme(legend.position="bottom", 
        legend.text.align = 0) + # omit plot title saying 'color'
  scale_fill_distiller(palette ="Greys", direction = 1)  +
  #geom_sf(data = st_sfc(st_point(zoom_to), crs = 4326),
  #        color = 'red') +
  coord_sf(xlim = lon_bounds, ylim = lat_bounds) 

# Zooming into map for Norway
library(sf)
#zoom_to <- c(13.38, 52.52)  # Berlin
zoom_to <- c(6.2, 60.8)  # Ticino
zoom_level <- 5
lon_span <- 360 / 2^zoom_level
lat_span <- 180 / 2^zoom_level

lon_bounds <- c(zoom_to[1] - lon_span / 2, zoom_to[1] + lon_span / 2)
lat_bounds <- c(zoom_to[2] - lat_span / 2, zoom_to[2] + lat_span / 2)

base_world +
  geom_point(data=metadat, 
             aes(x=as.numeric(long), y=as.numeric(lat), fill=Yearrange, size=Elevrange), 
             pch=21, alpha=I(0.9)) + 
  theme(legend.position="bottom", 
        legend.text.align = 0) + # omit plot title saying 'color'
  scale_fill_distiller(palette ="Greys", direction = 1)  +
  #geom_sf(data = st_sfc(st_point(zoom_to), crs = 4326),
  #        color = 'red') +
  coord_sf(xlim = lon_bounds, ylim = lat_bounds) 






