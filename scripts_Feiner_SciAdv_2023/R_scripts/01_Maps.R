#check here for example on how to plot average across landscape:  
#https://gadm.org/download_country.html

rm(list=ls())

#load required packages
library(ggplot2)
library(rgdal)
library(tmap)
library(dplyr)
library(sf)
library(raster)
library(rnaturalearth)

############### Fig. 1A - P. muralis distribution map

world <- ne_countries(scale = "large", returnclass = "sf") #change to large for final plot!

Pmur_dist <- st_read("Input_files/data_0.shp")
pdf("Output_plots/DistributionMap.pdf", useDingbats = F)
ggplot(data = world) + 
  geom_sf(data = world, fill= "white", color = "grey", size = 0.1) +
  geom_sf(data = Pmur_dist, fill = "darkgrey", color = NA) +
  coord_sf(xlim = c(-8.00, 34.00), ylim = c(35.00, 52.00), expand = T) +
  annotate("rect", xmin = 7.57, xmax = 13.73, ymin = 41.18, ymax = 44.93,color = "black", fill=NA) +
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.background = element_rect(fill = "#EBF9FA"), panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
dev.off()

Locs <- read.csv("Input_files/PopsLocs.csv")

ggplot(data = world) + 
  geom_sf(data = world, fill= "white", color = "grey", size = 0.1) +
  geom_point(data = Locs, aes(x = Longitude, y = Latitude, fill=lineage), size = 3, shape = 21) + 
  scale_fill_manual(values=c("#E59F24","#D26027","#5BB3E4")) +
  coord_sf(xlim = c(7.57, 13.73), ylim = c(41.18, 44.93), expand = T) +
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.background = element_rect(fill = "aliceblue"))

###############  Fig. 1E - Coastal cline pops

Locs_SAintro <- subset(Locs, abbpop == "LO" | abbpop == "NL" | abbpop == "VA" | 
                       abbpop == "ME" | abbpop == "GN" | abbpop == "RA" | abbpop == "SL" | abbpop == "LE" | 
                       abbpop == "ST" | abbpop == "MG" | abbpop == "VI" | abbpop == "CA" | abbpop == "CR" |
                       abbpop == "CN" | abbpop == "VE")
world_cropped <- st_crop(world, xmin = 8.206791, xmax = 11.16181,
                          ymin = 43.37148, ymax = 44.49316)

pdf("Output_plots/Intro_cline.pdf", useDingbats = F)
ggplot(data = world_cropped) + 
  geom_sf(data = world_cropped, fill= "white", color = "grey", size = 0.1) +
  coord_sf(xlim = c(8.206791, 11.16181), ylim = c(43.37148, 44.49316), expand = F, clip="on") +
  geom_line(data = Locs_SAintro, aes(x = Longitude, y = Latitude), size=0.2) +
  geom_point(data = Locs_SAintro, aes(x = Longitude, y = Latitude, fill=lineage), size = 3, shape = 21) + 
  scale_fill_manual(values=c("#E59F24","#D26027","#5BB3E4")) +
  xlab("Longitude") + ylab("Latitude") +
  theme(panel.background = element_rect(fill = "aliceblue"))
dev.off()

############### Fig. 1B
Data <- read.csv("Input_files/Phenotype_data_final_June2023.csv")
Data <- Data[!is.na(Data$GreenResolved),]
Data <- subset(Data, sex == "M")
observ <- Data %>% group_by(abbpop) %>% count()
Data_mean <- Data %>% group_by(abbpop) %>% summarize(Green = mean(GreenResolved))

Data_all <- merge(Data_mean,Locs, by = "abbpop")

# Load shapefile:
shapename <- rgdal::readOGR(dsn= "Input_files/gadm36_ITA_shp",layer ="gadm36_ITA_0")

#This code below is the treatment of the raster and our data to have everything in the correct coordinates
g <- spTransform(shapename, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
Borders <- ggplot() + geom_polygon(data=g,aes(x=long,y=lat,group=group),fill='white',color = "black")

coordinates(Data_all) <- ~Longitude+Latitude
proj4string(Data_all) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84" )
sp_mydata <- spTransform(Data_all,CRS(proj4string(g)))
cord.UTM <- spTransform(g, CRS("+init=epsg:7791"))
sp_mydata <- spTransform(sp_mydata,CRS("+init=epsg:7791"))

#extent(Data_all) #this gives the extreme coordinates of all populations, used in next step (plus margin)
grd <- expand.grid(x = seq(from = 7.57, to = 13.73, length.out = 10000),
                    y = seq(from = 41.18, to = 44.93,  length.out = 10000)) # change this for modifying the frame around data points
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object
proj4string(Data_all) <- proj4string(grd)
str(grd)

Krig = gstat::idw(Green~1,Data_all,newdata=grd, idp=4.0)

# Convert to raster object the interpolation and clip it to the Italian raster
r       <- raster(Krig)
r.m     <- mask(r,g)
#cbp1 <- c("#4d3e1d","#675326", "#816830","#336600","#339900","#66CC33")
#cbp1 <- c("#4d3e1d","#675326", "#816830","#66CC33","#339900","#336600")
cbp1 <- c("#816830","#675326","#4d3e1d","#336600","#339900","#66CC33")# Manual set of the palette

# Plot
green<- tm_shape(r.m, raster.downsample = FALSE) + #plot interpolation in raster
  tm_raster(n=10,palette = cbp1, interpolate=TRUE, #plot using the palet and interpolating across map
            legend.show = TRUE,title = "") + # add legend 
  tm_shape(Data_all) + #add the sampling locations
  tm_dots("abbpop", col ='black', size=0.3, shape=21) +
 # tm_dots("abbpop", fill = "white", size=0.3, shape=21) +
    tm_legend(position = c("left", "bottom"), frame = FALSE)+ #position of the legend
tm_layout(bg.color = "#EBF9FA", inner.margins = c(0, 0, 0, 0))+  #lay out of the graph
tm_legend(outside = F,
          legend.title.size = 1.5,
          legend.text.size = 1,
          legend.bg.alpha = 0.6)+
    tm_layout(frame = T) #With no frame, delete if color not working!

pdf("Output_plots/Map_Greenness.pdf", height = 4, useDingbats=FALSE)
green
dev.off()
