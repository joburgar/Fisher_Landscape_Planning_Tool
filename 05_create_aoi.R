# Copyright 2021 Province of British Columbia
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and limitations under the License.
#
#####################################################################################
# 05_create_aoi.R
# script to create canned example aoi for Boreal and Columbian fisher populations
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 01-Mar-2022
#####################################################################################
version$major
version$minor
R_version <- paste0("R-",version$major,".",version$minor)

.libPaths(paste0("C:/Program Files/R/",R_version,"/library")) # to ensure reading/writing libraries from C drive
tz = Sys.timezone() # specify timezone in BC

# Load Packages
list.of.packages <- c("tidyverse", "lubridate","chron","bcdata", "bcmaps","sf", "rgdal",
                      "Cairo","OpenStreetMap", "ggmap","PNWColors","units","nngeo","raster")

# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

#####################################################################################
# keep in mind that WGS84 lat/long espg = 4326; BC Albers espg = 3005; NAD83 / UTM zone 10N espg = 26910

###--- function to retrieve geodata from BCGW
droplevels.sfc = function(x, except, exclude, ...) x

retrieve_geodata_aoi <- function (ID=ID){
  aoi.geodata <- bcdc_query_geodata(ID) %>%
    filter(BBOX(st_bbox(aoi))) %>%
    collect()
  aoi.geodata <- aoi.geodata %>% st_intersection(aoi)
  aoi.geodata$Area_km2 <- st_area(aoi.geodata)*1e-6
  aoi.geodata <- drop_units(aoi.geodata)
  aoi.geodata <- droplevels.sfc(aoi.geodata, except = geometry)
  return(aoi.geodata)
}

#####################################################################################
###--- function to retrieve data from downloaded gdb or shapefile
# Read the feature class

retrieve_gdb_shp_aoi <- function (dsn=dsn, layer=layer){
  # if a gdb then the dsn should be the path all the way to the '.gdb'
  aoi.geodata <- st_read(dsn=dsn,layer=layer) %>%
    st_transform(crs=3005) %>% st_intersection(aoi)
  aoi.geodata$Area_km2 <- st_area(aoi.geodata)*1e-6
  aoi.geodata <- drop_units(aoi.geodata)
  aoi.geodata <- droplevels.sfc(aoi.geodata, except = geometry)
  return(aoi.geodata)
}


#####################################################################################
###--- function to retrieve FETA raster data for aoi

retrieve_raster_aoi <- function(rfile=rfile, raoi=raoi){
  # rfile=FETA.rasters[3]
  rfile <- raster(paste0("data/FETA_fromKyle/",rfile))
  rfile <- projectRaster(rfile,crs = 26910)
  # Crop FETA rasters by extent of raoi
  rfile.aoi <- crop(rfile, extent(raoi))
  # Sum FETA rasters for value of FETA within fisher territory
  # recall FETA rasters are 100x100 m while aoi is 5500x5500
  # need to sum by 55 cells in x and y direction to make fisher equivalent territory
  ra <- aggregate(rfile.aoi, fact=c(55,55), fun=sum, na.rm=TRUE)
  return(ra)
}

################################################################################
#- function to create grid cell from aoi (shapefile)
create_grid <- function (input=input, cellsize=cellsize){
  # input=traps.sf
  # cellsize=5000
  aoi_utm <- st_transform(input, crs=26910) # to have in metres for specifying grid cell size
  aoi_utm$TrapNum <- rownames(aoi_utm)
  aoi_grid <- st_make_grid(st_bbox(aoi_utm), what="polygons", cellsize=cellsize, square=TRUE) #  grid for entire AOI (rectangle)

  # subset grid to just cells of interest
  aoi_grid <- aoi_grid[aoi_utm]
  # To sf and add grid ID
  fishnet_grid_sf = st_sf(aoi_grid) %>%
    # add grid ID
    mutate(grid_id = 1:length(lengths(aoi_grid)))

  fishnet_grid_sf$Area_km2 <- st_area(fishnet_grid_sf)*1e-6
  fishnet_grid_sf <- drop_units(fishnet_grid_sf)
  fishnet_grid_sf %>% summarise(sum(Area_km2)) %>% st_drop_geometry() #714.5 km2 study area each cell is 21.65 km2

  tg.dist <- st_nn(aoi_utm, fishnet_grid_sf, k=1, returnDist = T)
  aoi_utm$Trap_Grp <- unlist(tg.dist$nn)

  return(list(aoi_utm=aoi_utm, fishnet_grid_sf=fishnet_grid_sf))
}

#####################################################################################
###--- CREATE AGGREGATED RASTER STACK OF FETA VALUES
#####################################################################################
# Created, and no need to re-run, just pull in appropriate raster files

# First create fishnet for entire area (aoi = buffer of 30 km2 around cutblocks and bounding box)
# Then bring in underlying habitat for aoi
# make raster stack for Fisher population (combined) extent
# aoi <- st_make_grid(aoi.Fpop, n=1)
# aoi <- st_as_sf(aoi)
# aoi <- st_transform(aoi, crs=26910)
#
# raoi <- raster(ext=extent(aoi), crs=26910, res=c(5500,5500))
# raoi <- rasterize(aoi, raoi)
#
# FETA.rasters <- list.files("data/FETA_fromKyle")
# # for now, just want to bring in fisher habitat requirement data
# FETA.rasters <- FETA.rasters[grepl("movement|rest|denning", FETA.rasters)]
#
# Fraster_list=list()
# for(i in 1:length(FETA.rasters)){
#   Fraster_list[[i]] <- retrieve_raster_aoi(rfile=FETA.rasters[i], raoi=raoi)
# }
# Fraster_stack = stack(Fraster_list)
# writeRaster(Fraster_stack, file=paste0("out/Fraster_stack_CBpops_aoi.tif"), bylayer=TRUE, overwrite=TRUE)
# plot(Fraster_stack)

# saved as individual rasters: "Fraster_stack_CBpops_aoi_#.tif" with number 1:5 for 5 habitat types

#####################################################################################
###--- CREATE CANNED SCENARIO AOI SHAPEFILE
#####################################################################################
# Bring in population polygons
aoi.Fpop <- st_read(dsn=paste0(getwd(),"/data"), layer="aoi.Fpop")
aoi.Fpop$Fpop <- case_when(aoi.Fpop$Fpop=="Central Interior" ~ "Columbian",
                           TRUE ~ as.character(aoi.Fpop$Fpop))

###--- Create an example polygon, with potential cutblocks
# Consolidated cutblock
# 1: Harvested Areas of BC (Consolidated Cutblocks) (multiple, fgdb, wms, kml, pdf)
# ID: b1b647a6-f271-42e0-9cd0-89ec24bce9f7
# Name: harvested-areas-of-bc-consolidated-cutblocks-
# aoi <- st_read(dsn=paste0(getwd(),"/data/ThompOkan_FGT/1_CUTBLOCKS"), layer="cutblocks_sample") %>% st_transform(crs=26910)
#
# aoi <- aoi %>% st_transform(crs=3005)
# aoi.CUT <- retrieve_geodata_aoi(ID = "2ebb35d8-c82f-4a17-9c96-612ac3532d55")
# aoi.CUT <- aoi.CUT %>% dplyr::select(id, FEATURE_ID, HARVEST_DATE, Area_km2)
# aoi.CUT$HARVEST_Yr <- year(aoi.CUT$HARVEST_DATE)
# aoi.CUT <- aoi.CUT %>% filter(!is.na(HARVEST_Yr))
# aoi.CUT %>% dplyr::count(HARVEST_Yr) %>% st_drop_geometry()
#
# aoi.CUT$HARVEST_DECADE <- as.numeric(substr(aoi.CUT$HARVEST_Yr,3,3))
# unique(aoi.CUT$HARVEST_DECADE)
# aoi.CUT$HARVEST_DECADE <- case_when(aoi.CUT$HARVEST_DECADE==1 ~2010,
#                                     aoi.CUT$HARVEST_DECADE==0 ~2000,
#                                     aoi.CUT$HARVEST_DECADE==9 ~1990,
#                                     aoi.CUT$HARVEST_DECADE==8 ~1980,
#                                     aoi.CUT$HARVEST_DECADE==7 ~1970,
#                                     aoi.CUT$HARVEST_DECADE==6 ~1960,
#                                     aoi.CUT$HARVEST_DECADE==5 ~1950)
# st_write(aoi.CUT, dsn=paste0(getwd(),"/data/aoi_CUT_example_TO.shp"), delete_layer=TRUE)

aoi.CUT_Bex <- st_read(dsn=paste0(getwd(),"/data"), layer="aoi_CUT_example") %>%
  rename(HARVEST_DATE = HARVEST_DA, HARVEST_YEAR = HARVEST_Y, HARVEST_DECADE = HARVEST_DE, Area_km2 = Are_km2)
aoi.CUT <- aoi.CUT_Bex

aoi.CUT <- st_read(dsn=paste0(getwd(),"/data"), layer="aoi_CUT_example_Skeena") %>%
  rename(HARVEST_DATE = HARVEST_DA, HARVEST_YEAR = HARVEST_Y, HARVEST_DECADE = HARVEST_DE, Area_km2 = Are_km2)

aoi.CUT %>% group_by(HARVEST_DECADE) %>% summarise(sum(Area_km2)) %>% st_drop_geometry()
aoi.CUT %>% dplyr::count(HARVEST_DECADE) %>%  st_drop_geometry()

ggplot()+
  geom_sf(data=aoi.CUT, aes(fill=HARVEST_DECADE))

# Run through canned examples with different configurations of cuts from 1970
# actual cutblocks, but long enough ago to potentially be classified as 'suitable' by our hypothetical VRI test
# first generate random numbers for each cutblock and select different sets per example
# 420 cutblocks covering 55.6 km2 in 1980, try 4 examples
# 1. all cutblocks harvested (420)
# 2. half cutblocks harvested (210) random selection
# 3. half cutblocks harvested (210), in 'lower quality' (i.e., unsuitable) habitat
# 4. half cutblocks harvested (210), in 'higher quality' (i.e., suitable) habitat

aoi <- st_make_grid(st_buffer(aoi.CUT %>% st_transform(crs=26910), dist=30000), n=1)
aoi <- st_as_sf(aoi)

# Baoi <- st_make_grid(st_buffer(aoi.CUT_Bex %>% st_transform(crs=26910), dist=30000), n=1)
# Baoi <- st_as_sf(Baoi)

pal = pnw_palette(name="Winter",n=2,type="discrete")

Cairo(file="out/BCex_Fpop_aoi_plot.PNG",type="png",width=2400,height=2400,pointsize=16,bg="white",dpi=300)
ggplot()+
  geom_sf(data=aoi.Fpop, aes(fill=Fpop))+
  geom_sf(data=Baoi)+
  geom_sf(data=aoi.CUT_Bex)+
  geom_sf(data=aoi)+
  geom_sf(data=aoi.CUT)+
  scale_fill_manual(values=rev(pal))+
  theme(legend.position="bottom")+
  theme(legend.title=element_blank())
dev.off()

aoi <- st_join(aoi, aoi.Fpop %>% st_transform(crs=26910),left=TRUE, largest=TRUE)

ggplot()+
  geom_sf(data=aoi.CUT %>% filter(HARVEST_DECADE %in% c(1980)), aes(fill=HARVEST_YEAR))+
  geom_sf(data=aoi, fill=NA)


### CREATNG SUITABLE HABITAT WITH FIRST ROUND OF MAHALANOBIS DISTANCE OUTPUTS
aoi.CUTex <- aoi.CUT %>% filter(HARVEST_DECADE==1980)
aoi.CUTex$rndmrnk <- rank(round(runif(nrow(aoi.CUTex), min=10000, max=99999)))

# maldist <- st_read(dsn=paste0(getwd(),"/data"), layer="Mahalanobis_predictions_2021_220304")
# summary(maldist %>% st_drop_geometry())

aoi <- aoi %>% st_transform(crs=3005)
aoi.MAL <- retrieve_gdb_shp_aoi(dsn=paste0(getwd(),"/data"), layer="Mahalanobis_predictions_2021_220304")
# aoi.MAL %>% count(D2) %>% st_drop_geometry()
summary(aoi.MAL$D2); nrow(aoi.MAL) # 439 FETA in the Boreal aoi, 115 in Columbian

aoi.MAL <- aoi.MAL %>% mutate(D2_grp = cut(D2, c(0, 10, 20, 30, Inf)))
levels(aoi.MAL$D2_grp)
aoi.MAL$D2_grp <- recode(aoi.MAL$D2_grp, "(0,10]"="Suitable Territory",
                  "(10,20]"="Suitable Movement",
                  "(20,30]"="Possible Movement",
                  "(30,Inf]"="Unsuitable")
glimpse(aoi.MAL)
aoi.MAL %>% group_by(D2_grp) %>% summarize(n=n()) %>% st_drop_geometry()
# Boreal
# 1 Suitable Territory   140
# 2 Suitable Movement    124
# 3 Possible Movement     68
# 4 Unsuitable           107

# Columbian
# 1 Suitable Territory     7
# 2 Suitable Movement      5
# 3 Possible Movement      4
# 4 Unsuitable            99

# indicate which FETA are targeted for harvesting
aoi.MAL <- st_join(aoi.MAL, aoi.CUTex %>% dplyr::select(rndmrnk),left=TRUE, largest=TRUE)
aoi.MAL %>% count(rndmrnk)
aoi.MAL$Harvest <- ifelse(is.na(aoi.MAL$rndmrnk), 0,1)
aoi.MAL %>% group_by(D2_grp) %>% count(Harvest) %>% st_drop_geometry()
aoi.MAL %>% count(Harvest) %>% st_drop_geometry() # only 54 cells targeted for harvesting
# 54/439 12% set for harvesting in 1980
# 23 or 5% in suitable territory (i.e., 16% of suitable territories to be removed)

# if want to make more dire, can change to all being "suitable" pre model run, but first try with real data
# aoi.MAL <- aoi.MAL %>% mutate(D2_grp2 = case_when(Harvest==1 ~ "Suitable Territory",
#                                                   TRUE ~ as.character(D2_grp)))
# aoi.MAL %>% group_by(D2_grp2) %>% count(Harvest) %>% st_drop_geometry() # now have all harvesting blocks in 'suitable habitat'

pal = pnw_palette(name="Cascades",n=4,type="discrete")

aoi.MAL.plot <- ggplot()+
  geom_sf(data=aoi.MAL, aes(fill=D2_grp))+
  geom_sf(data=aoi, fill=NA)+
  scale_fill_manual(values=rev(pal))+
  theme(legend.position="bottom")+
  theme(legend.title=element_blank())

Cairo(file="out/Cex_aoi_MAL_plot.PNG",type="png",width=2400,height=2400,pointsize=14,bg="white",dpi=300)
aoi.MAL.plot
dev.off()

################################################################################
###--- Set up 4 possible scenarios
# (1) 0% harvesting - run simulation with underlying landbase for baseline number
# (2) 25% harvesting of cutblocks in suitable habitat
# (3) 50% harvesting of cutblocks in suitable habitat
# (4) 75% harvesting of cutblocks in suitable habitat

# determine number of cutblocks in suitable territory (FETA)
###NEED TO REWORK THIS - NOT WORKING NOW
num.cutblocks.FETA <- aoi.MAL %>% filter(D2_grp=="Suitable Territory") %>% count(D2_grp) %>% st_drop_geometry()
num.cutblocks.FETA <- num.cutblocks.FETA$n

# Scenario 1 -0% harvesting - run simulation with underlying landbase for baseline number - Columbian - only 1 scenario because cutblocks not in suitable habitat
aoi.CUTCex <- st_join(aoi.CUTex %>% filter(HARVEST_DECADE==1980), aoi.MAL %>% dplyr::select(D2_grp),left=TRUE, largest=TRUE)
aoi.CUTCex <- aoi.CUTCex %>% arrange(D2_grp,rndmrnk)
aoi.CUTCex %>% count(D2_grp)

aoi.MAL.canCEx1.plot <- ggplot()+
  geom_sf(data=aoi.MAL, aes(fill=D2_grp))+
  geom_sf(data=aoi.CUTCex,fill="black",col="black",lwd=1)+
  geom_sf(data=aoi, fill=NA)+
  scale_fill_manual(values=rev(pal))+
  theme(legend.position="bottom")+
  theme(legend.title=element_blank())+
  ggtitle("Scenario 1\n 0% harvesting in Fisher Equivalent Territory Areas") +
  theme(plot.title = element_text(color="grey2", size=12, face="bold"))

Cairo(file="out/canCex_aoi_MAL_plot.PNG",type="png",width=2400,height=2400,pointsize=14,bg="white",dpi=300)
aoi.MAL.canCEx1.plot
dev.off()

aoi.MAL <- st_join(aoi.MAL, aoi.CUTCex %>% dplyr::select(HARVEST_DECADE),left=TRUE, largest=TRUE)
aoi.MAL$HARVEST_DECADE[is.na(aoi.MAL$HARVEST_DECADE)] <- 0
habitat <- ifelse(aoi.MAL$D2_grp=="Suitable Territory" & aoi.MAL$HARVEST_DECADE==0,1,0)
aoi.MAL$HARVEST_DECADE <- NULL

aoi.MAL$Habitat_canCex1 <- habitat

aoi.MAL %>% group_by(D2_grp) %>%
  summarise_at("Habitat_canCex1", sum, na.rm = TRUE) %>%
  st_drop_geometry()

glimpse(aoi.MAL)
# because multiple cutblocks within a FETA ended up with 130 FETAs in CanBex1, 121, 118, and 117 in CanBex2, CanBex3, and CanBex4, respectively

aoi <- aoi %>% st_transform(crs=26910)

canCex_raster <- list()
raoi <- raster(ext=extent(aoi), crs=26910, res=c(5500,5500))
raoi <- rasterize(aoi.MAL %>% st_transform(crs=26910), raoi, field="Habitat_canCex1", fun="max")
raoi[is.na(raoi[])] <- 0

plot(raoi)
#plot of the raster showing habitat as binary (green = suitable, white = unsuitable)
# sum(raoi@data@values) # 129
# round(sum(raoi@data@values) / length(raoi@data@values)*100) # 27% habitat
Cairo(file="out/IBM_raster_Habitat_canCex1.PNG", type="png", width=2200, height=2000,pointsize=15,bg="white",dpi=300)
plot(raoi,main=c(paste0("Fisher Equivalent Territory Areas (FETA) in the Area of Interest,"),
                 paste0(sum(raoi@data@values)," or ",round(sum(raoi@data@values) / length(raoi@data@values)*100),"% FETA in Scenario canCex1")),
                 cex.main=0.8)
dev.off()
canCex1_raster <- raoi

sum(canCex1_raster@data@values) # 7 suitable habitat cells
dim(canCex1_raster) # 156 cells so 7 chances for fisher

IBM_aoi <- list(aoi=aoi.MAL, canCex1_raster=canCex1_raster)
save(IBM_aoi, file="data/IBM_aoi_canCex.RData")


#################################################################################
# Boreal
# Scenario 1 -0% harvesting - run simulation with underlying landbase for baseline number
aoi.CUTBex <- st_join(aoi.CUTex %>% filter(HARVEST_DECADE==1980), aoi.MAL %>% dplyr::select(D2_grp),left=TRUE, largest=TRUE)
aoi.CUTBex <- aoi.CUTBex %>% arrange(D2_grp,rndmrnk)
aoi.CUTBex %>% count(D2_grp)

# Scenario 1 = 0% harvesting aka 0 cutblocks in suitable territory
canBEx1 <- aoi.CUTBex
canBEx1 <- canBEx1 %>% filter(D2_grp!="Suitable Territory")
canBEx1 %>% count(D2_grp)

# Scenario 2 = 25% harvesting in suitable habitat aka 59 cutblocks in suitable territory
canBEx2 <- aoi.CUTBex
canBEx2 <- canBEx2[round(num.cutblocks.FETA*.75):nrow(aoi.CUTBex),]
canBEx2 %>% dplyr::count(D2_grp)

# Scenario 3 = 50% harvesting in suitable habitat aka 94 cutblocks in suitable territory
canBEx3 <- aoi.CUTBex
canBEx3 <- canBEx3[round(num.cutblocks.FETA*.5):nrow(aoi.CUTBex),]
canBEx3 %>% dplyr::count(D2_grp)

# Scenario 4 = 75% harvesting in suitable habitat aka 129 cutblocks in suitable territory
canBEx4 <- aoi.CUTBex
canBEx4 <- canBEx4[round(num.cutblocks.FETA*.25):nrow(aoi.CUTBex),]
canBEx4 %>% dplyr::count(D2_grp)

aoi.MAL.canBEx1.plot <- ggplot()+
  geom_sf(data=aoi.MAL, aes(fill=D2_grp))+
  geom_sf(data=canBEx1,fill="black",col="black",lwd=1)+
  geom_sf(data=aoi, fill=NA)+
  scale_fill_manual(values=rev(pal))+
  theme(legend.position="bottom")+
  theme(legend.title=element_blank())+
  ggtitle("Scenario 1\n 0% harvesting in Fisher Equivalent Territory Areas") +
  theme(plot.title = element_text(color="grey2", size=12, face="bold"))

Cairo(file="out/canBex1_aoi_MAL_plot.PNG",type="png",width=2400,height=2400,pointsize=14,bg="white",dpi=300)
aoi.MAL.canBEx1.plot
dev.off()

aoi.MAL.canBEx2.plot <- ggplot()+
  geom_sf(data=aoi.MAL, aes(fill=D2_grp))+
  geom_sf(data=canBEx2,fill="black",col="black",lwd=1)+
  geom_sf(data=aoi, fill=NA)+
  scale_fill_manual(values=rev(pal))+
  theme(legend.position="bottom")+
  theme(legend.title=element_blank())+
  ggtitle("Scenario 1\n 25% harvesting in Fisher Equivalent Territory Areas") +
  theme(plot.title = element_text(color="grey2", size=12, face="bold"))

Cairo(file="out/canBex2_aoi_MAL_plot.PNG",type="png",width=2400,height=2400,pointsize=14,bg="white",dpi=300)
aoi.MAL.canBEx2.plot
dev.off()

aoi.MAL.canBEx3.plot <- ggplot()+
  geom_sf(data=aoi.MAL, aes(fill=D2_grp))+
  geom_sf(data=canBEx3,fill="black",col="black",lwd=1)+
  geom_sf(data=aoi, fill=NA)+
  scale_fill_manual(values=rev(pal))+
  theme(legend.position="bottom")+
  theme(legend.title=element_blank())+
  ggtitle("Scenario 1\n 50% harvesting in Fisher Equivalent Territory Areas") +
  theme(plot.title = element_text(color="grey2", size=12, face="bold"))

Cairo(file="out/canBex3_aoi_MAL_plot.PNG",type="png",width=2400,height=2400,pointsize=14,bg="white",dpi=300)
aoi.MAL.canBEx3.plot
dev.off()

aoi.MAL.canBEx4.plot <- ggplot()+
  geom_sf(data=aoi.MAL, aes(fill=D2_grp))+
  geom_sf(data=canBEx4,fill="black",col="black",lwd=1)+
  geom_sf(data=aoi, fill=NA)+
  scale_fill_manual(values=rev(pal))+
  theme(legend.position="bottom")+
  theme(legend.title=element_blank())+
  ggtitle("Scenario 1\n 75% harvesting in Fisher Equivalent Territory Areas") +
  theme(plot.title = element_text(color="grey2", size=12, face="bold"))

Cairo(file="out/canBex4_aoi_MAL_plot.PNG",type="png",width=2400,height=2400,pointsize=14,bg="white",dpi=300)
aoi.MAL.canBEx4.plot
dev.off()

################################################################################
# (1) transform aoi into fisher grid
# (2) join cells to identify as habitat based on D2_grp category
# (3) still binary "suitable" vs all other quality types for model

################################################################################
# for canned scenarios - i.e., with any harvesting

# (1) 0% harvesting - run simulation with underlying landbase for baseline number
# (2) 25% harvesting of cutblocks in suitable habitat
# (3) 50% harvesting of cutblocks in suitable habitat
# (4) 75% harvesting of cutblocks in suitable habitat

# bring in harvesting data
habitat <- list()
canexamples <- list(canBEx1, canBEx2, canBEx3, canBEx4)
for(i in 1:length(canexamples)){
aoi.MAL <- st_join(aoi.MAL, canexamples[[i]] %>% dplyr::select(HARVEST_DECADE),left=TRUE, largest=TRUE)
aoi.MAL$HARVEST_DECADE[is.na(aoi.MAL$HARVEST_DECADE)] <- 0
habitat[[i]] <- ifelse(aoi.MAL$D2_grp=="Suitable Territory" & aoi.MAL$HARVEST_DECADE==0,1,0)
aoi.MAL$HARVEST_DECADE <- NULL
}

aoi.MAL$Habitat_canBex1 <- habitat[[1]]
aoi.MAL$Habitat_canBex2 <- habitat[[2]]
aoi.MAL$Habitat_canBex3 <- habitat[[3]]
aoi.MAL$Habitat_canBex4 <- habitat[[4]]

aoi.MAL %>% group_by(D2_grp) %>%
  summarise_at(c("Habitat_canBex1","Habitat_canBex2","Habitat_canBex3","Habitat_canBex4"), sum, na.rm = TRUE) %>%
  st_drop_geometry()

glimpse(aoi.MAL)
# because multiple cutblocks within a FETA ended up with 130 FETAs in CanBex1, 121, 118, and 117 in CanBex2, CanBex3, and CanBex4, respectively

Habitat_canBex <- c("canBex1", "canBex2", "canBex3", "canBex4")
aoi <- aoi %>% st_transform(crs=26910)

canBex_raster <- list()
for(i in 1:length(Habitat_canBex)){
raoi <- raster(ext=extent(aoi), crs=26910, res=c(5500,5500))
raoi <- rasterize(aoi.MAL %>% st_transform(crs=26910), raoi, field=paste0("Habitat_",Habitat_canBex[[i]]), fun="max")
raoi[is.na(raoi[])] <- 0

#plot of the raster showing habitat as binary (green = suitable, white = unsuitable)
# sum(raoi@data@values) # 129
# round(sum(raoi@data@values) / length(raoi@data@values)*100) # 27% habitat
Cairo(file=paste0("out/IBM_raster_",Habitat_canBex[i],".PNG"), type="png", width=2200, height=2000,pointsize=15,bg="white",dpi=300)
plot(raoi,main=c(paste0("Fisher Equivalent Territory Areas (FETA) in the Area of Interest,"),
                 paste0(sum(raoi@data@values)," or ",round(sum(raoi@data@values) / length(raoi@data@values)*100),"% FETA in Scenario ",Habitat_canBex[i])),
     cex.main=0.8)
dev.off()
canBex_raster[[i]] <- raoi
}

sum(canBex_raster[[1]]@data@values) # 118 suitable habitat cells
sum(canBex_raster[[2]]@data@values) # 109 suitable habitat cells
sum(canBex_raster[[3]]@data@values) # 106 suitable habitat cells
sum(canBex_raster[[4]]@data@values) # 105 suitable habitat cells

IBM_aoi <- list(aoi=aoi.MAL, canBex_raster=canBex_raster)
save(IBM_aoi, file="data/IBM_aoi_canBex.RData")

save.image("05_create_aoi.RData")
# #####################################################################################

# ###--- view OSM map area of interest (need to figure out how to add polygons...)
# aoi.latlon <- aoi %>% st_transform(crs=4326)
# st_bbox(aoi)
#
# LAT1 = st_bbox(aoi.latlon)[2] ; LAT2 = st_bbox(aoi.latlon)[4]
# LON1 = st_bbox(aoi.latlon)[3] ; LON2 = st_bbox(aoi.latlon)[1]
#
# #our background map
# map <- openmap(c(LAT2,LON1), c(LAT1,LON2), zoom = NULL,
#                type = c("osm", "stamen-toner", "stamen-terrain","stamen-watercolor", "esri","esri-topo")[6],
#                mergeTiles = TRUE)
#
# ## OSM CRS :: "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs"
# map.latlon <- openproj(map, projection = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
#
# OpenStreetMap::autoplot.OpenStreetMap(map.latlon)

#####################################################################################
