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
# 01_load.R
# script to load spatial data for fisher (provincial scale)
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 09-Nov-2021
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

###--- Area of Interest (AOI) is larger than Study Area
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


################################################################################
###--- function to create study grid (i.e., female fisher sized pixels)
# note that this function creates a hexagon grid (use square=TRUE for squares)
# this function originally dissolved multipolygon into single polygon,
# assumes original shapefile is the entire aoi,
# for this exercise, commented out the st_union and st_combine code to keep the two Fisher populations distinct

create_study_grid <- function (dsn=dsn, layer=layer, square=TRUE, cellsize=cellsize, output=output){

  aoi <- st_read(dsn=dsn, layer=layer) %>% st_transform(crs = 3005) # ensures poly is in Albers
  aoi <- aoi %>%
    summarise(across(geometry, ~ st_union(.))) %>%
    summarise(across(geometry, ~ st_combine(.)))
  aoi_utm <- st_transform(aoi, crs=26910) # to have in metres for specifying grid cell size
  aoi_grid <- sa_grid <- st_make_grid(st_bbox(aoi_utm), cellsize=cellsize, square=square) #  grid for entire AOI (rectangle)
  aoi_grid <- aoi_grid %>% st_intersection(aoi_utm)

  sa_points <- st_point_on_surface(sa_grid)  # if using portion of aoi
  sa_points <- st_intersection(sa_points, aoi_utm)

  st_write(aoi_utm, paste0(getwd(),"/out/",output,"_utm.shp"), delete_layer = TRUE)
  st_write(aoi_grid %>% st_intersection(aoi_utm), paste0(getwd(),"/out/",output,"_grid.shp"), delete_layer = TRUE)
  st_write(sa_points, paste0(getwd(),"/out/",output,"_points.shp"), delete_layer = TRUE)

  return(list(aoi_utm, aoi_grid, sa_points))
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
###--- NO NEED TO RUN THIS NEXT SECTION
###--- USED THIS TO CREATE TMP AOI SHAPEFILE FOR EXAMPLE SCENARIOS

# Natural Region District
# bcdc_search("Natural Region District", res_format = "wms")
# 1: Natural Resource (NR) Districts (multiple, wms)
# ID: 0bc73892-e41f-41d0-8d8e-828c16139337
# Name: natural-resource-nr-district
# aoi <- aoi.Fpop
# aoi.NRD <- retrieve_geodata_aoi(ID = "0bc73892-e41f-41d0-8d8e-828c16139337")
# unique(aoi.NRD$DISTRICT_NAME) # 15 unique Natural Resource Districts

# remove the slivers of area, only map polygons with >10km2
# ggplot()+
#   geom_sf(data=aoi.NRD %>% filter(Area_km2>10), aes(fill=DISTRICT_NAME))
#
# aoi.PEACE <- aoi.NRD %>% filter(DISTRICT_NAME=="Peace Natural Resource District")
# ggplot()+
#   geom_sf(data=aoi.PEACE %>% filter(Area_km2>10), aes(fill=DISTRICT_NAME))
# aoi.PEACE$Area_km2
#
# aoi.QUESNEL <- aoi.NRD %>% filter(DISTRICT_NAME=="Quesnel Natural Resource District")
# ggplot()+
#   geom_sf(data=aoi.QUESNEL %>% filter(Area_km2>10), aes(fill=DISTRICT_NAME))
# aoi.QUESNEL$Area_km2

# want to create a fishnet using the aoi extent and have as 30 km2 cells, and then join covariates to each cell
# use 30km2 for now but recognize that some data may be at the smaller scale, and will need to be transformed to fit within the female territory

# tmp <- create_grid(input=aoi.QUESNEL, cellsize=5500)
# ggplot(tmp$fishnet_grid_sf)+
#   geom_sf(data=tmp$fishnet_grid_sf)+
#   geom_sf_label(aes(label=grid_id))

# tmp <- create_grid(input=aoi.PEACE, cellsize=5500)
# ggplot(tmp$fishnet_grid_sf)+
#   geom_sf(data=tmp$fishnet_grid_sf)+
#   geom_sf_label(aes(label=grid_id))
#
# tmp2 <- tmp$fishnet_grid_sf %>% filter(grid_id %in% seq_len(150))
#
# ggplot(tmp2)+
#   geom_sf(data=tmp2)+
#   geom_sf_label(aes(label=grid_id))
#
# tmp3 <- tmp$fishnet_grid_sf %>%
#   filter(grid_id %in% c(61:68,78:85,96:103,111:118,131:138))
#   # filter(grid_id %in% c(1:16,22:29,38:45,59:66))
#
# tmp3$grid_id <- rownames(tmp3)
#
# ggplot(tmp3)+
#   geom_sf(data=tmp3)+
#   geom_sf_label(aes(label=grid_id))
#
# aoi <- tmp3
# st_write(aoi, dsn=paste0(getwd(),"/data/aoi_PTSA_example2.shp"), delete_layer=TRUE)
# st_write(aoi, dsn=paste0(getwd(),"/data/aoi_QTSA_example2.shp"), delete_layer=TRUE)
#####################################################################################
###--- Create an example polygon, with potential cutblocks
# Consolidated cutblock
# 1: Harvested Areas of BC (Consolidated Cutblocks) (multiple, fgdb, wms, kml, pdf)
# ID: b1b647a6-f271-42e0-9cd0-89ec24bce9f7
# Name: harvested-areas-of-bc-consolidated-cutblocks-
# aoi <- st_read(dsn=paste0(getwd(),"/data"), layer="cutblocks_sample") %>% st_transform(crs=26910)
aoi <- st_make_grid(st_buffer(aoi %>% st_transform(crs=26910), dist=30000), n=1)
aoi <- st_as_sf(aoi)

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
#                                     aoi.CUT$HARVEST_DECADE==6 ~1960)
# st_write(aoi.CUT, dsn=paste0(getwd(),"/data/aoi_CUT_example.shp"), delete_layer=TRUE)

aoi.CUT <- st_read(dsn=paste0(getwd(),"/data"), layer="aoi_CUT_example")

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

# First create fishnet for entire area (aoi = buffer of 30 km2 around cutblocks and bounding box)
# Then bring in VRI for aoi
aoi.grid <- create_grid(input=aoi %>% st_transform(crs=26910), cellsize=5500)

aoi.utm <- aoi.grid$aoi_utm
aoi.grid <- aoi.grid$fishnet_grid_sf

ggplot()+
  geom_sf(data=aoi.grid)

# add in population field
aoi.Fpop <- st_read(dsn=paste0(getwd(),"/data"), layer="aoi.Fpop")
aoi.Fpop$Fpop <- case_when(aoi.Fpop$Fpop=="Central Interior" ~ "Columbian",
                           TRUE ~ as.character(aoi.Fpop$Fpop))

aoi.utm <- st_join(aoi.utm, aoi.Fpop %>% st_transform(crs=26910))
aoi.grid <- st_join(aoi.grid, aoi.Fpop %>% st_transform(crs=26910))

ggplot()+
  geom_sf(data=aoi.CUT %>% filter(HARVEST_DECADE==1970), aes(fill=HARVEST_Yr))+
  # geom_sf(data=aoi.CUT %>% filter(HARVEST_DECADE==1980), aes(fill=HARVEST_Yr))+
  geom_sf(data=aoi.utm, fill=NA)


### WILL HAVE TO CREATE "FAKE" SUITABLE HABITAT TO START - I.E., MAKE MORE SUITABLE TO SHOW TOOL USE

# Canned Example 1
aoi.CUT$rndmrnk <- rank(round(runif(nrow(aoi.CUT), min=10000, max=99999)))
# example 1 = all cutblocks
canBEx1 <- aoi.CUT %>% filter(HARVEST_DECADE==1970)
# example 2 = random half of cutblocks
canBEx2 <- aoi.CUT %>% filter(HARVEST_DECADE==1980) %>% arrange(rndmrnk)
canBEx2 <- canBEx2[1:210,]
# example 3 = "unsuitable quality" cuts
tmp <- st_join(aoi.CUT %>% filter(HARVEST_DECADE==1960), aoi.VRI %>% dplyr::select(TREE20_prop),left=TRUE, largest=TRUE)
tmp %>% dplyr::count(TREE20_prop)

# Choose



################################################################################
###--- UPDATE - 2022-02-23 ---###
###- PICK A SMALL AOI AS PROOF OF CONCEPT (4 x 4 fisher cells)
###- USE SQUARE CELLS TO EMULATE IBM
###- USE BASIC VRI CONDITION FOR GOOD HABITAT

# now reduce aoi to proof-of-concept study area (same as grid above - Prince George TSA)
# female fisher home range = 5.5 * 5.5 = 30 km2
# the CLUS model (Kyle Lochhead and Tyler Muhly) analyses at ha pixel size and sums up to fisher home range
# will want to go even smaller, but for now try the Prince George TSA as the study area (aoi)

# takes too long for such a huge study area
# for needs of a small test model, choose a portion of a RDA and go from there
# eventually create 100 m square pixels - keep in mind that this means 3025 square pixels in a female home range
# start with fisher sized square pixels

###--- LOAD COVARIATES
# Use cached data, downloaded through the bcgeodata warehouse website for big files that are updated annually
# download only needs to happen once per year or as notified of updates
# for this example, where possible, downloaded to aoi (Prince George TSA - Prince George TSA Block E) extent
# https://catalogue.data.gov.bc.ca/

# to see the downloaded gdb and shp files
# list.files("./data")

aoi <- aoi.utm %>% st_transform(crs=3005)

# vegetation data (VRI)
# bcdc_search("VRI", res_format = "wms")
aoi.VRI <- retrieve_geodata_aoi(ID = "2ebb35d8-c82f-4a17-9c96-612ac3532d55")

aoi.VRI %>% filter(!is.na(PROJ_HEIGHT_1)) %>%
  summarise(mean = mean(PROJ_HEIGHT_1), min = min(PROJ_HEIGHT_1), max=max(PROJ_HEIGHT_1), sd = sd(PROJ_HEIGHT_1)) %>% st_drop_geometry()
#mean   min   max    sd
#17.3   0.1  40.2  7.57

# proportion of VRI with projected height >= 20 m
# aoi.VRI$TREE20_prop <- NA
# aoi.VRI$TREE20_prop <- ifelse(aoi.VRI$PROJ_HEIGHT_1>=20, 1,0)
# aoi.VRI$TREE20_prop[is.na(aoi.VRI$TREE20_prop)] <- 0

sum(aoi.VRI$Area_km2)
# proportion of VRI with projected height >= 20 m
aoi$TREE20_prop <- NA
for(i in seq_len(nrow(aoi))){
  tmp <- aoi %>% st_transform(crs=26910) %>% filter(aoi$grid_id==i)
  tmp2 <- st_intersection(tmp, aoi.VRI %>% filter(PROJ_HEIGHT_1>=20) %>% st_transform(crs=26910))
  cov.area <- sum(drop_units(tmp2 %>% st_area()*1e-6))
  cov.prop <- cov.area/tmp$Area_km2
  aoi$TREE20_prop[i] <- cov.prop
}

ggplot()+
  geom_sf(data=aoi.VRI, aes(fill=TREE20_prop))
summary(aoi$TREE20_prop)

aoi$Habitat <- ifelse(aoi$TREE20_prop>0.2, 1, 0)

ggplot()+
  geom_sf(data=aoi, aes(fill=Habitat))

aoi <- aoi %>% st_transform(crs=26910)
raoi <- raster(ext=extent(aoi), crs=26910, res=c(5500,5500))
raoi <- rasterize(aoi, raoi, field="Habitat")

IBM_aoi <- list(aoi=aoi, raoi=raoi)
# sum(IBM_aoi$raoi@data@values)
# nrow(IBM_aoi$aoi)


#plot of the raster showing habitat as binary (green = suitable, white = unsuitable)
Cairo(file=paste0("out/IBM_raoi_Pex2.PNG"), type="png", width=2200, height=2000,pointsize=15,bg="white",dpi=300)
plot(IBM_aoi$raoi, legend=FALSE,main="Potential Fisher Equivalent Territories in the Area of Interest")
dev.off()

# plot of the underlying habitat as a gradation of suitable habitat
Cairo(file=paste0("out/IBM_haoi_Pex2.PNG"), type="png", width=2200, height=2000,pointsize=15,bg="white",dpi=300)
ggplot()+
  geom_sf(data=IBM_aoi$aoi, aes(fill=TREE20_prop))+
  guides(fill=guide_legend(title="Proportion Suitable Habitat"))+
  theme(legend.position="bottom")+
  ggtitle("Fisher Equivalent Territories within the Area of Interest")
dev.off()

save(IBM_aoi, file="data/IBM_aoi_Pex2.RData")



# rm(aoi.Fpop); rm(aoi.CUT) # for housekeeping, big file
save.image("01_load.RData")
# #####################################################################################

# ###--- view OSM map area of interest (need to figure out how to add polygons...)
aoi.latlon <- aoi %>% st_transform(crs=4326)
st_bbox(aoi)

LAT1 = st_bbox(aoi.latlon)[2] ; LAT2 = st_bbox(aoi.latlon)[4]
LON1 = st_bbox(aoi.latlon)[3] ; LON2 = st_bbox(aoi.latlon)[1]

#our background map
map <- openmap(c(LAT2,LON1), c(LAT1,LON2), zoom = NULL,
               type = c("osm", "stamen-toner", "stamen-terrain","stamen-watercolor", "esri","esri-topo")[6],
               mergeTiles = TRUE)

## OSM CRS :: "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs"
map.latlon <- openproj(map, projection = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

OpenStreetMap::autoplot.OpenStreetMap(map.latlon)
ggplot()+
  geom_sf(data=aoi.latlon)

#####################################################################################
