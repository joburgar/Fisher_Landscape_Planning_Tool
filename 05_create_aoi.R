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
list.of.packages <- c("tidyverse", "lubridate","chron","bcdata", "bcmaps","sf", "rgdal","viridis",
                      "Cairo","OpenStreetMap", "ggmap","PNWColors","units","nngeo","raster","qs")

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


se <- function(x) sqrt(var(x)/length(x))

#####################################################################################
###--- CREATE AGGREGATED RASTER STACK OF FETA VALUES
#####################################################################################
# Created, and no need to re-run, just pull in appropriate raster files

# First create fishnet for entire area (aoi = buffer of 30 km2 around cutblocks and bounding box)
# Then bring in underlying habitat for aoi
# make raster stack for Fisher population (combined) extent
# aoi.Fpop <- st_read(dsn=paste0(getwd(),"/data"), layer="aoi.Fpop_simpl")
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
###--- EXPLORE FOR LANDSCAPE CONFIGURATION EXPERIMENTAL DESIGN
#####################################################################################

Fpop <- st_read(dsn=paste0(getwd(),"/data"), layer="aoi.Fpop_simpl") %>% st_transform(crs=3005)
FHEzones <- st_read(dsn=paste0(getwd(),"/data"), layer="FHE_zones")
FHEzones <- st_simplify(FHEzones %>% st_transform(crs=26910), preserveTopology = FALSE, dTolerance = 1000) %>% st_transform(crs=3005)

ggplot()+
  geom_sf(data=FHEzones, aes(fill=Hab_zone)) +
    geom_sf(data=Fpop, aes(fill=Fpop))

aoi <- Fpop %>%
  summarise(across(geometry, ~ st_union(.))) %>%
  summarise(across(geometry, ~ st_combine(.)))

# Natural Resource Region
# 1: Natural Resource (NR) Regions (multiple, wms, oracle_sde)
# ID: dfc492c0-69c5-4c20-a6de-2c9bc999301f
# Name: natural-resource-nr-regions
# 2: Natural Resource (NR) Areas (multiple, wms, oracle_sde)
# ID: c1861ba4-abb8-4947-b3e5-7f7c4d7257d5
# Name: natural-resource-nr-areas
aoi.NRR <- retrieve_geodata_aoi(ID = "dfc492c0-69c5-4c20-a6de-2c9bc999301f")
unique(aoi.NRR$REGION_NAME) # 6 unique Natural Resource Regions


# bring in the different NR Districts
# Natural Region District
# bcdc_search("Natural Resource District", res_format = "wms")
# 1: Natural Resource (NR) Districts (multiple, wms, oracle_sde)
# ID: 0bc73892-e41f-41d0-8d8e-828c16139337
# Name: natural-resource-nr-district
aoi.NRD <- retrieve_geodata_aoi(ID = "0bc73892-e41f-41d0-8d8e-828c16139337")
unique(aoi.NRD$DISTRICT_NAME) # 15 unique Natural Resource Districts
aoi.NRD %>% group_by(REGION_ORG_UNIT_NAME, DISTRICT_NAME) %>% summarise(sum(Area_km2)) %>% st_drop_geometry()
aoi.NRD %>% group_by(DISTRICT_NAME,REGION_ORG_UNIT_NAME) %>% summarise(sum(Area_km2)) %>% st_drop_geometry()

# bring in the different Timber Supply Areas
# bcdc_search("Timber Supply Area", res_format = "wms")
# 1: FADM - Timber Supply Area (TSA) (multiple, wms, kml, oracle_sde)
# ID: 8daa29da-d7f4-401c-83ae-d962e3a28980
# Name: fadm-timber-supply-area-tsa
aoi.TSA <- retrieve_geodata_aoi(ID = "8daa29da-d7f4-401c-83ae-d962e3a28980")
# TSA is a BIG file, simplify to 100 m lines to reduce?
format(object.size(aoi.TSA), units="Mb")
aoi.TSA <- st_simplify(aoi.TSA %>% st_transform(crs=26910), preserveTopology = FALSE, dTolerance = 100) %>% st_transform(crs=3005)
format(object.size(aoi.TSA), units="Mb") # went from 24.7 MB to 1.9 MB
aoi.TSA$TSA_NUMBER <- as.numeric(aoi.TSA$TSA_NUMBER)

length(unique(aoi.TSA$TSA_NUMBER)) # 25 unique TSAs but 181 features...need to sum area per feature
options(scipen = 10)
as.data.frame(aoi.TSA %>% group_by(TSA_NUMBER, TSA_NUMBER_DESCRIPTION) %>% summarise(sum(Area_km2)) %>% st_drop_geometry())

# bring in the different Tree Farm Licenses
bcdc_search("Tree Farm License Current View", res_format = "wms")
# 1: FADM - Tree Farm License Current View (TFL) (wms, kml, multiple)
# ID: 454f2153-efbd-4a6e-8966-a6d9755da9a6
# Name: fadm-tree-farm-license-current-view-tfl-
aoi.TFL <- retrieve_geodata_aoi(ID = "454f2153-efbd-4a6e-8966-a6d9755da9a6") # not a ton of info here...


# bring in the different Land Use Plan boundaries
bcdc_search("Land Use Plan", res_format = "wms")
# 1: Strategic Land and Resource Plans - Current (multiple, wms, kml, other, oracle_sde)
# ID: 4b142d4c-83d6-4ecc-b66c-66601ae65992
# Name: strategic-land-and-resource-plans-current
aoi.LUP <- retrieve_geodata_aoi(ID = "4b142d4c-83d6-4ecc-b66c-66601ae65992")
length(unique(aoi.LUP$STRGC_LAND_RSRCE_PLAN_ID))


# to union each TSA number
aoi.TSA.union <-aoi.TSA %>% group_by(TSA_NUMBER,TSA_NUMBER_DESCRIPTION) %>%
  summarise(across(geometry, ~ st_union(.)), .groups = "keep") %>%
  summarise(across(geometry, ~ st_combine(.)))
aoi.TSA.union$Area_km2 <- st_area(aoi.TSA.union)*1e-6
aoi.TSA.union <- drop_units(aoi.TSA.union)
aoi.TSA.union <- droplevels.sfc(aoi.TSA.union, except = geometry)
aoi.TSA.union %>% arrange(Area_km2)

# Summary stats of area within fisher habitat per TSA, NRD, and LUP
summary(aoi.TSA.union$Area_km2)
summary(aoi.NRD$Area_km2)
summary(aoi.LUP$Area_km2)


# Histograms showing range of sizes per admin type
# TSA
aoi.TSA.union$TSA_Name <- aoi.TSA.union$TSA_NUMBER_DESCRIPTION %>% str_replace(" TSA","")
aoi.TSA.union_FHEzones <- st_distance(aoi.TSA.union %>% st_transform(26910), FHEzones %>% st_transform(26910))
colnames(aoi.TSA.union_FHEzones)<-FHEzones$Hab_zone
aoi.TSA.union_FHEzones <- drop_units(aoi.TSA.union_FHEzones)
aoi.TSA.union_FHEzones <- ifelse(aoi.TSA.union_FHEzones>0,0,1) # a 1 if that habitat zone is in the district
aoi.TSA.union <- cbind(aoi.TSA.union, aoi.TSA.union_FHEzones)

TSA_NRD <- st_distance(aoi.TSA.union,aoi.NRD)
rownames(TSA_NRD) <- aoi.TSA.union$TSA_NUMBER_DESCRIPTION
colnames(TSA_NRD) <- aoi.NRD$DISTRICT_NAME
TSA_NRD <- drop_units(TSA_NRD)
TSA_NRD <- ifelse(TSA_NRD>0,0,1) # a 1 if that habitat zone is in the district
TSA_NRD <- as.data.frame(TSA_NRD)
TSA_NRD$Fpop <- aoi.TSA.union$Fpop
TSA_NRD$Area_km2 <- aoi.TSA.union$Area_km2
TSA_NRD <- TSA_NRD %>% arrange(Fpop, desc(Area_km2))

write.csv(TSA_NRD, "out/TSA_NRD_distance.csv")

Cairo(file=paste0("out/TSA_map_plot.PNG"), type="png", width=2000, height=2200,pointsize=15,bg="white",dpi=300)
ggplot() +
  geom_sf(data=aoi.TSA.union, aes(fill=TSA_Name), lwd=0.2, col="black")+
  scale_fill_viridis_d(option="D")+
  theme(legend.position = "bottom", legend.title = element_blank())+
  ggtitle("Timber Supply Areas within Fisher Habitat")
dev.off()

aoi.TSA.union$Fpop <- as.factor(ifelse(aoi.TSA.union$Boreal==1,"Boreal","Columbian"))

pal=pnw_palette("Lake",2, type = "discrete")

Cairo(file=paste0("out/TSA_area_plot.PNG"), type="png", width=2200, height=1400,pointsize=15,bg="white",dpi=300)
ggplot(aoi.TSA.union %>% filter(Area_km2>3), aes(x = reorder(TSA_Name, -Area_km2), y = Area_km2, fill=Fpop)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values=pal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Timber Supply Areas within Fisher Habitat") +
  theme(axis.title.x=element_blank()) +
  ylab(bquote('Area '(km^2)))
dev.off()

TSA.area.summarised <- aoi.TSA.union %>% filter(Area_km2>3) %>% group_by(Fpop) %>%
  summarise(sum=sum(Area_km2),mean=mean(Area_km2),se=se(Area_km2),min=min(Area_km2), max=max(Area_km2)) %>% st_drop_geometry()

write.csv(TSA.area.summarised,"out/TSA.area.summarised.csv")


# NRD and NRR
aoi.NRD$Name <- aoi.NRD$DISTRICT_NAME %>% str_replace(" Natural Resource District","")
aoi.NRD$Region <- aoi.NRD$REGION_ORG_UNIT_NAME %>% str_replace(" Natural Resource Region","")
NRD_FHEzones <- st_distance(aoi.NRD %>% st_transform(26910), FHEzones %>% st_transform(26910))
colnames(NRD_FHEzones)<-FHEzones$Hab_zone
NRD_FHEzones <- drop_units(NRD_FHEzones)
NRD_FHEzones <- ifelse(NRD_FHEzones>0,0,1) # a 1 if that habitat zone is in the district
aoi.NRD <- cbind(aoi.NRD, NRD_FHEzones)

write.csv(as.data.frame(aoi.NRD %>% group_by(Name,Region) %>% summarise(sum(Area_km2)) %>% st_drop_geometry()),
          "out/Natural_Resource_Districts_in_Regions.csv", row.names = F)


pal=pnw_palette("Cascades",6, type = "discrete")
glimpse(aoi.NRD)

Cairo(file=paste0("out/NRD_area_plot.PNG"), type="png", width=2200, height=1400,pointsize=15,bg="white",dpi=300)
ggplot(aoi.NRD %>% filter(Area_km2>3), aes(x = reorder(Name, -Area_km2), y = Area_km2, fill=Region)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values=pal)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("Natural Resource Districts within Fisher Habitat") +
  theme(axis.title.x=element_blank()) +
  ylab(bquote('Area '(km^2)))
dev.off()

Cairo(file=paste0("out/NRD_map_plot.PNG"), type="png", width=2200, height=2000,pointsize=15,bg="white",dpi=300)
ggplot() +
  geom_sf(data=aoi.NRD %>% filter(Area_km2>3), aes(fill=Name), lwd=0.2, col="black")+
  scale_fill_viridis_d(option="D")+
  theme(legend.position = "right", legend.title = element_blank())+
  ggtitle("Natural Resource Districts within Fisher Habitat")
dev.off()

aoi.NRD %>% arrange(Area_km2)

## LUP
glimpse(aoi.LUP)
aoi.LUP %>% group_by(PLAN_TYPE)%>% count(STRGC_LAND_RSRCE_PLAN_NAME) %>% st_drop_geometry()
aoi.LUP %>% count(PLAN_TYPE)%>%st_drop_geometry()
aoi.LUP %>% filter(grepl("Cariboo",STRGC_LAND_RSRCE_PLAN_NAME))
as.data.frame(aoi.LUP %>% filter(PLAN_TYPE=="LRMP") %>% group_by(PLAN_STATUS) %>% count(STRGC_LAND_RSRCE_PLAN_NAME) %>% st_drop_geometry())

write.csv(aoi.LUP %>% st_drop_geometry(),"out/LUP_FHEzones.csv", row.names = F)

LUP.area.summarised <- aoi.LUP %>% filter(Area_km2>3) %>% filter(PLAN_TYPE==c("SRMP","LRMP")) %>% group_by(PLAN_TYPE, Fpop) %>%
  summarise(sum=sum(Area_km2),mean=mean(Area_km2),se=se(Area_km2),min=min(Area_km2), max=max(Area_km2)) %>% st_drop_geometry()

write.csv(LUP.area.summarised,"out/LUP.area.summarised.csv")

# use LRMP & SRMP only
Cairo(file=paste0("out/LRMP_SRMP_area_plot.PNG"), type="png", width=2200, height=1400,pointsize=15,bg="white",dpi=300)
ggplot(aoi.LUP %>% filter(PLAN_TYPE==c("LRMP","SRMP")), aes(x = reorder(STRGC_LAND_RSRCE_PLAN_ID, -Area_km2), y = Area_km2, fill=Fpop)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values=pal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle("LRMPs and SRMPs within Fisher Habitat") +
  xlab("Plan ID")+
  ylab(bquote('Area '(km^2)))
dev.off()

aoi.LRMP <- aoi.LUP %>% filter(PLAN_TYPE=="LRMP")

LRMP_FHEzones <- st_distance(aoi.LRMP %>% st_transform(26910), FHEzones %>% st_transform(26910))
colnames(LRMP_FHEzones)<-FHEzones$Hab_zone
LRMP_FHEzones <- drop_units(LRMP_FHEzones)
LRMP_FHEzones <- ifelse(LRMP_FHEzones>0,0,1) # a 1 if that habitat zone is in the district
aoi.LRMP <- cbind(aoi.LRMP, LRMP_FHEzones)

pal=pnw_palette("Lake",2, type = "discrete")
aoi.LRMP$Fpop <- as.factor(ifelse(aoi.LRMP$Boreal==1,"Boreal","Columbian"))
aoi.LRMP$Name <- aoi.LRMP$STRGC_LAND_RSRCE_PLAN_NAME %>% str_replace(" Land and Resource Management Plan","")
aoi.LRMP$Name <- aoi.LRMP$Name %>% str_replace(" Land Resource Management Plan","")
aoi.LRMP$Abb_Name <- aoi.LRMP$STRGC_LAND_RSRCE_ABRVN %>% str_replace("LRMP","")


Cairo(file=paste0("out/LRMP_area_plot.PNG"), type="png", width=2200, height=1400,pointsize=15,bg="white",dpi=300)
ggplot(aoi.LRMP, aes(x = reorder(Name, -Area_km2), y = Area_km2, fill=Fpop)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values=pal)+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+
  ggtitle("Land and Resource Management Plans within Fisher Habitat") +
  ylab(bquote('Area '(km^2)))
dev.off()

aoi.LRMP <- droplevels.sfc(aoi.LRMP)

Cairo(file=paste0("out/LRMP_map_plot.PNG"), type="png", width=2000, height=2200,pointsize=15,bg="white",dpi=300)
ggplot() +
  geom_sf(data=aoi.LRMP, aes(fill=as.factor(STRGC_LAND_RSRCE_PLAN_ID)), lwd=0.2, col="black")+
  scale_fill_viridis_d(option="D")+
  theme(legend.position = "right", legend.title = element_blank())+
  ggtitle("Land and Resource Management Plans \nwithin Fisher Habitat")
dev.off()


aoi.SRMP <- aoi.LUP %>% filter(PLAN_TYPE=="SRMP")
aoi.SRMP <- droplevels.sfc(aoi.SRMP)
Cairo(file=paste0("out/SRMP_map_plot.PNG"), type="png", width=2000, height=2200,pointsize=15,bg="white",dpi=300)
ggplot() +
  geom_sf(data=aoi.LUP)+
  geom_sf(data=aoi.SRMP,aes(fill=as.factor(STRGC_LAND_RSRCE_PLAN_ID)), lwd=0.2, col="black")+
  scale_fill_viridis_d(option="D")+
  theme(legend.position = "right", legend.title = element_blank())+
  ggtitle("Strategic Resource Management Plans \nwithin Fisher Habitat")
dev.off()
#####################################################################################
###--- CREATE CANNED SCENARIO AOI SHAPEFILE
#####################################################################################
# Bring in population polygons
# aoi.Fpop <- st_read(dsn=paste0(getwd(),"/data"), layer="aoi.Fpop")
# aoi.Fpop$Fpop <- case_when(aoi.Fpop$Fpop=="Central Interior" ~ "Columbian",
#                            TRUE ~ as.character(aoi.Fpop$Fpop))
# ggplot()+
#   geom_sf(data=aoi.Fpop)
#
# aoi.Fpop_simpl <- st_simplify(aoi.Fpop %>% st_transform(crs=26910), preserveTopology = FALSE, dTolerance = 1000) %>% st_transform(crs=3005)
# # preserveToplogy=TRUE ensures that no polygons or inner holes are reduced to lines/removed - not necessary because just boundary for other intersections
#
# ggplot()+
#   geom_sf(data=aoi.Fpop_simpl, aes(fill=Fpop))

# st_write(aoi.Fpop_simpl, paste0(getwd(),"/data/aoi.Fpop_simpl.shp"), delete_layer=TRUE)
# round(c(object.size(aoi.Fpop),object.size((aoi.Fpop_simpl)))/1024) # huge reduction in memory usage 100x magnitude less memory storage!

aoi.Fpop_simpl <- st_read(dsn="./data", layer="aoi.Fpop_simpl")

FHE_zones <- st_read(dsn=paste0(getwd(),"/data"), layer="FHE_zones")
FHE_zones <- st_simplify(FHE_zones %>% st_transform(crs=26910), preserveTopology = FALSE, dTolerance = 1000) %>% st_transform(crs=3005)

# bring in the different Natural Resource Districts
# Natural Region District
# bcdc_search("Natural Region District", res_format = "wms")
# 1: Natural Resource (NR) Districts (multiple, wms)
# ID: 0bc73892-e41f-41d0-8d8e-828c16139337
# Name: natural-resource-nr-district
aoi <- aoi.Fpop_simpl
aoi.NRD <- retrieve_geodata_aoi(ID = "0bc73892-e41f-41d0-8d8e-828c16139337")
unique(aoi.NRD$DISTRICT_NAME) # 15 unique Natural Resource Districts

# create example for each of the Regions, choosing one Natural Resource District
aoi.NRD %>% filter(grepl("Omineca",REGION_ORG_UNIT_NAME)) %>% group_by(Fpop,OBJECTID) %>% summarise(sum(Area_km2))
# small amount of Mackenzie as Boreal pop...should this be? Should it be all the way to end of Boreal Hab_zone?
# keep as is for now, to consider for later iterations of FLEX and shiny

# For Omineca use Prince George Natural Resource District
ggplot()+
  geom_sf(data=aoi.NRD %>% filter(grepl("Omineca",REGION_ORG_UNIT_NAME)), aes(fill=DISTRICT_NAME))
aoi.Omineca <- aoi.NRD %>% filter(grepl("Prince",DISTRICT_NAME))

# ggplot()+
#   geom_sf(data=aoi.Fpop_simpl)+
#   geom_sf(data = aoi.Omineca, col="blue", fill="blue")

# For Skeena use Nadina Natural Resource District
ggplot()+
  geom_sf(data=aoi.NRD %>% filter(grepl("Skeena",REGION_ORG_UNIT_NAME)), aes(fill=DISTRICT_NAME))
aoi.Skeena <- aoi.NRD %>% filter(grepl("Nadina",DISTRICT_NAME))

# ggplot()+
#   geom_sf(data=aoi.Fpop_simpl)+
#   geom_sf(data = aoi.Skeena, col="blue", fill="blue")

# For Northeast use Peace Natural Resource District
ggplot()+
  geom_sf(data=aoi.NRD %>% filter(grepl("Northeast",REGION_ORG_UNIT_NAME)), aes(fill=DISTRICT_NAME))
aoi.Northeast <- aoi.NRD %>% filter(grepl("Peace",DISTRICT_NAME))

# ggplot()+
#   geom_sf(data=aoi.Fpop_simpl)+
#   geom_sf(data = aoi.Northeast, col="blue", fill="blue")

# For Cariboo use Cariboo-Chilcotin Natural Resource District
ggplot()+
  geom_sf(data=aoi.NRD %>% filter(grepl("Cariboo",REGION_ORG_UNIT_NAME)), aes(fill=DISTRICT_NAME))
aoi.Cariboo <- aoi.NRD %>% filter(grepl("Cariboo",DISTRICT_NAME))

# ggplot()+
#   geom_sf(data=aoi.Fpop_simpl)+
#   geom_sf(data = aoi.Cariboo, col="blue", fill="blue")

###--- Create artificial example with exaggerated scenarios
# One for each region, may need to subset if too large
# 1. Create grid of FETAs
# 2. Identify Mahalanobis distance for each FETA
# 3. Keep as probability (consistent with Rich Weir's usage)
# 4. Generate random numbers for FETAs and change a portion of suitable to movement / unsuitable
# 1. all cutblocks harvested (420)
# 2. half cutblocks harvested (210) random selection
# 3. half cutblocks harvested (210), in 'lower quality' (i.e., unsuitable) habitat
# 4. half cutblocks harvested (210), in 'higher quality' (i.e., suitable) habitat

# Create examples for Omineca, Skeena, Northeast, and Cariboo

nrd_examples <- c("Omineca","Skeena","Northeast","Cariboo")
nrd_ex_Laoi <- list(aoi.Omineca, aoi.Skeena, aoi.Northeast, aoi.Cariboo)

for(i in 1:length(nrd_examples)){
  nrd_name <- nrd_examples[i]
  nrd_ex_aoi <- nrd_ex_Laoi[[i]]

  aoi <- st_make_grid(st_buffer(nrd_ex_aoi %>% st_transform(crs=26910), dist=30000), n=1)
  aoi <- st_as_sf(aoi)

  pal = pnw_palette(name="Winter",n=2,type="discrete")

  Cairo(file=paste0("out/Ex_",nrd_name,"_Fpop_aoi_plot.PNG"),type="png",width=2400,height=2400,pointsize=16,bg="white",dpi=300)
  ggplot()+
    geom_sf(data=aoi.Fpop_simpl, aes(fill=Fpop))+
    geom_sf(data=aoi)+
    geom_sf(data=nrd_ex_aoi)+
    scale_fill_manual(values=rev(pal))+
    theme(legend.position="bottom")+
    theme(legend.title=element_blank())
  dev.off()

  aoi <- st_join(aoi, aoi.Fpop_simpl %>% st_transform(crs=26910),left=TRUE, largest=TRUE)
  aoi <- aoi %>% st_transform(crs=3005)
  aoi.MAL <- retrieve_gdb_shp_aoi(dsn=paste0(getwd(),"/data"), layer="Mahalanobis_predictions_2021_220304")
  summary(aoi.MAL); nrow(aoi.MAL) # 2785 FETA in Omineca example; 1651 in Skeena; 3337 in Northeast; 2471 in Cariboo

  # aoi.MAL <- aoi.MAL %>% mutate(D2_grp = cut(D2, c(0, 6, 10, 20, Inf))) # this is the input for Mahalanobis threshold
  # levels(aoi.MAL$D2_grp)
  # aoi.MAL$D2_grp <- recode(aoi.MAL$D2_grp, "(0,6]"="Suitable Territory",
  #                   "(6,10]"="Possible Territory",
  #                   "(10,20]"="Suitable Movement",
  #                   "(20,Inf]"="Unsuitable")
  # glimpse(aoi.MAL)
  # aoi.MAL %>% group_by(D2_grp) %>% summarize(n=n(), min=min(p), max=max(p)) %>% st_drop_geometry()
  # Omineca example - Prince George
  # D2_grp                 n      min     max
  # 1 Suitable Territory    67 0.207    0.992
  # 2 Possible Territory   107 0.0445   0.303
  # 3 Suitable Movement    303 0.000503 0.0720
  # 4 Unsuitable          2307 0        0.00123
  # 5 NA                     1 0        0

  # Skeena example - Nadina
  # D2_grp                 n      min     max
  # 1 Suitable Territory    25 0.207    0.961
  # 2 Possible Territory    43 0.0432   0.291
  # 3 Suitable Movement    126 0.000540 0.0676
  # 4 Unsuitable          1451 0        0.00122

  # Northeast example - Peace
  # D2_grp                 n      min      max
  # 1 Suitable Territory   188 2.00e- 1 0.981
  # 2 Possible Territory   265 4.10e- 2 0.303
  # 3 Suitable Movement    642 5.03e- 4 0.0672
  # 4 Unsuitable          2239 3.49e-52 0.000835

  # Cariboo example - Cariboo-Chilcotin
  # D2_grp                 n      min     max
  # 1 Suitable Territory   101 0.207    0.992
  # 2 Possible Territory   159 0.0405   0.299
  # 3 Suitable Movement    461 0.000504 0.0720
  # 4 Unsuitable          1734 0        0.00123


  # pal = pnw_palette(name="Cascades",n=4,type="discrete")
  #
  # aoi.MAL.plot <- ggplot()+
  #   geom_sf(data=aoi.MAL %>% filter(!is.na(D2_grp)), aes(fill=D2_grp))+
  #   geom_sf(data=aoi, fill=NA)+
  #   scale_fill_manual(values=rev(pal))+
  #   theme(legend.position="bottom")+
  #   theme(legend.title=element_blank())+
  #   ggtitle(paste0("Fisher Equivalent Territory Areas in the ",nrd_name," Scenario"))
  #
  #
  # Cairo(file=paste0("out/Ex_",nrd_name,"_aoi_MAL_plot.PNG"),type="png",width=2400,height=2400,pointsize=18,bg="white",dpi=300)
  # aoi.MAL.plot
  # dev.off()

  ################################################################################
  ###--- DECISION 4 May 2022 - keep D2 values as raster and use actual value in IBM
  # No need to create individual scenarios, will do that as part of SpaDES parameterization

  aoi <- aoi %>% st_transform(crs=26910)
  raoi <- raster(ext=extent(aoi), crs=26910, res=c(5500,5500))
  rMahal <- rasterize(aoi.MAL %>% st_transform(crs=26910), raoi, field="D2", fun="min")
  # raoi[is.na(raoi[])] <- 0 # keep NAs as NAs (at least for now)

  tmpST <- rMahal<=6

  #plot of the raster showing habitat as binary (green = suitable, white = gray as unsuitable and white as NA [not fisher habitat])
  sum(tmpST@data@values, na.rm=TRUE) # 110
  # round(sum(tmpST@data@values) / length(raoi@data@values)*100) # 1% habitat
  Cairo(file=paste0("out/IBM_raster_",nrd_name,"_Scenario.PNG"), type="png", width=2200, height=2000,pointsize=15,bg="white",dpi=300)
  plot(tmpST,main=c(paste0("Fisher Equivalent Territory Areas (FETA) in the ",nrd_name," Scenario,"),
                    paste0(sum(tmpST@data@values,na.rm=T)," or ",round(sum(tmpST@data@values,na.rm=T) / length(rMahal@data@values)*100),"% potentially suitable FETAs")),
       cex.main=0.8, legend=FALSE)
  dev.off()

  aoi.MAL$Fpop_num <- ifelse(aoi.MAL$Fpop=="Boreal",1,2)
  rFpop <- rasterize(aoi.MAL %>% st_transform(crs=26910), raoi, field="Fpop_num", fun="min")

  FHE_zones$Hab_zone_num <- ifelse(FHE_zones$Hab_zone=="Boreal",1,
                                   ifelse(FHE_zones$Hab_zone=="Sub-Boreal moist",2,
                                          ifelse(FHE_zones$Hab_zone=="Sub-Boreal dry",3,
                                                 ifelse(FHE_zones$Hab_zone=="Dry Forest",4,NA))))
  rFHzone <- rasterize(FHE_zones %>% st_transform(crs=26910), raoi, field="Hab_zone_num", fun="min")

  # dynamic raster = raster brick with rasters at 5 year intervals; static raster = single raster layer each; pixel size = 5.5*5.5 or 1 FETA
  # for now dummy variables in dynamic rasters
  rMove <- rMahal <- stack(rMahal, rMahal)

  values(rMove[[1]]) = runif(ncell(rMove[[1]]), 0, 1)
  values(rMove[[2]]) = runif(ncell(rMove[[1]]), 0, 1)

  # r_dynamic <- stack(rMahal, rMove)  # raster stack for dynamic values = rMahal for mahalanobis distances; rMove = movement values
  # r_static <-  stack(rFpop, rFHzone) # raster stack for constant values = Fpop for population, rFHzone for Mahalanobis dist value

  r_list <- list(rMahal=rMahal, rMove=rMove, rFpop=rFpop, rFHzone=rFHzone)

  myfile <- paste0(getwd(),"/data/EX_",nrd_name,"_IBM_aoi.qs")
  qsave(r_list, myfile)
}

# ################################################################################
# ###--- Set up 4 possible scenarios
# # (1) BAU - keep FETAs as currently mapped; can only establish in 'Suitable'
# # (2) Increase FETAs - can establish in 'suitable' and 'possible'
# # (3) Increase movement - can establish in 'suitable' and move in all FETAs
# # (4) Increase FETAs & movement - can establish in 'suitable' and 'possible' and move in all FETAs
#
# aoi.MAL$Scenario1 <- ifelse(aoi.MAL$D2_grp=="Suitable Territory",2,
#                             ifelse(aoi.MAL$D2_grp=="Possible Territory" | aoi.MAL$D2_grp=="Suitable Movement",1,0))
#
# aoi.MAL$Scenario2 <- ifelse(aoi.MAL$D2_grp=="Suitable Territory" | aoi.MAL$D2_grp=="Possible Territory",2,
#                             ifelse(aoi.MAL$D2_grp=="Suitable Movement",1,0))
#
# aoi.MAL$Scenario3 <- ifelse(aoi.MAL$D2_grp=="Suitable Territory",2,1)
#
# aoi.MAL$Scenario4 <- ifelse(aoi.MAL$D2_grp=="Suitable Territory" | aoi.MAL$D2_grp=="Possible Territory",2,1)
#
# aoi.MAL <- aoi.MAL[complete.cases(aoi.MAL$D2_grp),]
# aoi.MAL %>% count(Scenario1) %>% st_drop_geometry() # 67 suitable FETAs, 410 movement FETAs - Omineca example
# aoi.MAL %>% count(Scenario2) %>% st_drop_geometry() # 174 suitable FETAs, 303 movement FETAs
# aoi.MAL %>% count(Scenario3) %>% st_drop_geometry() # 67 suitable FETAs, 2717 movement FETAs
# aoi.MAL %>% count(Scenario4) %>% st_drop_geometry() # 174 suitable FETAs, 2610 movement FETAs
#
# aoi <- aoi %>% st_transform(crs=26910)
# Ex_raster <- list()
# for(i in 1:4){
#   raoi <- raster(ext=extent(aoi), crs=26910, res=c(5500,5500))
#   raoi <- rasterize(aoi.MAL %>% st_transform(crs=26910), raoi, field=paste0("Scenario",i), fun="max")
#   raoi[is.na(raoi[])] <- 0
#
#   tmpST <- raoi==2
#   #plot of the raster showing habitat as binary (green = suitable, white = unsuitable)
#   # sum(tmpST@data@values) # 62
#   # round(sum(tmpST@data@values) / length(raoi@data@values)*100) # 1% habitat
#   Cairo(file=paste0("out/IBM_raster_",nrd_name,"_Scenario",i,".PNG"), type="png", width=2200, height=2000,pointsize=15,bg="white",dpi=300)
#   plot(raoi,main=c(paste0("Fisher Equivalent Territory Areas (FETA) in the Area of Interest,"),
#                    paste0(sum(tmpST@data@values)," or ",round(sum(tmpST@data@values) / length(raoi@data@values)*100),"% FETA in ",nrd_name," Scenario ",i)),
#        cex.main=0.8)
#   dev.off()
#   Ex_raster[[i]] <- raoi
# }
#
# # plot(Ex_raster[[4]])
#
# IBM_aoi <- list(aoi=aoi.MAL, Ex_raster=Ex_raster)
#
# myfile <- paste0(getwd(),"/data/EX_",nrd_name,"_IBM_aoi.qs")
# qsave(IBM_aoi, myfile)
