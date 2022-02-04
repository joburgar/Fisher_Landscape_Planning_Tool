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
                      "Cairo","OpenStreetMap", "ggmap","PNWColors","units","nngeo")

# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
#####################################################################################

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
###--- NO NEED TO RUN THIS NEXT SECTION, INSTEAD JUST IMPORT aoi.Fpop shapefile
# aoi <- st_read(dsn=paste0(getwd(),"/data"), layer="FHE_two_pops")
#
# ggplot()+
#   geom_sf(data=aoi, aes(fill=Hab_zone))
#
# aoi <- aoi %>%
#   summarise(across(geometry, ~ st_union(.))) %>%
#   summarise(across(geometry, ~ st_combine(.)))
#
# ggplot()+
#   geom_sf(data=aoi, fill="black")

###--- Use TSA data to differentiate the two fisher populations
# May need to break the Mackenzie TSA into 2?

# Timber Supply Area
# bcdc_search("Timber Supply Area", res_format = "wms")
# 1: FADM - Timber Supply Area (TSA) (multiple, wms, kml)
# ID: 8daa29da-d7f4-401c-83ae-d962e3a28980
# Name: fadm-timber-supply-area-tsa
# aoi.TSA <- retrieve_geodata_aoi(ID = "8daa29da-d7f4-401c-83ae-d962e3a28980")
# ggplot()+
#   geom_sf(data=aoi.TSA, aes(fill=TSA_NUMBER_DESCRIPTION))
# unique(aoi.TSA$TSA_NUMBER_DESCRIPTION) # 25 unique TSAs
#
# ggplot()+
#   geom_sf(data=aoi.TSA %>% filter(TSA_NUMBER_DESCRIPTION=="MacKenzie TSA"), aes(fill=TSA_NUMBER_DESCRIPTION))
#
# ggplot()+
#   geom_sf(data=aoi.TSA %>% filter(TSA_NUMBER_DESCRIPTION %in% c("Cassiar TSA",
#                                                                 "Fort St. John TSA",
#                                                                 "Fort Nelson TSA",
#                                                                 "Dawson Creek TSA")),aes(fill=TSA_NUMBER_DESCRIPTION))
#
# aoi.TSA$Fpop <- "Central Interior"
# aoi.TSA$Fpop <- as.factor(case_when(aoi.TSA$TSA_NUMBER_DESCRIPTION %in% c("Cassiar TSA",
#                                                                 "Fort St. John TSA",
#                                                                 "Fort Nelson TSA",
#                                                                 "Dawson Creek TSA") ~ "Boreal",TRUE ~ as.character(aoi.TSA$Fpop)))
#
# ggplot()+
#   geom_sf(data=aoi.TSA, aes(fill=Fpop))
#
# aoi.Fpop <- aoi.TSA %>%
#   group_by(Fpop) %>%
#   summarize(geometry = st_union(geometry))
#
# aoi.Fpop
# ggplot()+
#   geom_sf(data=aoi.Fpop, aes(fill=Fpop))

aoi <- aoi.Fpop <- st_read(dsn=paste0(getwd(),"/data"), layer="aoi.Fpop")

ggplot()+
  geom_sf(data=aoi, aes(fill=Fpop))


# Natural Region District
# bcdc_search("Natural Region District", res_format = "wms")
# 1: Natural Resource (NR) Districts (multiple, wms)
# ID: 0bc73892-e41f-41d0-8d8e-828c16139337
# Name: natural-resource-nr-district
aoi.NRD <- retrieve_geodata_aoi(ID = "0bc73892-e41f-41d0-8d8e-828c16139337")
unique(aoi.NRD$DISTRICT_NAME) # 15 unique Natural Resource Districts

# remove the slivers of area, only map polygons with >10km2
ggplot()+
  geom_sf(data=aoi.NRD %>% filter(Area_km2>10), aes(fill=DISTRICT_NAME))

# want to create a fishnet using the aoi extent and have as 30 km2 hexagons, and then join covariates to each cell
# use 30km2 for now but recognize that some data may be at the smaller scale, and will need to be transformed to fit within the female territory

###--- function to create study grid (i.e., female fisher sized pixels)
# note that this function creates a hexagon grid (use square=TRUE for squares)
# this function originally dissolved multipolygon into single polygon,
# assumes original shapefile is the entire aoi,
# for this exercise, commented out the st_union and st_combine code to keep the two Fisher populations distinct

create_study_grid <- function (dsn=dsn, layer=layer, square=FALSE, cellsize=cellsize, output=output){

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

#####################################################################################

###--- Area of Interest (AOI) is larger than Study Area
# keep in mind that WGS84 lat/long espg = 4326; BC Albers espg = 3005; NAD83 / UTM zone 10N espg = 26910

aoi$Fpop <- factor(aoi$Fpop,levels = c("Boreal", "Central Interior"))
aoi_grid <- st_read(dsn="./out", layer="aoi_grid") %>% st_transform(crs = 3005)

Fpop_PGTSA <- ggplot()+
  geom_sf(data=aoi.Fpop,aes(fill=Fpop))+
  scale_fill_manual(values = rev(pnw_palette("Lake",2)))+
  theme(legend.title = element_blank()) +
  geom_sf(data=aoi_grid, col="slategrey")

Cairo(file="out/Fpop_PGTSA.PNG",
      type="png",
      width=3600,
      height=2400,
      pointsize=10,
      bg="white",
      dpi=300)
Fpop_PGTSA
dev.off()

Cairo(file="out/Fpop.PNG",type="png",width=3600,height=2400,pointsize=10,bg="white",dpi=300)
ggplot()+
  geom_sf(data=aoi.Fpop,aes(fill=Fpop))+
  scale_fill_manual(values = rev(pnw_palette("Lake",2)))+
  theme(legend.title = element_blank())
dev.off()
# now reduce aoi to proof-of-concept study area (same as grid above - Prince George TSA)
# female fisher home range = 5.5 * 5.5 = 30 km2
# the CLUS model (Kyle Lochhead and Tyler Muhly) analyses at ha pixel size and sums up to fisher home range
# will want to go even smaller, but for now try the Prince George TSA as the study area (aoi)

# takes too long for such a huge study area
# for needs of a small test model, choose a portion of a TSA and go from there
# using the Prince George TSA - Block E
# create 100 m square pixels - keep in mind that this means 3025 square pixels in a female home range
# start with fisher sized hexagon pixels

list.files("./data")
PGTSA_BLCKE <- create_study_grid(dsn="./data", layer="aoi_PGTSA_BLCKE", cellsize=100, square=TRUE, output="PGTSA_BLCKE") # female fisher hexagon
aoi <- PGTSA_BLCKE[[1]]
aoi_grid <- PGTSA_BLCKE[[2]]
sa_points <- PGTSA_BLCKE[[3]]

# housekeeping - remove large files
rm(PGTSA_BLCKE, aoi.NRD, aoi.Fpop)

###--- LOAD COVARIATES
# Use cached data, downloaded through the bcgeodata warehouse website for big files that are updated annually
# download only needs to happen once per year or as notified of updates
# for this example, where possible, downloaded to aoi (Prince George TSA - Prince George TSA Block E) extent
# https://catalogue.data.gov.bc.ca/

# to see the downloaded gdb and shp files
list.files("./data")

aoi <- aoi %>% st_transform(crs=3005)

# Timber Supply Area
# cached from online as too slow through bcdata package
# downloaded 2022-01-20
aoi.TSA <- retrieve_gdb_shp_aoi(dsn="./data/FADM_TSA", layer="FADM_TSA_polygon")
aoi.TSA %>% group_by(TSNMBRDSCR) %>% summarise(sum(Area_km2)) %>% st_drop_geometry()

# Timber Harvesting Landbase
fgdb = "./data/tsa_thlb.gdb"
# List all feature classes in a file geodatabase
# Read the feature class
st_layers(fgdb)
aoi.THBL <- retrieve_gdb_shp_aoi(dsn=fgdb, layer="tsa_thlb")

# VRI
# https://www2.gov.bc.ca/assets/gov/farming-natural-resources-and-industry/forestry/stewardship/forest-analysis-inventory/data-management/standards/vegcomp_toc_data_dictionaryv5_2019.pdf
# https://www2.gov.bc.ca/assets/gov/farming-natural-resources-and-industry/forestry/stewardship/forest-analysis-inventory/data-management/standards/vegcomp_poly_rank1_data_dictionaryv5_2019.pdf
# TOO BIG! DOWNLOAD THROUGH BC GEODATA AND CACHE FOR LARGER AREAS
# will need to either use smaller aoi or grab this once and save - will take WAY TOO LONG on an as need basis unless going to cell size

# fgdb = "./data/VEG_COMP_POLY_AND_LAYER_2020.gdb"
# List all feature classes in a file geodatabase

# vegetation data (VRI)
# bcdc_search("VRI", res_format = "wms")
aoi.VRI <- retrieve_geodata_aoi(ID = "2ebb35d8-c82f-4a17-9c96-612ac3532d55")

head(aoi.VRI)
aoi.VRI$PROJ_HEIGHT_1_cat <- as.factor(ifelse(aoi.VRI$PROJ_HEIGHT_1 < 10, "H0-10",
                                              ifelse(aoi.VRI$PROJ_HEIGHT_1 < 20, "H10-20",
                                                     ifelse(aoi.VRI$PROJ_HEIGHT_1 < 30, "H20-30",
                                                            ifelse(aoi.VRI$PROJ_HEIGHT_1 < 40, "H30-40",
                                                                   ifelse(aoi.VRI$PROJ_HEIGHT_1 < 50, "H40-50", "H50+"))))))# remove NAs
aoi.VRI <- aoi.VRI[complete.cases(aoi.VRI$PROJ_HEIGHT_1),]


# plot to check - clipped the Anderson release area


Cairo(file="out/PGTSA_BLCKE_VRI.PNG",
      type="png",
      width=2200,
      height=2400,
      pointsize=10,
      bg="white",
      dpi=300)
ggplot()+
  geom_sf(data = aoi , lwd=0.6, col="black", fill="azure3")+
  geom_sf(data = aoi.VRI, aes(fill=PROJ_HEIGHT_1_cat, col=NA)) +
  scale_fill_brewer(palette="Greens") +
  scale_color_brewer(palette="Greens") +
  geom_sf(data = aoi , lwd=0.6, col="black", fill=NA)+
  theme(legend.title=element_blank())
dev.off()

# Biogeoclimatic zones
# bcdc_search("Biogeoclimatic zone", res_format = "wms")
# 3: BEC Map (other, wms, kml)
# ID: f358a53b-ffde-4830-a325-a5a03ff672c3
# Name: bec-map
aoi.BEC <- retrieve_geodata_aoi(ID = "f358a53b-ffde-4830-a325-a5a03ff672c3")
# remove the slivers of area, only map polygons with >10km2
aoi.BEC$Zone_subzone <- paste0(aoi.BEC$ZONE, aoi.BEC$SUBZONE)
ggplot()+
  geom_sf(data=aoi.BEC %>% filter(Area_km2>10), aes(fill=Zone_subzone))

# fisher.bec.zones <- c("CWH","ICH","IDF","BWBS","MS","SBPS","SBS")
# fisher.bec.subzones <- c("BWBSdk","BWBSmw","BWBSwk","SBSwk","ICHwc","SBSmc","SBSmh","SBSmk","SBSmm","SBSmw","SBSdh","SBSdk","SBSdw","SBPSxc","SBPSmc","SBPSdc","SBPSmk","IDFdk","IDFmw","IDFdw","IDFww",
# "MSxc","MSxk","MSdv","MSdm","MSdk","MSdc","ICHmk","ICHmw")
#
# sum(unique(aoi.BEC$Zone_subzone) %in% fisher.bec.subzones, na.rm=TRUE) # 11 of the fisher bec subzones, as per the fisher website, in the area
# sum(unique(aoi.BEC$ZONE) %in% fisher.bec.zones, na.rm=TRUE) # 5 of the 7 BEC zones, as per the fisher website
# not sure if that matters...


# Cutblocks
# bcdc_search("cutblock", res_format = "wms")
# 1: Harvested Areas of BC (Consolidated Cutblocks) (multiple, fgdb, wms, kml, pdf)
# ID: b1b647a6-f271-42e0-9cd0-89ec24bce9f7
# Name: harvested-areas-of-bc-consolidated-cutblocks-
aoi.CTBLK <- retrieve_geodata_aoi(ID = "b1b647a6-f271-42e0-9cd0-89ec24bce9f7")
summary(aoi.CTBLK)
aoi.CTBLK %>% count(VEG_CONSOLIDATED_CUT_BLOCK_ID)

Cairo(file="out/PGTSA_BLCKE_CTBLK.PNG",
      type="png",
      width=2200,
      height=2400,
      pointsize=10,
      bg="white",
      dpi=300)
ggplot()+
  geom_sf(data = aoi , lwd=0.6, col="black", fill="azure3")+
  geom_sf(data=aoi.CTBLK, aes(fill=HARVEST_YEAR))+
  scale_fill_gradient(low="darkgreen", high="white")+
  geom_sf(data = aoi , lwd=0.6, col="black", fill=NA)
dev.off()

save.image("01_load.RData")


###--- workaround to find suitable habitat based on ecoprovince and BEC
# first find the general aoi (will refine)
ecoprov <- ecoprovinces() %>% filter(ECOPROVINCE_CODE %in% c("BOP","SBI","CEI","NBM","SAL"))

# create aoi of dissolved area
# aoi <- ecoprov %>% st_transform(crs = 3005) # ensures poly is in Albers
# sf::st_crs(aoi) = 3005


# plot the ecoprovinces, ordering from CEI to BOP
ecoprov$ECOPROVINCE_NAME <- factor(ecoprov$ECOPROVINCE_NAME,
                                   levels = c("CENTRAL INTERIOR", "SUB-BOREAL INTERIOR", "BOREAL PLAINS"))

ggplot()+
  geom_sf(data = ecoprov, aes(fill=ECOPROVINCE_NAME))+
  geom_sf(data = aoi, fill=NA, color="black", lwd=1)+
  scale_fill_manual(values = pnw_palette("Lake",5))


# Import BEC, Ecoprovinces, Ecoregions and Ecosections
#SBSwk in Peace Region only


aoi.BEC <- bec() %>% st_intersection(aoi) %>% filter(ZONE %in% fisher.bec.zones) %>% filter(grepl(fisher.bec.subzones, MAP_LABEL))
aoi.BEC %>% count(ZONE) %>% st_drop_geometry()

aoi.BEC$Area_km2 <- st_area(aoi.BEC)*1e-6
aoi.BEC <- drop_units(aoi.BEC)
aoi.BEC %>% group_by(ZONE) %>% summarise(sum(Area_km2)) %>% st_drop_geometry()


ggplot()+
  geom_sf(data = aoi, fill=NA, color="black", lwd=1)+
  geom_sf(data = ecoprov, aes(fill=ECOPROVINCE_NAME))+
  geom_sf(data = aoi.BEC)+
  scale_fill_manual(values = pnw_palette("Lake",3))


#- add distance to spatial join
aoi.ecoprov.dist <- st_nn(aoi.BEC, ecoprov, k=1, returnDist = T)
# sum(as.numeric(aoi.ecoprov.dist$dist)) # check to make sure all = 0
aoi.BEC$ecoprov <- unlist(aoi.ecoprov.dist$nn)
aoi.BEC$ecoprov <- ecoprov$ECOPROVINCE_NAME[match(aoi.BEC$ecoprov,rownames(ecoprov))]
summary(aoi.BEC$ecoprov)

aoi.sub <- aoi.BEC %>% filter(str_detect(MAP_LABEL, fisher.bec.subzones))
aoi.sub$ecoprov <- as.character(aoi.sub$ecoprov)

aoi.sub$habitat <- ifelse(aoi.sub$ecoprov=="BOREAL PLAINS" & grepl("BWBSdk|BWBSmw|BWBSwk|SBSwk|ICHwc", aoi.sub$MAP_LABEL), 1,
                          ifelse(aoi.sub$ecoprov=="SUB-BOREAL INTERIOR" & grepl("CWH|ICH|IDF", aoi.sub$MAP_LABEL),1,
                                 ifelse(aoi.sub$ecoprov=="SUB-BOREAL INTERIOR" & grepl("SBSwk|SBSmc|SBSmh|SBSmk|SBSmm|SBSmw|SBSdh|SBSdk|SBSdw", aoi.sub$MAP_LABEL),1,
                                        ifelse(aoi.sub$ecoprov=="CENTRAL INTERIOR" & grepl("SBPSxc|SBPSmc|SBPSdc|SBPSmk|IDFdk|IDFmw|IDFdw|IDFww|MSxc|MSxk|MSdv|MSdm|MSdk|MSdc|ICHmk|ICHmw|ICHmk|SBSdw|SBSmc", aoi.sub$MAP_LABEL),
                                               1,0))))


aoi.fisher.habitat <- aoi.sub %>% filter(habitat==1) %>%
  summarise(across(geometry, ~ st_combine(.))) %>%
  summarise(across(geometry, ~ st_union(.)))


aoi2 <- aoi.BEC %>% st_intersection(aoi.fisher.habitat)

ggplot()+
  geom_sf(data = aoi, color="blue")+
  geom_sf(data=aoi2, color="red")

# not really what I'm wanting....probably best to just see if I can get the fisher distribution shapefile from someone else
# aoi.fisher.habitat <- aoi.sub %>% filter(habitat==1) %>% st_combine() %>% st_union(by_feature=FALSE, is_coverage=TRUE)
# aoi.fisher.habitat <- aoi.fisher.habitat %>% st_zm(drop=TRUE, what="ZM")
#
# ggplot()+
#   # geom_sf(data = aoi)+
#   geom_sf(data=aoi.fisher.habitat)
#
# glimpse(aoi.fisher.habitat)
# st_write(aoi.fisher.habitat, dsn=paste0(getwd(),"/out/aoi.fisher.habitat.kml"), delete_layer=TRUE)

rm(aoi)
save.image("01_load.RData")
# load("01_load.RData")
