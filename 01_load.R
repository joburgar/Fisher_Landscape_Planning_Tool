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

#####################################################################################
aoi <- st_read(dsn=paste0(getwd(),"/data"), layer="FHE_two_pops")

ggplot()+
  geom_sf(data=aoi, aes(fill=Hab_zone))

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
aoi.TSA <- retrieve_geodata_aoi(ID = "8daa29da-d7f4-401c-83ae-d962e3a28980")
ggplot()+
  geom_sf(data=aoi.TSA, aes(fill=TSA_NUMBER_DESCRIPTION))
unique(aoi.TSA$TSA_NUMBER_DESCRIPTION) # 25 unique TSAs

ggplot()+
  geom_sf(data=aoi.TSA %>% filter(TSA_NUMBER_DESCRIPTION=="MacKenzie TSA"), aes(fill=TSA_NUMBER_DESCRIPTION))

ggplot()+
  geom_sf(data=aoi.TSA %>% filter(TSA_NUMBER_DESCRIPTION %in% c("Cassiar TSA",
                                                                "Fort St. John TSA",
                                                                "Fort Nelson TSA",
                                                                "Dawson Creek TSA")),aes(fill=TSA_NUMBER_DESCRIPTION))

aoi.TSA$Fpop <- "Central Interior"
aoi.TSA$Fpop <- as.factor(case_when(aoi.TSA$TSA_NUMBER_DESCRIPTION %in% c("Cassiar TSA",
                                                                "Fort St. John TSA",
                                                                "Fort Nelson TSA",
                                                                "Dawson Creek TSA") ~ "Boreal",TRUE ~ as.character(aoi.TSA$Fpop)))

ggplot()+
  geom_sf(data=aoi.TSA, aes(fill=Fpop))

aoi.Fpop <- aoi.TSA %>%
  group_by(Fpop) %>%
  summarize(geometry = st_union(geometry))

aoi.Fpop
ggplot()+
  geom_sf(data=aoi.Fpop, aes(fill=Fpop))

st_write(aoi.Fpop, paste0(getwd(),"/data/aoi.Fpop.shp"), delete_layer = TRUE)
st_write(aoi.TSA, paste0(getwd(),"/data/aoi.TSA.shp"), delete_layer = TRUE)

# want to create a fishnet using the aoi extent and have as 30 km2 hexagons, and then join covariates to each cell
# use 30km2 for now but recognize that some data may be at the smaller scale, and will need to be transformed to fit within the female territory

###--- function to create study grid (i.e., female fisher sized pixels)
# note that this function creates a hexagon grid (use square=TRUE for squares)
# this function originally dissolved multipolygon into single polygon,
# assumes original shapefile is the entire aoi,
# for this exercise, commented out the st_union and st_combine code to keep the two Fisher populations distinct

create_study_grid <- function (dsn=dsn, layer=layer, Fpop="Boreal", cellsize=cellsize, output=output){

  aoi <- st_read(dsn=dsn, layer=layer) %>% st_transform(crs = 3005) # ensures poly is in Albers
  aoi <- aoi %>% filter(Fpop==Fpop)
    # summarise(across(geometry, ~ st_union(.))) %>% # want to keep as distinct Fisher populations (already did this to create aoi)
    # summarise(across(geometry, ~ st_combine(.)))
  aoi_utm <- st_transform(aoi, crs=26910) # to have in metres for specifying grid cell size
  aoi_grid <- sa_grid <- st_make_grid(st_bbox(aoi_utm), cellsize=cellsize, square=FALSE) #  grid for entire AOI (rectangle)

  sa_points <- st_point_on_surface(sa_grid)  # if using portion of aoi
  sa_points <- st_intersection(sa_points, aoi_utm)

  st_write(aoi_utm, paste0(getwd(),"/out/",output,"_utm.shp"), delete_layer = TRUE)
  st_write(aoi_grid %>% st_intersection(aoi_utm), paste0(getwd(),"/out/",output,"_point.shp"), delete_layer = TRUE)
  st_write(sa_points, paste0(getwd(),"/out/",output,"_grid.shp"), delete_layer = TRUE)

  return(list(aoi_utm, aoi_grid, sa_points))
}

# takes too long for such a huge study area
# FFH <- create_study_grid(dsn= "./data", layer = "aoi.Fpop", Fpop="Boreal", cellsize = 5500, output = "FFH") # female fisher hexagon
# not run
# aoi <- FFH[[1]]
# aoi_utm <- FFH[[2]]
# sa_points <- FFH[[3]]
# for needs of a small test model, choose a few TSAs and go from there


#####################################################################################

###--- Area of Interest (AOI) is larger than Study Area
# keep in mind that WGS84 lat/long espg = 4326; BC Albers espg = 3005; NAD83 / UTM zone 10N espg = 26910

# cellsize <- 5500 # female fisher home range = 5.5 * 5.5 = 30 km2
# unique(aoi.TSA$TSA_NUMBER_DESCRIPTION)
#
# aoi <- aoi.TSA %>% filter(TSA_NUMBER_DESCRIPTION %in% c("Prince George TSA")) %>% st_transform(crs = 3005) # ensures poly is in Albers
# aoi <- aoi %>%
#   summarise(across(geometry, ~ st_union(.))) %>% # to create dissolved aoi
#   summarise(across(geometry, ~ st_combine(.)))
# aoi_utm <- st_transform(aoi, crs=26910) # to have in metres for specifying grid cell size
# aoi_grid <- sa_grid <- st_make_grid(st_bbox(aoi_utm), cellsize=5500, square=FALSE) #  grid for entire AOI (rectangle)
# aoi_grid <- st_intersection(aoi_grid, aoi_utm)
#
# st_write(aoi_utm, paste0(getwd(),"/out/aoi_utm.shp"), delete_layer = TRUE)
# st_write(aoi_grid, paste0(getwd(),"/out/aoi_grid.shp"), delete_layer = TRUE)


###--- will want user specified options for downloading layers, using the bcdata search function and then inputting the ID
# but also will have some pre-loaded options for the ones most likely to be used
# is it worth downloading on a nightly/weekly/monthly basis to save time when the user is interacting with the app? think so, need to check on good update interval


aoi_grid <- st_read(dsn="./out", layer="aoi_grid") %>% st_transform(crs = 3005)
aoi <- st_read(dsn="./out", layer="aoi_utm") %>% st_transform(crs = 3005)
aoi

ggplot()+
  geom_sf(data=aoi, lwd=1, fill=NA)+
  geom_sf(data=aoi_grid)

# Natural Region District
# bcdc_search("Natural Region District", res_format = "wms")
# 1: Natural Resource (NR) Districts (multiple, wms)
# ID: 0bc73892-e41f-41d0-8d8e-828c16139337
# Name: natural-resource-nr-district
aoi.NRD <- retrieve_geodata_aoi(ID = "0bc73892-e41f-41d0-8d8e-828c16139337")
unique(aoi.NRD$DISTRICT_NAME) # 15 unique Natural Resource Districts
aoi.NRD %>% group_by(ORG_UNIT) %>% summarise(sum=Area_km2) %>% st_drop_geometry()

# remove the slivers of area, only map polygons with >10km2
ggplot()+
  geom_sf(data=aoi.NRD %>% filter(Area_km2>10), aes(fill=DISTRICT_NAME))

# Timber Supply Area
# bcdc_search("Timber Supply Area", res_format = "wms")
# 1: FADM - Timber Supply Area (TSA) (multiple, wms, kml)
# ID: 8daa29da-d7f4-401c-83ae-d962e3a28980
# Name: fadm-timber-supply-area-tsa
aoi.TSA <- retrieve_geodata_aoi(ID = "8daa29da-d7f4-401c-83ae-d962e3a28980")
ggplot()+
  geom_sf(data=aoi.TSA %>% filter(Area_km2>10), aes(fill=TSA_NUMBER_DESCRIPTION))
# slivers of neighbouring TSAs but really only the one (should be as aoi came from one TSA)

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

fisher.bec.zones <- c("CWH","ICH","IDF","BWBS","MS","SBPS","SBS")
fisher.bec.subzones <- c("BWBSdk","BWBSmw","BWBSwk","SBSwk","ICHwc","SBSmc","SBSmh","SBSmk","SBSmm","SBSmw","SBSdh","SBSdk","SBSdw","SBPSxc","SBPSmc","SBPSdc","SBPSmk","IDFdk","IDFmw","IDFdw","IDFww",
"MSxc","MSxk","MSdv","MSdm","MSdk","MSdc","ICHmk","ICHmw")

sum(unique(aoi.BEC$Zone_subzone) %in% fisher.bec.subzones, na.rm=TRUE) # 11 of the fisher bec subzones, as per the fisher website, in the area
sum(unique(aoi.BEC$ZONE) %in% fisher.bec.zones, na.rm=TRUE) # 5 of the 7 BEC zones, as per the fisher website
# not sure if that matters...

# VRI
# bcdc_search("vri", res_format = "wms")
# 1: VRI - 2020 - Forest Vegetation Composite Rank 1 Layer (R1) (fgdb, multiple, wms, kml)
# ID: 2ebb35d8-c82f-4a17-9c96-612ac3532d55
# Name: vri-2020-forest-vegetation-composite-rank-1-layer-r1-
aoi.VRI <- retrieve_geodata_aoi(ID = "2ebb35d8-c82f-4a17-9c96-612ac3532d55")
# will need to either use smaller aoi or grab this once and save - will take WAY TOO LONG on an as need basis unless going to cell size

# Cutblocks
# bcdc_search("cutblock", res_format = "wms")
# 1: Harvested Areas of BC (Consolidated Cutblocks) (multiple, fgdb, wms, kml, pdf)
# ID: b1b647a6-f271-42e0-9cd0-89ec24bce9f7
# Name: harvested-areas-of-bc-consolidated-cutblocks-
aoi.CTBLK <- retrieve_geodata_aoi(ID = "b1b647a6-f271-42e0-9cd0-89ec24bce9f7")
# will need to grab this once and save - will take awhile to process (142959 records and 143 pages) on an as need basis unless going to cell size

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
