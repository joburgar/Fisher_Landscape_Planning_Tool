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
# 06_landscape_configuration.R
# script to create example data to run landscape configuration study
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 23-June-2022
#####################################################################################
version$major
version$minor
R_version <- paste0("R-",version$major,".",version$minor)

.libPaths(paste0("C:/Program Files/R/",R_version,"/library")) # to ensure reading/writing libraries from C drive
tz = Sys.timezone() # specify timezone in BC

# Load Packages
list.of.packages <- c("tidyverse", "lubridate","bcdata", "bcmaps","sf", "rgdal","nngeo","raster",
                      "viridis","Cairo","PNWColors","units")

# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

#####################################################################################
# keep in mind that WGS84 lat/long espg = 4326; BC Albers espg = 3005; NAD83 / UTM zone 10N espg = 26910
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
  # ra <- aggregate(rfile.aoi, fact=c(55,55), fun=sum, na.rm=TRUE)
  return(ra)
}

###--- uncertainty function
se <- function(x) sqrt(var(x)/length(x))


#####################################################################################
###--- EXPLORE FOR LANDSCAPE CONFIGURATION EXPERIMENTAL DESIGN
#####################################################################################

Fpop <- st_read(dsn=paste0(getwd(),"/data"), layer="aoi.Fpop_simpl") %>% st_transform(crs=3005)
FHEzones <- st_read(dsn=paste0(getwd(),"/data"), layer="FHE_zones")
FHEzones <- st_simplify(FHEzones %>% st_transform(crs=26910), preserveTopology = FALSE, dTolerance = 1000) %>% st_transform(crs=3005)

ggplot()+
  geom_sf(data=FHEzones, aes(fill=Hab_zone))+
  scale_fill_viridis_d(option="D", alpha=0.7) # alpha does the transparency


aoi <- Fpop %>%
  summarise(across(geometry, ~ st_union(.))) %>%
  summarise(across(geometry, ~ st_combine(.)))

#####################################################################################
###--- BRING IN FEMALE FISHER HOME RANGES
#####################################################################################

GISDir <- "//Sfp.idir.bcgov/s140/S40203/Ecosystems/Conservation Science/Species Conservation Science/Mesocarnivores"
fgdb = "FisherStudyDatabases-do not edit/Female_HRs_compiled.gdb"

# List all feature classes in a file geodatabase
# st_layers(paste(GISDir,fgdb, sep="/"))

Female_HRs_compiled <- retrieve_gdb_shp_aoi(dsn=paste(GISDir,fgdb, sep="/"),layer="Female_HRs_compiled")
# Think the other two layers in the gdb are also in the compiled layer but contain more info (accidental duplication?)
# WFI_F04_HR <- retrieve_gdb_shp_aoi(dsn=paste(GISDir,fgdb, sep="/"),layer="WFI_F04_HR")
# WFI_F02 <- retrieve_gdb_shp_aoi(dsn=paste(GISDir,fgdb, sep="/"),layer="WFI_F02")

ggplot()+
  geom_sf(data=FHEzones, aes(fill=Hab_zone))+
  scale_fill_viridis_d(option="D", alpha=0.7)+ # alpha does the transparency
  geom_sf(data=Female_HRs_compiled, lwd=2, col="black")

Female_HRs_compiled <- st_intersection(Female_HRs_compiled, Fpop)
Female_HRs_compiled <- st_join(Female_HRs_compiled, FHEzones %>% dplyr::select(Hab_zone), left=TRUE, largest=TRUE)

ggplot()+
  geom_sf(data=Female_HRs_compiled, aes(fill=SA_Fisher_ID))+
  scale_fill_viridis_d(option="D", alpha=0.7) # alpha does the transparency

Female_HRs_compiled %>% group_by(Fpop, Hab_zone) %>% summarise_at(vars(Area_km2), list(Min=min, Mean=mean, Max=max, SE=se)) %>% st_drop_geometry()
Female_HRs_compiled %>% count(Fpop, Hab_zone) %>% st_drop_geometry()

ggplot()+
  geom_sf(data=Female_HRs_compiled, aes(fill=Hab_zone))+
  scale_fill_viridis_d(option="D", alpha=0.7) # alpha does the transparency
# For proof of concept, try for Boreal first

Female_HRs_compiled <- Female_HRs_compiled %>% st_transform(crs=26910)
Female_HRs_compiled$HR <- 1

Female_HRs_compiled$SA_Fisher_ID

KPF_FHR <- Female_HRs_compiled %>% filter(grepl("KPF",SA_Fisher_ID))
KPF_FHR <- KPF_FHR %>% arrange(SA_Fisher_ID)

# home ranges overlap each other - need to create rasters for each home range and stack them
rBHR <- raster(ext=extent(KPF_FHR), crs=26910, res=c(100,100))

rBHR_list=list()
for(i in 1:nrow(KPF_FHR)){
  rBHR_list[[i]] <- rasterize(KPF_FHR[i,], rBHR, field="HR", background=0) # interim work around until terra and new raster package uploaded
}
rBHR_list
rBHR_stack = stack(rBHR_list)

names(rBHR_stack) <- KPF_FHR$SA_Fisher_ID
rBHR_stack

plot(rBHR_stack)
plot(rBHR_stack,2)
# quick visual check that individual rasters look the same as vector file
ggplot()+
  geom_sf(data=KPF_FHR, aes(fill=SA_Fisher_ID))+
  scale_fill_viridis_d(option="D", alpha=0.7) # alpha does the transparency

#####################################################################################
###--- CREATE AGGREGATED RASTER STACK OF FETA VALUES
#####################################################################################
# Created, and no need to re-run, just pull in appropriate raster files

FETA.rasters <- list.files("data/FETA_fromKyle")
# for now, just want to bring in fisher habitat requirement data
FETA.rasters <- FETA.rasters[grepl("movement|rest|denning", FETA.rasters)]

Fraster_list=list()
for(i in 1:length(FETA.rasters)){
  Fraster_list[[i]] <- raster(paste0("data/FETA_fromKyle/",FETA.rasters[i]))
}
Fraster_stack = stack(Fraster_list[[1]],Fraster_list[[2]],Fraster_list[[4]],Fraster_list[[5]])

# convert home range raster extent to same as FETA rasters to save time on the crop
FETAproj <- "+proj=aea +lat_0=45 +lon_0=-126 +lat_1=50 +lat_2=58.5 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs"
tmpprj <- projectExtent(rBHR, crs=FETAproj)
res(tmpprj) <- 100

rBHR_FETAproj <- projectRaster(rBHR,tmpprj)
# crop FETA raster stack to same extent as home range rasters
Fraster_tmp <- crop(Fraster_stack, extent(rBHR_FETAproj))
# reproject back into UTM, same as home range rasters
Fraster_BHR <- projectRaster(Fraster_tmp, crs=26910, res=100)

extent(Fraster_BHR) <- alignExtent(Fraster_BHR, rBHR_stack)
Fraster_BHR <- crop(Fraster_BHR, extent(rBHR_stack),snap='out')

#plot to check it worked
plot(Fraster_BHR) # no cavity_resting habitat (highly correlated so not necessary)
Fraster_BHR[is.na(Fraster_BHR[])] <- 0

plot(Fraster_BHR[[2]])
plot(rBHR_stack[[2]])

# jnk=layerStats(Fraster_BHR, 'pearson', na.rm=T)
# corr_matrix=jnk$'pearson correlation coefficient'

# want to extract the denning raster values for the home range
# x = FETA raster brick (can do multiple layers at a time)
# home range raster is the
values(Fraster_BHR[[1]])

# extract habitat data for each home range
rBHR_hab_stack <- list()
for(i in 1:dim(rBHR_stack)[3]){
  rBHR_hab_stack[[i]] <- raster::mask(x=Fraster_BHR, mask=rBHR_stack[[i]], maskvalue=0)
}

# print area for each home range (at 1 ha pixel scale)
for(i in 1:dim(rBHR_stack)[3]){
  print(sum(rBHR_stack[[i]]@data@values))
}

# quick check = seems to match up with polygons, considering now pixelated
KPF_FHR %>% dplyr::select(Area_km2) %>% st_drop_geometry()

# print area for each home range (at 1 ha pixel scale)
for(i in 1:length(rBHR_hab_stack)){
  print(sum(rBHR_hab_stack[[i]]$denning@data@values, na.rm=T))
}
