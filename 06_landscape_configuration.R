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

# retrieve_geodata_aoi <- function (ID=ID){
#   aoi.geodata <- bcdc_query_geodata(ID) %>%
#     filter(BBOX(st_bbox(aoi))) %>%
#     collect()
#   aoi.geodata <- aoi.geodata %>% st_intersection(aoi)
#   aoi.geodata$Area_km2 <- st_area(aoi.geodata)*1e-6
#   aoi.geodata <- drop_units(aoi.geodata)
#   aoi.geodata <- droplevels.sfc(aoi.geodata, except = geometry)
#   return(aoi.geodata)
# }

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

# retrieve_raster_aoi <- function(rfile=rfile, raoi=raoi){
#   # rfile=FETA.rasters[3]
#   rfile <- raster(paste0("data/FETA_fromKyle/",rfile))
#   rfile <- projectRaster(rfile,crs = 26910)
#   # Crop FETA rasters by extent of raoi
#   rfile.aoi <- crop(rfile, extent(raoi))
#   # Sum FETA rasters for value of FETA within fisher territory
#   # recall FETA rasters are 100x100 m while aoi is 5500x5500
#   # need to sum by 55 cells in x and y direction to make fisher equivalent territory
#   # ra <- aggregate(rfile.aoi, fact=c(55,55), fun=sum, na.rm=TRUE)
#   return(ra)
# }

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

unique(Female_HRs_compiled$Hab_zone)

### Now create Fisher Habitat Extension Zone specific home range raster stacks
# create function and loop it over Hab_zone

HR_rasterstack_function <- function(HRsfobj=HRsfobj, FHEzone=FHEzone){
  # FHEzone = "Boreal"
  FHR <- HRsfobj %>% filter(Hab_zone==FHEzone)
  FHR <- FHR %>% st_transform(crs=26910)
  FHR <- FHR %>% arrange(SA_Fisher_ID)
  FHR$HR <- 1

  # home ranges overlap each other - need to create rasters for each home range and stack them
  rHR <- raster(ext=extent(FHR), crs=26910, res=c(100,100))

  rHR_list=list()
  for(i in 1:nrow(FHR)){
    rHR_list[[i]] <- rasterize(FHR[i,], rHR, field="HR")
  }

  rHR_list
  rHR_stack = stack(rHR_list)
  names(rHR_stack) <- FHR$SA_Fisher_ID

  return(rHR_stack)
}

Hab_zone <- unique(Female_HRs_compiled$Hab_zone)
Hab_HR_rasters <- list()
for(i in 1:length(Hab_zone)){
  Hab_HR_rasters[[i]] <- HR_rasterstack_function(HRsfobj=Female_HRs_compiled, FHEzone = Hab_zone[i])
}

plot(Hab_HR_rasters[[1]])
plot(Hab_HR_rasters[[2]])
plot(Hab_HR_rasters[[3]])
plot(Hab_HR_rasters[[4]])

#####################################################################################
###--- CREATE AGGREGATED RASTER STACK OF FETA VALUES
#####################################################################################
# currently just habitat 1 or 0 if it meets fisher criteria, should be a 0-1 value?
###--- grab data from Rich's gdb

# GISDir <- "//spatialfiles.bcgov/work/srm/vic/tib/Workarea/RWeir/Fisher Habitat Extension/"
# list.files(GISDir)
# Boreal_models_210215.gdb
# Dry_forest_models_210215.gdb
# Sub_Boreal_dry_models_210215.gdb
# Sub_Boreal_moist_models_210215.gdb

# fgdb = "Dry_forest_models_210215.gdb"
# st_layers(paste(GISDir,fgdb, sep="/"))

# massive files, R is not great at geo-processing with vector data
# work with Rich to have completed in ArcGIS and export rasters for easy upload

# DRY_Denning <- st_read(dsn=paste(GISDir,fgdb, sep="/"),layer="Denning_primary") %>% st_transform(crs=26910)
# DRY_Denning$Primary <- 1
# rDRY_Denning <- raster(ext=extent(DRY_Denning), crs=26910, res=c(100,100))
# rDRY_Denning <- rasterize(DRY_Denning, rDRY_Denning, field="Primary", background=0)
# plot(rDRY_Denning)

# Vector of the 2020 would be in same folder, but different gdb for each FHE zone:
# [FHE_zone]_models_210215.gdb (e.g., Sub_Boreal_moist_models_210215.gdb would have the Denning, Branch_resting, Cavity_resting, CWD_resting primary stands).

# List all feature classes in a file geodatabase
#####################################################################################
###--- Go back to using Kyle's FETA until 2021 VRI rasters created

FETA.rasters <- list.files("data/FETA_fromKyle")
# bring in fisher habitat data - necessary for Mahalanobis
FETA.rasters <- FETA.rasters[grepl("movement|rest|denning", FETA.rasters)]

Fraster_list=list()
for(i in 1:length(FETA.rasters)){
  Fraster_list[[i]] <- raster(paste0("data/FETA_fromKyle/",FETA.rasters[i]))
}
Fraster_stack = stack(Fraster_list)
plot(Fraster_stack)
extent(Fraster_stack)
projection(Fraster_stack)

# sum(Fraster_list[[3]]$rest_cavity@data@values)

### Now create Fisher Habitat Extension Zone specific FETA raster stacks
# convert home range raster extent to same as FETA rasters to save time on the crop
# create function and loop it over Hab_zone_HR_raster hab zones


FHEzone_rasterstack_function <- function(HRsfobj=HRsfobj, FHEzone=FHEzone, Hab_raster=Fraster_stack){

  FHR <- HRsfobj %>% filter(Hab_zone==FHEzone)
  FHR <- FHR %>% st_transform(crs=projection(Hab_raster))

  Hab_raster_cropped <- crop(Hab_raster, FHR, snap='out')

  # Fraster_HR[is.na(Fraster_HR[])] <- 0
  return(Hab_raster_cropped)
}

Hab_rasters <- list()
for(i in 1:length(Hab_zone)){
  Hab_rasters[[i]] <- FHEzone_rasterstack_function(HRsfobj=Female_HRs_compiled,
                                               FHEzone = Hab_zone[i],
                                               Hab_raster=Fraster_stack)
}

plot(Hab_rasters[[1]])
plot(Hab_rasters[[2]])
plot(Hab_rasters[[3]])
plot(Hab_rasters[[4]])

# curiosity - how correlated are the underlying habitat layers?
# jnk=layerStats(Fraster_BHR, 'pearson', na.rm=T)
# corr_matrix=jnk$'pearson correlation coefficient'

################################################################################

################################################################################
###--- First bring in Rich's vector D2 data
MAHALDir <- "//Sfp.idir.bcgov/s140/S40203/Ecosystems/Conservation Science/Species Conservation Science/Mesocarnivores/Projects/Fisher Critical Habitat/Analysis/FHE_hexes/"
list.files(MAHALDir)

# Boreal_FHE_home_ranges.Rda
# Dry_forest_FHE_home_ranges.Rda
# Sub_Boreal_dry_FHE_home_ranges.Rda
# Sub_Boreal_moist_FHE_home_ranges.Rda

load(paste0(MAHALDir,'Boreal_FHE_home_ranges.Rda'))
Boreal_FHE_home_ranges

load(paste0(MAHALDir,'Dry_forest_FHE_home_ranges.Rda'))
Dry_forest_FHE_home_ranges

load(paste0(MAHALDir,'Sub_Boreal_dry_FHE_home_ranges.Rda'))
Sub_Boreal_dry_FHE_home_ranges

load(paste0(MAHALDir,'Sub_Boreal_moist_FHE_home_ranges.Rda'))
Sub_Boreal_moist_FHE_home_ranges

###################################################################################
###--- Try Mahalanobis for an entire region, pixel by pixel
plot(Fraster_HR) # no cavity_resting habitat (highly correlated so not necessary)

###--- NEED TO FIX THIS FUNCTION UP, ALIGN EXTENTS AND THEN RUN MAHALANOBIS FHEZones, THEN CLIP TO HR RASTERS

###--- Mahalanobis distance - how to make this work for rasters
options(scipen = 10)

add_Mahal_raster_layer_function <- function(Hab_raster=Hab_rasters[[1]],HR_raster=Hab_HR_rasters[[1]]){

  tmp.mat <- matrix(NA,nrow=dim(FHEzone_rasters)[1]*dim(FHEzone_rasters)[2],ncol=dim(FHEzone_rasters)[3])
  dim(tmp.mat)
  for(i in 1:dim(tmp.mat)[2]){
    tmp.mat[,i] <- values(FHEzone_rasters[[i]])
  }

  # remove the non HR habitat but subsetting to non-NA values
  # tmp.hab <- tmp.mat[complete.cases(tmp.mat),]
  tmp.hab <- tmp.mat
  # tmp.hab[is.na(tmp.hab)]<-0
  colSums(tmp.hab)

  Sx <- cov(tmp.hab)
  D2 <- mahalanobis(tmp.hab, colMeans(tmp.hab),Sx)

  # summary(D2)
  # hist(D2)
  # plot(density(D2, bw = 0.5),
  #      main="Squared Mahalanobis distances, n=100, p=3") ; rug(D2)
  # qqplot(qchisq(ppoints(100), df = 3), D2,
  #        main = expression("Q-Q plot of Mahalanobis" * ~D^2 *
  #                            " vs. quantiles of" * ~ chi[3]^2))
  # abline(0, 1, col = 'gray')

  D2raster <- raster(ext=extent(FHEzone_rasters), crs=26910, res=c(100,100))
  values(D2raster) <- D2

  # add D2 raster layer to raster stack and add name to layer
  FHEzone_rasters <- stack(FHEzone_rasters,D2raster)
  names(FHEzone_rasters)[5] <- "D2"

  # want to extract the denning raster values for the home range
  # x = FETA raster brick (can do multiple layers at a time)
  # home range raster is the
  values(FHEzone_rasters[[1]])

  # extract habitat data for each home range
  rHR_hab_stack <- list()
  for(i in 1:dim(FHEzone_HR_rasters)[3]){
    rHR_hab_stack[[i]] <- raster::mask(x=FHEzone_rasters, mask=FHEzone_HR_rasters[[i]], maskvalue=0)
  }

  # print area for each home range (at 1 ha pixel scale)
  for(i in 1:dim(rBHR_stack)[3]){
    print(sum(rBHR_stack[[i]]@data@values))
  }

  # quick check = seems to match up with polygons, considering now pixelated
  KPF_FHR %>% dplyr::select(Area_km2) %>% st_drop_geometry()

  # print area for each home range (at 1 ha pixel scale)
  # nrow = number of HR polygons (i.e., HR raster stack)
  # ncol = number of habitat rasters + 1 for total habitat
  HR_hab_df <- as.data.frame(matrix(NA,nrow=dim(rBHR_stack)[3],ncol=dim(Fraster_BHR)[3]+1))
  colnames(HR_hab_df) <- c("Total_area","Denning","Movement","Rest_cwd","Rest_rust","D2_sum")
  for(i in 1:length(rBHR_hab_stack)){
    HR_hab_df$Total_area[i] <- sum(rBHR_stack[[i]]@data@values, na.rm=T)
    HR_hab_df$Denning[i] <- sum(rBHR_hab_stack[[i]]$denning@data@values, na.rm=T)
    HR_hab_df$Movement[i] <- sum(rBHR_hab_stack[[i]]$movement@data@values, na.rm=T)
    HR_hab_df$Rest_cwd[i] <- sum(rBHR_hab_stack[[i]]$rest_cwd@data@values, na.rm=T)
    HR_hab_df$Rest_rust[i] <- sum(rBHR_hab_stack[[i]]$rest_rust@data@values, na.rm=T)
    HR_hab_df$D2_sum[i] <- sum(rBHR_hab_stack[[i]]$D2@data@values, na.rm=T)
  }

  HR_hab_df$Prop_Denning <- HR_hab_df$Denning / HR_hab_df$Total_area
  HR_hab_df$Prop_Movement <- HR_hab_df$Movement / HR_hab_df$Total_area
  HR_hab_df$Prop_Rest_cwd <- HR_hab_df$Rest_cwd / HR_hab_df$Total_area
  HR_hab_df$Prop_Rest_rust <- HR_hab_df$Rest_rust / HR_hab_df$Total_area
  HR_hab_df$HR_D2 <- HR_hab_df$D2_sum / HR_hab_df$Total_area
  HR_hab_df$SA_Fisher_ID <- names(rBHR_stack)

  HR_hab_df <- cbind(HR_hab_df, Boreal_FHE_home_ranges[5:17,])
  glimpse(HR_hab_df)
  cor(HR_hab_df$HR_D2, HR_hab_df$D2)

  write.csv(HR_hab_df, "out/KPF_HR_hab_df.csv")
  plot(HR_hab_df$HR_D2 ~HR_hab_df$Total_area)
  plot(HR_hab_df$D2 ~HR_hab_df$Total_area)

  # plot the D2 values for each home range
  D2_stack <- stack(rBHR_hab_stack[[1]]$D2)
  for(i in 2:length(rBHR_hab_stack)){
    D2_stack <- stack(D2_stack,rBHR_hab_stack[[i]]$D2)
  }
  names(D2_stack) <- names(rBHR_stack)
  plot(D2_stack)
}
