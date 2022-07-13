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

Female_HRs_compiled %>% filter(!grepl("SiteC",SA_Fisher_ID)) %>% group_by(Fpop, Hab_zone) %>% summarise_at(vars(Area_km2), list(Min=min, Mean=mean, Max=max, SE=se)) %>% st_drop_geometry()

Female_HRs_compiled %>% count(Fpop, Hab_zone) %>% st_drop_geometry()
Female_HRs_compiled %>% dplyr::select(Hab_zone, SA_Fisher_ID) %>% st_drop_geometry


Female_HRs_compiled %>% filter(Hab_zone=="Boreal") %>% filter(grepl("SiteC",SA_Fisher_ID)) %>% st_drop_geometry()


Female_HRs_compiled %>% filter(Hab_zone=="Boreal") %>% count(SA_Fisher_ID) %>% st_drop_geometry()
Female_HRs_compiled <- Female_HRs_compiled %>% mutate(Telem_years = case_when(Hab_zone=="Boreal" ~ "2005-2009",
                                                                              Hab_zone=="Sub-Boreal moist" ~ "1996-2000",
                                                                              Hab_zone=="Sub-Boreal dry" ~ "1990-1993",
                                                                              Hab_zone=="Dry Forest" ~ "2005-2009"))
# to clip VRI, export HR shapefiles by Hab_zone
# st_write(Female_HRs_compiled %>% filter(Hab_zone=="Boreal"), paste0(getwd(),"/data/Female_HRs_Boreal.shp"), delete_layer=TRUE)
# st_write(Female_HRs_compiled %>% filter(Hab_zone=="Sub-Boreal moist"), paste0(getwd(),"/data/Female_HRs_SBM.shp"), delete_layer=TRUE)
# st_write(Female_HRs_compiled %>% filter(Hab_zone=="Sub-Boreal dry"), paste0(getwd(),"/data/Female_HRs_SBD.shp"), delete_layer=TRUE)
# st_write(Female_HRs_compiled %>% filter(Hab_zone=="Dry Forest"), paste0(getwd(),"/data/Female_HRs_Dry.shp"), delete_layer=TRUE)

Hab_zone <- unique(Female_HRs_compiled$Hab_zone)


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

  Hab_raster_cropped <- crop(Hab_raster, extent(FHR)+2000, snap='out') # add 1 km buffer area around extent)

  return(Hab_raster_cropped)
}

# need to check if normally distributed and then log if not - or just check in with Rich on

Hab_rasters <- list()
for(i in 1:length(Hab_zone)){
  Hab_rasters[[i]] <- FHEzone_rasterstack_function(HRsfobj=Female_HRs_compiled,
                                               FHEzone = Hab_zone[i],
                                               Hab_raster=Fraster_stack)
}

plot(Hab_rasters[[1]]) # no rest_cavity
plot(Hab_rasters[[2]])
plot(Hab_rasters[[3]])
plot(Hab_rasters[[4]]) # no rest_cavity

# curiosity - how correlated are the underlying habitat layers?
# jnk=layerStats(Fraster_BHR, 'pearson', na.rm=T)
# corr_matrix=jnk$'pearson correlation coefficient'


###################################################################################
###--- Mahalanobis for an entire region, pixel by pixel

add_Mahal_raster_layer_function <- function(Hab_raster=Hab_raster){

  # Hab_raster=Hab_rasters[[1]]
  # extract values from each raster and put into matrix (formatting for Mahalanobis)
  tmp.mat <- matrix(NA,nrow=ncell(Hab_raster),ncol=nlayers(Hab_raster))

  for(i in 1:ncol(tmp.mat)){
    tmp.mat[,i] <- values(Hab_raster[[i]])
  }

  # remove the raster data that was completely NA (i.e., rest_cavity)
  # remove first from df used in Mahalanobis
  tmp.hab <- as.data.frame(tmp.mat)
  colnames(tmp.hab) <- names(Hab_raster)


  tmp.hab[is.na(tmp.hab)]<-0
  all.hab <- tmp.hab[colSums(tmp.hab)!=0]

  # then drop from raster stack
  layers_to_drop <- setdiff(names(Hab_raster),colnames(all.hab))
  Hab_raster <- dropLayer(Hab_raster, layers_to_drop)

  # run mahalanobis using all data from raster extent
  D2 <- mahalanobis(all.hab, colMeans(all.hab),cov(all.hab))

  # summary(D2)
  # hist(D2)
  # length(D2)
  # ncell(Hab_raster)
#   hist(D2)
#   plot(density(D2, bw = 0.5),
#        main="Squared Mahalanobis distances, n=100, p=3") ; rug(D2)
#   qqplot(qchisq(ppoints(100), df = 3), D2,
#          main = expression("Q-Q plot of Mahalanobis" * ~D^2 *
#                              " vs. quantiles of" * ~ chi[3]^2))
#   abline(0, 1, col = 'gray')

  D2raster <- raster(ext=extent(Hab_raster), crs=projection(Hab_raster), res=c(100,100))
  values(D2raster) <- D2

  # plot(D2raster<10)

  # add D2 raster layer to raster stack and add name to layer
  Hab_raster_use <- addLayer(Hab_raster, D2raster)
  names(Hab_raster_use)[nlayers(Hab_raster_use)] <- "D2"

  return(Hab_raster_use)

}

D2_rasters <- list()
for(i in 1:length(Hab_rasters)){
  D2_rasters[[i]] <- add_Mahal_raster_layer_function(Hab_raster=Hab_rasters[[i]])
}


options(scipen = 1)
plot(D2_rasters[[3]]) # because using all the data, looks like the higher D2 is reflective or what might be good habitat..i.e., what comprises most of the area (is this a flip from how it would be if using home ranges to quantify)


### Now create Fisher Habitat Extension Zone specific home range raster stacks
# create function and loop it over Hab_zone

HR_D2_function <- function(HRsfobj=HRsfobj, FHEzone=FHEzone, D2_raster=D2_raster){
  # HRsfobj = Female_HRs_compiled
  # FHEzone = Hab_zone[2]
  # D2_raster = D2_rasters[[2]]

  # home ranges overlap each other - need to create rasters for each home range and stack them

  FHR <- HRsfobj %>% filter(Hab_zone==FHEzone)
  FHR <- FHR %>% st_transform(crs=projection(D2_raster))
  FHR <- FHR %>% arrange(SA_Fisher_ID)

  FHR <- FHR %>% st_transform(crs=projection(D2_raster))

  D2_raster_cropped <- list()
  for(i in 1:nrow(FHR)){
    D2_raster_cropped[[i]] <- raster::mask(D2_raster, FHR[i,], maskvalue=TRUE) # add 1 km buffer area around extent)
  }

  # habitat and D2 values for each HR
  D2_tmp_list <- list()
  for(i in 1:length(D2_raster_cropped)){
    D2_tmp_list[[i]] <- D2_raster_cropped[[i]]$D2
  }

  # create raterstack of just D2 for all HRs, easier for plotting
  D2_HR_raster <- stack(D2_tmp_list)
  names(D2_HR_raster) <- paste0(FHR$SA_Fisher_ID,"_D2")

  # plot(D2_HR_raster)

  # HR_Hab_D2 <- as.data.frame(matrix(NA,nrow=nrow(FHR),ncol=ncol(rValues)+2))
  # bring in habitat values from raster

  ###- try this two ways as a test
  # re-run Mahalanobis for each HR with just the data within the HR
  # delete the rows with NAs for each habitat feature?

  # D2_raster_cropped[[1]]
  # summary(D2_raster_cropped[[1]])

  HR_Hab_D2 <- as.data.frame(matrix(NA,nrow=length(D2_raster_cropped),ncol=2*nlayers(D2_raster_cropped[[1]])+3))
  if(ncol(HR_Hab_D2)==13){
  colnames(HR_Hab_D2) <- c("denning_prop","movement_prop", "rest_cwd_prop","rest_rust_prop","D2_mean","D2_subset_mean",
                           "denning_sum","movement_sum", "rest_cwd_sum","rest_rust_sum","D2_sum","D2_subset_sum","Total_Area_ha")
  } else {
    colnames(HR_Hab_D2) <- c("denning_prop","movement_prop", "rest_cavity_prop","rest_cwd_prop","rest_rust_prop","D2_mean","D2_subset_mean",
                             "denning_sum","movement_sum", "rest_cavity_sum","rest_cwd_sum","rest_rust_sum","D2_sum","D2_subset_sum","Total_Area_ha")
  }

  for(i in 1:length(D2_raster_cropped)){
    tmp.raster <- D2_raster_cropped[[i]]

    tmp.mat <- matrix(NA,nrow=ncell(tmp.raster),ncol=nlayers(tmp.raster))

    for(y in 1:ncol(tmp.mat)){
      tmp.mat[,y] <- values(tmp.raster[[y]])
    }

    tmp.hab <- as.data.frame(tmp.mat)
    colnames(tmp.hab) <- names(tmp.raster)

    tmp.hab <- tmp.hab[complete.cases(tmp.hab$D2),]
    tmp.hab[is.na(tmp.hab)]<-0
    D2.hab <- tmp.hab[colSums(tmp.hab)!=0]

    D2_subset <- mahalanobis(D2.hab[,1:ncol(D2.hab)-1],
                             colMeans(D2.hab[,1:ncol(D2.hab)-1]),
                             cov(D2.hab[,1:ncol(D2.hab)-1]))

    tmp.hab$D2_subset <- D2_subset
    prop.data <- t(as.data.frame(colSums(tmp.hab)/nrow(tmp.hab)))
    sum.data <- t(as.data.frame(colSums(tmp.hab)))

    HR_Hab_D2[i,] <- c(prop.data, sum.data,nrow(tmp.hab))
  }


  HR_Hab_D2$SA_Fisher_ID <- FHR$SA_Fisher_ID
  HR_Hab_D2$FISHER_ID <- str_remove(HR_Hab_D2$SA_Fisher_ID,".*_")

  return(list(HR_Hab_D2=HR_Hab_D2, D2_HR_raster=D2_HR_raster,D2_raster_cropped=D2_raster_cropped))

}

HR_D2_perzone <- list()
for(i in 1:length(Hab_zone)){
HR_D2_perzone[[i]] <- HR_D2_function(HRsfobj=Female_HRs_compiled, FHEzone=Hab_zone[i], D2_raster=D2_rasters[[i]])
}

saveRDS(HR_D2_perzone,"HR_D2_perzone.Rda")

# readRDS("HR_D2_perzone.Rda")
################################################################################
###--- Now compare to Rich's vector D2 data
MAHALDir <- "//Sfp.idir.bcgov/s140/S40203/Ecosystems/Conservation Science/Species Conservation Science/Mesocarnivores/Projects/Fisher Critical Habitat/Analysis/FHE_hexes/"

# list.files(MAHALDir)
# Boreal_FHE_home_ranges.Rda
# Sub_Boreal_moist_FHE_home_ranges.Rda
# Sub_Boreal_dry_FHE_home_ranges.Rda
# Dry_forest_FHE_home_ranges.Rda

Hab_zone

load(paste0(MAHALDir,'Boreal_FHE_home_ranges.Rda'))
colnames(Boreal_FHE_home_ranges) <- c("FISHER_ID","denning","rest_rust","rest_cwd","movement","D2","p")
Boreal_FHE_home_ranges <- left_join(Boreal_FHE_home_ranges, HR_D2_perzone[[1]]$HR_Hab_D2)
Boreal_FHE_home_ranges$FHE <- Hab_zone[1]

load(paste0(MAHALDir,'Sub_Boreal_moist_FHE_home_ranges.Rda'))
colnames(Sub_Boreal_moist_FHE_home_ranges) <- c("FISHER_ID","denning","rest_rust","rest_cavity","rest_cwd","movement","D2","p")
Sub_Boreal_moist_FHE_home_ranges <- left_join(Sub_Boreal_moist_FHE_home_ranges, HR_D2_perzone[[2]]$HR_Hab_D2)
Sub_Boreal_moist_FHE_home_ranges$FHE <- Hab_zone[2]

load(paste0(MAHALDir,'Sub_Boreal_dry_FHE_home_ranges.Rda'))
colnames(Sub_Boreal_dry_FHE_home_ranges) <- c("FISHER_ID","denning","rest_rust","rest_cavity","rest_cwd","movement","D2","p")
Sub_Boreal_dry_FHE_home_ranges <- left_join(Sub_Boreal_dry_FHE_home_ranges, HR_D2_perzone[[3]]$HR_Hab_D2)
Sub_Boreal_dry_FHE_home_ranges$FHE <- Hab_zone[3]

load(paste0(MAHALDir,'Dry_forest_FHE_home_ranges.Rda'))
colnames(Dry_forest_FHE_home_ranges) <- c("FISHER_ID","rest_rust","denning","rest_cwd","movement","D2","p")
Dry_forest_FHE_home_ranges <- left_join(Dry_forest_FHE_home_ranges, HR_D2_perzone[[4]]$HR_Hab_D2)
Dry_forest_FHE_home_ranges$FHE <- Hab_zone[4]

names(Boreal_FHE_home_ranges)
names(Sub_Boreal_moist_FHE_home_ranges)
FHE_home_ranges_hab_D2 <- bind_rows(Boreal_FHE_home_ranges,Sub_Boreal_moist_FHE_home_ranges,Sub_Boreal_dry_FHE_home_ranges,Dry_forest_FHE_home_ranges)

# now try running Mahal for the summed data for each home range
head(FHE_home_ranges_hab_D2)
FHE_home_ranges_hab_D2 %>% group_by(FHE) %>% summarise(sum(rest_cavity))

BoR_D2_data_rProp <- FHE_home_ranges_hab_D2 %>% filter(FHE=="Boreal") %>% dplyr::select(denning_prop, movement_prop, rest_cwd_prop, rest_rust_prop)
SBM_D2_data_rProp <- FHE_home_ranges_hab_D2 %>% filter(FHE=="Sub-Boreal moist") %>% dplyr::select(denning_prop, movement_prop, rest_cwd_prop, rest_rust_prop, rest_cavity_prop)
SBD_D2_data_rProp <- FHE_home_ranges_hab_D2 %>% filter(FHE=="Sub-Boreal dry") %>% dplyr::select(denning_prop, movement_prop, rest_cwd_prop, rest_rust_prop, rest_cavity_prop)
DRY_D2_data_rProp <- FHE_home_ranges_hab_D2 %>% filter(FHE=="Dry Forest") %>% dplyr::select(denning_prop, movement_prop, rest_cwd_prop, rest_rust_prop)


FHE_home_ranges_hab_D2$D2_rProp <- c(mahalanobis(BoR_D2_data_rProp,colMeans(BoR_D2_data_rProp), cov(BoR_D2_data_rProp)),
                                     mahalanobis(SBM_D2_data_rProp,colMeans(SBM_D2_data_rProp), cov(SBM_D2_data_rProp)),
                                     mahalanobis(SBD_D2_data_rProp,colMeans(SBD_D2_data_rProp), cov(SBD_D2_data_rProp)),
                                     mahalanobis(DRY_D2_data_rProp,colMeans(DRY_D2_data_rProp), cov(DRY_D2_data_rProp)))


BoR_D2_data_rSum <- FHE_home_ranges_hab_D2 %>% filter(FHE=="Boreal") %>% dplyr::select(denning_sum, movement_sum, rest_cwd_sum, rest_rust_sum)
SBM_D2_data_rSum <- FHE_home_ranges_hab_D2 %>% filter(FHE=="Sub-Boreal moist") %>% dplyr::select(denning_sum, movement_sum, rest_cwd_sum, rest_rust_sum, rest_cavity_sum)
SBD_D2_data_rSum <- FHE_home_ranges_hab_D2 %>% filter(FHE=="Sub-Boreal dry") %>% dplyr::select(denning_sum, movement_sum, rest_cwd_sum, rest_rust_sum, rest_cavity_sum)
DRY_D2_data_rSum <- FHE_home_ranges_hab_D2 %>% filter(FHE=="Dry Forest") %>% dplyr::select(denning_sum, movement_sum, rest_cwd_sum, rest_rust_sum)

FHE_home_ranges_hab_D2$D2_rSum <- c(mahalanobis(BoR_D2_data_rSum,colMeans(BoR_D2_data_rSum), cov(BoR_D2_data_rSum)),
                                     mahalanobis(SBM_D2_data_rSum,colMeans(SBM_D2_data_rSum), cov(SBM_D2_data_rSum)),
                                     mahalanobis(SBD_D2_data_rSum,colMeans(SBD_D2_data_rSum), cov(SBD_D2_data_rSum)),
                                     mahalanobis(DRY_D2_data_rSum,colMeans(DRY_D2_data_rSum), cov(DRY_D2_data_rSum)))

D2_columns <- colnames(FHE_home_ranges_hab_D2)[grepl("D2",colnames(FHE_home_ranges_hab_D2))]
cor(FHE_home_ranges_hab_D2 %>% dplyr::select(D2_columns))

#                         D2     D2_mean D2_subset_mean      D2_sum D2_subset_sum    D2_rProp   D2_rSum
# D2              1.00000000  0.01956824     0.03158791 -0.03381094   -0.04311010  0.53405969 0.1789417
# D2_mean         0.01956824  1.00000000     0.28830361  0.90664197   -0.07556379  0.21019559 0.1844352
# D2_subset_mean  0.03158791  0.28830361     1.00000000  0.33998842    0.36599634  0.19242593 0.2090766
# D2_sum         -0.03381094  0.90664197     0.33998842  1.00000000    0.25786407  0.19592786 0.3758040
# D2_subset_sum  -0.04311010 -0.07556379     0.36599634  0.25786407    1.00000000 -0.03656169 0.4121138
# D2_rProp        0.53405969  0.21019559     0.19242593  0.19592786   -0.03656169  1.00000000 0.6352578
# D2_rSum         0.17894172  0.18443521     0.20907659  0.37580400    0.41211378  0.63525785 1.0000000

write.csv(FHE_home_ranges_hab_D2,"out/FHE_home_ranges_hab_D2.csv")

head(FHE_home_ranges_hab_D2)

FHE_HR_Hab_D2_longer <- FHE_home_ranges_hab_D2 %>% dplyr::select(-SA_Fisher_ID,-D2_subset_sum,-D2_sum) %>%
  pivot_longer(!c(FISHER_ID, FHE), names_to = "Variable", values_to = "Values")

unique(FHE_HR_Hab_D2_longer$Variable)
FHE_HR_Hab_D2_longer$Group <- "Vector_area"
FHE_HR_Hab_D2_longer$Group <- case_when(grepl("D2",FHE_HR_Hab_D2_longer$Variable) ~ "Mahalanobis",
                                        FHE_HR_Hab_D2_longer$Variable=="p" ~ "Mahalanobis",
                                        grepl("prop",FHE_HR_Hab_D2_longer$Variable) ~ "Raster_prop",
                                        grepl("sum",FHE_HR_Hab_D2_longer$Variable) ~ "Raster_area",
                                        FHE_HR_Hab_D2_longer$Variable=="Total_Area_ha" ~ "Raster_area",
                                        TRUE ~ as.character(FHE_HR_Hab_D2_longer$Group))

# as.data.frame(FHE_HR_Hab_D2_longer %>% group_by(Group) %>% dplyr::count(Variable)) # check to make sure it worked


pal = pnw_palette(name="Cascades",n=4,type="discrete")
unique(FHE_HR_Hab_D2_longer$FHE)

FHE_HR_Hab_D2_longer$FHE <- fct_relevel(FHE_HR_Hab_D2_longer$FHE,  "Boreal","Sub-Boreal moist", "Sub-Boreal dry","Dry Forest")
FHE_HR_Hab_D2_longer$Group <- fct_relevel(FHE_HR_Hab_D2_longer$Group,  "Mahalanobis","Vector_area", "Raster_prop","Raster_area")

unique(FHE_HR_Hab_D2_longer$Variable)
FHE_HR_Hab_D2_plot <- ggplot(FHE_HR_Hab_D2_longer %>% filter(Variable!="p"),
                         aes(x=Variable, y=Values, fill=as.factor(FHE))) +
  geom_boxplot(notch=FALSE) +
  # facet_wrap(~Variable, scale="free",nrow=2)+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank()) +
  theme(legend.position="bottom")+
  theme(legend.title=element_blank())+
  scale_fill_manual(values=pal)+
  facet_wrap(~Group, scale="free",ncol=1)

Cairo(file="out/FHE_HR_Hab_D2_plot.PNG",type="png",width=3200,height=4000,pointsize=20,bg="white",dpi=300)
FHE_HR_Hab_D2_plot
dev.off()

Hab_zone
Cairo(file="out/FHE_Boreal_rD2_plot.PNG",type="png",width=4000,height=4000,pointsize=20,bg="white",dpi=300)
plot(HR_D2_perzone[[1]]$D2_HR_raster)
dev.off()

Cairo(file="out/FHE_SBM_rD2_plot.PNG",type="png",width=4000,height=4000,pointsize=20,bg="white",dpi=300)
plot(HR_D2_perzone[[2]]$D2_HR_raster)
dev.off()

Cairo(file="out/FHE_SBD_rD2_plot.PNG",type="png",width=4000,height=4000,pointsize=20,bg="white",dpi=300)
plot(HR_D2_perzone[[3]]$D2_HR_raster)
dev.off()

Cairo(file="out/FHE_Dry_rD2_plot.PNG",type="png",width=4000,height=4000,pointsize=20,bg="white",dpi=300)
plot(HR_D2_perzone[[4]]$D2_HR_raster)
dev.off()
