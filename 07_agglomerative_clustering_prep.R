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
# 07_agglomerative_clustering_prep.R
# script to create trial data to run and test agglomerative clustering
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 07-July-2022
#####################################################################################
version$major
version$minor
R_version <- paste0("R-",version$major,".",version$minor)

.libPaths(paste0("C:/Program Files/R/",R_version,"/library")) # to ensure reading/writing libraries from C drive
tz = Sys.timezone() # specify timezone in BC

# Load Packages
list.of.packages <- c("tidyverse", "lubridate","bcdata", "sf", "rgdal","nngeo","raster",
                      "viridis","Cairo","PNWColors","units","rstatix")

# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

#####################################################################################
# keep in mind that WGS84 lat/long espg = 4326; BC Albers espg = 3005; NAD83 / UTM zone 10N espg = 26910
#####################################################################################

#####################################################################################
###--- BRING IN FEMALE FISHER HOME RANGE AOIs
#####################################################################################

Female_HRs_Boreal <- st_read(dsn=paste0(getwd(),"/data"),layer="Female_HRs_Boreal")
Female_HRs_Boreal <- Female_HRs_Boreal %>% filter(grepl("KPF",SA_F_ID))
Female_HRs_SBM <- st_read(dsn=paste0(getwd(),"/data"),layer="Female_HRs_SBM")
Female_HRs_SBD <- st_read(dsn=paste0(getwd(),"/data"),layer="Female_HRs_SBD")
Female_HRs_Dry <- st_read(dsn=paste0(getwd(),"/data"),layer="Female_HRs_Dry")

ggplot()+
  geom_sf(data = Female_HRs_Boreal, fill="black") +
  geom_sf(data = Female_HRs_SBM, fill="darkblue") +
  geom_sf(data = Female_HRs_SBD, fill="darkgreen") +
  geom_sf(data = Female_HRs_Dry, fill="darkred")

#####################################################################################
###--- BRING IN HISTORIC VRI TO MATCH WITH HR AOIs
#####################################################################################

VRI2008_Dir <- "//spatialfiles.bcgov/Work/wlap/sry/Workarea/jburgar/Fisher_rasters/2008/"
tif.files <- list.files(VRI2008_Dir, pattern="\\.tif$")

# bring in fisher habitat data - necessary for Mahalanobis
VRI_aoi_raster <- function(FHEzone=FHEzone, Raster_Dir=Raster_Dir){
  # FHEzone = "Dry"
  # Raster_Dir = VRI2008_Dir

  rFiles = list.files(Raster_Dir, pattern="\\.tif$")
  rFHEZ <- rFiles[grepl(FHEzone, rFiles)]
  rFHEZ_HR <- rFHEZ[grepl("HR", rFHEZ)]

  rFHEZ_list=list()
  for(i in 1:length(rFHEZ_HR)){
    rFHEZ_list[[i]] <- raster(paste0(Raster_Dir,rFHEZ_HR[i]))
    # rTMP[is.na(rTMP[])] <- 0
    # rFHEZ_list[[i]] <- rTMP
  }

  FHEZraster_stack = stack(rFHEZ_list)

  stack.names.to.use <- word(rFHEZ_HR,3,sep="_")
  names(FHEZraster_stack) <- stack.names.to.use
  # plot(FHEZraster_stack)

  return(FHEZraster_stack)
}

BOR_aoi_raster <- VRI_aoi_raster(FHEzone="Boreal", Raster_Dir=VRI2008_Dir)
DRY_aoi_raster <- VRI_aoi_raster(FHEzone="Dry", Raster_Dir=VRI2008_Dir)

HR_D2_function <- function(HRsfobj=HRsfobj, FHEzone=FHEzone, rStack_aoi=rStack_aoi){
  # home ranges overlap each other - need to create rasters for each home range and stack them

  FHR <- HRsfobj %>% filter(Hab_zon==FHEzone)
  FHR <- FHR %>% st_transform(crs=projection(rStack_aoi))
  FHR <- FHR %>% arrange(SA_F_ID)
  FHR$Fisher <- 1

  rFHR <- list()
  for(i in 1:nrow(FHR)){
    rFHR[[i]] <- rasterize(FHR[i,], rStack_aoi, field="Fisher", background=0)
  }
  rFHR <- stack(rFHR)
  names(rFHR) <- FHR$SA_F_ID

  ### Create a list with raster stacks of habitat values for each HR
  # length of list is number of HR
  # for each HR, raster stack of habitat rasters, specific to HR but same extent as aoi
  raster_cropped <- list()
  for(i in 1:nrow(FHR)){
    rTMP <- list()
    for(t in 1:nlayers(rStack_aoi)){
      rTMP[[t]] <- mask(rFHR[[i]], rStack_aoi[[t]])
    }
    rTMP_stack <- stack(rTMP)
    rTMP_stack <- crop(rTMP_stack, FHR[i,])
    names(rTMP_stack) <- paste(names(rFHR[[i]]),names(rStack_aoi),sep="_")
    raster_cropped[[i]] <- rTMP_stack
  }


  ### Run Mahalanobis distance for each HR
  # delete the rows with NAs for each habitat feature
  # include the transformations to deal with non-normal data (although note it doesn't totally fix normality)

  HR_Hab_D2 <- as.data.frame(matrix(NA,nrow=length(raster_cropped),ncol=2*nlayers(raster_cropped[[1]])+1))
  if(ncol(HR_Hab_D2)==9){
    colnames(HR_Hab_D2) <- c("Branch_sum","CWD_sum", "Denning_sum","Movement_sum","Total_Area_ha",
                             "Branch_prop","CWD_prop", "Denning_prop","Movement_prop")
  } else {
    colnames(HR_Hab_D2) <- c("Branch_sum","CWD_sum", "Denning_sum","Movement_sum","Cavity_sum","Total_Area_ha",
                             "Branch_prop","CWD_prop", "Denning_prop","Movement_prop","Cavity_prop")
  }

  HR_Hab_D2$HR <- FHR$SA_F_ID

  for(i in 1:nrow(HR_Hab_D2)){
    hab_sums <- colSums(raster_cropped[[i]]@data@values, na.rm=T)
    total_area <- sum(rFHR[[i]]@data@values)
    hab_prop <- hab_sums/total_area

    HR_Hab_D2[i,1:4] <- hab_sums
    HR_Hab_D2[i,5] <- total_area
    HR_Hab_D2[i,6:9] <- hab_prop
  }

  D2_inputs <- HR_Hab_D2 %>% dplyr::select(ends_with(("prop")))

  if(FHEzone=="Dry Forest"|FHEzone=="Boreal"){
    D2_inputs$Branch_prop <- log(D2_inputs$Branch_prop+1)
  } else { if(FHEzone=="Sub-Boreal moist"){
    D2_inputs$Denning_prop <- log(D2_inputs$Denning_prop+1)
    D2_inputs$Cavity_prop <- log(D2_inputs$Cavity+1)
  } else { if(FHEzone=="Sub-Boreal dry"){
    D2_inputs$Denning_prop <- log(D2_inputs$Denning_prop+1)
  }
  }
  }

  # check for normality on post-transformed data
  # D2_inputs_names <- names(D2_inputs)
  # D2_inputs %>% shapiro_test(D2_inputs_names)

  D2 <- mahalanobis(D2_inputs, colMeans(D2_inputs),cov(D2_inputs))
  HR_Hab_D2$D2 <- D2

  return(list(rFHR=rFHR, raster_cropped=raster_cropped, HR_Hab_D2=HR_Hab_D2))

}


HR_BOR_raster_D2 <- HR_D2_function(HRsfobj=Female_HRs_Boreal, FHEzone="Boreal", rStack_aoi=BOR_aoi_raster)
plot(HR_BOR_raster_D2$rFHR)
HR_BOR_raster_D2$HR_Hab_D2

HR_DRY_raster_D2 <- HR_D2_function(HRsfobj=Female_HRs_Dry, FHEzone="Dry Forest", rStack_aoi=DRY_aoi_raster)
HR_DRY_raster_D2$HR_Hab_D2



# try  using covariance matrix (Mahal) as weighting tool, conditional on raster layer input


###############################################################################
###--- Now build home range 'clusters' and run them with Mahal
# https://www.r-bloggers.com/2021/03/clustering-similar-spatial-patterns/

# library(motif)
# library(stars)

plot(D2_rasters[[1]]) # 8800 km2
8800/30 # potential for ~290 fisher territories if ave size = 30 km2
test <- st_as_stars(D2_rasters[[1]][[1]],D2_rasters[[1]][[2]],D2_rasters[[1]][[3]],D2_rasters[[1]][[4]])

# ?lsp_signature
test_signature = lsp_signature(test,
                               type = "coma", #co-occurence matrix
                               window = 10)

test_dist = lsp_to_dist(test_signature, dist_fun="jensen-shannon")

test_hclust = hclust(test_dist, method = "ward.D2")
plot(test_hclust)


clusters = cutree(test_hclust, k = 4)
test_grid_sf = lsp_add_clusters(test_signature,clusters)
test_grid_sf %>% count(clust)

ggplot()+
  geom_sf(data=test_grid_sf, aes(fill=clust))

# https://www.r-bloggers.com/2021/02/finding-similar-spatial-patterns/
# maybe better to use the home range polygons to then find the spatial patterns in the landscape
#

