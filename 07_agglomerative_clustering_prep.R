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
                      "viridis","Cairo","PNWColors","units","rstatix","stars")

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
hist(Female_HRs_Boreal$Are_km2)
sort(Female_HRs_Boreal$Are_km2)
summary(Female_HRs_Boreal[Female_HRs_Boreal$Are_km2<80,]$Are_km2)
# rpois(n=1, lambda=29)


Female_HRs_SBM <- st_read(dsn=paste0(getwd(),"/data"),layer="Female_HRs_SBM")
hist(Female_HRs_SBM$Are_km2)
sort(Female_HRs_SBM$Are_km2)
summary(Female_HRs_SBM[Female_HRs_SBM$Are_km2<80,]$Are_km2)

Female_HRs_SBD <- st_read(dsn=paste0(getwd(),"/data"),layer="Female_HRs_SBD")
hist(Female_HRs_SBD$Are_km2, breaks=10)
sort(Female_HRs_SBD$Are_km2)
summary(Female_HRs_SBD[Female_HRs_SBD$Are_km2<80,]$Are_km2)
sd(Female_HRs_SBD[Female_HRs_SBD$Are_km2<80,]$Are_km2)

Female_HRs_Dry <- st_read(dsn=paste0(getwd(),"/data"),layer="Female_HRs_Dry")
hist(Female_HRs_Dry$Are_km2, breaks=10)
sort(Female_HRs_Dry$Are_km2)
summary(Female_HRs_Dry[Female_HRs_Dry$Are_km2<80,]$Are_km2)
sd(Female_HRs_Dry$Are_km2)

ggplot()+
  geom_sf(data = Female_HRs_Boreal, fill="black") +
  geom_sf(data = Female_HRs_SBM, fill="darkblue") +
  geom_sf(data = Female_HRs_SBD, fill="darkgreen") +
  geom_sf(data = Female_HRs_Dry, fill="darkred")

#####################################################################################
###--- BRING IN HISTORIC VRI TO MATCH WITH HR AOIs
#####################################################################################

VRI2002_Dir <- "//spatialfiles.bcgov/Work/wlap/sry/Workarea/jburgar/Fisher_rasters/2002/"
VRI2008_Dir <- "//spatialfiles.bcgov/Work/wlap/sry/Workarea/jburgar/Fisher_rasters/2008/"
# tif.files <- list.files(VRI2008_Dir, pattern="\\.tif$")

# bring in fisher habitat data - necessary for Mahalanobis
VRI_aoi_raster <- function(FHEzone=FHEzone, Raster_Dir=Raster_Dir){
  # FHEzone = "SBD"
  # Raster_Dir = VRI2002_Dir

  rFiles = list.files(Raster_Dir, pattern="\\.tif$")
  rFHEZ <- rFiles[grepl(FHEzone, rFiles)]
  rFHEZ_HR <- rFHEZ[grepl("HR", rFHEZ)]

  rFHEZ_list=list()
  for(i in 1:length(rFHEZ_HR)){
    rFHEZ_list[[i]] <- raster(paste0(Raster_Dir,rFHEZ_HR[i]))
  }

  FHEZraster_stack = stack(rFHEZ_list)

  # creating additive raster layer (additive)
  # rCostadd # additive layer, equal weight
  rCostadd <- sum(FHEZraster_stack, na.rm = T)
  # plot(rCostadd)

  # rCostmult # multiplicative layer, resting hab all gets a 1, denning and movement are 3 and 2, respectively
  # need to add in if statement to account for cavity in SBD and SBM
  # to deal with the 4 vs 5 habitat characteristics per FHEzone, need to add if statement on number of raster layers
  if(length(rFHEZ_list)==4){
    denning <- rFHEZ_list[[3]]
  } else {
    denning <- rFHEZ_list[[4]]
  }
  denning[denning ==1] <- 3

  if(length(rFHEZ_list)==4){
    movement <- rFHEZ_list[[4]]
  } else {
    movement <- rFHEZ_list[[5]]
  }
  movement[movement==1] <- 2

  if(length(rFHEZ_list)==4){
    tmpstack <- stack(rFHEZ_list[[1]],rFHEZ_list[[2]], denning, movement)
  } else {
    tmpstack <- stack(rFHEZ_list[[1]],rFHEZ_list[[2]], rFHEZ_list[[3]], denning, movement)
  }

  rCostmult <-prod(tmpstack, na.rm=T)
  # plot(rCostmult)

  FHEZraster_stack <- stack(FHEZraster_stack, rCostadd, rCostmult)

  stack.names.to.use <- word(rFHEZ_HR,3,sep="_")
  names(FHEZraster_stack) <- c(stack.names.to.use,"costAdd","costMult")
  # plot(FHEZraster_stack[[6:7]])

  return(FHEZraster_stack)
}

BOR_aoi_raster <- VRI_aoi_raster(FHEzone="Boreal", Raster_Dir=VRI2008_Dir)
writeRaster(BOR_aoi_raster, file="out/BOR_aoi_raster.tif", bylayer=TRUE, overwrite=TRUE)
plot(raster("out/BOR_aoi_raster_5.tif"))

DRY_aoi_raster <- VRI_aoi_raster(FHEzone="Dry", Raster_Dir=VRI2008_Dir)
writeRaster(DRY_aoi_raster, file="out/DRY_aoi_raster.tif", bylayer=TRUE, overwrite=TRUE)
plot(raster("out/DRY_aoi_raster_5.tif"))

SBD_aoi_raster <- VRI_aoi_raster(FHEzone="SBD", Raster_Dir=VRI2002_Dir)
writeRaster(SBD_aoi_raster, file="out/SBD_aoi_raster.tif", bylayer=TRUE, overwrite=TRUE)
plot(raster("out/SBD_aoi_raster_6.tif"))
plot(raster("out/SBD_aoi_raster_7.tif"))

SBM_aoi_raster <- VRI_aoi_raster(FHEzone="SBM", Raster_Dir=VRI2002_Dir)
writeRaster(SBM_aoi_raster, file="out/SBM_aoi_raster.tif", bylayer=TRUE, overwrite=TRUE)
plot(raster("out/SBM_aoi_raster_6.tif"))
plot(raster("out/SBM_aoi_raster_7.tif"))
colSums(SBM_aoi_raster, na.rm=T)

# to plot, convert raster to df
SBM_aoi_raster_df <- as.data.frame(SBM_aoi_raster, xy = TRUE)
head(SBM_aoi_raster_df)
options(scipen = 10)
colSums(SBM_aoi_raster_df, na.rm=T)
SBM_aoi_raster_df %>%
  # filter(costAdd!=0) %>%
  count(costAdd, costMult)

(g_BOR_cost_map <- ggplot(data = SBM_aoi_raster_df) +
    geom_raster(aes(x = x, y = y, fill = `cost`)) +
    geom_sf(data = Female_HRs_Boreal, fill=NA, lwd=1, color="white") +
    scale_fill_viridis_c() +
    theme_void() +
    theme(legend.position = "bottom"))

test <- raster("out/test.tif")
plot(test)
# to plot, convert raster to df
test_df <- as.data.frame(test, xy = TRUE)
head(test_df)

(g_BOR_cost_map <- ggplot(data = test_df) +
    geom_raster(aes(x = x, y = y, fill = `test`)) +
    geom_sf(data = Female_HRs_Boreal, fill=NA, lwd=1, color="white") +
    scale_fill_viridis_c() +
    theme_void() +
    theme(legend.position = "bottom"))

test_multi <- raster("out/dry_test.tif")
plot(test_multi)
# to plot, convert raster to df
test_multi_df <- as.data.frame(test_multi, xy = TRUE)
head(test_multi_df)

(g_DRY_cost_map <- ggplot(data = test_multi_df) +
    geom_raster(aes(x = x, y = y, fill = `dry_test`)) +
    geom_sf(data = Female_HRs_Dry, fill=NA, lwd=1, color="white") +
    scale_fill_viridis_c() +
    # theme_void() +
    theme(legend.position = "bottom"))

(g_BOR_cost_map <- ggplot(data = test_multi_df) +
    geom_raster(aes(x = x, y = y, fill = `test_multi_dt`)) +
    geom_sf(data = Female_HRs_Boreal, fill=NA, lwd=1, color="white") +
    scale_fill_viridis_c() +
    theme_void() +
    theme(legend.position = "bottom"))

test_multi <- read_stars("out/dry_test.tif", package="stars")

multi_poly <- st_as_sf(test_multi, as_points=FALSE, merge=TRUE) %>% st_transform(3005)
colnames(multi_poly)[1] <- "Polygon_ID"
multi_poly$Area_km2 <- st_area(multi_poly)*1e-6
multi_poly <- drop_units(multi_poly)
summary(multi_poly$Area_km2)

# multi_poly$HR_id <- as.numeric(rownames(multi_poly))
multi_poly$SA_F_ID <- paste0("HRId_",rownames(multi_poly))
multi_poly$Hab_zon <- "Dry Forest"
glimpse(multi_poly)

ggplot()+
  geom_sf(data=multi_poly, aes(fill=Polygon_ID))


######################################################################################
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
  if(nlayers(rStack_aoi)==6){
    rStack_aoi <- dropLayer(rStack_aoi,c(5,6))
  } else {
    rStack_aoi <- dropLayer(rStack_aoi,c(6,7))
  }

  raster_cropped <- list()
  for(i in 1:nrow(FHR)){
    rTMP <- list()
    for(t in 1:nlayers(rStack_aoi)){
      rTMP[[t]] <- mask(rFHR[[i]], rStack_aoi[[t]])
    }

    rTMP_stack <- stack(rTMP)
    rTMP_stack <- crop(rTMP_stack, FHR[i,])

    if(nlayers(rTMP_stack)==4){
      rRest <- stackApply(rTMP_stack[[1:2]], indices=1, fun=min)
      rDM <- stackApply(rTMP_stack[[3:4]], indices=1, fun=sum)
    } else {
      rRest <- stackApply(rTMP_stack[[1:3]], indices=1, fun=min)
      rDM <- stackApply(rTMP_stack[[4:5]], indices=1, fun=sum)
    }

    rCost <- stack(rRest, rDM)
    rCost <- stackApply(rCost, indices=1, fun=sum)
    rTMP_stack <- stack(rTMP_stack, rRest, rCost)

    stack.names.to.use <- c(names(rStack_aoi),"Rest","cost")
    names(rTMP_stack) <- stack.names.to.use

    names(rTMP_stack) <- paste(names(rFHR[[i]]), stack.names.to.use, sep="_")
    raster_cropped[[i]] <- rTMP_stack

  }

  ### Run Mahalanobis distance for each HR
  # delete the rows with NAs for each habitat feature
  # include the transformations to deal with non-normal data (although note it doesn't totally fix normality)

  HR_Hab_D2 <- as.data.frame(matrix(NA,nrow=length(raster_cropped),ncol=2*nlayers(raster_cropped[[1]])+1))
  if(ncol(HR_Hab_D2)==13){
    colnames(HR_Hab_D2) <- c("Branch_sum","CWD_sum", "Denning_sum","Movement_sum","Rest_sum","Hab_sum","Total_Area_ha",
                             "Branch_prop","CWD_prop", "Denning_prop","Movement_prop","Rest_prop","Hab_prop")
  } else {
    colnames(HR_Hab_D2) <- c("Branch_sum","Cavity_sum","CWD_sum", "Denning_sum","Movement_sum","Rest_sum","Hab_sum","Total_Area_ha",
                             "Branch_prop","Cavity_prop","CWD_prop", "Denning_prop","Movement_prop","Rest_prop","Hab_prop")
  }

  HR_Hab_D2$HR <- FHR$SA_F_ID


  for(i in 1:nrow(HR_Hab_D2)){
    cost_layer <- nlayers(raster_cropped[[i]])
    values(raster_cropped[[i]][[cost_layer]])[values(raster_cropped[[i]][[cost_layer]]) > 0] = 1

    hab_sums <- cellStats(raster_cropped[[i]], stat="sum", na.rm=T)
    total_area <- cellStats(rFHR[[i]], stat="sum", na.rm=T)
    hab_prop <- hab_sums/total_area

    if(ncol(HR_Hab_D2)==14){
      HR_Hab_D2[i,1:6] <- hab_sums[1:6]
      HR_Hab_D2[i,7] <- total_area
      HR_Hab_D2[i,8:13] <- hab_prop[1:6]
    } else {
      HR_Hab_D2[i,1:7] <- hab_sums[1:7]
      HR_Hab_D2[i,8] <- total_area
      HR_Hab_D2[i,9:15] <- hab_prop[1:7]
    }
  }

  D2_inputs <- HR_Hab_D2 %>% dplyr::select(ends_with("prop")) %>% dplyr::select(-c("Hab_prop","Rest_prop"))

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

Female_HRs_Boreal$FHEzone = Female_HRs_Boreal$Hab_zon
HR_BOR_D2 <- HR_D2_function(HRsfobj=Female_HRs_Boreal, FHEzone="Boreal", rStack_aoi=BOR_aoi_raster)
HR_BOR_D2$HR_Hab_D2$FHE <- "Boreal"

Female_HRs_Dry$FHEzone = Female_HRs_Dry$Hab_zon
HR_DRY_D2 <- HR_D2_function(HRsfobj=Female_HRs_Dry, FHEzone="Dry Forest", rStack_aoi=DRY_aoi_raster)
HR_DRY_D2$HR_Hab_D2$FHE <- "Dry Forest"

Female_HRs_SBM$FHEzone = Female_HRs_SBM$Hab_zon
HR_SBM_D2 <- HR_D2_function(HRsfobj=Female_HRs_SBM, FHEzone="Sub-Boreal moist", rStack_aoi=SBM_aoi_raster)
HR_SBM_D2$HR_Hab_D2$FHE <- "Sub-Boreal moist"

Female_HRs_SBD$FHEzone = Female_HRs_SBD$Hab_zon
HR_SBD_D2 <- HR_D2_function(HRsfobj=Female_HRs_SBD, FHEzone="Sub-Boreal dry", rStack_aoi=SBD_aoi_raster)
HR_SBD_D2$HR_Hab_D2$FHE <- "Sub-Boreal dry"


HR_df <- bind_rows(HR_BOR_D2$HR_Hab_D2, HR_DRY_D2$HR_Hab_D2, HR_SBM_D2$HR_Hab_D2, HR_SBD_D2$HR_Hab_D2)
HR_df[is.na(HR_df)] <- 0 # change NAs to 0 for cavity resting (Boreal and Dry)
HR_df$Fpop <- ifelse(HR_df$FHE=="Boreal","Boreal","Columbian")

write.csv(HR_df, "out/FHE_home_ranges_hab_D2.csv")


######################################################################################
HR_DRY_multi_D2 <- HR_D2_function(HRsfobj=multi_poly %>% filter(Area_km2>10), FHEzone="Dry Forest", rStack_aoi=DRY_aoi_raster)
summary(HR_DRY_multi_D2$HR_Hab_D2$D2)

HR_BOR_multi_D2 <- HR_D2_function(HRsfobj=multi_poly %>% filter(Area_km2>10), FHEzone="Boreal", rStack_aoi=BOR_aoi_raster)
write.csv(HR_BOR_multi_D2$HR_Hab_D2, "out/HR_BOR_multi_D2.csv")
summary(HR_BOR_multi_D2$HR_Hab_D2$D2)
HR_BOR_multi_hab_D2 <- read.csv("out/HR_BOR_multi_D2.csv", row.names=1)

# HR_BOR_multi_hab_D2 <- HR_BOR_multi_hab_D2 %>% filter(Total_Area_ha>1000)
# D2_inputs <- HR_BOR_multi_hab_D2 %>% dplyr::select(ends_with(("prop")))
# D2_inputs$Branch_prop <- log(D2_inputs$Branch_prop+1)
# D2 <- mahalanobis(D2_inputs, colMeans(D2_inputs),cov(D2_inputs))
# HR_BOR_multi_hab_D2$D2_new <- D2

multi_poly$D2 <- HR_BOR_multi_D2$HR_Hab_D2$D2[match(multi_poly$SA_F_ID, HR_BOR_multi_D2$HR_Hab_D2$HR)]
multi_poly$D2 <- HR_DRY_multi_D2$HR_Hab_D2$D2[match(multi_poly$SA_F_ID, HR_DRY_multi_D2$HR_Hab_D2$HR)]

summary(multi_poly$D2)
hist(multi_poly$Area_km2)
multi_poly %>% filter(Area_km2>20)
multi_poly <- multi_poly %>% st_transform(crs=3005)

ggplot()+
  geom_sf(data=multi_poly %>% filter(D2<10), fill="lightblue") +
  geom_sf(data=multi_poly %>% filter(D2<5), fill="darkblue")

hist(HR_BOR_multi_D2$HR_Hab_D2$D2)
hist(HR_DRY_multi_D2$HR_Hab_D2$D2)

ggplot()+
  geom_sf(data=multi_poly %>% filter(!is.na(D2)), aes(fill=D2)) # not projecting properly...
  # geom_sf(data = Female_HRs_Boreal, fill=NA, lwd=1, color="white")
  # geom_sf(data = Female_HRs_Dry, fill=NA, lwd=1, color="white")


sort(D2)
summary(D2)
length(D2)
nrow(HR_BOR_multi_hab_D2)
multi_poly


HR_BOR_raster_D2 <- HR_D2_function(HRsfobj=Female_HRs_Boreal, FHEzone="Boreal", rStack_aoi=BOR_aoi_raster)

plot(HR_BOR_raster_D2$raster_cropped[[1]])
HR_BOR_raster_D2$HR_Hab_D2
write.csv(HR_BOR_raster_D2$HR_Hab_D2, "out/HR_BOR_Hab_D2.csv")


# to plot, convert raster to df
BOR_aoi_raster_df <- as.data.frame(BOR_aoi_raster, xy = TRUE)
head(BOR_aoi_raster_df)

(g_BOR_cost_map <- ggplot(data = BOR_aoi_raster_df) +
    geom_raster(aes(x = x, y = y, fill = `cost`)) +
    geom_sf(data = Female_HRs_Boreal, fill=NA, lwd=1, color="white") +
    scale_fill_viridis_c() +
    theme_void() +
    theme(legend.position = "bottom"))

# check for normality on post-transformed data
# D2_inputs_names <- names(D2_inputs)
# D2_inputs %>% shapiro_test(D2_inputs_names)

D2 <- mahalanobis(D2_inputs, colMeans(D2_inputs),cov(D2_inputs))
HR_Hab_D2$D2 <- D2

HR_DRY_raster_D2 <- HR_D2_function(HRsfobj=Female_HRs_Dry, FHEzone="Dry Forest", rStack_aoi=DRY_aoi_raster)
HR_DRY_raster_D2$HR_Hab_D2
write.csv(HR_DRY_raster_D2$HR_Hab_D2, "out/HR_DRY_Hab_D2.csv")

# to plot, convert raster to df
DRY_aoi_raster_df <- as.data.frame(DRY_aoi_raster, xy = TRUE)
head(DRY_aoi_raster_df)

(g_DRY_cost_map <- ggplot(data = DRY_aoi_raster_df) +
    geom_raster(aes(x = x, y = y, fill = `cost`)) +
    geom_sf(data = Female_HRs_Dry, fill=NA, lwd=1, color="white") +
    scale_fill_viridis_c() +
    theme_void() +
    theme(legend.position = "bottom"))


# checking out options for cost layer...very little overlap so perhaps no reason to do anything but additive
# try  using covariance matrix (Mahal) as weighting tool, conditional on raster layer input
D2_inputs <- HR_BOR_raster_D2$HR_Hab_D2 %>% dplyr::select(ends_with(("prop")))
D2_inputs$Branch_prop <- log(D2_inputs$Branch_prop+1)
cor(D2_inputs)


# habitat and D2 values for each HR
D2_tmp_values <- matrix(NA,nrow=length(HR_BOR_raster_D2$raster_cropped[[1]]),ncol=nlayers(HR_BOR_raster_D2$raster_cropped[[1]]))
for(i in 1:nlayers(HR_BOR_raster_D2$raster_cropped[[1]])){
  D2_tmp_values[,i] <- HR_BOR_raster_D2$raster_cropped[[1]][[i]]@data@values
}

# subset to where at least movement (most number of pixels) has an input value
D2_tmp2 <- D2_tmp_values[complete.cases(D2_tmp_values[,4]),]
D2_tmp2[is.na(D2_tmp2)] <-0 # change all NAs to 0s

sum(rowSums(D2_tmp2)==4)
sum(rowSums(D2_tmp2)==3)
sum(rowSums(D2_tmp2)==2)
sum(rowSums(D2_tmp2)==1)
length(D2_tmp2)

#################################################################################
###--- look at proportions of habitat for input parameter specification
HR_df <- read.csv("out/FHE_home_ranges_hab_D2.csv", row.names=1)
HR_df %>% as_tibble()

HR_df %>% group_by(Fpop) %>% summarise(across(Hab_prop, list(min=min, mean=mean, max=max, sd=sd)))
HR_df %>% group_by(Fpop) %>% summarise(across(Denning_prop, list(min=min, mean=mean, max=max, sd=sd)))
HR_df %>% group_by(Fpop) %>% summarise(across(Movement_prop, list(min=min, mean=mean, max=max, sd=sd)))
HR_df %>% group_by(Fpop) %>% summarise(across(Rest_prop, list(min=min, mean=mean, max=max, sd=sd)))

HR_df %>% group_by(Fpop) %>% summarise(across(Denning_sum, list(min=min, mean=mean, max=max, sd=sd)))
HR_df %>% group_by(Fpop) %>% summarise(across(Movement_sum, list(min=min, mean=mean, max=max, sd=sd)))
HR_df %>% group_by(Fpop) %>% summarise(across(Rest_sum, list(min=min, mean=mean, max=max, sd=sd)))

# plot(HR_df[HR_df$FHE=="Boreal",]$denning_prop ~ HR_df[HR_df$FHE=="Boreal",]$Total_Area_ha)
# plot(HR_df[HR_df$FHE=="Sub-Boreal moist",]$denning_prop ~ HR_df[HR_df$FHE=="Sub-Boreal moist",]$Total_Area_ha)
# plot(HR_df[HR_df$FHE=="Sub-Boreal dry",]$denning_prop ~ HR_df[HR_df$FHE=="Sub-Boreal dry",]$Total_Area_ha)
# plot(HR_df[HR_df$FHE=="Dry Forest",]$denning_prop ~ HR_df[HR_df$FHE=="Dry Forest",]$Total_Area_ha)

names(HR_df)
HR_df %>% dplyr::select(-starts_with("D2"),-p)

HR_sum_longer <- HR_df %>% dplyr::select(ends_with("sum"), FHE, SA_Fisher_ID, Total_Area_ha) %>%
  pivot_longer(cols=ends_with("sum"), names_to="hab_char", values_to="sum_values")
HR_sum_longer$hab_char <- str_replace(HR_sum_longer$hab_char, "_sum", "")

HR_prop_longer <- HR_df %>% dplyr::select(ends_with("prop"),SA_Fisher_ID) %>%
  pivot_longer(cols=ends_with("prop"), names_to="hab_char", values_to="prop_values")
HR_prop_longer$hab_char <- str_replace(HR_prop_longer$hab_char, "_prop", "")

HR_longer <- left_join(HR_sum_longer %>% filter(!grepl("D2",  hab_char)),
                       HR_prop_longer %>% filter(!grepl("D2", hab_char)),
                       by=c("SA_Fisher_ID", "hab_char"))


Cairo(file=paste0("out/hab_vs_area_sum_plot.PNG"), type="png", width=4000, height=3000,pointsize=15,bg="white",dpi=300)
ggplot(HR_longer %>% filter(Total_Area_ha<14000), aes(x=Total_Area_ha, y=sum_values, color=FHE))+
  geom_point() +
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE)+
  facet_wrap(~hab_char)+
  # theme_classic()+
  theme(legend.position = "bottom", legend.title =element_blank())
dev.off()

Cairo(file=paste0("out/hab_vs_area_prop_plot.PNG"), type="png", width=4000, height=3000,pointsize=15,bg="white",dpi=300)
ggplot(HR_longer %>% filter(Total_Area_ha<14000), aes(x=Total_Area_ha, y=prop_values, color=FHE))+
  geom_point() +
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE)+
  facet_wrap(~hab_char)+
  theme(legend.position = "bottom", legend.title =element_blank())
dev.off()
