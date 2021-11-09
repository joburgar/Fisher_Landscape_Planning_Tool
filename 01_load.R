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

.libPaths("C:/Program Files/R/R-4.1.1/library") # to ensure reading/writing libraries from C drive
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

retrieve_geodata_aoi <- function (ID=ID){
  aoi.geodata <- bcdc_query_geodata(ID) %>%
    filter(BBOX(st_bbox(aoi))) %>%
    collect()
  aoi.geodata <- aoi.geodata %>% st_intersection(aoi)
  aoi.geodata$Area_km2 <- st_area(aoi.geodata)*1e-6
  aoi.geodata <- drop_units(aoi.geodata)
  return(aoi.geodata)
}

#####################################################################################

# first find the general aoi (will refine)
ecoprov <- ecoprovinces() %>% filter(ECOPROVINCE_CODE %in% c("BOP","SBI","CEI"))

# create aoi of dissolved area
aoi <- ecoprov %>% st_transform(crs = 3005) # ensures poly is in Albers
sf::st_crs(aoi) = 3005

aoi <- aoi %>%
  summarise(across(geometry, ~ st_union(.))) %>%
  summarise(across(geometry, ~ st_combine(.)))


# plot the ecoprovinces, ordering from CEI to BOP
ecoprov$ECOPROVINCE_NAME <- factor(ecoprov$ECOPROVINCE_NAME,
                                   levels = c("CENTRAL INTERIOR", "SUB-BOREAL INTERIOR", "BOREAL PLAINS"))

ggplot()+
  geom_sf(data = aoi, fill=NA, color="black", lwd=1)+
  geom_sf(data = ecoprov, aes(fill=ECOPROVINCE_NAME))+
  scale_fill_manual(values = pnw_palette("Lake",3))


# Import BEC, Ecoprovinces, Ecoregions and Ecosections
#SBSwk in Peace Region only

fisher.bec.zones <- c("CWH","ICH","IDF","BWBS","MS","SBPS","SBS")
fisher.bec.subzones <- c("BWBSdk|BWBSmw|BWBSwk|SBSwk|ICHwc|SBSmc|SBSmh|SBSmk|SBSmm|SBSmw|SBSdh|SBSdk|SBSdw|SBPSxc|SBPSmc|SBPSdc|SBPSmk|IDFdk|IDFmw|IDFdw|IDFww|
MSxc|MSxk|MSdv|MSdm|MSdk|MSdc|ICHmk|ICHmw")


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
aoi.sub %>% filter(ecoprov=="BOREAL PLAINS") %>% filter(str_detect(MAP_LABEL, "BWBSdk|BWBSmw|BWBSwk|SBSwk|ICHwc"))

glimpse(aoi.sub)
aoi.sub$habitat <- ifelse(aoi.sub$ecoprov=="BOREAL PLAINS" & grepl("BWBSdk|BWBSmw|BWBSwk|SBSwk|ICHwc", aoi.sub$MAP_LABEL), 1,
                          ifelse(aoi.sub$ecoprov=="SUB-BOREAL INTERIOR" & grepl("CWH|ICH|IDF", aoi.sub$MAP_LABEL),1,
                                 ifelse(aoi.sub$ecoprov=="SUB-BOREAL INTERIOR" & grepl("SBSwk|SBSmc|SBSmh|SBSmk|SBSmm|SBSmw|SBSdh|SBSdk|SBSdw", aoi.sub$MAP_LABEL),1,
                                        ifelse(aoi.sub$ecoprov=="CENTRAL INTERIOR" & grepl("SBPSxc|SBPSmc|SBPSdc|SBPSmk|IDFdk|IDFmw|IDFdw|IDFww|MSxc|MSxk|MSdv|MSdm|MSdk|MSdc|ICHmk|ICHmw|ICHmk|SBSdw|SBSmc", aoi.sub$MAP_LABEL),
                                               1,0))))


aoi.fisher.habitat <- aoi.sub %>% filter(habitat==1) %>%
  summarise(across(geometry, ~ st_combine(.))) %>%
  summarise(across(geometry, ~ st_union(.)))


ggplot()+
  # geom_sf(data = aoi)+
  geom_sf(data=aoi.fisher.habitat)

glimpse(aoi.fisher.habitat)
st_write(aoi.sub, dsn=paste0(getwd(),"/out/aoi.fisher.habitat.kml"))

