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

#####################################################################################
# 04_output.R
# script to produce outputs of Individual Based Models (IBMs) for fisher
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 25-Jan-2022
#####################################################################################
version$major
version$minor
R_version <- paste0("R-",version$major,".",version$minor)

.libPaths(paste0("C:/Program Files/R/",R_version,"/library")) # to ensure reading/writing libraries from C drive
tz = Sys.timezone() # specify timezone in BC

# Load Packages
list.of.packages <- c("tidyverse", "NetLogoR","nnls","lcmix","MASS","Cairo","PNWColors", "ggplot2",
                      "sf","raster","rgdal","fasterize")
# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

source("00_IBM_functions.R")
#####################################################################################
# Create 3 sets of 100 simulations - vary the proportion of habitat and survival
# Low, medium and high habitat = 0.5, 0.6, and 0.7 (same world set up, get actual values)
# Low, medium and high survival = 0.7, 0.8, 0.9

load("out/CI_noescape_rfsurv.RData")
w1 <- CI_noescape_rfsurv[[1]]; w1$actual.prop.hab # 0.52
w2 <- CI_noescape_rfsurv[[2]]; w2$actual.prop.hab # 0.67
w3 <- CI_noescape_rfsurv[[3]]; w3$actual.prop.hab # 0.77

load("out/Boreal_noescape_rfsurv.RData")

###--- plot the simulated landbases
Cairo(file="out/BCI_noescape_w1.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
plot(w1$land, main="Simulated Landbase\n52% Suitable Habitat")
points(w1$t0, pch = w1$t0$shape, col = of(agents = w1$t0, var = "color"))
dev.off()

Cairo(file="out/BCI_noescape_w2.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
plot(w2$land, main="Simulated Landbase\n67% Suitable Habitat")
points(w2$t0, pch = w2$t0$shape, col = of(agents = w2$t0, var = "color"))
dev.off()

Cairo(file="out/BCI_noescape_w3.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
plot(w3$land, main="Simulated Landbase\n77% Suitable Habitat")
points(w3$t0, pch = w3$t0$shape, col = of(agents = w3$t0, var = "color"))
dev.off()

# Run 100 simulations for each, save as objects
# Calculate mean # of animals per cell at 10 years for each simulation to produce a heat map
# Create a figure with mean number of animals (+/- SE) for each time step and graph for each simulation

sim_output <- function(sim_out=sim_out, sim_order=sim_order, numsims=numsims){
  ABM.df <- as.data.frame(array(NA,c(200,23)))
  colnames(ABM.df) <- paste0("TimeStep_",str_pad(seq_len(23),2,pad="0"))
  for(i in 1:numsims){
    ABM.df[i,] <- unlist(lapply(lapply(sim_out[[sim_order]][[i]], as.array), ncol)) # if 14 then know that has at least one fisher
    ABM.df[i+numsims,] <- unlist(lapply(lapply(sim_out[[sim_order]][[i]], as.array), nrow)) # if 14 then know that has at least one fisher
  }

  ABM.df$Type <- rep(c("Pfisher","Count"), each=numsims)
  ABM.df$Run <- rep(seq_len(numsims), times=2)

  ABM.df <- ABM.df %>% pivot_longer(cols = TimeStep_01:TimeStep_23,names_to = "TimeStep",values_to = "Value" )
  ABM.df <- ABM.df %>% pivot_wider(names_from = Type, values_from = Value)

  ABM.df$NewCount <- as.numeric(ABM.df$Count)
  ABM.df$NewCount <- case_when(is.na(ABM.df$Pfisher) ~ 0,
                                        TRUE ~ ABM.df$NewCount)

  ABM.df <- ABM.df %>% dplyr::select(Run, TimeStep, NewCount)
  ABM.df$Sim <- paste0("Sim",str_pad(sim_order,2,pad="0"))
  return(ABM.df)
}


# Now format all of the simulated output from lists into one df with number of fisher per time step
# create the dataframe
B_ABM.df <- CI_ABM.df <- as.data.frame(array(NA,c(3*2300,4)))
colnames(B_ABM.df) <- colnames(CI_ABM.df) <- c("Run","TimeStep","NewCount","Sim")

# loop to put in all of the values
# Central Interior
# starting point of data frame
a=1
b=2300
for(i in 4:6){
  CI_ABM.df[a:b,] <- sim_output(sim_out=CI_noescape_rfsurv, sim_order=i, numsims=100)
  a=a+2300
  b=b+2300
}

# Boreal
# starting point of data frame
a=1
b=2300
for(i in 4:6){
  B_ABM.df[a:b,] <- sim_output(sim_out=Boreal_noescape_rfsurv, sim_order=i, numsims=100)
  a=a+2300
  b=b+2300
}

CI_ABM.df$Pop <- "Central Interior"
B_ABM.df$Pop <- "Boreal"

###---
# IBM.w1.rfsurv.sim100 # Sim04
# IBM.w2.rfsurv.sim100 # Sim05
# IBM.w3.rfsurv.sim100 # Sim06

ABM.df <- rbind(CI_ABM.df, B_ABM.df)

ABM.df <- ABM.df %>% mutate(Prophab = case_when(Sim %in% c("Sim04") ~ w1$actual.prop.hab,
                                                Sim %in% c("Sim05") ~ w2$actual.prop.hab,
                                                Sim %in% c("Sim06") ~ w3$actual.prop.hab))

ABM.TS.mean <- ABM.df %>% dplyr::select(-Run) %>% pivot_wider(names_from=TimeStep, values_from=NewCount, values_fn=mean)
ABM.TS.mean$Param <- "Mean"

se <- function(x) sqrt(var(x)/length(x))
ABM.TS.se <- ABM.df %>% dplyr::select(-Run) %>% pivot_wider(names_from=TimeStep, values_from=NewCount, values_fn=se)
ABM.TS.se$Param <- "SE"

ABM.TS <- rbind(ABM.TS.mean, ABM.TS.se)

ABM.TS.df <- ABM.TS %>% pivot_longer(cols = TimeStep_01:TimeStep_23,names_to = "TimeStep",values_to = "Value" )
ABM.TS.df <- ABM.TS.df %>% pivot_wider(names_from = Param, values_from = Value)


ABM.TS.use <- ABM.TS.df %>% filter(!TimeStep %in% c("TimeStep_01", "TimeStep_02"))

sim.TS.plot <- ggplot(data = ABM.TS.use) +
  theme_bw() + theme(strip.background = element_rect(fill = "white", colour = "white")) +
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = "TimeStep_11", col="grey", lty=4) +
  geom_point(aes(x = TimeStep, y = Mean), size=2) +
  geom_errorbar(aes(x = TimeStep, y = Mean, ymin=Mean-SE, ymax= Mean+SE),
                width=.2, position=position_dodge(0.05)) +
  theme(axis.text.x = element_blank()) +
  xlab("Time Step in 6 Month Intervals over 10 years") +
  ylab("Number of Fishers Alive (Mean \u00B1 1 SE)")+ # \u00B1 is ± in unicode
  ggtitle("Simulations of Fisher Populations (100 Runs)\nBy Population and Proportion of Suitable Habitat")+
  facet_wrap(~Pop+Prophab)

sim.TS.plot

#- Plot

Cairo(file="out/BCI_sim_noescape.TS.plot.PNG",
      type="png",
      width=3000,
      height=2200,
      pointsize=15,
      bg="white",
      dpi=300)
sim.TS.plot
dev.off()

################################################################################

### Create heatmaps for the w3 outputs
# keep in mind that WGS84 lat/long espg = 4326; BC Albers espg = 3005; NAD83 / UTM zone 10N espg = 26910

Nozero.runs <- ABM.df %>% filter(TimeStep=="TimeStep_23") %>%
  group_by(Sim) %>%
  filter(NewCount!=0)

nozerosims <- function(sim=Sim, pop=Pop){

  nozerosims <- Nozero.runs %>% filter(Sim==sim & Pop==pop) %>% dplyr::select(Run)
  nozerosims <- unique(nozerosims$Run)

  return(nozerosims)
}

Bph52_nozero <- nozerosims(sim="Sim04", pop="Boreal")
Bph67_nozero <- nozerosims(sim="Sim05", pop="Boreal")
Bph77_nozero <- nozerosims(sim="Sim06", pop="Boreal")

CIph52_nozero <- nozerosims(sim="Sim04", pop="Central Interior")
CIph67_nozero <- nozerosims(sim="Sim05", pop="Central Interior")
CIph77_nozero <- nozerosims(sim="Sim06", pop="Central Interior")


length(Bph52_nozero); length(Bph67_nozero); length(Bph77_nozero) # 52, 78, 83 reached 10 years
length(CIph52_nozero); length(CIph67_nozero); length(CIph77_nozero) # none reached 10 years


# find coordinates for each fisher at 11 year mark
# create a heat map based on number of times fisher is on pixel
# need to consider mean # of fishers vs fisher present/absent

# convert the worlds to rasters (for plotting habitat...but doesn't matter for extent)
rw1 <- world2raster(w1$land)
rw2 <- world2raster(w2$land)
rw3 <- world2raster(w3$land)

extent(rw1) # extents of three worlds the same so can use the same raster as base
plot(rw1)


raster_output <- function(sim_out=sim_out, sim_order=sim_order, sim_use=sim_use, land=land,
                         TS=TS, rExtent=rExtent, rFun=rFun, sFun=sFun){

  # sim_out=IBM_noescape
  # sim_order=6
  # sim_use=Sim06_nozero
  # land=w1$land
  # TS=23
  # rExtent=rw1
  # rFun="sum"
  # sFun="sum"

  r <- raster()
  r <- setExtent(r, rExtent, keepres=TRUE)

  r_list=list()

  # for simulations where at least one fisher survived
  for(i in 1:length(sim_use)){
    ftmp <- as.data.frame(patchHere(land, sim_out[[sim_order]][[sim_use[i]]][[TS]]))
    ftmp$Fisher <- 1
    ftmp.sf <- st_as_sf(ftmp, coords = c("pxcor", "pycor"))
    ftmp.sfp <- st_buffer(ftmp.sf, dist=.1)

    # r_list[[i]] <- fasterize(ftmp.sfp, r, field="Fisher", fun=rFun, background=0)
    r_list[[i]] <- rasterize(ftmp.sfp, r, field="Fisher", fun=rFun, background=0) # interim work around until terra and new raster package uploaded
  }

  r_zeroes <- raster()
  r_zeroes <- setExtent(r_zeroes, rExtent, keepres=TRUE)
  values(r_zeroes) <- 0

  r_zeroes_list=list()

  if(length(sim_use)!=100){
    for(i in 1:(100-length(sim_use))){
      r_zeroes_list[[i]] <- r_zeroes
    }
  }

  r_stack = stack(r_list, r_zeroes_list)
  r_stackApply <- stackApply(r_stack, indices=1, fun=sFun)
  writeRaster(r_stackApply, file=paste0("out/rSim",str_pad(sim_order,2,pad="0"),".tif"), bylayer=TRUE, overwrite=TRUE)

  Fisher_Nmean <- mean(r_stackApply@data@values)
  Fisher_Nse <- se(r_stackApply@data@values)

  return(list(raster=r_stackApply, Fisher_Nmean=Fisher_Nmean, Fisher_Nse=Fisher_Nse))

}

###--- For Sim04 in Boreal
# Bph52_nozero
rBph52 <- raster_output(sim_out=Boreal_noescape_rfsurv, sim_order=4, sim_use=Bph52_nozero,land=w1$land,
                        TS=23, rExtent=rw1, rFun="sum",sFun="sum")

rBph52
plot(rBph52$raster)

length(Bph52_nozero) # 52
w1$t0
Cairo(file="out/rBph52_title.PNG", type="png", width=2200, height=2000,pointsize=15,bg="white",dpi=300)
plot(rBph52$raster,
     main="Estimated Fisher Abundance over 100 Simulations",
     sub="Starting with 20 fishers and 52% suitable habitat\npredicted 27.9 \u00B1 2.6 (mean \u00B1 1 SE) fishers after 10 years.")
dev.off()

###--- For Sim05 in Boreal
# Bph67_nozero
rBph67 <- raster_output(sim_out=Boreal_noescape_rfsurv, sim_order=5, sim_use=Bph67_nozero,land=w2$land,
                        TS=23, rExtent=rw2, rFun="sum",sFun="sum")

rBph67
plot(rBph67$raster)

length(Bph67_nozero) # 78
w2$t0
Cairo(file="out/rBph67_title.PNG", type="png", width=2200, height=2000,pointsize=15,bg="white",dpi=300)
plot(rBph67$raster,
     main="Estimated Fisher Abundance over 100 Simulations",
     sub="Starting with 20 fishers and 67% suitable habitat\npredicted 54.3 \u00B1 6.5 (mean \u00B1 1 SE) fishers after 10 years.")
dev.off()

###--- For Sim06 in Boreal
# Bph77_nozero
rBph77 <- raster_output(sim_out=Boreal_noescape_rfsurv, sim_order=6, sim_use=Bph77_nozero,land=w3$land,
                        TS=23, rExtent=rw3, rFun="sum",sFun="sum")

rBph77
plot(rBph77$raster)

length(Bph77_nozero) # 83
w3$t0
Cairo(file="out/rBph77_title.PNG", type="png", width=2200, height=2000,pointsize=15,bg="white",dpi=300)
plot(rBph77$raster,
     main="Estimated Fisher Abundance over 100 Simulations",
     sub="Starting with 20 fishers and 77% suitable habitat\npredicted 69.7 \u00B1 6.6 (mean \u00B1 1 SE) fishers after 10 years.")
dev.off()
