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

load("out/IBM_noescape.RData")
w1 <- IBM_noescape[[1]]; w1$actual.prop.hab # 0.45
w2 <- IBM_noescape[[2]]; w2$actual.prop.hab # 0.56
w3 <- IBM_noescape[[3]]; w3$actual.prop.hab # 0.71

###--- plot the simulated landbases
Cairo(file="out/IBM_noescape_w1.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
plot(w1$land, main="Simulated Landbase\n45% Suitable Habitat")
points(w1$t0, pch = w1$t0$shape, col = of(agents = w1$t0, var = "color"))
dev.off()

Cairo(file="out/IBM_noescape_w2.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
plot(w2$land, main="Simulated Landbase\n56% Suitable Habitat")
points(w2$t0, pch = w2$t0$shape, col = of(agents = w2$t0, var = "color"))
dev.off()

Cairo(file="out/IBM_noescape_w3.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
plot(w3$land, main="Simulated Landbase\n71% Suitable Habitat")
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
ABM.df <- as.data.frame(array(NA,c(9*2300,4)))
colnames(ABM.df) <- c("Run","TimeStep","NewCount","Sim")

# starting point of data frame
a=1
b=2300

# loop to put in all of the values
for(i in 4:12){
  ABM.df[a:b,] <- sim_output(sim_out=IBM_noescape, sim_order=i, numsims=100)
  a=a+2300
  b=b+2300
}


###---
# IBM.w1.surv7.sim100 # Sim04
# IBM.w1.surv8.sim100 # Sim05
# IBM.w1.surv9.sim100 # Sim06
# IBM.w2.surv7.sim100 # Sim07
# IBM.w2.surv8.sim100 # Sim08
# IBM.w2.surv9.sim100 # Sim09
# IBM.w3.surv7.sim100 # Sim10
# IBM.w3.surv8.sim100 # Sim11
# IBM.w3.surv9.sim100 # Sim12

ABM.df <- ABM.df %>% mutate(Prophab = case_when(Sim %in% c("Sim04", "Sim05", "Sim06") ~ w1$actual.prop.hab,
                                                Sim %in% c("Sim07", "Sim08", "Sim09") ~ w2$actual.prop.hab,
                                                Sim %in% c("Sim10", "Sim11", "Sim12") ~ w3$actual.prop.hab))

ABM.df <- ABM.df %>% mutate(Survival = case_when(Sim %in% c("Sim04", "Sim07", "Sim10") ~ 0.7,
                                                Sim %in% c("Sim05", "Sim08", "Sim11") ~ 0.8,
                                                Sim %in% c("Sim06", "Sim09", "Sim12") ~ 0.9))

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
  ggtitle("Simulations of Fisher Populations (100 Runs)\nBy Proportion of Suitable Habitat and Survival Rate")+
  facet_wrap(~Prophab+Survival)

sim.TS.plot

#- Plot

Cairo(file="out/sim_noescape.TS.plot.PNG",
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

# IBM.w1.surv9.sim100 # Sim06
# IBM.w2.surv9.sim100 # Sim09
# IBM.w3.surv9.sim100 # Sim12

Nozero.runs <- ABM.df %>% filter(TimeStep=="TimeStep_23") %>%
  filter(Sim %in% c("Sim06","Sim09","Sim12")) %>%
  group_by(Sim) %>%
  filter(NewCount!=0)

Sim06_nozero <- Nozero.runs %>% filter(Sim=="Sim06") %>% dplyr::select(Run)
Sim06_nozero <- unique(Sim06_nozero$Run)

Sim09_nozero <- Nozero.runs %>% filter(Sim=="Sim09") %>% dplyr::select(Run)
Sim09_nozero <- unique(Sim09_nozero$Run)

Sim12_nozero <- Nozero.runs %>% filter(Sim=="Sim12") %>% dplyr::select(Run)
Sim12_nozero <- unique(Sim12_nozero$Run)

length(Sim06_nozero); length(Sim09_nozero); length(Sim12_nozero)

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

    r_list[[i]] <- fasterize(ftmp.sfp, r, field="Fisher", fun=rFun, background=0)
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

###--- For Sim06 - # IBM.w1.surv9.sim100 # Sim06

rSim06 <- raster_output(sim_out=IBM_noescape, sim_order=6, sim_use=Sim06_nozero,land=w1$land,
                        TS=23, rExtent=rw1, rFun="sum",sFun="sum")

rSim06
plot(rSim06$raster)

length(Sim06_nozero)
w1$t0
Cairo(file="out/rSim06_title.PNG", type="png", width=2200, height=2000,pointsize=15,bg="white",dpi=300)
plot(rSim06$raster,
     main="Estimated Fisher Abundance over 100 Simulations",
     sub="Starting with 20 fishers and 45% suitable habitat\npredicted 7.8 \u00B1 1.2 (mean \u00B1 1 SE) fishers after 10 years.")
dev.off()

###--- For Sim09 - # IBM.w2.surv9.sim100 # Sim09

rSim09 <- raster_output(sim_out=IBM_noescape, sim_order=9, sim_use=Sim09_nozero,land=w2$land,
                        TS=23, rExtent=rw2, rFun="sum",sFun="sum")

rSim09
plot(rSim09$raster)
w2$actual.prop.hab
length(Sim09_nozero)

Cairo(file="out/rSim09_title.PNG", type="png", width=2200, height=2000,pointsize=15,bg="white",dpi=300)
plot(rSim09$raster,
     main="Estimated Fisher Abundance over 100 Simulations",
     sub="Starting with 20 fishers and 56% suitable habitat\npredicted 18.8 \u00B1 2.3 (mean \u00B1 1 SE) fishers after 10 years.")
dev.off()

###--- For Sim12 - # IBM.w3.surv9.sim100 # Sim12

rSim12 <- raster_output(sim_out=IBM_noescape, sim_order=12, sim_use=Sim12_nozero,land=w3$land,
                        TS=23, rExtent=rw3, rFun="sum",sFun="sum")

rSim12
plot(rSim12$raster)
w3$actual.prop.hab
length(Sim12_nozero)

Cairo(file="out/rSim12_title.PNG", type="png", width=2200, height=2000,pointsize=15,bg="white",dpi=300)
plot(rSim12$raster,
     main="Estimated Fisher Abundance over 100 Simulations",
     sub="Starting with 20 fishers and 71% suitable habitat\nyielded 8.2 \u00B1 1.5 (mean \u00B1 1 SE) fishers after 10 years.")
dev.off()
