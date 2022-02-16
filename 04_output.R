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
# version$major
# version$minor
# R_version <- paste0("R-",version$major,".",version$minor)
# 
# .libPaths(paste0("C:/Program Files/R/",R_version,"/library")) # to ensure reading/writing libraries from C drive
# tz = Sys.timezone() # specify timezone in BC
# 
# # Load Packages
# list.of.packages <- c("tidyverse", "NetLogoR","nnls","lcmix","MASS","Cairo","PNWColors", "ggplot2",
#                       "sf","raster","rgdal")
# # Check you have them and load them
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)
# lapply(list.of.packages, require, character.only = TRUE)
# 
# source("00b_Female_IBM_functions.R")
#####################################################################################
# Create 3 sets of 100 simulations - vary the proportion of habitat and survival
# Low, medium and high habitat = 0.5, 0.6, and 0.7 (same world set up, get actual values)
# Low, medium and high survival = 0.7, 0.8, 0.9


load("out/Boreal_escape_FEMALE_bern.RData")
w1 <- Boreal_escape_FEMALE_bern[[1]]; w1$actual.prop.hab # 0.52
w2 <- Boreal_escape_FEMALE_bern[[2]]; w2$actual.prop.hab # 0.65
w3 <- Boreal_escape_FEMALE_bern[[3]]; w3$actual.prop.hab # 0.69

load("out/Columbian_escape_FEMALE_bern.RData")

###--- plot the simulated landbases
Cairo(file="outputs/BCI_Fescape_w1.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
plot(w1$land, main=c(paste0("Simulated Landbase"),paste0(w1$actual.prop.hab*100,"% Suitable Habitat")))
points(w1$t0, pch = w1$t0$shape, col = of(agents = w1$t0, var = "color"))
dev.off()

Cairo(file="outputs/BCI_Fescape_w2.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
plot(w2$land, main=c(paste0("Simulated Landbase"),paste0(w2$actual.prop.hab*100,"% Suitable Habitat")))
points(w2$t0, pch = w2$t0$shape, col = of(agents = w2$t0, var = "color"))
dev.off()

Cairo(file="outputs/BCI_Fescape_w3.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
plot(w3$land, main=c(paste0("Simulated Landbase"),paste0(w3$actual.prop.hab*100,"% Suitable Habitat")))
points(w3$t0, pch = w3$t0$shape, col = of(agents = w3$t0, var = "color"))
dev.off()

Cairo(file="outputs/BCI_Fescape_w3_nofisher.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
plot(w3$land, main=c(paste0("Simulated Landbase"),paste0(w3$actual.prop.hab*100,"% Suitable Habitat")))
dev.off()


# Run 100 simulations for each, save as objects
# Calculate mean # of animals per cell at 10 years for each simulation to produce a heat map
# Create a figure with mean number of animals (+/- SE) for each time step and graph for each simulation

sim_output <- function(sim_out=sim_out, sim_order=sim_order, numsims=numsims, yrs_sim=yrs_sim){
  num.runs <- yrs_sim + 2
  ABM.df <- as.data.frame(array(NA,c(2*numsims,num.runs)))

  colnames(ABM.df) <- paste0("TimeStep_",str_pad(seq_len(num.runs),2,pad="0"))
  for(i in 1:numsims){
    ABM.df[i,] <- unlist(lapply(lapply(sim_out[[sim_order]][[i]], as.array), ncol)) # if 14 then know that has at least one fisher
    ABM.df[i+numsims,] <- unlist(lapply(lapply(sim_out[[sim_order]][[i]], as.array), nrow)) # if 14 then know that has at least one fisher
  }

  ABM.df$Type <- rep(c("Pfisher","Count"), each=numsims)
  ABM.df$Run <- rep(seq_len(numsims), times=2)

  ABM.df <- ABM.df %>% pivot_longer(cols = TimeStep_01:TimeStep_12,names_to = "TimeStep",values_to = "Value" )
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
B_ABM.df <- C_ABM.df <- as.data.frame(array(NA,c(3600,4)))
colnames(B_ABM.df) <- colnames(C_ABM.df) <- c("Run","TimeStep","NewCount","Sim")

# loop to put in all of the values
# Columbian
# starting point of data frame
a=1
b=1200
for(i in 4:6){
  C_ABM.df[a:b,] <- sim_output(sim_out=Columbian_escape_FEMALE_bern, sim_order=i, numsims=100, yrs_sim=10)
  a=a+1200
  b=b+1200
}

# Boreal
# starting point of data frame
a=1
b=1200
for(i in 4:6){
  B_ABM.df[a:b,] <- sim_output(sim_out=Boreal_escape_FEMALE_bern, sim_order=i, numsims=100, yrs_sim=10)
  a=a+1200
  b=b+1200
}

C_ABM.df$Pop <- "Columbian"
B_ABM.df$Pop <- "Boreal"

###---
ABM.df <- rbind(C_ABM.df, B_ABM.df)

ABM.df <- ABM.df %>% mutate(Prophab = case_when(Sim %in% c("Sim04") ~ w1$actual.prop.hab,
                                                Sim %in% c("Sim05") ~ w2$actual.prop.hab,
                                                Sim %in% c("Sim06") ~ w3$actual.prop.hab))

ABM.TS.mean <- ABM.df %>% dplyr::select(-Run) %>% pivot_wider(names_from=TimeStep, values_from=NewCount, values_fn=mean)
ABM.TS.mean$Param <- "Mean"

se <- function(x) sqrt(var(x)/length(x))
ABM.TS.se <- ABM.df %>% dplyr::select(-Run) %>% pivot_wider(names_from=TimeStep, values_from=NewCount, values_fn=se)
ABM.TS.se$Param <- "SE"

LCL <- function(x) quantile(x, probs=0.05)
ABM.TS.LCL <- ABM.df %>% dplyr::select(-Run) %>% pivot_wider(names_from=TimeStep, values_from=NewCount, values_fn=LCL)
ABM.TS.LCL$Param <- "LCL"

UCL <- function(x) quantile(x, probs=0.95)
ABM.TS.UCL <- ABM.df %>% dplyr::select(-Run) %>% pivot_wider(names_from=TimeStep, values_from=NewCount, values_fn=UCL)
ABM.TS.UCL$Param <- "UCL"

ABM.TS <- rbind(ABM.TS.mean, ABM.TS.se, ABM.TS.LCL, ABM.TS.UCL)

ABM.TS.df <- ABM.TS %>% pivot_longer(cols = TimeStep_01:TimeStep_12,names_to = "TimeStep",values_to = "Value" )
ABM.TS.df <- ABM.TS.df %>% pivot_wider(names_from = Param, values_from = Value)

ABM.TS.use <- ABM.TS.df %>% filter(!TimeStep %in% c("TimeStep_01", "TimeStep_02"))

ABM.TS.use$TimeStepNum <- as.numeric(substr(ABM.TS.use$TimeStep,10,11))

pal_col <- pnw_palette(name="Starfish",n=7,type="discrete")

sim.TS.plot <- ggplot(data = ABM.TS.use) +
  theme_bw() + theme(strip.background = element_rect(fill = "white", colour = "white")) +
  theme(panel.grid = element_blank())+
  geom_ribbon(aes(x = TimeStepNum, ymin = LCL, ymax = UCL), fill = "#2c6184") +
  geom_vline(xintercept = 7, col="darkgrey", lty=4) +
  geom_line(aes(x = TimeStepNum, y = Mean)) +
  # geom_errorbar(aes(x = TimeStepNum, y = Mean, ymin=Mean-SE, ymax= Mean+SE),
  #               width=.2, position=position_dodge(0.05)) +
  theme(axis.text.x = element_blank()) +
  xlab("Time Step in 6 Month Intervals over 10 years") +
  ylab("Number of Fishers Alive (Mean + 95% Confidence Intervals)")+
  ggtitle("Simulations of Fisher Populations (100 Runs)\nBy Population and Proportion of Suitable Habitat")+
  facet_wrap(~Pop+Prophab)

sim.TS.plot

#- Plot

Cairo(file="out/BCI_sim_escape_FEMALE_bern.TS.plot_CL.PNG",
      type="png",
      width=3000,
      height=2200,
      pointsize=15,
      bg="white",
      dpi=300)
sim.TS.plot
dev.off()

sim.TS.plot_se <- ggplot(data = ABM.TS.use) +
  theme_bw() + theme(strip.background = element_rect(fill = "white", colour = "white")) +
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = "TimeStep_07", col="grey", lty=4) +
  geom_point(aes(x = TimeStep, y = Mean), size=2) +
  geom_errorbar(aes(x = TimeStep, y = Mean, ymin=Mean-SE, ymax= Mean+SE),
                width=.2, position=position_dodge(0.05)) +
  theme(axis.text.x = element_blank()) +
  xlab("Time Step in 6 Month Intervals over 10 years") +
  ylab("Number of Fishers Alive (Mean \u00B1 1 SE)")+ # \u00B1 is ± in unicode
  ggtitle("Simulations of Fisher Populations (100 Runs)\nBy Population and Proportion of Suitable Habitat")+
  facet_wrap(~Pop+Prophab)

sim.TS.plot_se

#- Plot
Cairo(file="out/BCI_sim_escape_FEMALE_bern.TS.plot_SE.PNG",type="png",width=3000,height=2200,pointsize=15,bg="white",dpi=300)
sim.TS.plot_se
dev.off()

################################################################################

### Create heatmaps for the w3 outputs
# keep in mind that WGS84 lat/long espg = 4326; BC Albers espg = 3005; NAD83 / UTM zone 10N espg = 26910

Nozero.runs <- ABM.df %>% filter(TimeStep=="TimeStep_12") %>%
  group_by(Sim) %>%
  filter(NewCount!=0)

nozerosims <- function(sim=Sim, pop=Pop){

  nozerosims <- Nozero.runs %>% filter(Sim==sim & Pop==pop) %>% dplyr::select(Run)
  nozerosims <- unique(nozerosims$Run)

  return(nozerosims)
}

BL_nozero <- nozerosims(sim="Sim04", pop="Boreal")
BM_nozero <- nozerosims(sim="Sim05", pop="Boreal")
BH_nozero <- nozerosims(sim="Sim06", pop="Boreal")

CL_nozero <- nozerosims(sim="Sim04", pop="Columbian")
CM_nozero <- nozerosims(sim="Sim05", pop="Columbian")
CH_nozero <- nozerosims(sim="Sim06", pop="Columbian")


length(BL_nozero); length(BM_nozero); length(BH_nozero) # all reach 10 years
length(CL_nozero); length(CM_nozero); length(CH_nozero) # 1,2,1 reached 10 years


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
  writeRaster(r_stackApply, file=paste0("outputs/rSim",str_pad(sim_order,2,pad="0"),".tif"), bylayer=TRUE, overwrite=TRUE)

  Fisher_Nmean <- mean(r_stackApply@data@values)
  Fisher_Nse <- se(r_stackApply@data@values)

  return(list(raster=r_stackApply, Fisher_Nmean=Fisher_Nmean, Fisher_Nse=Fisher_Nse))

}

###--- For Sim04 in Boreal
# lowest suitable habitat
rBL <- raster_output(sim_out=Boreal_escape_FEMALE_bern, sim_order=4, sim_use=BL_nozero,land=w1$land,
                        TS=12, rExtent=rw1, rFun="sum",sFun="sum")

rBL
plot(rBL$raster)

length(BL_nozero) #
# w1$t0
plot(rBL$raster, oma=c(2, 3, 5, 2))
mytitle = "Estimated Fisher Abundance over 100 Simulations"
mysubtitle1 = paste0("Starting with 20 fishers and ",w1$actual.prop.hab*100,"% suitable habitat")
mysubtitle2 = paste0("predicted ",round(rBL$Fisher_Nmean)," \u00B1 ",round(rBL$Fisher_Nse)," (mean \u00B1 1 SE) fishers after 10 years.")
mtext(side=3, line=3, at=-0.07, adj=0, cex=1, mytitle)
mtext(side=3, line=2, at=-0.07, adj=0, cex=0.8, mysubtitle1)
mtext(side=3, line=1, at=-0.07, adj=0, cex=0.8, mysubtitle2)
Cairo(file="out/rB_escape_F_bern_L.PNG", type="png", width=2200, height=2000,pointsize=15,bg="white",dpi=300)
dev.off()

###--- For Sim05 in Boreal
rBM <- raster_output(sim_out=Boreal_escape_FEMALE_bern, sim_order=5, sim_use=BM_nozero,land=w2$land,
                        TS=12, rExtent=rw2, rFun="sum",sFun="sum")

rBM
plot(rBM$raster)

length(BM_nozero) # 55

plot(rBM$raster, oma=c(2, 3, 5, 2))
mytitle = "Estimated Fisher Abundance over 100 Simulations"
mysubtitle1 = paste0("Starting with 20 fishers and ",w2$actual.prop.hab*100,"% suitable habitat")
mysubtitle2 = paste0("predicted ",round(rBM$Fisher_Nmean)," \u00B1 ",round(rBM$Fisher_Nse)," (mean \u00B1 1 SE) fishers after 10 years.")
mtext(side=3, line=3, at=-0.07, adj=0, cex=1, mytitle)
mtext(side=3, line=2, at=-0.07, adj=0, cex=0.8, mysubtitle1)
mtext(side=3, line=1, at=-0.07, adj=0, cex=0.8, mysubtitle2)
Cairo(file="out/rB_escape_F_bern_M.PNG", type="png", width=2200, height=2000,pointsize=15,bg="white",dpi=300)
dev.off()

###--- For Sim06 in Boreal
rBH <- raster_output(sim_out=Boreal_escape_FEMALE_bern, sim_order=6, sim_use=BH_nozero,land=w3$land,
                        TS=12, rExtent=rw3, rFun="sum",sFun="sum")

rBH
plot(rBH$raster)

length(BH_nozero) # 84
# w3$t0
Cairo(file="out/rB_escape_F_bern_H.PNG", type="png", width=2200, height=2000,pointsize=15,bg="white",dpi=300)
plot(rBH$raster, oma=c(2, 3, 5, 2))
mytitle = "Estimated Fisher Abundance over 100 Simulations"
mysubtitle1 = paste0("Starting with 20 fishers and ",w2$actual.prop.hab*100,"% suitable habitat")
mysubtitle2 = paste0("predicted ",round(rBH$Fisher_Nmean)," \u00B1 ",round(rBH$Fisher_Nse)," (mean \u00B1 1 SE) fishers after 10 years.")
mtext(side=3, line=3, at=-0.07, adj=0, cex=1, mytitle)
mtext(side=3, line=2, at=-0.07, adj=0, cex=0.8, mysubtitle1)
mtext(side=3, line=1, at=-0.07, adj=0, cex=0.8, mysubtitle2)
dev.off()
