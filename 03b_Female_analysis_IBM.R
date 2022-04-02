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
# 03b_Female_analysis_IBM.R
# script to simulate Individual Based Model (IBM) spatial data for female fisher only
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 11-Jan-2022
#####################################################################################
version$major
version$minor
R_version <- paste0("R-",version$major,".",version$minor)

.libPaths(paste0("C:/Program Files/R/",R_version,"/library")) # to ensure reading/writing libraries from C drive
tz = Sys.timezone() # specify timezone in BC

# Load Packages
# install.packages("R.methodsS3")
# install.packages("lcmix", repos="http://R-Forge.R-project.org")

list.of.packages <- c("tidyverse", "NetLogoR","nnls","lcmix","MASS","SpaDES.core","SpaDES.tools",
                      "Cairo","PNWColors")

# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# for(i in 1:length(list.of.packages)){
#   print(packageDescription(list.of.packages[i], fields=c("Package","Version")))
# }

# Package: tidyverse
# Version: 1.3.1
#
# Package: NetLogoR
# Version: 0.3.9
#
# Package: nnls
# Version: 1.4
#
# Package: lcmix
# Version: 0.3
#
# Package: MASS
# Version: 7.3-54
#
# Package: SpaDES.core
# Version: 1.0.10
#
# Package: SpaDES.tools
# Version: 0.3.9
#
# Package: Cairo
# Version: 1.5-14
#
# Package: PNWColors
# Version: 0.1.0

source("00b_Female_IBM_functions.R")
#####################################################################################

# Start with a  simple example - habitat patch is either suitable or not suitable
# Die by 2 if haven’t found a territory
# Females max age length = 9 (as per discussion with team - see Rory's spreadsheet)
# Female territory = 30 km2 (1 pixel / cell)
# Female transient animal = 70 km2 Euclidian distance within a couple months
# Finding a mate is also part of denning rate
# When a landscape is full of animals an animal has to move quite far – so if territories are all full, then the animal has to leave the study area (one went to Alberta)
# When the landscape is poor and isn’t able to support fishers, then the distance moved is also far

################################################################################

# t0	April	Kits are born
# 	Number of kits born conditional on denning rate and litter size
# t0 <- repro(fishers=t0, repro_estimates=repro.CI, Fpop="C")

# t1	October	Kits are kicked out of natal territory
# 	Age 1 year (t0 and t1)

# tApril	April	Establish / maintain territory & reproduce & scent territory
# 	Females ≥ 0.5 years without established territories disperse (up to 30 FETAs / 6 months of movement)
# 	Females > 2 years give birth and with established territories give birth; number of kits born conditional on denning rate and litter size
# 	Age 0.5 years

# tOctober	October	Establish / maintain territory & survive
# 	Females ≥ 0.5 years without established territories disperse
# 	Fishers subject to mortality; number who die is conditional on survival estimates, maximum lifespan (for adults and dispersing juveniles)
# 	Age 0.5 years

################################################################################
# data to read in; already prepped / formatted in 00_surv_repro_estimates_prep.R

# reproductive rates from manuscript
repro.CI <- read.csv("data/repro.CI.csv", header=TRUE)

# survival probability estimates
# taken from Rory's updated survival, trapping mortality excluded
rf_surv_estimates <- read.csv("data/rf_surv_estimates.csv", header=TRUE)


################################################################################

# create function to loop through functions, allow sub-function specification
# now that the function is using the cohort survival data, have the survival run on an annual basis, not per time step
FEMALE_IBM_simulation_same_world <- function(land=land, t0=t0,                        # import world
                                             repro_estimates=repro.CI, Fpop=Fpop,      # reproduction
                                             surv_estimates=rf_surv_estimates,        # survive
                                             maxAgeFemale=maxAgeFemale,               # survive
                                             dist_mov=1.0, out=TRUE, torus=TRUE,      # disperse
                                             yrs.to.run=10){                          # number of years to run simulations ()

  IBM.sim.out <- vector('list', yrs.to.run+1)

  # land = canBex.FEMALE.world[[1]]$land
  # t0 = canBex.FEMALE.world[[1]]$t0

  # t0	April	Kits are born
  # 	Number of kits born conditional on denning rate and litter size
  # t0 <- repro(fishers=t0, repro_estimates=repro.CI, Fpop="C")

  # t0=tmp$t0; land=tmp$land
  t0 <- repro_FEMALE(fishers=t0, repro_estimates=repro_estimates, Fpop=Fpop)

  # t1	October	Kits are kicked out of natal territory
  # 	Age 1 year (t0 and t1)

  age.val <- of(agents=t0, var=c("age"))+1
  t1 <- NLset(turtles = t0, agents=turtle(t0, who=t0$who),var="age", val=age.val)

  # plot(land)
  # points(t1, pch = t1$shape, col = of(agents = t1, var = "color"))

  print(NLcount(t1))
  IBM.sim.out[[1]] <- t1 # time step ends at April; gives the number of fishers at the start

  # tmp <- patchHere(land, t1[t1$disperse=="E"])
  # duplicated(tmp)

  tApr <- t1

  for(tcount in 2:(yrs.to.run+1)){

    # tApril	April	Establish / maintain territory & reproduce & scent territory
    # 	Females ≥ 0.5 years without established territories disperse (up to 30 FETAs / 6 months of movement)
    # 	Females > 2 years give birth and with established territories give birth; number of kits born conditional on denning rate and litter size
    # 	Age 0.5 years

    if(NLcount(tApr)!=0){

      tApr <- repro_FEMALE(fishers=tApr, repro_estimates=repro_estimates, Fpop=Fpop)

      for(i in 1:30){
        tApr <- disperse_FEMALE(land=land, fishers=tApr, dist_mov=dist_mov, out=out, torus=torus)
      }

      # tmp <- patchHere(land, tApr[tApr$disperse=="E"])
      # duplicated(tmp)

      age.val <- of(agents=tApr, var=c("age"))+0.5
      tApr <- NLset(turtles = tApr, agents=turtle(tApr, who=tApr$who),var="age", val=age.val)

    # plot(land)
    # points(tApr, pch = tApr$shape, col = of(agents = tApr, var = "color"))
    }

    # tOctober	October	Establish / maintain territory & survive
    # 	Females ≥ 0.5 years without established territories disperse
    # 	Fishers subject to mortality; number who die is conditional on survival estimates, maximum lifespan (for adults and dispersing juveniles)
    # 	Age 0.5 years
    tOct <- tApr

    if(NLcount(tOct)!=0){

      for(i in 1:30){
        tOct <- disperse_FEMALE(land=land, fishers=tOct, dist_mov=dist_mov, out=out, torus=torus)
      }

      # tmp <- patchHere(land, tOct)
      # duplicated(tmp)

      age.val <- of(agents=tOct, var=c("age"))+0.5
      tOct <- NLset(turtles = tOct, agents=turtle(tOct, who=tOct$who),var="age", val=age.val)

      breed.val <- as.data.frame(of(agents=tOct, var=c("breed","age")))
      breed.val$breed <- case_when(breed.val$age>2 ~ "adult",
                                   TRUE ~ as.character(breed.val$breed))

      tOct <- NLset(turtles = tOct, agents=turtle(tOct, who=tOct$who),var="breed", val=breed.val$breed)

      tOct <- survive_FEMALE(fishers=tOct, surv_estimates=surv_estimates, Fpop=Fpop, maxAgeFemale=maxAgeFemale)

      IBM.sim.out[[tcount]] <- tOct
    # plot(land)
    # points(tOct, pch = tOct$shape, col = of(agents = tOct, var = "color"))

    } else {
        IBM.sim.out[[tcount]] <- 0 }

    print(NLcount(tOct))
    tApr <- tOct
  }
   return(IBM.sim.out)
}


######################################################################################
# Run a canned example with real world habitat
# Ex 1 = Columbian pop, 4x4 grid, 16 potential fisher equivalent territories (Quesnel)
# Ex 2 = Columbian pop, 5x8 grid, 40 potential fisher equivalent territories (Quesnel)
# Run 100 simulations for each, save as objects
######################################################################################

# Read in aoi shapefile (with attribute info) and aoi raster (binary 0/1 for habitat)
# load("data/IBM_aoi_Pex2.RData")
load("data/IBM_aoi_canBex.RData")
# load("data/IBM_aoi_canCex.RData")

# glimpse(IBM_aoi)
Fpop <- unique(substr(IBM_aoi$aoi$Fpop,1,1))
Fpop = "C" # for using the Boreal example data with Columbian demographic rates

# print(sum(IBM_aoi$canCex1_raster@data@values)) # number of equivalent territories with suitable habitat
# 7
# length(IBM_aoi$canCex1_raster@data@values)
# 156

for(i in 1:length(IBM_aoi$canBex_raster)){
print(sum(IBM_aoi$canBex_raster[[i]]@data@values)) # number of equivalent territories with suitable habitat
}
length(IBM_aoi$canBex_raster[[1]]@data@values)

# go with ~1/3 for actual female territories to start with, rounded to nearest 5
nFemales = plyr::round_any(sum(IBM_aoi$canBex_raster[[1]]@data@values)*0.3,5); nFemales
nFemales = round(sum(IBM_aoi$canCex1_raster@data@values)*0.3); nFemales
#35 for Boreal; 2 for Columbian

canBex.FEMALE.world <- list()
for(i in 1:length(IBM_aoi$canBex_raster)){
  canBex.FEMALE.world[[i]] <- set_up_REAL_world_FEMALE(nFemales=nFemales, maxAgeFemale=9, raoi=IBM_aoi$canBex_raster[[i]])
}

# Scenario1
start_time <- Sys.time()
canBex1.FEMALE.sim100 <- vector('list',100)
  for(i in 1:100){
    canBex1.FEMALE.sim100[[i]] <- FEMALE_IBM_simulation_same_world(land=canBex.FEMALE.world[[1]]$land,
                                                                  t0=canBex.FEMALE.world[[1]]$t0,             # import world
                                                                  repro_estimates=repro.CI, Fpop=Fpop,   # reproduction
                                                                  surv_estimates=rf_surv_estimates,      # survive
                                                                  maxAgeFemale=9,                        # survive
                                                                  dist_mov=1.0, out=TRUE, torus=TRUE,    # disperse
                                                                  yrs.to.run=10)                         # number of years to run simulation post set up
}
end_time <- Sys.time(); print(end_time - start_time)

canBex1.FEMALE <- list(canBex.FEMALE.world[[1]], canBex1.FEMALE.sim100)

# tmp <- canBex1.FEMALE.sim100[[1]][[11]]
# canBex1.FEMALE.sim100[[1]]
# plot(canBex.FEMALE.world[[1]]$land)
# points(tmp[tmp$breed=="adult"], pch = tmp[tmp$breed=="adult"]$shape, col = "black")
# points(tmp[tmp$breed=="juvenile"], pch = tmp[tmp$breed=="adult"]$shape, col = "blue")
#
# plot(canBex.FEMALE.world[[1]]$land)
# points(tmp[tmp$disperse=="E"], pch = tmp[tmp$disperse=="E"]$shape, col = "black")
# points(tmp[tmp$disperse=="D"], pch = tmp[tmp$disperse=="D"]$shape, col = "blue")
#
# tmp <- as.data.frame(patchHere(canBex.FEMALE.world[[1]]$land, tmp[tmp$disperse=="E"]))
# tmp %>% dplyr::arrange(pxcor)
# duplicated(tmp)

# Scenario2
start_time <- Sys.time()
canBex2.FEMALE.sim100 <- vector('list',100)
for(i in 1:100){
  canBex2.FEMALE.sim100[[i]] <- FEMALE_IBM_simulation_same_world(land=canBex.FEMALE.world[[2]]$land,
                                                                 t0=canBex.FEMALE.world[[2]]$t0,             # import world
                                                                 repro_estimates=repro.CI, Fpop=Fpop,   # reproduction
                                                                 surv_estimates=rf_surv_estimates,      # survive
                                                                 maxAgeFemale=9,                        # survive
                                                                 dist_mov=1.0, out=TRUE, torus=TRUE,    # disperse
                                                                 yrs.to.run=10)                         # number of years to run simulation post set up
}
end_time <- Sys.time(); print(end_time - start_time)

canBex2.FEMALE <- list(canBex.FEMALE.world[[2]], canBex2.FEMALE.sim100)

# Scenario3
start_time <- Sys.time()
canBex3.FEMALE.sim100 <- vector('list',100)
for(i in 1:100){
  canBex3.FEMALE.sim100[[i]] <- FEMALE_IBM_simulation_same_world(land=canBex.FEMALE.world[[3]]$land,
                                                                 t0=canBex.FEMALE.world[[3]]$t0,             # import world
                                                                 repro_estimates=repro.CI, Fpop=Fpop,   # reproduction
                                                                 surv_estimates=rf_surv_estimates,      # survive
                                                                 maxAgeFemale=9,                        # survive
                                                                 dist_mov=1.0, out=TRUE, torus=TRUE,    # disperse
                                                                 yrs.to.run=10)                         # number of years to run simulation post set up
}
end_time <- Sys.time(); print(end_time - start_time)

canBex3.FEMALE <- list(canBex.FEMALE.world[[3]], canBex3.FEMALE.sim100)

# Scenario4
start_time <- Sys.time()
canBex4.FEMALE.sim100 <- vector('list',100)
for(i in 1:100){
  canBex4.FEMALE.sim100[[i]] <- FEMALE_IBM_simulation_same_world(land=canBex.FEMALE.world[[4]]$land,
                                                                 t0=canBex.FEMALE.world[[4]]$t0,             # import world
                                                                 repro_estimates=repro.CI, Fpop=Fpop,   # reproduction
                                                                 surv_estimates=rf_surv_estimates,      # survive
                                                                 maxAgeFemale=9,                        # survive
                                                                 dist_mov=1.0, out=TRUE, torus=TRUE,    # disperse
                                                                 yrs.to.run=10)                         # number of years to run simulation post set up
}
end_time <- Sys.time(); print(end_time - start_time)

canBex4.FEMALE <- list(canBex.FEMALE.world[[4]], canBex4.FEMALE.sim100)

canBex.FEMALE<- list(canBex1.FEMALE,canBex2.FEMALE,canBex3.FEMALE,canBex4.FEMALE)
save(canBex.FEMALE, file="out/canBex.FEMALE.RData")


#################################################################################
# For columbian
canCex.FEMALE.world <- set_up_REAL_world_FEMALE(nFemales=nFemales, maxAgeFemale=9, raoi=IBM_aoi$canCex1_raster)

# Scenario1
start_time <- Sys.time()
canCex1.FEMALE.sim100 <- vector('list',100)
for(i in 1:100){
  canCex1.FEMALE.sim100[[i]] <- FEMALE_IBM_simulation_same_world(land=canCex.FEMALE.world$land,
                                                                 t0=canCex.FEMALE.world$t0,             # import world
                                                                 repro_estimates=repro.CI, Fpop=Fpop,   # reproduction
                                                                 surv_estimates=rf_surv_estimates,      # survive
                                                                 maxAgeFemale=9,                        # survive
                                                                 dist_mov=1.0, out=TRUE, torus=TRUE,    # disperse
                                                                 yrs.to.run=10)                         # number of years to run simulation post set up
}
end_time <- Sys.time(); print(end_time - start_time)

canCex1.FEMALE <- list(canCex.FEMALE.world, canCex1.FEMALE.sim100)
save(canCex1.FEMALE, file="out/canCex.FEMALE_actual.RData") # all went to 0 well before 10 years


### For illustration purposes, use the Boreal landscape and Columbian demographics
canCex.FEMALE.world <- list()
for(i in 1:length(IBM_aoi$canBex_raster)){
  canCex.FEMALE.world[[i]] <- set_up_REAL_world_FEMALE(nFemales=nFemales, maxAgeFemale=9, raoi=IBM_aoi$canBex_raster[[i]])
}

# Scenario1
start_time <- Sys.time()
canCex1.FEMALE.sim100 <- vector('list',100)
for(i in 1:100){
  canCex1.FEMALE.sim100[[i]] <- FEMALE_IBM_simulation_same_world(land=canCex.FEMALE.world[[1]]$land,
                                                                 t0=canCex.FEMALE.world[[1]]$t0,             # import world
                                                                 repro_estimates=repro.CI, Fpop="C",   # reproduction
                                                                 surv_estimates=rf_surv_estimates,      # survive
                                                                 maxAgeFemale=9,                        # survive
                                                                 dist_mov=1.0, out=TRUE, torus=TRUE,    # disperse
                                                                 yrs.to.run=10)                         # number of years to run simulation post set up
}
end_time <- Sys.time(); print(end_time - start_time)
canCex1.FEMALE <- list(canCex.FEMALE.world[[1]], canCex1.FEMALE.sim100)


# Scenario2
start_time <- Sys.time()
canCex2.FEMALE.sim100 <- vector('list',100)
for(i in 1:100){
  canCex2.FEMALE.sim100[[i]] <- FEMALE_IBM_simulation_same_world(land=canCex.FEMALE.world[[2]]$land,
                                                                 t0=canCex.FEMALE.world[[2]]$t0,             # import world
                                                                 repro_estimates=repro.CI, Fpop="C",   # reproduction
                                                                 surv_estimates=rf_surv_estimates,      # survive
                                                                 maxAgeFemale=9,                        # survive
                                                                 dist_mov=1.0, out=TRUE, torus=TRUE,    # disperse
                                                                 yrs.to.run=10)                         # number of years to run simulation post set up
}
end_time <- Sys.time(); print(end_time - start_time)

canCex2.FEMALE <- list(canCex.FEMALE.world[[2]], canCex2.FEMALE.sim100)

# Scenario3
start_time <- Sys.time()
canCex3.FEMALE.sim100 <- vector('list',100)
for(i in 1:100){
  canCex3.FEMALE.sim100[[i]] <- FEMALE_IBM_simulation_same_world(land=canCex.FEMALE.world[[3]]$land,
                                                                 t0=canCex.FEMALE.world[[3]]$t0,             # import world
                                                                 repro_estimates=repro.CI, Fpop="C",   # reproduction
                                                                 surv_estimates=rf_surv_estimates,      # survive
                                                                 maxAgeFemale=9,                        # survive
                                                                 dist_mov=1.0, out=TRUE, torus=TRUE,    # disperse
                                                                 yrs.to.run=10)                         # number of years to run simulation post set up
}
end_time <- Sys.time(); print(end_time - start_time)

canCex3.FEMALE <- list(canCex.FEMALE.world[[3]], canCex3.FEMALE.sim100)

# Scenario4
start_time <- Sys.time()
canCex4.FEMALE.sim100 <- vector('list',100)
for(i in 1:100){
  canCex4.FEMALE.sim100[[i]] <- FEMALE_IBM_simulation_same_world(land=canCex.FEMALE.world[[4]]$land,
                                                                 t0=canCex.FEMALE.world[[4]]$t0,             # import world
                                                                 repro_estimates=repro.CI, Fpop="C",   # reproduction
                                                                 surv_estimates=rf_surv_estimates,      # survive
                                                                 maxAgeFemale=9,                        # survive
                                                                 dist_mov=1.0, out=TRUE, torus=TRUE,    # disperse
                                                                 yrs.to.run=10)                         # number of years to run simulation post set up
}
end_time <- Sys.time(); print(end_time - start_time)

canCex4.FEMALE <- list(canCex.FEMALE.world[[4]], canCex4.FEMALE.sim100)

canCex.FEMALE<- list(canCex1.FEMALE,canCex2.FEMALE,canCex3.FEMALE,canCex4.FEMALE)
save(canCex.FEMALE, file="out/canCex.FEMALE.RData")

