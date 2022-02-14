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


source("00_IBM_functions.R")
#####################################################################################

# Start with a very simple example - habitat patch is either suitable or not suitable
# Die by 2 if haven’t found a territory
# Females max age length = 9 (as per discussion with team - see Rory's spreadsheet)
# Female territory = 30 km2 (1 pixel / cell)
# Female transient animal = 70 km2 Euclidian distance within a couple months
# Finding a mate is also part of denning rate
# When a landscape is full of animals an animal has to move quite far – so if territories are all full, then the animal has to leave the study area (one went to Alberta)
# When the landscape is poor and isn’t able to support fishers, then the distance moved is also far

################################################################################

# *** Step 1. START ***
# The assumption is that there is 100% survival during the first year (i.e., the set up), at the first time step no fishers die
# •	t0 = October to April = kits are born; need to run through the reproduce functions
# i.	t0 <- repro(fishers=t0, repro_estimates=repro.CI, Fpop="C")
# •	all fishers age 0.5 years
#
# *** Step 2. AGE ***
# •	The assumption is that there is 100% survival during the first year (i.e., the set up), at the second time step no fishers die
# •	t1 = April to October = kits kicked out of natal territory
# •	all fishers age 0.5 years
#
# *** Step 3. ESTABLISH / MAINTAIN TERRITORY & SCENT TERRITORY (MATE) & SURVIVE ***
# •	 t2 = October to April = females with established territory find mate
# •	 3a = the first step is for individuals without territories to disperse; run through the disperse function up to 30 times to allow 6 months of movement
# i.	t2 <- disperse(land=land, fishers=t2, dist_mov=dist_mov, out=FALSE)
# •	all fishers age 0.5 years
# •	at the end of this time step, all fishers subject to mortality; run through the survive function
# i.	t2 <- survive(fishers=t2, surv_estimates=rf_surv_estimates, Fpop=Fpop, maxAgeFemale=maxAgeFemale)
#
# *** Step 4.  ESTABLISH / MAINTAIN TERRITORY ***
# •	t3 = April to October = keep surviving
# •	4a = the first step is for individuals without territories to disperse; run through the disperse function up to 30 times to allow 6 months of movement
# i.	TOct <- disperse(land=land, fishers=tOct, dist_mov=dist_mov, out=FALSE)
# •	all fishers age 0.5 years
# •	update the fisher table to change juveniles to adults as they age out of (i.e., age > 2)
#
# *** Step 5. ESTABLISH / MAINTAIN TERRITORY & REPRODUCE & SCENT TERRITORY (MATE) & SURVIVE ***
# •	t4 = October to April = females with established territory produce kits
# •	5a = the first step is for pregnant female fishers to reproduce; run through the reproduce denning and kits_produced functions
# i.	tApr <- repro(fishers=t0, repro_estimates=repro.CI, Fpop="C")
# •	5b = the second step is for juvenile fishers without established territories to move; loop through the disperse function up to 30 times
# i.	tApr <- disperse(land=land, fishers=tApr, dist_mov=dist_mov, out=out)
# •	all fishers age 0.5 years
# •	at the end of this time step, all fishers subject to mortality; run through the survive function
# i.	tApr <- survive(fishers=tApr, surv_estimates=rf_surv_estimates, Fpop=Fpop, maxAgeFemale=maxAgeFemale)


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
FEMALE_IBM_simulation_same_world <- function(land=land, t0=t0,                                # import world
                                             repro_estimates=repro.CI, Fpop="C",      # reproduction
                                             surv_estimates=rf_surv_estimates,                # survive
                                             maxAgeFemale=maxAgeFemale,                       # survive
                                             dist_mov=1.0, out=TRUE, torus=TRUE,              # disperse
                                             yrs.to.run=10){                                             # number of years to run simulations ()

  # 2 times steps per year so yrs.to.run*2 plus the initial 3 time steps (start in Apr=t0, Oct=t1, Apr=t2)
  IBM.sim.out <- vector('list', (yrs.to.run*2)+3)

  # *** Step 1. START ***
  # The assumption is that there is 100% survival during the first year (i.e., the set up), at the first time step no fishers die
  # •	t0 = October to April = kits are born; need to run through the reproduce functions
  # i.	t0 <- denning(fishers=t0, denLCI=denLCI, denUCI=denUCI)
  # ii.	t0 <- kits_produced(fishers=t0, ltrM=ltrM, ltrSD=ltrSD)
  # •	all fishers age 0.5 years

  # t0=w1$t0
  t0 <- repro_FEMALE(fishers=t0, repro_estimates=repro_estimates, Fpop=Fpop)

  print(NLcount(t0))
  IBM.sim.out[[1]] <- t0 # time step ends at April

  # *** Step 2. AGE ***
  # •	The assumption is that there is 100% survival during the first year (i.e., the set up), at the second time step no fishers die
  # •	t1 = April to October = kits kicked out of natal territory
  # •	all fishers age 0.5 years

  age.val <- of(agents=t0, var=c("age"))+0.5
  t1 <- NLset(turtles = t0, agents=turtle(t0, who=t0$who),var="age", val=age.val)

  print(NLcount(t1))
  IBM.sim.out[[2]] <- t1 # time step ends at October

  # *** Step 3. ESTABLISH / MAINTAIN TERRITORY & SCENT TERRITORY (MATE) & SURVIVE ***
  # •	 t2 = October to April = females with established territory find mate
  # •	 3a = the first step is for individuals without territories to disperse; run through the disperse function up to 30 times to allow 6 months of movement
  # i.	t2 <- disperse(land=land, fishers=t2, dist_mov=dist_mov, out=FALSE)
  # •	all fishers age 0.5 years
  # •	at the end of this time step, all fishers subject to mortality; run through the survive function
  # i.	t2 <- survive(fishers=t2, surv_estimates=rf_surv_estimates, Fpop=Fpop, maxAgeFemale=maxAgeFemale)

  # land=w1$land; dist_mov=1; out=TRUE; torus=TRUE

  t2 <- t1
  for(i in 1:30){
    t2 <- disperse_FEMALE(land=land, fishers=t2, dist_mov=dist_mov, out=out, torus=torus)
  }

  age.val <- of(agents=t2, var=c("age"))+0.5
  t2 <- NLset(turtles = t2, agents=turtle(t2, who=t2$who),var="age", val=age.val)

  # Fpop="C"; maxAgeFemale=9; surv_estimates=rf_surv_estimates
  t2 <- survive_FEMALE(t2, surv_estimates=surv_estimates, Fpop=Fpop, maxAgeFemale=maxAgeFemale)


  print(NLcount(t2))
  IBM.sim.out[[3]] <- t2 # time step ends at April


  ################################################################################

  tOct <- t2

  for(tcount in 4:(yrs.to.run*2+3)){

    #    # *** Step 4.  ESTABLISH / MAINTAIN TERRITORY ***
    # •	t3 = April to October = keep surviving
    # •	4a = the first step is for individuals without territories to disperse; run through the disperse function up to 30 times to allow 6 months of movement
    # i.	TOct <- disperse(land=land, fishers=tOct, dist_mov=dist_mov, out=FALSE)
    # •	all fishers age 0.5 years
    # •	update the fisher table to change juveniles to adults as they age out of (i.e., age > 2)


    if(NLcount(tOct)!=0){

      for(i in 1:30){
        tOct <- disperse_FEMALE(land=land, fishers=tOct, dist_mov=dist_mov, out=out, torus=torus)
      }

      age.val <- of(agents=tOct, var=c("age"))+0.5
      tOct <- NLset(turtles = tOct, agents=turtle(tOct, who=tOct$who),var="age", val=age.val)

      breed.val <- as.data.frame(of(agents=tOct, var=c("breed","age")))
      breed.val$breed <- case_when(breed.val$age>2 ~ "adult",
                                   TRUE ~ as.character(breed.val$breed))

      tOct <- NLset(turtles = tOct, agents=turtle(tOct, who=tOct$who),var="breed", val=breed.val$breed)

      print(NLcount(tOct))
      IBM.sim.out[[tcount]] <- tOct

      # *** Step 5. ESTABLISH / MAINTAIN TERRITORY & REPRODUCE & SCENT TERRITORY (MATE) & SURVIVE ***
      # •	t4 = October to April = females with established territory produce kits
      # •	5a = the first step is for pregnant female fishers to reproduce; run through the reproduce denning and kits_produced functions
      # i.	tApr <- denning(fishers=tOct, denLCI=denLCI, denUCI=denUCI)
      # ii.	tApr <- kits_produced(fishers=tApr, ltrM=ltrM, ltrSD=ltrSD)
      # •	5b = the second step is for juvenile fishers without established territories to move; loop through the disperse function up to 30 times
      # i.	tApr <- disperse(land=land, fishers=tApr, dist_mov=dist_mov, out=out)
      # •	all fishers age 0.5 years
      # •	at the end of this time step, all fishers subject to mortality; run through the survive function
      # i.	tApr <- survive(fishers=tApr, surv_estimates=rf_surv_estimates, Fpop=Fpop, maxAgeFemale=maxAgeFemale)
      # *** Step 4.  ESTABLISH / MAINTAIN TERRITORY ***
      # t3 = April to October = keep surviving
      # 4a. function DISPERSE - run through DISPERSE function for individuals without territories, up to 30 times to allow 6 months of movement


      # repro_estimates=repro.CI
      tApr <- repro_FEMALE(fishers=tOct, repro_estimates=repro_estimates, Fpop=Fpop)

      for(i in 1:30){
        tApr <- disperse_FEMALE(land=land, fishers=tApr, dist_mov=dist_mov, out=out, torus=torus)
      }

      age.val <- of(agents=tApr, var=c("age"))+0.5
      tApr <- NLset(turtles = tApr, agents=turtle(tApr, who=tApr$who),var="age", val=age.val)

      tApr <- survive_FEMALE(tApr, surv_estimates=surv_estimates, Fpop=Fpop, maxAgeFemale=maxAgeFemale)

      print(NLcount(tApr))

      tcount <- tcount+1
      IBM.sim.out[[tcount]] <- tApr

      tOct <- tApr

    } else {
      IBM.sim.out[[tcount]] <- 0 }
  }

  return(IBM.sim.out)

}


w1 <- set_up_world_FEMALE(nFemales=20, maxAgeFemale=9,xlim=c(1,10), ylim=c(1,10), prophab=0.5)
w1

test <- FEMALE_IBM_simulation_same_world(land=w1$land, t0=w1$t0,                 # import world
                                         repro_estimates=repro.CI, Fpop="C",    # reproduction
                                         surv_estimates=rf_surv_estimates,      # survive
                                         maxAgeFemale=9,                        # survive
                                         dist_mov=1.0, out=TRUE, torus=TRUE,    # disperse
                                         yrs.to.run=10)


################################################################################
# Create 3 sets of 100 simulations - vary the proportion of habitat and survival
# Low, medium and high habitat = 0.5, 0.6, and 0.7 (same world set up, get actual values)
# Low, medium and high survival = 0.7, 0.8, 0.9
# Run 100 simulations for each, save as objects
# Calculate mean # of animals per cell at 10 years for each simulation to produce a heat map

# Create a figure with mean number of animals (+/- SE) for each time step and graph for each simulation
# test.surv <- rf_surv_estimates
# test.surv$L95CL <- test.surv$U95CL <- test.surv$Surv
# test.surv$L95CL <- test.surv$U95CL <- 0.9

rf_surv_estimates

################################################################################
###--- RUN FOR BOREAL
###--- Run with low habitat (prop hab ~ 0.5)
w1 <- set_up_world(nMales=7, nFemales=13, maxAgeMale=6, maxAgeFemale=9,
                   xlim=c(1,10), ylim=c(1,10), prophab=0.5)

B.w1.FEMALE.sim100 <- vector('list',100)
for(i in 1:100){
  B.w1.FEMALE.sim100[[i]] <- FEMALE_IBM_simulation_same_world(land=w1$land, t0=w1$t0,                # import world
                                                              repro_estimates=repro.CI, Fpop="B",    # reproduction
                                                              surv_estimates=rf_surv_estimates,      # survive
                                                              maxAgeFemale=9,                        # survive
                                                              dist_mov=1.0, out=TRUE, torus=TRUE,    # disperse
                                                              yrs.to.run=10)                                             # number of years to run simulation post set up
}

###--- Run with medium habitat (prop hab ~ 0.6)
w2 <- set_up_world(nMales=7, nFemales=13, maxAgeMale=6, maxAgeFemale=9,
                   xlim=c(1,10), ylim=c(1,10), prophab=0.6)

B.w2.FEMALE.sim100 <- vector('list',100)
for(i in 1:100){
  B.w2.FEMALE.sim100[[i]] <- FEMALE_IBM_simulation_same_world(land=w2$land, t0=w2$t0,                # import world
                                                              repro_estimates=repro.CI, Fpop="B",    # reproduction
                                                              surv_estimates=rf_surv_estimates,      # survive
                                                              maxAgeFemale=9,                        # survive
                                                              dist_mov=1.0, out=TRUE, torus=TRUE,    # disperse
                                                              yrs.to.run=10)                                             # number of years to run simulation post set up
}



###--- Run with high habitat (prop hab ~ 0.7)
w3 <- set_up_world(nMales=7, nFemales=13, maxAgeMale=6, maxAgeFemale=9,
                   xlim=c(1,10), ylim=c(1,10), prophab=0.7)

start_time <- Sys.time()
B.w3.FEMALE.sim100 <- vector('list',100)
for(i in 1:100){
  B.w3.FEMALE.sim100[[i]] <- FEMALE_IBM_simulation_same_world(land=w3$land, t0=w3$t0,                # import world
                                                              repro_estimates=repro.CI, Fpop="B",    # reproduction
                                                              surv_estimates=rf_surv_estimates,      # survive
                                                              maxAgeFemale=9,                        # survive
                                                              dist_mov=1.0, out=TRUE, torus=TRUE,    # disperse
                                                              yrs.to.run=10)                                             # number of years to run simulation post set up
}

end_time <- Sys.time()

# end_time - start_time

Boreal_escape_FEMALE <- list(w1, w2, w3,B.w1.FEMALE.sim100,B.w2.FEMALE.sim100,B.w3.FEMALE.sim100)

save(Boreal_escape_FEMALE, file="out/Boreal_escape_FEMALE.RData")

################################################################################
###--- RUN FOR CENTRAL INTERIOR
# use the same worlds as with Boreal

C.w1.FEMALE.sim100 <- vector('list',100)
for(i in 1:100){
  C.w1.FEMALE.sim100[[i]] <- FEMALE_IBM_simulation_same_world(land=w1$land, t0=w1$t0,                # import world
                                                              repro_estimates=repro.CI, Fpop="C",    # reproduction
                                                              surv_estimates=rf_surv_estimates,      # survive
                                                              maxAgeFemale=9,                        # survive
                                                              dist_mov=1.0, out=TRUE, torus=TRUE,    # disperse
                                                              yrs.to.run=10)                                             # number of years to run simulation post set up
}

###--- Run with medium habitat (prop hab ~ 0.6)
C.w2.FEMALE.sim100 <- vector('list',100)
for(i in 1:100){
  C.w2.FEMALE.sim100[[i]] <- FEMALE_IBM_simulation_same_world(land=w2$land, t0=w2$t0,                # import world
                                                              repro_estimates=repro.CI, Fpop="C",    # reproduction
                                                              surv_estimates=rf_surv_estimates,      # survive
                                                              maxAgeFemale=9,                        # survive
                                                              dist_mov=1.0, out=TRUE, torus=TRUE,    # disperse
                                                              yrs.to.run=10)                                             # number of years to run simulation post set up
}


###--- Run with high habitat (prop hab ~ 0.7)
C.w3.FEMALE.sim100 <- vector('list',100)
for(i in 1:100){
  C.w3.FEMALE.sim100[[i]] <- FEMALE_IBM_simulation_same_world(land=w3$land, t0=w3$t0,                # import world
                                                              repro_estimates=repro.CI, Fpop="C",    # reproduction
                                                              surv_estimates=rf_surv_estimates,      # survive
                                                              maxAgeFemale=9,                        # survive
                                                              dist_mov=1.0, out=TRUE, torus=TRUE,    # disperse
                                                              yrs.to.run=10)                                             # number of years to run simulation post set up
}


Columbian_escape_FEMALE <- list(w1, w2, w3,C.w1.FEMALE.sim100,C.w2.FEMALE.sim100,C.w3.FEMALE.sim100)

save(Columbian_escape_FEMALE, file="out/Columbian_escape_FEMALE.RData")
