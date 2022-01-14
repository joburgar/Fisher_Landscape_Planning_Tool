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
# 03_analysis_IBM.R
# script to simulate Individual Based Model (IBM) spatial data for fisher
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 09-Dec-2021
#####################################################################################
version$major
version$minor
R_version <- paste0("R-",version$major,".",version$minor)

.libPaths(paste0("C:/Program Files/R/",R_version,"/library")) # to ensure reading/writing libraries from C drive
tz = Sys.timezone() # specify timezone in BC

# Load Packages
list.of.packages <- c("tidyverse", "NetLogoR","nnls","lcmix","MASS","SpaDES.core","SpaDES.tools",
                      "Cairo","PNWColors")
# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

source("00_IBM_functions.R")
#####################################################################################

# Start with a very simple example - habitat patch is either suitable or not suitable
# Die by 2 if haven’t found a territory
# Females max age length = 8 (might instead go with survival analysis data)
# Males max age length = 4 (may also go with survival analysis data)
# Female territory = 30 km2 (1 pixel / cell)
# Male territory = 5*female territory (have 3 female territories) 150-200 km2
# Female transient animal = 70 km2 Euclidian distance within a couple months
# Male transient animal = moves much farther
# When a landscape is full of animals an animal has to move quite far – so if territories are all full, then the animal has to leave the study area (one went to Alberta)
# When the landscape is poor and isn’t able to support fishers, then the distance moved is also far

################################################################################
# *** Step 1. START ***
# t0 = April
# function FIND_MATE - assign status of nearby mates to females (based on distance of nearby males)
# function DENNING - assign 1 or 0 if female dens (based on nearby mates and published probabilities)
# function KITS_PRODUCED - kits produced (based on two preceding functions and published probabilities)
# simple scenario is to start with equal numbers of adult males and females, all on quality habitat with established territories

# *** Step 2. SURVIVE ***
# t1 = October
# function SURVIVE - add 0.5 to all fishers, kill off individuals who do not survive through t1

# *** Step 3. ESTABLISH / MAINTAIN TERRITORY & SCENT TERRITORY (MATE) & SURVIVE ***
# t2 = April
# 3a. function DISPERSE - run through DISPERSE function for individuals without territories, up to 30 times to allow 6 months of movement
# 3b. function FIND_MATE - for female fishers with ESTABLISHED territory, if male is within 2 cells in either direction or 8 adjacent cells plus same cell, assign mated status (i.e., if male is in same cell or ± 1 cell either via xlim and/or ylim)
# 3c. function SURVIVE - add 0.5 to all fishers, kill off individuals who do not survive this 6 month time step

# *** Step 4.  ESTABLISH / MAINTAIN TERRITORY & SURVIVE ***
# t3 = October
# 4a. function DISPERSE - run through DISPERSE function for individuals without territories, up to 30 times to allow 6 months of movement
# 4b. function SURVIVE - add 0.5 to all fishers, kill off individuals who do not survive through this 6 month time step

# *** Step 5. ESTABLISH / MAINTAIN TERRITORY & REPRODUCE & SCENT TERRITORY (MATE) & SURVIVE ***
# t4 = April
# 5a. function DENNING - assign 1 or 0 if female dens (based on nearby mates and published probabilities)
# 5b. function KITS_PRODUCED - kits produced (based on denning rate and published probabilities)

# 5c. function DISPERSE - run through DISPERSE function for individuals without territories, up to 30 times to allow 6 months of movement
# 5d. function FIND_MATE - for female fishers with ESTABLISHED territory, if male is within 2 cells in either direction or 8 adjacent cells plus same cell, assign mated status (i.e., if male is in same cell or ± 1 cell either via xlim and/or ylim)
# 5e. function SURVIVE - add 0.5 to all fishers, kill off individuals who do not survive this 6 month time step

# *** Step 6. LOOP THROUGH Steps 4 and 5 for X number of years ***
# 6a. Print or save each tn to keep details of population over time - create a list and have it populated by each output (normal simulation stuff)

################################################################################
# data to read in; already prepped / formatted in 00_surv_repro_estimates_prep.R

# reproductive rates in Rich and Eric's paper -
repro.CI <- read.csv("data/repro.CI.csv", header=TRUE, row.names = 1)

# survival probability estimates
km_surv_estimates <- read.csv("data/km_surv_estimates.csv", header=TRUE)

################################################################################
# *** Step 1. START ***
# t0 = April
# function FIND_MATE - assign status of nearby mates to females (based on distance of nearby males)
# function DENNING - assign 1 or 0 if female dens (based on nearby mates and published probabilities)
# function KITS_PRODUCED - kits produced (based on two preceding functions and published probabilities)
# simple scenario is to start with equal numbers of adult males and females, all on quality habitat with established territories

###--- AGENTS
# There are two types of 'agents' in NetLogoR: patches and turtles
# Patches cannot move (i.e., landbase) while turtles can (i.e., fishers)

# Create a square landscape of 20 by 20 cells (400 cells total)
# Each cell is assumed to be the same size as one fisher territory
# Cell values are randomly chosen either 1 or 2
# Assume 0 = habitat unsuitable for fisher territory; 1 = suitable fisher habitat
# Create the patches
land <- createWorld(minPxcor = 1, maxPxcor = 20,
                    minPycor = 1, maxPycor = 20,
                    sample(c(0, 1), 400, replace = TRUE))
plot(land) # visualize the landscape

# randomly select 40 "good" habitat cells for adult female fishers
rtmp <- world2raster(land)
fishers_start <- as.data.frame(sampleStratified(rtmp, size=60, xy=TRUE)) %>% filter(layer==1)
fishers_start <- as.matrix(fishers_start[c("x","y")])

# Start with a landscape of reproductively capable females and males
# Create 40 adult fishers, all with established territories
# Place the fishers on "good" habitat
nfishers = 60
t0 <- createTurtles(n = nfishers, coords=fishers_start, breed="adult")

# assign each turtle as a female with an established territory
t0 <- turtlesOwn(turtles = t0, tVar = c("sex"), tVal = rep(c("F","M"), each=nfishers/2))
t0 <- turtlesOwn(turtles = t0, tVar = c("shape"), tVal = rep(c(16,15), each=nfishers/2)) # females are circles, males are squares
t0 <- turtlesOwn(turtles = t0, tVar = c("disperse"), tVal = c(rep("E", each=nfishers)))
t0 <- turtlesOwn(turtles = t0, tVar = c("mate_avail"), tVal = c(rep("NA", each=nfishers)))
t0 <- turtlesOwn(turtles = t0, tVar = c("repro"), tVal = 0)

# create a random age for the fishers
# randomly assign females with ages from 2.5-8 and males with ages from 2.5-4
# keep in mind that time steps are 6 months so have ages in 6 month increments
# the oldest a female fisher can be is 8 or 16 time steps
# the youngest time step for an adult is 2.5 or 5 time steps (juvenile = up to 2 years of 4 time steps)
yrs.adult <- c(sample(5:16, nfishers/2, replace=TRUE), sample(5:8, nfishers/2, replace=TRUE))
t0 <- turtlesOwn(turtles=t0, tVar = c("age"), tVal = yrs.adult/2)
t0

# Visualize the turtles on the landscape with their respective color
plot(land)
points(t0, pch = t0$shape, col = of(agents = t0, var = "color"))

###--- REPRODUCE
# check if mates are available for females
t0 <- find_mate(t0, dx=c(-4:4), dy=c(-4:4)) # give unrealistic mating distance to start off
t0[t0$mate_avail=="Y"] # number of 30 females able to find a mate with unrealistic distance

# for Central interior population
t1 <- denning(fishers=t0, denLCI=repro.CI$drC[3], denUCI=repro.CI$drC[4]); t1
t1 <- kits_produced(fishers=t1, ltrM=repro.CI$lsC[1], ltrSD=repro.CI$lsC[2]); t1

plot(land)
points(t1, pch = t1$shape, col = of(agents = t1, var = "color")) # looks the same because kit are at same location as their moms

################################################################################
# *** Step 2. SURVIVE *** #
# t1 = October
# function SURVIVE - add 0.5 to all fishers, kill off individuals who do not survive through t1

age.val <- of(agents=t1, var=c("age"))+0.5
t2 <- NLset(turtles = t1, agents=turtle(t1, who=t1$who),var="age", val=age.val)
NLcount(t2)

t2 <- survive(t2)
NLcount(t2)

plot(land)
points(t2, pch = t2$shape, col = of(agents = t2, var = "color")) # looks the same because kit are at same location as their moms

################################################################################
# *** Step 3. ESTABLISH / MAINTAIN TERRITORY & SCENT TERRITORY (MATE) & SURVIVE ***
# t2 = April
# 3a. function DISPERSE - run through DISPERSE function for individuals without territories, up to 30 times to allow 6 months of movement
# 3b. function FIND_MATE - for female fishers with ESTABLISHED territory, if male is within 2 cells in either direction or 8 adjacent cells plus same cell, assign mated status (i.e., if male is in same cell or ± 1 cell either via xlim and/or ylim)
# 3c. function SURVIVE - add 0.5 to all fishers, kill off individuals who do not survive this 6 month time step

t3 <- t2
for(i in 1:30){
  t3 <- disperse(land=land, fishers=t3, dist_mov=1.0)
}

t3 <- find_mate(t3)

age.val <- of(agents=t3, var=c("age"))+0.5
t3 <- NLset(turtles = t3, agents=turtle(t3, who=t3$who),var="age", val=age.val)
NLcount(t3)

t3 <- survive(t3)
NLcount(t3)

plot(land)
points(t3, pch = t3$shape, col = of(agents = t3, var = "color")) # looks the same because kit are at same location as their moms

################################################################################
# *** Step 4.  ESTABLISH / MAINTAIN TERRITORY & SURVIVE ***
# t3 = October
# 4a. function DISPERSE - run through DISPERSE function for individuals without territories, up to 30 times to allow 6 months of movement
# 4b. function SURVIVE - add 0.5 to all fishers, kill off individuals who do not survive through this 6 month time step

t4 <- t3
for(i in 1:30){
  t4 <- disperse(land=land, fishers=t4, dist_mov=1.0)
}

age.val <- of(agents=t4, var=c("age"))+0.5
t4 <- NLset(turtles = t4, agents=turtle(t4, who=t4$who),var="age", val=age.val)
NLcount(t4)

t4 <- survive(t4)
NLcount(t4)

plot(land)
points(t4, pch = t4$shape, col = of(agents = t4, var = "color")) # looks the same because kit are at same location as their moms

################################################################################
# *** Step 5. ESTABLISH / MAINTAIN TERRITORY & REPRODUCE & SCENT TERRITORY (MATE) & SURVIVE ***
# t4 = April
# 5a. function DENNING - assign 1 or 0 if female dens (based on nearby mates and published probabilities)
# 5b. function KITS_PRODUCED - kits produced (based on denning rate and published probabilities)

# 5c. function DISPERSE - run through DISPERSE function for individuals without territories, up to 30 times to allow 6 months of movement
# 5d. function FIND_MATE - for female fishers with ESTABLISHED territory, if male is within 2 cells in either direction or 8 adjacent cells plus same cell, assign mated status (i.e., if male is in same cell or ± 1 cell either via xlim and/or ylim)
# 5e. function SURVIVE - add 0.5 to all fishers, kill off individuals who do not survive this 6 month time step


# for Central interior population
t5 <- denning(fishers=t4, denLCI=repro.CI$drC[3], denUCI=repro.CI$drC[4])
t5 <- kits_produced(fishers=t5, ltrM=repro.CI$lsC[1], ltrSD=repro.CI$lsC[2])

for(i in 1:30){
  t5 <- disperse(land=land, fishers=t5, dist_mov=1.0)
}

t5 <- find_mate(t5)

age.val <- of(agents=t5, var=c("age"))+0.5
t5 <- NLset(turtles = t5, agents=turtle(t5, who=t5$who),var="age", val=age.val)
NLcount(t5)

t5 <- survive(t5)
NLcount(t5)

plot(land)
points(t5, pch = t5$shape, col = of(agents = t5, var = "color"))

################################################################################
