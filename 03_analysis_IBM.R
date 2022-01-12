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

# Start with a very simple example
# Die by 2 if haven’t found a territory.
# Females mean = 5 years; 8 = max
# Males = 4 max
# Female territory = 30 km2
# Male territory = 5*female territory (have 3 female territories) 150-200 km2
# Female transient animal = 70 km2 Euclidian distance within a couple months
# Male transient animal = moves much farther
# When a landscape is full of animals an animal has to move quite far – so if territories are all full, then the animal has to leave the study area (one went to Alberta)
# When the landscape is poor and isn’t able to support fishers, then the distance moved is also far

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

# randomly select 10 "good" habitat cells for adult female fishers
rtmp <- world2raster(land)
fishers_start <- as.data.frame(sampleStratified(rtmp, size=10, xy=TRUE)) %>% filter(layer==1)
fishers_start <- as.matrix(fishers_start[c("x","y")])

# FEMALE ONLY MODEL
# Create ten female fishers and have them produce kits
# place the females on "good" habitat
nfishers = 10
t1 <- createTurtles(n = nfishers, coords=fishers_start, breed="adult")

# Visualize the turtles on the landscape with their respective color
plot(land)
points(t1, pch = 16, col = of(agents = t1, var = "color"))

# assign each turtle as a female with an established territory
t1 <- turtlesOwn(turtles = t1, tVar = c("sex"), tVal = c(rep("F", each=nfishers)))
t1 <- turtlesOwn(turtles = t1, tVar = c("disperse"), tVal = c(rep("E", each=nfishers)))

# create a random age for the fishers
# keep in mind that time steps are 6 months so have ages in 6 month increments
# the oldest a female fisher can be is 8 or 16 time steps
# the youngest time step for an adult is 5 (juvenile = up to 2 years of 4 time steps)
####  QUESTION - SHOULD AGES BE IN 0.5 INCREMENTS TO REFLECT YEARS? OR IS IT NOT TOO CONFUSING FOR AN AGE OF 2 TO MEAN 1 YEAR?
yrs.adult <- sample(5:16, nfishers, replace=TRUE)

t1 <- turtlesOwn(turtles=t1, tVar = c("age"), tVal = yrs.adult/2)
t1

# read in csv created from reproductive rates in Rich and Eric's paper - 00_surv_repro_estimates_prep.R
repro.CI <- read.csv("data/repro.CI.csv", header=TRUE, row.names = 1)

# for Central interior population
t2 <- reproduce(fishers=t1, denLCI=repro.CI$drC[3], denUCI=repro.CI$drC[4], ltrM=repro.CI$lsC[1], ltrSD=repro.CI$lsC[2])
t2



# should have first round of survival in here up to the first year
# so go through two time steps and have kits who survive have an age of 2

###--- SURVIVE
# have all fisher progress 1 time step (i.e., age 6 months)
valt2 <- of(agents=t2, var=c("age"))+0.5
t2 <- NLset(turtles = t2, agents=turtle(t2, who=t2$who),var="age", val=valt2)
t2

# load Eric Lofroth's survival data (recieved Dec 2021)
# load("./data/fisher_survival.RData")
# or read the csv of the already processed / formatted survival probability estimates
km_surv_estimates <- read.csv("data/km_surv_estimates.csv", header=TRUE)

# subset to estimates needed for survival function
km_surv_estimates <- km_surv_estimates %>% filter(Use==1 & age<8.5) %>% dplyr::select(-Use)
glimpse(km_surv_estimates)

# data check - to make sure it makes sense for each age class
km_surv_estimates %>% group_by(Cohort) %>% summarise(max(age))
km_surv_estimates %>% filter(grepl("J", Cohort))
# km_surv_estimates %>% filter(grepl("A", Cohort))

# now run function for up to 30 times for one season

# now need to go through 1 time step before kits are kicked out of natal territory and 1 time step associated with dispersal
# ran through the dispersal for one "season" or 6 month period so need to update the age for each fisher

# For a female:
#   Start - Female is born (t0 - Apr 1)
#
# Step 1. Female is kicked out of natal territory (t1 - Oct 1)
# *** Probability of survival to t1 - either 0.41 or 0.50 depending on population
# *** Assume female fisher can move ~35 km in a month, and if each pixel is 5.5 km in length or 7.8 km in diameter than a female fisher can move between 5-6 pixels per month or 30-36 pixels in each time step. (Will need to think of a movement model to use - random walk? Need to code in that bearing can change within timestep???)
# *** If the fisher encounter a vacant territory then move to Step 2, otherwise go back to Step 1.
# *** Assume that first available territory beyond some base threshold will be taken if vacant (later can add in increased mortality risk if territory quality is lower, but sufficient; to start have territories as 1 = 1 suitable and 0 = 0 unsuitable).
# *** Can only survive until 2 without a territory - this means that if no territory by t4 (or 3 loops) then fisher dies.
# *** Cannot breed unless in vacant territory - will need to code this in (if pixel occupied, can travel through but not stay / breed)
# Step 2. Establishes / maintains territory & scents territory (t2 - Apr 1)


for(i in 1:30){
  t3 <- disperse(land=land, fishers=t3, dist_mov=1.0)
}


# should think about adding in a "break" to stop it once all "D" turn to "E"
# currently just keeps running through

t3 # all have now established territory
plot(land)
points(t1)
points(t3, pch = 16, col = of(agents = t3, var = "color"))

# next step is to add in survival probability for each fisher to survive a full year
# and then to repeat the reproducing and dispersing functions


