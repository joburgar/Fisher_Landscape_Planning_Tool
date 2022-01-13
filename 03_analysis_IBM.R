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
# Females max age length = 8 (might instead go with survival analysis data)
# Males max age length = 4 (may also go with survival analysis data)
# Female territory = 30 km2
# Male territory = 5*female territory (have 3 female territories) 150-200 km2
# Female transient animal = 70 km2 Euclidian distance within a couple months
# Male transient animal = moves much farther
# When a landscape is full of animals an animal has to move quite far – so if territories are all full, then the animal has to leave the study area (one went to Alberta)
# When the landscape is poor and isn’t able to support fishers, then the distance moved is also far

###--- Need to run the survival function every 6 months to start the next time step with appropriate number of individuals

# Start - Pregnant females on the landscape
# the start of the scenario, set it up with pregnant females in established territories

# *** Step 1. START ***
# t0 = April
# function REPRODUCE - all kits are 0 and random age of adults to start scenario
# simple assumption is to start with 10 adult males and 10 females, all on quality habitat

# *** Step 2. SURVIVE ***
# t1 = October
# function SURVIVE - add 0.5 to all fishers, kill off individuals who do not survive through t1

# *** Step 3. ESTABLISH / MAINTAIN TERRITORY & SCENT TERRITORY (MATE) & SURVIVE ***
# t2 = April
# 3a. function DISPERSE - run through DISPERSE function for individuals without territories, up to 30 times to allow 6 months of movement
# 3b. function MATE - for female fishers with ESTABLISHED territory, if male is within 2 cells, assign mated status (i.e., if male is in same cell or ± 2 cell either via xlim and/or ylim)
# 3c. function SURVIVE - add add 0.5 to all fishers, kill off individuals who do not survive through t2

# *** Step 4.  ESTABLISH / MAINTAIN TERRITORY & SURVIVE ***
# t3 = October
# 4a. function DISPERSE - run through DISPERSE function for individuals without territories, up to 30 times to allow 6 months of movement
# 4b. function SURVIVE - add 0.5 to all fishers, kill off individuals who do not survive through t3

# *** Step 5. REPRODUCE ***
# t4 = April
# function REPRODUCE

# Step 1. Female is kicked out of natal territory (t1 - Oct 1)
# *** Probability of survival can be found in km_surv_estimates (population, age and sex dependent)
# *** Assume female fisher can move ~35 km in a month, and if each pixel is 5.5 km in length or 7.8 km in diameter than a female fisher can move between 5-6 pixels per month or 30-36 pixels in each time step. (Will need to think of a movement model to use - random walk? Need to code in that bearing can change within timestep???)
# *** If the dispersing female encounters a vacant territory then move to Step 2, otherwise go back to Step 1.
# *** Assume that first available territory beyond some base threshold will be taken if vacant (later can add in increased mortality risk if territory quality is lower, but sufficient; to start have territories as 1 = 1 suitable and 0 = 0 unsuitable).
# *** Can only survive until 2 without a territory - this means that if no territory by t4 (or 3 loops) then fisher dies.
# *** Cannot breed unless in vacant territory - will need to code this in (if pixel occupied, can travel through but not stay / breed)
#
# Step 2. Establishes / maintains territory & scents territory (t2 - Apr 1)
# *** If male within 2 cells of female (i.e., if male is in same cell or ± 1 cell either via xlim and/or ylim) , then female mates (or use denning rate if all female model)
# *** Estimate of denning rate: either 0.54 or 0.75 depending on population
#
# Step 3. Survived in territory (t3 - Oct)
# *** Necessary step for time steps to work properly
# *** Perhaps here we can add in adult survival or maybe at each time step using the probability of survival within the bounds of mean-max lifespan?  (Probability of survival in adulthood - 0.79 or 0.86 depending on population)
#
# Step 4. Kits are born (t4 - Apr 1)
# *** If on the most optimistic timeline, then female fisher is 2 years old when having first litter
# *** Number of kits based on mean litter size = 1.7 or 2.6 depending on population, round number of kits to whole number for realism
# *** Use Bernoulli distribution for female/male sex in kits (i.e., 50/50 probability)
# *** Loop back to Step 2 and allow female fisher to repeat steps 2-4 another 3-5 times (mean female fisher lifespan = 5 years, max = 8 years)
#
# For a male:
#   Start - Male is born (t0 - Apr 1)
#
# Step 1. Male is kicked out of natal territory.
# *** Probability of survival to Year 1 - either 0.86 or 1.0 depending on population
# *** Assume male fisher can move ~70 km in a month and if each pixel is 5.5 km in length or 7.8 km in diameter than a male fisher can move between 9-13 pixels per month or 45-78 pixels in each time step. (Will need to think of a movement model to use - random walk? Need to code in that bearing can change within timestep???)
# *** If there are no other males within 2 pixels (i.e., ± 1 cell either via xlim and/or ylim), consider the territory as 'vacant' and then move to Step 2, otherwise go back to Step 1.
# *** Can only survive until 2 without a territory - this means that if no territory by t4 (or 3 loops) then fisher dies.
#
# Step 2. Establishes / maintains territory & finds female(s) (t2 - Apr 1)
# *** If female within 2 cells of male (i.e., if female is in same cell or ± 1 cell either via xlim and/or ylim) ,
# *** This is really just a female step. What matters is that if a male establishes a territory he can live up to 4 years, otherwise he dies after 2. Can add in some sort of probability associated with survival for each time step from establishing territory (t2 at earliest, t4 at latest) till death at  end of t8. Rosiest scenario has male breeding for 4 time steps (t2, t4, t6, t8). Make mortality probabilistic based on male adulthood survival - 0.9 or 0.33 depending on population.


################################################################################
# *** Step 1. START ***
# t0 = April
# function DENNING & KITS_PRODUCED - all kits are 0 and random age of adults to start scenario

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

# check if mates are available for females
t0 <- find_mate(t0, dx=c(-4:4), dy=c(-4:4)) # give unrealistic mating distance to start off
t0[t0$mate_avail=="Y"] # number of 30 females able to find a mate with unrealistic distance

# read in csv created from reproductive rates in Rich and Eric's paper - 00_surv_repro_estimates_prep.R

repro.CI <- read.csv("data/repro.CI.csv", header=TRUE, row.names = 1)


# for Central interior population
t1 <- denning(fishers=t0, denLCI=repro.CI$drC[3], denUCI=repro.CI$drC[4]); t1
t1 <- kits_produced(fishers=t1, ltrM=repro.CI$lsC[1], ltrSD=repro.CI$lsC[2]); t1
points(t1, pch = 16, col = of(agents = t1, var = "color")) # looks the same because kit are at same location as their moms


################################################################################
# *** Step 2. SURVIVE *** #

# should have first round of survival in here up to the first year
# so go through two time steps and have kits who survive have an age of 2

# have all fisher progress 1 time step (i.e., age 6 months)
age.val <- of(agents=t1, var=c("age"))+0.5
t2 <- NLset(turtles = t1, agents=turtle(t1, who=t1$who),var="age", val=age.val)
t2

# load Eric Lofroth's survival data (recieved Dec 2021)
# load("./data/fisher_survival.RData")
# or read the csv of the already processed / formatted survival probability estimates
km_surv_estimates <- read.csv("data/km_surv_estimates.csv", header=TRUE)

# # subset to estimates needed for survival function
# km_surv_estimates <- km_surv_estimates %>% filter(Use==1 & age<8.5) %>% dplyr::select(-Use)
# glimpse(km_surv_estimates)
#
# # data check - to make sure it makes sense for each age class
# km_surv_estimates %>% group_by(Cohort) %>% summarise(max(age))
# km_surv_estimates %>% filter(grepl("J", Cohort))
# # km_surv_estimates %>% filter(grepl("A", Cohort))

# now run function for up to 30 times for one season
t2 <- survive(t2)
t2

plot(land)
points(t2, pch = 16, col = of(agents = t2, var = "color")) # looks the same because kit are at same location as their moms
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

t3 <- t2
for(i in 1:30){
  t3 <- disperse(land=land, fishers=t3, dist_mov=1.0)
}


# should think about adding in a "break" to stop it once all "D" turn to "E"
# currently just keeps running through

t3 # all have now established territory
plot(land)
points(t2)
points(t3, pch = 16, col = of(agents = t3, var = "color")) # looks the same because kit are at same location as their moms

# next step is to add in survival probability for each fisher to survive a full year
# and then to repeat the reproducing and dispersing functions


