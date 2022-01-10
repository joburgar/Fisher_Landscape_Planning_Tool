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
yrs.adult <- 5:16
t1 <- turtlesOwn(turtles=t1, tVar = c("age"), tVal = sample(yrs.adult, nfishers, replace=TRUE))


###--- MODEL FUNCTIONS
# Empirical / expert data as inputs
# use the denning rate and mean litter size to determine number of kits produced
# for Central Interior - denning rate mean = 0.54; sd = 0.41
# for Central Interior - litter size mean = 1.7; sd = 0.73
# create a reproduce function

###--- REPRODUCE
reproduce <- function(fishers, denM=0.54, denSD=0.41, ltrM=1.7, ltrSD=0.73) {

  # Random selection for which adult females reproduce, based on denning mean and SD (Central Interior)
  whoFishers <- of(agents = fishers, var = c("who","breed","sex")) # "who" of the fishers before they reproduce
  whoAFFishers <- whoFishers %>% filter(breed=="adult" & sex=="F") %>% dplyr::select(who)

  repro <- rnorm(n = nrow(whoAFFishers), mean=denM, sd=denSD) > 0.5
  whoAFFishers$repro <- repro

  reproWho <- whoAFFishers %>% filter(repro==TRUE) # "who" of fishers which reproduce
  reproInd <- turtle(fishers, who = reproWho$who) # fishers which reproduce

  # if there is at least one fisher reproducing
  # have those fishers have offspring, based on the mean and sd of litter size for Central Interior
  if (NLcount(reproInd) != 0) {
    # energyTurtles <- of(agents = reproInd, var = "energy") # might want to consider this for quality habitat
    # # Divide the energy between the parent and offspring
    # turtles <- NLset(turtles = turtles, agents = reproInd, var = "energy",
    #                  val = energyTurtles / 2)
    fishers <- hatch(turtles = fishers, who = reproWho$who, n=round(rnorm(n=1, mean=ltrM, sd=ltrSD)),breed="juvenile") # litter size based on Central Interior data

    whoNewFishers <- of(agents = fishers, var = "who") # "who" of the turtles after they reproduced
    whoOffspring <- fishers[fishers$breed=="juvenile",]$who # "who" of offspring
    offspring <- turtle(turtles = fishers, who = whoOffspring)

    # assign 50/50 male/female offspring, assign them all as dispersing
    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=offspring$who), var = "sex", val = sample(c("F","M"),NLcount(offspring),replace=TRUE))
    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=offspring$who), var = "disperse", val = "D")
    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=offspring$who), var = "age", val = "0") # just born so time step 0
    
  }

  return(fishers)
}

t2 <- reproduce(fishers=t1)

###--- SURVIVE
# should have first round of survival in here up to the first year
# so go through two time steps and have kits who survive have an age of 2

###--- DISPERSE
# Have the female fisher move 30 cells within dispersal season and the male fisher move 60 cells
# If she finds a good habitat cell without another female, she can take it, otherwise she keeps dispersing
# The disperse function movement distance can be changed, the default is 1 (i.e., moves 1 cell)
# Kits can move up to 30 times in one dispersal season so need to consider this when working into other loops
# Recall that 1 time step = 6 months or 30 potential moves
# "D" = disperse; "E" = establish territory

disperse <- function(land=land, fishers=fishers, dist_mov=1.0) {
  # Only want fishers without established territories to move
  # Assume female fisher can move ~35 km in a month, and that each pixel is 5.5 km in length or 7.8 km in diameter
  # Assume a male fisher can move ~70 km in a month
  # For ease of calculations, assume a dist_mov of 1.0 is one pixel
  # This means that a female fisher can move between 5-6 pixels per month or 30-36 pixels in each time step
  # And for simplicity, a male fisher can move 2*dist_mov (twice the distance of a female)
  # dist_mov relates to the number of cells (not quite right if fisher moving diagonally across a cell but works for our purposes)
  
  # fishers=t2
  whoDFishers <- fishers[fishers$disperse=="D",]$who
  disperseInd <- turtle(fishers, who = whoDFishers) # fishers who are dispersing (i.e., kits)

  # Have each fisher move 1 step in random heading
  # The landscape is not wrapped (torus = FALSE)
  # and the fishers can disperse outside of the landscape (out=TRUE)
  disperseInd <- right(disperseInd, angle = runif(n = NLcount(disperseInd), min = 0, max = 360))
  disperseInd <- fd(turtle(fishers, who=whoDFishers),dist=dist_mov, land, torus = FALSE, out = TRUE) # all kits move 1 cell
  
  ###--- for dispersing FEMALES
  # determine patch information for dispersing females
  # only run the loop if there are dispersing females
  if(NLcount(disperseInd[disperseInd$sex=="F"])!=0){
    disperseIndF <- turtle(disperseInd, who = disperseInd[disperseInd$sex=="F",]$who) # identify those dispersing (i.e., female kits)
    disperseHabitatF <- of(land, agents=patchHere(land, disperseIndF)) # identify habitat quality of current location
    dispersePatchF <- patchHere(land, disperseIndF) # the coordinates for the cells 
    
    # run loop to determine if females are dispersing
    # if the female kit finds a good quality cell (1) that is NOT occupied by another female (occupancy==1) can stay,
    # otherwise (if habitat = 0 OR occupancy>1) kit keeps moving
    # "D" = disperse; "E" = establish territory
    
    for(k in 1:nrow(disperseIndF)){
      # check how many fishers are currently on the cell
      disperseInd.patch.occ <- turtlesOn(world = land, turtles = disperseIndF[k],
                                       agents = patch(land, dispersePatchF[k,1], dispersePatchF[k,2]))
      
      if(disperseHabitatF[k]==1 & nrow(disperseInd.patch.occ)==1){ # if the habitat is good and there is only one turtle on the patch
        disperseIndF <- NLset(turtles = disperseIndF, agents = turtle(disperseIndF, who = disperseIndF[k]$who), var = "disperse", val = "E") # then establish
        } else {
          disperseIndF <- NLset(turtles = disperseIndF, agents = turtle(disperseIndF, who = disperseIndF[k]$who), var = "disperse", val = "D") # otherwise keep dispersing
        }
      }
  
  # now have updated object with kits dispersing or establishing
  # add the new values to the existing fishers 'turtle' object
  valdisperseIndF <- of(agents=disperseIndF, var=c("heading","xcor","ycor", "prevX","prevY","disperse"))
  fishers <- NLset(turtles = fishers, agents=turtle(fishers, who=disperseIndF$who),var=c("heading","xcor","ycor","prevX","prevY","disperse"), val=valdisperseIndF)
  }
  
  ###--- for dispersing MALES
  if(NLcount(disperseInd[disperseInd$sex=="M"])!=0){
    # determine patch information for dispersing males
    disperseIndM <- turtle(disperseInd, who = disperseInd[disperseInd$sex=="M",]$who) 
    disperseHabitatM <- of(land, agents=patchHere(land, disperseIndM)) 
    
    # if the male kit finds a good quality cell (1) and there aren't other established males within 2 cells (x and y direction) can stay,
    # otherwise (if habitat = 0 OR nearby male fishers>0) kit keeps moving
    for(k in 1:nrow(disperseIndM)){
      
      # check if any established male fisher is nearby (within 4 cells east:west and 4 cells north:south to consider within male territory)
      # k=2
      dispserseInd.neighbour <- turtlesAt(land, turtles = fishers[fishers$sex=="M" & fishers$disperse=="E"], agents = turtle(disperseIndM[k], who = disperseIndM[k]$who),
                                          dx=c(-2:2), dy=c(-2:2), torus = FALSE)
      
      if(disperseHabitatM[k]==1 & NLcount(dispserseInd.neighbour)==0){ # if the habitat is good and there are no other established male territories nearby
        disperseIndM <- NLset(turtles = disperseIndM, agents = turtle(disperseIndM, who = disperseIndM[k]$who), var = "disperse", val = "E") # then establish
        } else {
          disperseIndM <- NLset(turtles = disperseIndM, agents = turtle(disperseIndM, who = disperseIndM[k]$who), var = "disperse", val = "D") # otherwise keep dispersing
          
          # move dispersing male one more cell (can move 2 cells for every 1 cell female moves); keep the male moving in the same direction, moving one more dist.mov distance
          disperseIndTMP <- fd(turtle(disperseIndM, who=disperseIndM[k]$who),dist=dist_mov, land, torus = FALSE, out = TRUE)
          # patchHere(land, disperseIndTMP) # the coordinates for the cells

          disperseHabitatTMP <- of(land, agents=patchHere(land, disperseIndTMP))
          dispserseInd.neighbourTMP <- turtlesAt(land, fishers[fishers$sex=="M" & fishers$disperse=="E"], agents = turtle(disperseIndM, who = disperseIndTMP$who),
                                                dx=c(-2:2), dy=c(-2:2), torus = FALSE)
          if(disperseHabitatTMP==1 & is.na(NLcount(dispserseInd.neighbourTMP))==0){ # if the habitat is good and there are no other established male territories nearby
            disperseIndM <- NLset(turtles = disperseIndM, agents = turtle(disperseIndM, who = disperseIndTMP$who), var = "disperse", val = "E") # then establish
            } else {
              disperseIndM <- NLset(turtles = disperseIndM, agents = turtle(disperseIndM, who = disperseIndTMP$who), var = "disperse", val = "D") # otherwise keep dispersing
            }
        }
      }
    
    # now have updated object with kits dispersing or establishing
    # add the new values to the existing fishers 'turtle' object
    valdisperseIndM <- of(agents=disperseIndM, var=c("heading","xcor","ycor", "prevX","prevY","disperse"))
    fishers <- NLset(turtles = fishers, agents=turtle(fishers, who=disperseIndM$who),var=c("heading","xcor","ycor","prevX","prevY","disperse"), val=valdisperseIndM)
    }
  
  return(fishers)
  
  }



# now run function for up to 30 times for one season
t3 <- t2
t3 

for(i in 1:10){
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


###--- SURVIVE
# Have the fisher survive one time step depending on their age and cohort
# Use the survival function output from Eric's latest survival analysis
# Cohorts are broken down by population (Boreal / Central Interior), sex (M/F), and ageclass (A/J)
# Write the function with defaults to Central Interior but ability to update as necessary

###--- SURVIVE
# load Eric Lofroth's survival data (recieved Dec 2021)
# load("./data/fisher_survival.RData")
# library(survival)
# 
# # summary(fishersurvival)
# # Call: survfit(formula = fishersurv ~ Population + Sex + Ageclass, data = fisher1)
# # library(ggfortify)
# cox_full <- coxph(fishersurv ~ Population + Sex + Ageclass, data=fisher1)
# cox_nosex <- coxph(fishersurv ~ Population + Ageclass, data=fisher1)
# cox_nopop <- coxph(fishersurv ~ Sex + Ageclass, data=fisher1)
# cox_noac <- coxph(fishersurv ~ Population + Sex, data=fisher1)
# cox_pop <- coxph(fishersurv ~ Population, data=fisher1)
# cox_sex <- coxph(fishersurv ~ Sex, data=fisher1)
# cox_ac <- coxph(fishersurv ~ Ageclass, data=fisher1)
# anova(cox_full, cox_nosex, cox_noac, cox_ac, cox_pop, cox_nopop, cox_sex)
# 
# summary(cox_full)
# summary(cox_nosex)
# cox_fit <- survfit(cox_full)
# autoplot(cox_fit)

# summary(fishersurv)
# fisher1 %>% group_by(Population, Sex, Ageclass) %>% summarise(mean(DaysMonitored))
# fisher1$rownum <- as.numeric(rownames(fisher1))
# fisher_adult <- fisher1 %>% filter(Ageclass=="Adult")
# fisher_adult$rownum
# glimpse(fisher_adult)
# fishersurv_adult <- fishersurv[fisher_adult$rownum,]
# 
# km_full <- survfit(formula = fishersurv ~ Population + Sex + Ageclass, data=fisher1)
# km_adult <- survfit(formula = fishersurv_adult ~ Population + Sex, data=fisher_adult)
# 
# # no difference if running survfit with all data or just adult, so run with all data
# # ignore subadult output after 2 years as not biologically meaningful 
# # extract survival (mean and SE) for each cohort (Pop + Sex + Ageclass) in 6 month increments (365/2=182.5)
# # use these values as the probabilty of a fisher surviving in the IBM survival function
# surv_times <- c(182, 365, 548, 730, 912, 1095, 1278, 1460, 1642, 1825, 2008, 2190, 2372, 2555, 2738, 2920)
# 
# km_full_summary <- summary(km_full, times=surv_times, extend=TRUE)
# 
# str(km_full_summary)
# length(surv_times) # 16 times (every 6 months for 8 years)
# length(km_full_summary$strata) # 8 strata, repeating 16 values

# [1] "Population=Boreal, Sex=Female, Ageclass=Adult   "           "Population=Boreal, Sex=Female, Ageclass=Subadult"           "Population=Boreal, Sex=Male  , Ageclass=Adult   "
# [4] "Population=Boreal, Sex=Male  , Ageclass=Subadult"           "Population=Central Interior, Sex=Female, Ageclass=Adult   " "Population=Central Interior, Sex=Female, Ageclass=Subadult"
# [7] "Population=Central Interior, Sex=Male  , Ageclass=Adult   " "Population=Central Interior, Sex=Male  , Ageclass=Subadult"

# km_surv_estimates <- as.data.frame(cbind(km_full_summary$surv, km_full_summary$std.err, km_full_summary$lower, km_full_summary$upper))
# colnames(km_surv_estimates) <- c("Surv","SE","L95CL","U95CL")
# km_surv_estimates$Cohort <- rep(c("BFA","BFJ","BMA","BMJ","CFA","CFJ","CMA","CMJ"),each=16) # add in the group, e.g., BFA = Boreal Population, Female, Adult)
# km_surv_estimates$Time <- rep(surv_times, times=8) # add in the time interval for the estimates
# km_surv_estimates$Time_step <- rep(seq_len(16),times=8) # add in the time step (as per the IBM - every six months = 1 time step)
# km_surv_estimates$Use <- 1 # default to 1 (use) to set up case_when function
# km_surv_estimates <- km_surv_estimates %>% mutate(Use = case_when(grepl("J", Cohort) & Time_step > 4 ~ 0, TRUE ~ Use)) # only consider juvenile time steps for 2 years
# km_surv_estimates$Cohort <- as.factor(km_surv_estimates$Cohort)
# write.csv(km_surv_estimates, "km_surv_estimates.csv", row.names = FALSE)

# no need to create new survival estimates, load previously created unless new data from Eric
km_surv_estimates <- read.csv("km_surv_estimates.csv", header=TRUE)
glimpse(km_surv_estimates)

km_surv_estimates %>% group_by(Use,Cohort) %>% dplyr::summarise(max(Time_step))
km_surv_estimates %>% filter(Use==1 & grepl("J", Cohort)) 

# create a function that runs each time step to determine the probabilty of a fisher surviving to the next time step
survive <- function(fishers, vsex="F", vbreed="adult", dieFisher=0.5, surv_estimates=km_surv_estimates) {

  # Create
  whoFishers <- of(agents = fishers, var = c("who","breed","sex")) # "who" of the fishers before they reproduce
  whoSFishers <- whoFishers %>% filter(breed==vbreed & sex==vsex) # the select list of "who" to estimate survival

  if (nrow(whoSFishers[whoSFishers$breed=="adult" & whoSFishers$sex=="F",])!=0){  # adult female
    survival <- runif(n = nrow(whoSFishers), min=0.61, max=0.97)
    whoSFishers$survival <- survival

  } else if (nrow(whoSFishers[whoSFishers$breed=="adult" & whoSFishers$sex=="M",])!=0){ # adult male
    survival <- runif(n = nrow(whoSFishers), min=0.72, max=1.0)
    whoSFishers$survival <- survival

    } else if (nrow(whoSFishers[whoSFishers$breed=="juvenile" & whoSFishers$sex=="F",])!=0){ # subadult female
      survival <- runif(n = nrow(whoSFishers), min=0.06, max=0.76)
      whoSFishers$survival <- survival

      } else { # subadult male
        survival <- runif(n = nrow(whoSFishers), min=0.61, max=1.0)
        whoSFishers$survival <- survival
      }

  whoSFishers$die <- whoSFishers$survival^2 < dieFisher # this is the part that needs some thought!!!!!
  dieWho <- whoSFishers %>% filter(die==TRUE) # "who" of fishers which die

  fishers <- die(fishers, who=dieInd$who)
  return(fishers)
}

