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
t1 <- createTurtles(n = 10, coords=fishers_start, breed="adult")

# Visualize the turtles on the landscape with their respective color
plot(land)
points(t1, pch = 16, col = of(agents = t1, var = "color"))

# assign each turtle as a female with an established territory
t1 <- turtlesOwn(turtles = t1, tVar = c("sex"), tVal = c(rep("F", each=10)))
t1 <- turtlesOwn(turtles = t1, tVar = c("disperse"), tVal = c(rep("E", each=10)))


###--- MODEL FUNCTIONS
# Empirical / expert data as inputs
# use the denning rate and mean litter size to determine number of kits produced
# for Central Interior - denning rate mean = 0.54; sd = 0.41
# for Central Inerior - litter size mean = 1.7; sd = 0.73
# create a reproduce function

###--- REPRODUCE
reproduce <- function(fishers, denM=0.54, denSD=0.41, ltrM=1.7, ltrSD=0.73) {

  # Randomly selection for which adult females reproduce, based on denning mean and SD (Central Interior)
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

  }

  return(fishers)
}

t2 <- reproduce(fishers=t1)

###--- DISPERSE
# Have the female fisher move 30 times within dispersal season
# If she finds a good habitat cell without another female, she can take it
# Otherwise she keeps dispersing

disperse <- function(land, fishers, dist_mov=1.0) {
  # Only want fishers without established territories to move
  whoDFishers <- fishers[fishers$disperse=="D",]$who
  disperseInd <- turtle(fishers, who = whoDFishers) # fishers which reproduce

  # Have each fisher move 1 step in random heading
  # The landscape is not a torus (torus = FALSE)
  # and the fishers can disperse outside of the landscape (out=TRUE)
  disperseInd <- right(disperseInd, angle = runif(n = NLcount(disperseInd), min = 0, max = 360))
  disperseInd <- fd(disperseInd, dist=dist_mov, land, torus = FALSE, out = TRUE)

  cellTurtle <- patchHere(land, disperseInd)
  patchHere(land, fishers)
  disperseHabitat <- of(land, agents=patchHere(land, disperseInd))
  dispersePatch <- patchHere(land, tcount)

# the code we'll need for determining if there is a male within 2 cells of another dispersing male (establish territory)
# and also if there is a female within 2 cells (able to mate)
    # turtlesAt(land, fishers, agents=turtle(fishers,who=10), dx=c(0:5), dy=c(0:5), torus = FALSE)
# if the kit finds a good quality unoccupied cell, can stay, otherwise keeps moving
# "D" = disperse; "E" = establish territory

for(k in 1:nrow(tcount)){
  tcount.patch.occ <- turtlesOn(world = land, turtles = tcount[k],
                                agents = patch(land, tcount.patch[k,1], tcount.patch[k,2]))
  if(tcount.habitat[k]==1 & nrow(tcount.patch.occ)==1){
    tcount <- NLset(turtles = tcount, agents = turtle(tcount, who = tcount[k]$who), var = "Disperse", val = "E")
  } else {
    tcount <- NLset(turtles = tcount, agents = turtle(tcount, who = tcount[k]$who), var = "Disperse", val = "D")
  }
}
tcount <- right(turtles = tcount, angle = angleInd)
tcount.D <- of(agents=tcount, var="Disperse")
D.value <- which(tcount.D=="D")
tcount1 <- fd(turtles = tcount[tcount$Disperse=="D",], dist = distMoveRan[D.value], world = land, torus = FALSE, out = TRUE)
valtcount1 <- of(agents=tcount1, var=c("heading","xcor","ycor"))
tcount <- NLset(turtles=tcount, agents=turtle(tcount, who=tcount[tcount$Disperse=="D"]$who),
                var=c("heading","xcor","ycor"), val=valtcount1)


return(fishers)

}



# And the values of these cells (good quality habitat, where born)

#     # find who for female and male offspring
#     whoOffspringF <- fishers[fishers$breed=="juvenile" & fishers$sex=="F",]$who # "who" of female offspring
#     offspringF <- turtle(turtles = fishers, who = whoOffspringF)
#
#     whoOffspringM <- fishers[fishers$breed=="juvenile" & fishers$sex=="M",]$who # "who" of male offspring
#     offspringM <- turtle(turtles = fishers, who = whoOffspringM)
#
#     # Move the female offspring by up to 6 steps, stopping once on good habitat
#     offspringFMoved <- right(turtles = offspringF,
#                             angle = runif(n = NLcount(offspringF), min = 0, max = 360))
#     offspringFMoved <- fd(world = land, turtles = offspringF, dist = 1, torus=FALSE, out=TRUE)
#
#     offspringFhabitat <- of(world = land, agents = patchHere(world=land, turtles=offspringFMoved))
#     offspringFpatch <- patchHere(land, offspringFMoved)
#     offspringFpatchOcc <- turtlesOn(world=land, turtles=offspringFMoved, agents=patch(land, offspringFpatch[1], offspringFpatch[2]))
#
#     for(k in 1:nrow(tcount)){
#       tcount.patch.occ <- turtlesOn(world = land, turtles = tcount[k],
#                                     agents = patch(land, tcount.patch[k,1], tcount.patch[k,2]))
#       if(tcount.habitat[k]==1 & nrow(tcount.patch.occ)==1){
#         tcount <- NLset(turtles = tcount, agents = turtle(tcount, who = tcount[k]$who), var = "Disperse", val = "E")
#       } else {
#         tcount <- NLset(turtles = tcount, agents = turtle(tcount, who = tcount[k]$who), var = "Disperse", val = "D")
#       }
#     }
#
#     # Update the headings and coordinates of the offsprings inside the turtles
#     valOffspringF <- of(agents = offspringFMoved, var = c("heading", "xcor", "ycor"))
#     turtles <- NLset(turtles = turtles, agents = offspring, var = c("heading", "xcor", "ycor"),
#                      val = valOffspring)
#   }
#
#   return(turtles)
# }


###--- SURVIVE
# Have the fisher live for certain number of time steps and then die - maybe need to run this in a loop?
# what criteria do we want for fisher to die? right now written so that fisher will die if survival in 2 years (survival squared) falls below dieFisher (0.5)

###--- SURVIVE
survive <- function(fishers, vsex="F", vbreed="adult", dieFisher=0.5) {

  # Randomly selection for which adult females reproduce, based on denning mean and SD (Central Interior)
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

