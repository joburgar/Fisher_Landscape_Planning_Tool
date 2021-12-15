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
#
#####################################################################################
# 00_NetLogoR_tutorial.R
# script to run through NetLogoR tutorials
# collated by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 08-Dec-2021
#####################################################################################
rversion <- R.Version()
libpath_path <- paste0("C:/Program Files/R/R-",rversion$major,".",rversion$minor,"/library")
.libPaths(libpath_path) # to ensure reading/writing libraries from C drive

# Bauduin, S., McIntire, E.J.B., Chubaty, A.M.
# NetLogoR: A package to build and run spatially explicit agentbased models in R
# Ecography
# R script example of a spatially explicit agent-based model using NetLogoR

# Load Packages
list.of.packages <- c("NetLogoR","nnls","lcmix","MASS","SpaDES.core","SpaDES.tools","tidyverse")

# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages, new.packages) # for housekeeping

#install.packages("lcmix", repos="http://R-Forge.R-project.org")
# AGENTS


# Create a square landscape of 20 by 20 cells (400 cells total)
# Cell values are randomly chosen either 0 or 1 (bad or good habitat)
land <- createWorld(minPxcor = 1, maxPxcor = 20,
                    minPycor = 1, maxPycor = 20,
                    sample(c(0, 1), 400, replace = TRUE))
plot(land) # visualize the landscape
rtmp <- world2raster(land)
turtles_start <- as.data.frame(sampleStratified(rtmp, size=10, xy=TRUE)) %>% filter(layer==1)
turtles_start <- as.matrix(turtles_start[c("x","y")])

# FEMALE ONLY MODEL
# Create ten female fishers and have them produce kits
# place the females on "good" habitat
t1 <- createTurtles(n = 10, coords=turtles_start, breed="adult")

# Visualize the turtles on the landscape with their respective color
plot(land)
points(t1, pch = 16, col = of(agents = t1, var = "color"))

# use the denning rate and mean litter size to determine number of kits produced
# for Central Interior - denning rate mean = 0.54; sd = 0.41
# for Central Inerior - litter size mean = 1.7; sd = 0.73
# rnorm(n=20, mean=5, sd=10)
who_breeds <- of(agents=t1, var="who")
for(i in 1:length(who_breeds)){
  if(rnorm(n=1, mean=0.54, sd=0.41)>0.5){
    # if >0.5 then female reproduces
    t1 <- hatch(t1, who=who_breeds[i], n=round(rnorm(n=1, mean=1.7, sd=0.73)),breed="juvenile")
  } else {
    NULL
  }
}

# assign sex to each of the young - 50/50 male/female
t1 <- turtlesOwn (turtles = t1, tVar = "sex",
                  tVal = c(rep("F",10),
                           sample(c("F", "M"), nrow(t1[t1$breed=="juvenile",]), replace = TRUE)))

NLcount(t1) # fishers
t1

# have the female juveniles disperse
t2 <- t1[t1$breed=="juvenile" & t1$sex=="F",]
NLcount(t2)
# for when changing juveniles to adults, NLset may be helpful
# t2 <- NLset(turtles = t2, agents = turtle(t2, who = 0), var = "breed", val = "wolf")


# MODEL

# number <- 1:10
#
# for (val in number)  {
#   if (val == 7)  {
#     print(paste("Coming out from for loop Where i =  ", val))
#     break
#   }
#   print(paste("Values are :  ", val))
# }

# Have the female fisher move 30 times within dispersal season
# If she finds a good habitat cell without another female, she can take it
# Otherwise she keeps dispersing

distRate = 1.0

# Create a new variable for kits to establish territory or keep dispersing
tcount <- turtlesOwn(turtles = t2, tVar = "Disperse", tVal=c(rep("D",nrow(t2))))

# if (eventType == "init") {
#     sim <- scheduleEvent(sim, start(sim), "WolfSheepPredation", "event")
#
# } else if (eventType == "plot")


for(i in 1:30){
  if(NLcount(tcount[tcount$Disperse=="E"])==NLcount(tcount)) {
     t3 <- tcount
     }
     else

  # run the model 30 times - based on assumption that female fisher can move ~35 km per month
  # Identify the cells the turtles are on
  cellTurtle <- patchHere(land, tcount)
    # And the values of these cells (good quality habitat, where born)
     distMove <- of(land, cellTurtle)
  # A turtle moves with a mean of 1-cell distance
  # at the time (distMove), drawn from a multivariate gamma
  # distribution to show that all turtles move similar
  # distances, i.e., affected by unmeasured conditions
     distShape <- distMove * distRate
     rho <- matrix(rep(0.8, length = nrow(tcount)*nrow(tcount)), ncol=nrow(tcount))
     diag(rho) <- 1
     distMoveRan <- rmvgamma(2, distShape, distRate, rho)[1, ] # vector
  # The fishers tcount move with a step length of distMoveRan (one value each)
  # The landscape is not a torus (torus = FALSE)
  # and the fishers can disperse outside of the landscape (out=TRUE)
     tcount <- fd(turtles = tcount, dist = distMoveRan,world = land, torus = FALSE, out = TRUE)

  # if the kit finds a good quality unoccupied cell, can stay, otherwise keeps moving
  # "D" = disperse; "E" = establish territory
     tcount.habitat <- of(world = land, agents = patchHere(world=land, turtles=tcount))
     tcount.patch <- patchHere(land, tcount)

     for(k in 1:nrow(tcount)){
       tcount.patch.occ <- turtlesOn(world = land, turtles = tcount[k],
                                  agents = patch(land, tcount.patch[k,1], tcount.patch[k,2]))
       if(tcount.habitat[k]==1 & nrow(tcount.patch.occ)==1){
         tcount <- NLset(turtles = tcount, agents = turtle(tcount, who = tcount[k]$who), var = "Disperse", val = "E")
         } else {
           tcount <- NLset(turtles = tcount, agents = turtle(tcount, who = tcount[k]$who), var = "Disperse", val = "D")
         }
       }

  # If continuing on their dispersal, then
  # The fishers rotate with a multivariate normal turn angle,
  # based on the mean of the group, correlated at 0.8
     meanHeading <- mean(of(agents = tcount, var = "heading"))
     Sigma <- matrix(rep(0.8 * meanHeading, length = nrow(tcount)*nrow(tcount)), ncol = nrow(tcount))
     diag(Sigma) <- meanHeading
     angleInd = mvrnorm(n = 1, mu = of(agents = tcount, var = "heading"), Sigma = Sigma)
  # Turtles rotate to the right if angleInd > 0
  # or to the left if angleInd < 0
     tcount <- right(turtles = tcount, angle = angleInd)
     tcount.D <- of(agents=tcount, var="Disperse")
     D.value <- which(tcount.D=="D")
     tcount1 <- fd(turtles = tcount[tcount$Disperse=="D",], dist = distMoveRan[D.value], world = land, torus = FALSE, out = TRUE)
     valtcount1 <- of(agents=tcount1, var=c("heading","xcor","ycor"))
     tcount <- NLset(turtles=tcount, agents=turtle(tcount, who=tcount[tcount$Disperse=="D"]$who),
                  var=c("heading","xcor","ycor"), val=valtcount1)
     }



# Visualize the turtles' new position
plot(land)
points(t1[t1$breed=="adult"], pch = 16, col = of(agents = t1[t1$breed=="adult"], var = "color"))
points(t2, pch = 16, col = of(agents = t2, var = "color"))
points(t3, pch = 16, col = of(agents = t3, var = "color"))
###########################################################################################
# https://rdrr.io/cran/NetLogoR/f/vignettes/ProgrammingGuide.Rmd


# Create a world according to a given extent
w1 <- createWorld(minPxcor = 0, maxPxcor = 10, minPycor = 0, maxPycor = 10)

# Report the distance between the patch [pxcor = 0, pycor = 0] and the patch [pxcor = 1, pycor = 1]
pDist <- NLdist(agents = cbind(pxcor = 0, pycor = 0),
                agents2 = cbind(pxcor = 1, pycor = 1), world = w1, torus = TRUE)

# Create 10 turtles in the world w1
t1 <- createTurtles(n = 10, world = w1)

# Move all the turtles by a distance of 1
t1 <- fd(world = w1, turtles = t1, dist = 1)

# For all patches, assign a random value between 0 and 1
pQuality <- createWorld(minPxcor = 0, maxPxcor = 9, minPycor = 0, maxPycor = 9, data = runif(n = 100, min = 0, max = 1))

# Now each turtle in t1 has a "sex" variable
t1 <- turtlesOwn (turtles = t1, tVar = "sex",
                  tVal = c("M", "M", "M", "M", "M", "F", "F", "F", "F", "F"))

# 5 sheep and 5 wolves
t2 <- createTurtles(world = w1, n = 10, breed = c(rep("sheep", 5), rep("wolf", 5)))

# Or
sheep <- createTurtles(world = w1, n = 5, breed = "sheep") # 5 sheep
wolves <- createTurtles(world = w1, n = 5, breed = "wolf") # 5 wolves

# Turtle 0 which was "sheep" becomes "wolf"
t2 <- NLset(turtles = t2, agents = turtle(t2, who = 0), var = "breed", val = "wolf")

# Reports the pQuality value of the patches:
# [pxcor = 0, pycor = 0], [pxcor = 0, pycor = 1], and [pxcor = 0, pycor = 2]
of(world = pQuality, agents = patch(pQuality, c(0,0,0), c(0,1,2)))

# Define a variable
distRate <- 0.5

# MODEL
for(i in 1:10){ # run the model 10 times
  # Identify the cells the turtles are on
  cellTurtle <- patchHere(world = w1, turtles = t1)
  # And the values of these cells
  distMove <- of(world = w1, agents = cellTurtle)
  # A turtle moves with a mean of 1 or 2-cell distance
  # at the time (distMove), drawn from a multivariate gamma
  # distribution to show that all turtles move similar
  # distances, i.e., part of a social group or affected by
  # unmeasured conditions
  distShape <- distMove * distRate
  rho <- matrix(rep(0.8, length = nrow(t1) * nrow(t1)), ncol =
                  nrow(t1))
  diag(rho) <- 1
  distMoveRan <- rmvgamma(2, distShape, distRate, rho)[1, ] # vector
  # The turtles t1 move with a step length of distMoveRan (one value each)
  # The landscape is not a torus (torus = FALSE)
  # and the turtles cannot move outside of the landscape (out =FALSE)
  t1 <- fd(turtles = t1, dist = distMoveRan,
           world = w1, torus = FALSE, out = FALSE)
  # Then the turtles rotate with a multivariate normal turn angle,
  # based on the mean of the group, correlated at 0.8
  meanHeading <- mean(of(agents = t1, var = "heading"))
  Sigma <- matrix(rep(0.8 * meanHeading, length = nrow(t1) * nrow(t1)), ncol = nrow(t1))
  diag(Sigma) <- meanHeading
  angleInd = mvrnorm(n = 1, mu = rep(meanHeading, nrow(t1)), Sigma = Sigma)
  # Turtles rotate to the right if angleInd > 0
  # or to the left if angleInd < 0
  t1 <- right(turtles = t1, angle = angleInd)
  # Visualize the turtles' new position
  plot(w1)
  points(t1, pch = 16, col = of(agents = t1, var = "color"))
}

##################################################################################
# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "WolfSheepPredation",
  description = paste("Translation into R using the NetLogoR and SpaDES packages",
                      "of the Wolf-Sheep-Predation NetLogo model created by Wilensky (1997)"),
  keywords = c("NetLogo", "NetLogoR", "SpaDES", "wolf", "predator", "sheep", "prey", "predation"),
  authors = c(person("Sarah", "Bauduin", email = "sarahbauduin@hotmail.fr",
                     role = c("aut", "cre"))),
  childModules = character(),
  version = numeric_version("1.2.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "day", # e.g., "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "WolfSheepPredation.Rmd"),
  reqdPkgs = list("NetLogoR", "quickPlot", "SpaDES.core", "SpaDES.tools"),
  parameters = rbind(
    defineParameter(".plotInitialTime", "numeric", 0, NA, NA,
                    "This describes the simulation time at which the first plot event
                    should occur"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA,
                    "This describes the simulation interval time at which plot events
                    should occur"),
    defineParameter(".saveInitialTime", "numeric", 0, NA, NA,
                    "This describes the simulation time at which the first save event
                    should occur"),
    defineParameter(".saveInterval", "numeric", 1, NA, NA,
                    "This describes the simulation interval time at which save events
                    should occur"),
    defineParameter("grassOn", "logical", TRUE, NA, NA,
                    "TRUE to include grass in the model, FALSE to only include wolves and sheep"),
    defineParameter("grassTGrowth", "numeric", 30, 0, NA,
                    "How long it takes for grass to regrow once it is eaten"),
    defineParameter("nSheep", "numeric", 100, 0, NA, "Initial sheep population size"),
    defineParameter("gainFoodSheep", "numeric", 4, 0, NA,
                    "Amount of energy sheep get for every grass patch eaten"),
    defineParameter("reproSheep", "numeric", 4, 0, 100,
                    "Probability in % of a sheep reproducing at each time step"),
    defineParameter("nWolf", "numeric", 50, 0, NA, "Initial wolf population size"),
    defineParameter("gainFoodWolf", "numeric", 20, 0, NA,
                    "Amount of energy wolves get for every sheep eaten"),
    defineParameter("reproWolf", "numeric", 5, 0, 100,
                    "Probability in % of a wolf reproducing at each time step")
  ),
  inputObjects = data.frame(
    objectName = NA_character_,
    objectClass = NA_character_,
    sourceURL = "",
    other = NA_character_,
    stringsAsFactors = FALSE
  ),
  outputObjects = data.frame(
    objectName = NA_character_,
    objectClass = NA_character_,
    other = NA_character_,
    stringsAsFactors = FALSE
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.WolfSheepPredation <- function(sim, eventTime, eventType, debug = FALSE) {
  if (eventType == "init") {
    ### check for more detailed object dependencies:
    ### (use `checkObject` or similar)

    # Create the world, the sheep and the wolves
    sim <- sim$WolfSheepPredationInit(sim)

    # schedule future event(s)
    sim <- scheduleEvent(sim, params(sim)$WolfSheepPredation$.plotInitialTime,
                         "WolfSheepPredation", "plot")
    sim <- scheduleEvent(sim, params(sim)$WolfSheepPredation$.saveInitialTime,
                         "WolfSheepPredation", "save")
    sim <- scheduleEvent(sim, start(sim), "WolfSheepPredation", "event")

  } else if (eventType == "plot") {

    dev(4)
    sim <- sim$WolfSheepPredationPosition(sim)
    dev(5)
    sim <- sim$WolfSheepPredationPopSize(sim)

    sim <- scheduleEvent(sim, time(sim) + params(sim)$WolfSheepPredation$.plotInterval,
                         "WolfSheepPredation", "plot")

  } else if (eventType == "save") {

    sim <- sim$WolfSheepPredationSave(sim)
    sim <- scheduleEvent(sim, time(sim) + params(sim)$WolfSheepPredation$.saveInterval,
                         "WolfSheepPredation", "save")

  } else if (eventType == "event") {

    sim <- sim$WolfSheepPredationEvent(sim)
    sim <- scheduleEvent(sim, time(sim) + 1, "WolfSheepPredation", "event")

  } else {
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  }
  return(invisible(sim))
}

## event functions
#   - follow the naming convention `modulenameEventtype()`;
#   - `modulenameInit()` function is required for initiliazation;
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
WolfSheepPredationInit <- function(sim) {
  # Create the world
  grass <- createWorld(minPxcor = -25, maxPxcor = 25, minPycor = -25, maxPycor = 25)
  if (params(sim)$WolfSheepPredation$grassOn == FALSE) {
    grass <- NLset(world = grass, agents = patches(grass), val = 1) # cannot plot an empty world
  }
  # If grassOn is TRUE, assign grass and countdown values to patches
  # Because there are multiple patches variables, a worldArray is needed
  # If grassOn is TRUE, the grass grows and the sheep eat it, if FALSE, the sheep don't need to eat
  if (params(sim)$WolfSheepPredation$grassOn == TRUE) {
    # Initialize patch values (grass and countdown) at random
    # 0 or 1 (i.e., green or brown in the NetLogo model)
    grassVal <- sample(c(0, 1), size = NLcount(patches(grass)), replace = TRUE)
    grass <- NLset(world = grass, agents = patches(grass), val = grassVal)
    countdown <- grass # countdown is a new world with the same extent as grass
    # Grass grow clock
    countdownVal <- runif(n = NLcount(patches(grass)),
                          min = 0, max = params(sim)$WolfSheepPredation$grassTGrowth)
    countdown <- NLset(world = countdown, agents = patches(countdown), val = countdownVal)
    sim$field <- stackWorlds(grass, countdown)
  }
  # When no patches values are used,
  # using grass, countdown or field as the world argument required by a function
  # does not change anything because they all have the same extent and number of patches.
  # When patches values are used (e.g., when the sheep eat the grass),
  # use only field as the world argument for the functions
  # which update and retrieve the patches values.
  # When field is updated, the values on the individual world grass and countdown are not updated,
  # only the layers in field are.

  # Assign the created world to the sim object
  sim$grass <- grass
  if (params(sim)$WolfSheepPredation$grassOn == TRUE) {
    sim$field <- sim$field
  }

  # Create the sheep
  sheep <- createTurtles(
    n = params(sim)$WolfSheepPredation$nSheep,
    coords = randomXYcor(world = grass, n = params(sim)$WolfSheepPredation$nSheep),
    breed = "aSheep", color = rep("red", params(sim)$WolfSheepPredation$nSheep)
  )
  # Add the energy variable
  sheep <- turtlesOwn(turtles = sheep, tVar = "energy",
                      tVal = runif(n = params(sim)$WolfSheepPredation$nSheep, min = 0,
                                   max = 2 * params(sim)$WolfSheepPredation$gainFoodSheep))
  sim$sheep <- sheep # assign the created sheep to the sim object

  # Create the wolves
  wolves <- createTurtles(n = params(sim)$WolfSheepPredation$nWolf,
                          coords = randomXYcor(world = grass,
                                               n = params(sim)$WolfSheepPredation$nWolf),
                          breed = "wolf",
                          color = rep("black", params(sim)$WolfSheepPredation$nWolf))
  # Add the energy variable
  wolves <- turtlesOwn(turtles = wolves, tVar = "energy",
                       tVal = runif(n = params(sim)$WolfSheepPredation$nWolf, min = 0,
                                    max = 2 * params(sim)$WolfSheepPredation$gainFoodWolf))
  sim$wolves <- wolves # assign the created wolves to the sim object

  sim$numSheep <- numeric() # keep track of how many sheep there are
  sim$numWolves <- numeric() # keep track of how many wolves there are

  # Initialize the count of grass if grassOn == TRUE
  if (params(sim)$WolfSheepPredation$grassOn == TRUE) {
    sim$numGreen <- numeric() # keep track of how much grass there is
  }

  return(invisible(sim))
}

### template for save events
WolfSheepPredationSave <- function(sim) {
  sim$numSheep <- c(sim$numSheep, NLcount(sim$sheep)) # add the new number of sheep
  sim$numWolves <- c(sim$numWolves, NLcount(sim$wolves)) # add the new numbr of wolves

  if (params(sim)$WolfSheepPredation$grassOn == TRUE) {
    # patches equal to 1 (green)
    pGreen <- NLwith(world = sim$field, var = "grass", agents = patches(sim$field), val = 1)
    sim$numGreen <- c(sim$numGreen, NLcount(pGreen)) # add the new number of green patches
  }

  return(invisible(sim))
}

### template for plot events
# Plot the positions
WolfSheepPredationPosition <- function(sim) {
  if (time(sim) == start(sim)) clearPlot()

  if (params(sim)$WolfSheepPredation$grassOn == TRUE) {
    grassRas <- sim$field[["grass"]]
    Plot(grassRas, na.color = "white")
  } else {
    grassRas <- sim$grass
    Plot(grassRas, col = "green")
  }

  if (NLcount(sim$sheep) > 0)
    Plot(sim$sheep, addTo = "grassRas", cols = "blue")
  if (NLcount(sim$wolves) > 0)
    Plot(sim$wolves, addTo = "grassRas", cols = "red")

  return(invisible(sim))
}

# Plot the population sizes
WolfSheepPredationPopSize <- function(sim) {
  if (time(sim) == params(sim)$WolfSheepPredation$.plotInitialTime) {
    clearPlot()
    plot(time(sim), NLcount(sim$wolves), xlim = c(start(sim), end(sim)),
         col = "red", pch = 19, cex = 0.5,
         ylim = c(0, params(sim)$WolfSheepPredation$nSheep * 6))
    points(time(sim), NLcount(sim$sheep),
           col = "blue", pch = 19, cex = 0.5)
  } else {
    points(time(sim), NLcount(sim$wolves),
           col = "red", pch = 19, cex = 0.5)
    points(time(sim), NLcount(sim$sheep),
           col = "blue", pch = 19, cex = 0.5)
    if (params(sim)$WolfSheepPredation$grassOn == TRUE) {
      points(time(sim), sim$numGreen[time(sim)] / 4,
             col = "green", pch = 19, cex = 0.5)
    }
  }

  return(invisible(sim))
}

### template for the main event using the different functions defined under
WolfSheepPredationEvent <- function(sim){
  if (NLany(sim$sheep) | NLany(sim$wolves)) {
    # Ask sheep
    if (NLcount(sim$sheep) != 0) {
      moveSheep(sim)
      if (params(sim)$WolfSheepPredation$grassOn == TRUE) {
        energySheep <- of(agents = sim$sheep, var = "energy")
        sim$sheep <- NLset(turtles = sim$sheep, agents = sim$sheep, var = "energy",
                           val = energySheep - 1)
        eatGrass(sim)
      }
      dieSheep(sim)
      if (NLcount(sim$sheep) != 0) {
        reproduceSheep(sim)
      }
    }

    # Ask wolves
    if (NLcount(sim$wolves) != 0) {
      moveWolves(sim)
      energyWolves <- of(agents = sim$wolves, var = "energy")
      sim$wolves <- NLset(turtles = sim$wolves, agents = sim$wolves,
                          var = "energy", val = energyWolves - 1)
      catchSheep(sim)
      dieWolves(sim)
      if (NLcount(sim$wolves) != 0) {
        reproduceWolves(sim)
      }
    }

    # Ask grass
    if (params(sim)$WolfSheepPredation$grassOn == TRUE) {
      growGrass(sim)
    }

  }

  return(invisible(sim))
}

### template for moveSheep
moveSheep <- function(sim) {
  sim$sheep <- move(sim$sheep)

  return(invisible(sim))
}

### template for moveWolves
moveWolves <- function(sim) {
  sim$wolves <- move(sim$wolves)

  return(invisible(sim))
}

### template for dieSheep
dieSheep <- function(sim) {

  sim$sheep <- death(sim$sheep)

  return(invisible(sim))
}

### template for dieWolves
dieWolves <- function(sim) {
  sim$wolves <- death(sim$wolves)

  return(invisible(sim))
}

### template for reproduceSheep
reproduceSheep <- function(sim) {
  sim$sheep <- reproduce(sim$sheep, params(sim)$WolfSheepPredation$reproSheep)

  return(invisible(sim))
}

### template for reproduceWolves
reproduceWolves <- function(sim) {
  sim$wolves <- reproduce(sim$wolves, params(sim)$WolfSheepPredation$reproWolf)

  return(invisible(sim))
}

#### Sheep and Wolves procedures

move <- function(turtles) {
  # In NetLogo, two functions are used to give a random heading
  # by rotating the turtles to the right and then to the left.
  # With NetLogoR, it can be replaced by only one function,
  # as a negative value to turn right will turn left:
  turtles <- right(turtles, angle = runif(n = NLcount(turtles), min = -50, max = 50))
  turtles <- fd(world = grass, turtles = turtles, dist = 1, torus = TRUE)
  return(turtles)
}

death <- function(turtles) {
  # When energy dips below 0, die
  whoEnergy <- of(agents = turtles, var = c("who", "energy"))
  # "who" numbers of the turtles with their energy value below 0
  who0 <- whoEnergy[which(whoEnergy[, "energy"] < 0), "who"]

  if (length(who0) != 0) {
    turtles <- die(turtles = turtles, who = who0)
  }
  return(turtles)
}

reproduce(t1)
reproduce <- function(turtles, reproTurtles) {
  # Throw dice to see if the turtles will reproduce
  repro <- runif(n = NLcount(turtles), min = 0, max = 100) < reproTurtles
  whoTurtles <- of(agents = turtles, var = "who") # "who" of the turtles before they reproduce
  reproWho <- whoTurtles[repro] # "who" of turtles which reproduce
  reproInd <- turtle(turtles, who = reproWho) # turtles which reproduce

  # if there is at least one turtle reproducing
  if (NLcount(reproInd) != 0) {
    energyTurtles <- of(agents = reproInd, var = "energy")
    # Divide the energy between the parent and offspring
    turtles <- NLset(turtles = turtles, agents = reproInd, var = "energy",
                     val = energyTurtles / 2)
    turtles <- hatch(turtles = turtles, who = reproWho, n = 1) # hatch one offspring per parent

    # Move the offspring by 1 step
    whoNewTurtles <- of(agents = turtles, var = "who") # "who" of the turtles after they reproduced
    whoOffspring <- which(!whoNewTurtles %in% whoTurtles) # "who" of offspring
    offspring <- turtle(turtles = turtles, who = whoOffspring)
    offspringMoved <- right(turtles = offspring,
                            angle = runif(n = NLcount(offspring), min = 0, max = 360))
    offspringMoved <- fd(world = grass, turtles = offspring, dist = 1, torus = TRUE)
    # Update the headings and coordinates of the offsprings inside the turtles
    valOffspring <- of(agents = offspringMoved, var = c("heading", "xcor", "ycor"))
    turtles <- NLset(turtles = turtles, agents = offspring, var = c("heading", "xcor", "ycor"),
                     val = valOffspring)
  }

  return(turtles)
}


### template for eatGrass
eatGrass <- function(sim) {
  pGreen <- NLwith(world = sim$field, var = "grass", agents = patches(sim$field),
                   val = 1) # patches with grass equal to 1 (green)
  sheepOnGreen <- turtlesOn(world = sim$field, turtles = sim$sheep,
                            agents = pGreen) # sheep on green patches

  if (NLcount(sheepOnGreen) != 0) {
    # These sheep gain energy by eating
    energySheep <- of(agents = sheepOnGreen, var = "energy") # energy before eating
    # Update energy
    sim$sheep <- NLset(turtles = sim$sheep, agents = sheepOnGreen, var = "energy",
                       val = energySheep + params(sim)$WolfSheepPredation$gainFoodSheep)

    # If a sheep is on a green patch (value equal to 1),
    # it eats the grass and turns it to brown (value to 0).
    pHere <- patchHere(world = sim$field, turtles = sheepOnGreen)
    sim$field <- NLset(world = sim$field, agents = pHere, var = "grass", val = 0)

  }

  return(invisible(sim))
}

### template for catchSheep
catchSheep <- function(sim) {
  # "who" numbers of sheep that are on the same patches as the wolves
  sheepWolves <- turtlesOn(world = sim$grass, turtles = sim$sheep,
                           agents = sim$wolves, simplify = FALSE)
  if (nrow(sheepWolves) != 0) {
    # sheepWolves[,"whoTurtles"] are the "who" numbers of sheep
    # sheepWolves[,"id"] represent the rank/order of the individual wolf in the wolves
    # (! not the "who" numbers of the wolves)
    sheepGrabbed <- oneOf(agents = sheepWolves) # grab one random sheep

    sim$sheep <- die(turtles = sim$sheep, who = sheepGrabbed) # kill the grabbed sheep
    whoWolves <- of(agents = sim$wolves, var = "who")
    whoGrabbingWolves <- whoWolves[unique(sheepWolves[, "id"])]
    grabbingWolves <- turtle(turtles = sim$wolves, who = whoGrabbingWolves)
    energyGrabbingWolves <- of(agents = grabbingWolves, var = "energy")
    # Get energy from eating for the wolves who grabbed sheep
    sim$wolves <- NLset(turtles = sim$wolves, agents = grabbingWolves, var = "energy",
                        val = energyGrabbingWolves + params(sim)$WolfSheepPredation$gainFoodWolf)
  }

  return(invisible(sim))
}

### template for growGrass
growGrass <- function(sim) {
  # Identify patches with grass equal to 0 (brown) and countdown less or equal to 0
  pBrown <- NLwith(world = sim$field, var = "grass", agents = patches(sim$field), val = 0)
  # Countdown values for the patches equal to 0 (brown)
  pBrownCountdown <- of(world = sim$field, var = "countdown", agents = pBrown)

  pBrownCountdown0 <- which(pBrownCountdown <= 0) # patches with a countdown <= 0
  if (length(pBrownCountdown0) != 0) {
    # Patches with grass equal to 0 (brown) and countdown <= 0
    pGrow <- pBrown[pBrownCountdown0, , drop = FALSE]
    # Grow some grass on these patches and reset the countdown
    sim$field <- NLset(
      world = sim$field, var = c("grass", "countdown"),
      agents = pGrow,
      val = cbind(
        grass = rep(1, NLcount(pGrow)),
        countdown = rep(params(sim)$WolfSheepPredation$grassTGrowth, NLcount(pGrow)))
    )
  }

  pBrownCountdown1 <- which(!pBrownCountdown <= 0) # patches with a countdown > 0
  if (length(pBrownCountdown1) != 0) {
    # patches with grass equal to 0 (brown) and countdown > 0
    pWait <- pBrown[pBrownCountdown1, , drop = FALSE]
    # Decrease the countdown for the patches which wait
    sim$field <- NLset(world = sim$field, var = "countdown", agents = pWait,
                       val = pBrownCountdown[pBrownCountdown1] - 1)
  }

  return(invisible(sim))
}


# Define the parameters
wolfSheepParams <- list(.plotInitialTime = NA, .plotInterval = NA,
                        .saveInitialTime = 0, .saveInterval = 1,
                        grassOn = TRUE, grassTGrowth = 30,
                        nSheep = 100, gainFoodSheep = 4, reproSheep = 4,
                        nWolf = 50, gainFoodWolf = 20, reproWolf = 5)

# Model init
wolfSheepSim <- simInit(
  times = list(start = 0, end = 500),
  params = list(WolfSheepPredation = wolfSheepParams),
  modules = list("WolfSheepPredation"),
  paths = list(modulePath = paste0(libpath_path,"/NetLogoR/examples/Wolf-Sheep-Predation"))
)

# Run the model
wolfSheepRun <- spades(wolfSheepSim)
