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
# 00_IBM_functions.R
# script for Individual Based Model functions
# reproduce, survive, and disperse
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 02-Dec-2021
#####################################################################################


###--- SET-UP WORLD
set_up_world <- function(nMales=nMales, maxAgeMale=maxAgeMale, nFemales=nFemales, maxAgeFemale=maxAgeFemale,
                         xlim=xlim, ylim=ylim, prophab=prophab){
  # nMales = 7
  # maxAgeMale = 6 # Rory's spreadsheet notes that fishers
  # nFemales = 13 # Rory's VORTEX table suggests ~65% female and 35% male in population
  # maxAgeFemale = 9
  # xlim = ylim = c(1,10)
  # prophab = 0.7
  # There are two types of 'agents' in NetLogoR: patches and turtles
  # Patches cannot move (i.e., landbase) while turtles can (i.e., fishers)

  # Create a landscape with coordinates equal to xlim and ylim (default is 20*20 square landscape)
  # Each cell is assumed to be the same size as one fisher territory
  # Cell values are randomly chosen either 0 or 1
  # Assume 0 = habitat unsuitable for fisher territory; 1 = suitable fisher habitat
  # with the proportion of suitable "good" habitat (1) given in function (default = 0.5)
  # Create the patches
  numcells <- (xlim[2]-xlim[1]+1) * (ylim[2]-ylim[1]+1)

  land <- createWorld(minPxcor = xlim[1], maxPxcor = xlim[2],
                      minPycor = ylim[1], maxPycor = ylim[2],
                      rbinom(numcells, 1, prophab))

  # randomly select "good" habitat cells for each fishers
  rtmp <- world2raster(land)
  numhabitatcells <- sum(rtmp@data@values) # number of "good" habitat cells
  actual.prop.hab <- numhabitatcells/numcells

  nfishers = nMales + nFemales
  fishers_start <- as.data.frame(sampleStratified(rtmp, size=nfishers, xy=TRUE, )) %>% dplyr::filter(layer==1)
  fishers_start <- as.matrix(fishers_start[c("x","y")])

  # Start with a landscape of adult females and males, all on "good" habitat

  t0 <- createTurtles(n = nfishers, coords=fishers_start, breed="adult")

  # assign 50:50 sex ratio, all females have an established territory, males do if >2 cells from another male
  t0 <- turtlesOwn(turtles = t0, tVar = c("sex"), tVal = c(rep("F", time=nFemales),rep("M", time=nMales)))
  t0 <- turtlesOwn(turtles = t0, tVar = c("shape"), tVal = c(rep(16, time=nFemales),rep(15, time=nMales))) # females are circles, males are squares
  t0 <- turtlesOwn(turtles = t0, tVar = c("disperse"), tVal = c(rep("E", each=nfishers)))
  t0 <- turtlesOwn(turtles = t0, tVar = c("mate_avail"), tVal = c(rep("NA", each=nfishers)))
  t0 <- turtlesOwn(turtles = t0, tVar = c("repro"), tVal = 0)

  # if a male is within 2 cells (i.e., same one or 8 adjacent, must be dispersing "D", otherwise can stay)
  male.fishers <- t0[t0$sex=="M" & t0$disperse=="E"]
  for(k in 1:NLcount(male.fishers)){
    #k=1; rm(k)
    nearby.male <- turtlesAt(land, turtles = t0,
                             agents = turtle(male.fishers, who = male.fishers[k]$who),
                             dx=c(-1:1), dy=c(-1:1), torus=torus)

    if(NLcount(nearby.male)==1){ # if there are no established male territories nearby then can stay
      t0 <- NLset(turtles = t0, agents = turtle(male.fishers, who = male.fishers[k]$who), var = "disperse", val = "E") # able to stay
    } else {
      t0 <- NLset(turtles = t0, agents = turtle(male.fishers, who = male.fishers[k]$who), var = "disperse", val = "D") # not able to stay

    }
  }

  # male.fishers
  # create a random age for the fishers
  # randomly assign females with ages from 2.5-8 and males with ages from 2.5-5
  # assign ages to one year less lifespan to remove potential of fishers dying off during set up (first year)
  # keep in mind that time steps are 6 months so have ages in 6 month increments
  # the oldest a female fisher can be is 9 or 16 time steps
  # the youngest time step for an adult is 2.5 or 5 time steps (juvenile = up to 2 years of 4 time steps)
  yrs.male <- sample(5:((maxAgeMale-1)*2), nMales, replace=TRUE)
  yrs.female <- sample(5:((maxAgeFemale-1)*2), nFemales, replace=TRUE)
  yrs.adult <- c((yrs.female/2),(yrs.male/2))

  t0 <- turtlesOwn(turtles=t0, tVar = c("age"), tVal = yrs.adult)

  # Visualize the turtles on the landscape with their respective color
  plot(land)
  points(t0, pch = t0$shape, col = of(agents = t0, var = "color"))

  return(list("land"=land, "t0"=t0, "actual.prop.hab"=actual.prop.hab))

}


###--- SET-UP WORLD for female only IBM
set_up_world_FEMALE <- function(nFemales=nFemales, maxAgeFemale=maxAgeFemale,xlim=xlim, ylim=ylim, prophab=prophab){
  # nFemales = 20
  # maxAgeFemale = 9
  # xlim = ylim = c(1,10)
  # prophab = 0.7
  # There are two types of 'agents' in NetLogoR: patches and turtles
  # Patches cannot move (i.e., landbase) while turtles can (i.e., fishers)

  # Create a landscape with coordinates equal to xlim and ylim (default is 20*20 square landscape)
  # Each cell is assumed to be the same size as one fisher territory
  # Cell values are randomly chosen either 0 or 1
  # Assume 0 = habitat unsuitable for fisher territory; 1 = suitable fisher habitat
  # with the proportion of suitable "good" habitat (1) given in function (default = 0.5)
  # Create the patches
  numcells <- (xlim[2]-xlim[1]+1) * (ylim[2]-ylim[1]+1)

  land <- createWorld(minPxcor = xlim[1], maxPxcor = xlim[2],
                      minPycor = ylim[1], maxPycor = ylim[2],
                      rbinom(numcells, 1, prophab))

  # randomly select "good" habitat cells for each fishers
  rtmp <- world2raster(land)
  numhabitatcells <- sum(rtmp@data@values) # number of "good" habitat cells
  actual.prop.hab <- numhabitatcells/numcells

  nfishers = nFemales
  fishers_start <- as.data.frame(sampleStratified(rtmp, size=nfishers, xy=TRUE, )) %>% dplyr::filter(layer==1)
  fishers_start <- as.matrix(fishers_start[c("x","y")])

  # Start with a landscape of adult females and males, all on "good" habitat

  t0 <- createTurtles(n = nfishers, coords=fishers_start, breed="adult")

  # assign 50:50 sex ratio, all females have an established territory, males do if >2 cells from another male
  t0 <- turtlesOwn(turtles = t0, tVar = c("shape"), tVal =16) # females are circles, males are squares
  t0 <- turtlesOwn(turtles = t0, tVar = c("disperse"), tVal = c(rep("E", each=nfishers)))
  t0 <- turtlesOwn(turtles = t0, tVar = c("repro"), tVal = 0)

  yrs.adult <- (sample(5:((maxAgeFemale-1)*2), nfishers, replace=TRUE))/2

  t0 <- turtlesOwn(turtles=t0, tVar = c("age"), tVal = yrs.adult)

  # Visualize the turtles on the landscape with their respective color
  plot(land)
  points(t0, pch = t0$shape, col = of(agents = t0, var = "color"))

  return(list("land"=land, "t0"=t0, "actual.prop.hab"=actual.prop.hab))

}

###--- REPRODUCE
find_mate <- function(land=land, fishers=fishers, fmdx=c(-4:4), fmdy=c(-4:4)){

  whoMFishers <- fishers[fishers$sex=="F" & fishers$age>1 & fishers$disperse=="E"]$who

  if(length(whoMFishers)!=0){

    # if the female finds finds an established males within 2 cells (x and y direction) can mate,
    # otherwise (if nearby male fishers==0) no chance for mating
    # k=1; rm(k)
    for(k in 1:length(whoMFishers)){
      nearby.male <- turtlesAt(land, turtles = fishers[fishers$sex=="M" & fishers$age>0.5], agents = turtle(fishers, who = whoMFishers[k]),
                               dx=fmdx, dy=fmdy, torus = FALSE)

      if(NLcount(nearby.male)>0){ # if there are established male territories nearby (i.e., potential mate(s))
        fishers <- NLset(turtles = fishers, agents = turtle(fishers, who = whoMFishers[k]), var = "mate_avail", val = "Y") # able to mate
      } else {
        fishers <- NLset(turtles = fishers, agents = turtle(fishers, who = whoMFishers[k]), var = "mate_avail", val = "N") # otherwise no potential to mate

      }
    }
  }
  return(fishers)
}


denning <- function(fishers=fishers, denLCI=denLCI, denUCI=denUCI) {

  # Random selection for which adult females reproduce, based on denning mean and SD (Central Interior)
  # fishers=t0; rm(fishers)
  whoFishers <- of(agents = fishers, var = c("who","breed","sex","mate_avail")) # "who" of the fishers before they reproduce
  whoAFFishers <- whoFishers %>% filter(breed=="adult" & sex=="F" & mate_avail=="Y") %>% dplyr::select(who)


  # denMC=repro.CI$drC[1], densdC=repro.CI$drC[2]
  # repro <- runif(n = nrow(whoAFFishers), min=denLCI, max=denUCI) # after discussion on 3-Feb-2022, decided to try runif # issue is that not binary so have to make decision on what threshold to use.... might just up the chance of a female finding a mate instead
  # repro <- rnorm(n = nrow(whoAFFishers), mean=denMC, sd=densdC) # after discussion on 3-Feb-2022, decided to try rnorm with mean and standard deviation # this way have some probabilities >1 so could use a higher threshold...still thinking upping the distance for females to find a mate might be the best approach
  repro <- rbinom(n = nrow(whoAFFishers), size=1, prob=denLCI:denUCI) # prob can be a range - use confidence intervals
  whoAFFishers$repro <- repro

  fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=whoAFFishers$who), var = "repro", val = whoAFFishers$repro)

  return(fishers)
}


kits_produced <- function(fishers=fishers, ltrM=ltrM, ltrSD=ltrSD) {

  # Random selection for which adult females reproduce, based on denning mean and SD (Central Interior)
  # fishers=t1; rm(fishers)
  whoFishers <- as.data.frame(of(agents = fishers, var = c("who","repro"))) # "who" of the fishers before they reproduce
  reproWho <- whoFishers[whoFishers$repro==1,]$who # "who" of fishers which reproduce
  reproInd <- turtle(fishers, who = reproWho) # fishers which reproduce

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
    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=fishers[fishers$sex=="F",]$who), var = "shape", val = 16) # assign circle shape to females
    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=fishers[fishers$sex=="M",]$who), var = "shape", val = 15) # assign square shape to males
    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=offspring$who), var = "disperse", val = "D")
    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=offspring$who), var = "age", val = 0) # just born so time step 0
    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=offspring$who), var = "mate_avail", val = "NA") # just born so time step 0
    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=offspring$who), var = "repro", val = 0) # just born so time step 0

  }

  return(fishers)
}

repro_FEMALE <- function(fishers=fishers, repro_estimates=repro.CI, Fpop="C") {

  # Random selection for which adult females reproduce, based on denning mean and SD (Central Interior)
  # fishers=t0; rm(fishers)
  whoFishers <- of(agents = fishers, var = c("who","breed")) # "who" of the fishers before they reproduce
  whoAFFishers <- whoFishers[whoFishers$breed=="adult",]$who

  denLCI=repro.CI[repro.CI$Pop==Fpop & repro.CI$Param=="L95CI",]$dr
  denUCI=repro.CI[repro.CI$Pop==Fpop & repro.CI$Param=="U95CI",]$dr

  repro <- rbinom(n = length(whoAFFishers), size=1, prob=denLCI:denUCI) # prob can be a range - use confidence intervals
  fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=whoAFFishers), var = "repro", val = repro)

  # Random selection for which adult females reproduce, based on denning mean and SD (Central Interior)
  # fishers=t1; rm(fishers)
  whoFishers <- as.data.frame(of(agents = fishers, var = c("who","repro"))) # "who" of the fishers before they reproduce
  reproWho <- whoFishers[whoFishers$repro==1,]$who # "who" of fishers which reproduce
  reproInd <- turtle(fishers, who = reproWho) # fishers which reproduce

  ltrM=repro.CI[repro.CI$Pop==Fpop & repro.CI$Param=="mean",]$ls
  ltrSD=repro.CI[repro.CI$Pop==Fpop & repro.CI$Param=="sd",]$ls

  # if there is at least one fisher reproducing
  # have those fishers have offspring, based on the mean and sd of empirical data
  if (NLcount(reproInd) != 0) {
    fishers <- hatch(turtles = fishers, who = reproWho, n=round(rnorm(n=1, mean=ltrM, sd=ltrSD)/2),breed="juvenile") # litter size based on empirical data (divided by 2 for female only model)

    # assign all of the offsprig as dispersing, change repro and age values to reflect newborn kits rather than their moms
    allFishers <- of(agents=fishers, var="who")
    offspring <- allFishers[!(allFishers %in% whoFishers$who)]

    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=offspring), var = "disperse", val = "D")
    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=offspring), var = "age", val = 0) # just born so time step 0
    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=offspring), var = "repro", val = 0) # just born so time step 0
  }

  return(fishers)
}

###--- SURVIVE
# Have the fisher survive one time step depending on their age and cohort
# Use the survival function output from Eric's latest survival analysis
# Cohorts are broken down by population (Boreal / Central Interior), sex (M/F), and ageclass (A/J)
# create a function that runs each time step to determine the probability of a fisher surviving to the next time step
# also need to kill off any fishers that are over 8 years (female) and 4 years (male)
# *** UPDATE - not enough fishers were surviving when using age survival probabilities, changed to cohort level probabilities

survive <- function(fishers=fishers, surv_estimates=rf_surv_estimates, Fpop="C", maxAgeMale=6, maxAgeFemale=9) {

  # fishers=t1
  survFishers <- of(agents = fishers, var = c("who","breed","sex","disperse","age")) # "who" of the fishers before they reproduce
  survFishers$Cohort <- toupper(paste0(rep(Fpop,times=nrow(survFishers)),survFishers$sex,substr(survFishers$breed,1,1)))

  survFishers <- as.data.frame(left_join(survFishers,surv_estimates,by=c("Cohort")))

  survFishers[is.na(survFishers)] <- 0
  survFishers$live <- NA

  for(i in 1:nrow(survFishers)){
    # i=1
    if(survFishers[i,]$age!=0){ # can't kill off juveniles that haven't reached 6 months
    survFishers[i,]$live <- rbinom(n=1, size=1, prob=survFishers[i,]$L95CL:survFishers[i,]$U95CL)
    }
  }

  dieWho <- survFishers %>% filter(live!=TRUE) # "who" of fishers which die, based on probability
  oldM <- survFishers %>% filter(sex=="M" & age>maxAgeMale) # "who" of male fishers who die of 'old age' (i.e., > 4 yrs)
  oldF <- survFishers %>% filter(sex=="F" & age>maxAgeFemale) # "who" of female fishers who die of 'old age' (i.e., > 8 yrs)
  dispersing <- survFishers %>% filter(disperse=="D" & age>2) # "who" of dispersing fishers over 2

  fishers <- die(fishers, who=c(dieWho$who, oldM$who, oldF$who, dispersing$who))

  return(fishers)
}


survive_FEMALE <- function(fishers=fishers, surv_estimates=rf_surv_estimates, Fpop="C", maxAgeFemale=9) {

  # fishers=t1; fishers=tApr
  survFishers <- of(agents = fishers, var = c("who","breed","disperse","age")) # "who" of the fishers before they reproduce
  survFishers$Cohort <- toupper(paste0(rep(Fpop,times=nrow(survFishers)),rep("F",times=nrow(survFishers)),survFishers$sex,substr(survFishers$breed,1,1)))

  survFishers <- as.data.frame(left_join(survFishers,surv_estimates,by=c("Cohort")))

  survFishers[is.na(survFishers)] <- 0
  survFishers$live <- NA

  for(i in 1:nrow(survFishers)){
    if(survFishers[i,]$age!=0){ # can't kill off juveniles that haven't reached 6 months
      survFishers[i,]$live <- rbinom(n=1, size=1, prob=survFishers[i,]$L95CL:survFishers[i,]$U95CL)
    }
  }

  dieWho <- survFishers %>% filter(live==0) # "who" of fishers which die, based on probability
  oldF <- survFishers %>% filter(age>maxAgeFemale) # "who" of female fishers who die of 'old age' (i.e., > 8 yrs)
  dispersing <- survFishers %>% filter(disperse=="D" & age>2) # "who" of dispersing fishers over 2

  fishers <- die(fishers, who=c(dieWho$who, oldF$who, dispersing$who))
  return(fishers)
}


###--- DISPERSE
# Have the female fisher move 30 cells within dispersal season and the male fisher move 60 cells
# If she finds a good habitat cell without another female, she can take it, otherwise she keeps dispersing
# The disperse function movement distance can be changed, the default is 1 (i.e., moves 1 cell)
# Fishers CANNOT move out of "world"
# Kits can move up to 30 times in one dispersal season so need to consider this when working into other loops
# Recall that 1 time step = 6 months or 30 potential moves
# "D" = disperse; "E" = establish territory

disperse <- function(land=land, fishers=fishers, dist_mov=1.0, out=TRUE, torus=TRUE) {
  # Only want fishers without established territories to move
  # Assume female fisher can move ~35 km in a month, and that each pixel is 5.5 km in length or 7.8 km in diameter
  # Assume a male fisher can move ~70 km in a month
  # For ease of calculations, assume a dist_mov of 1.0 is one pixel
  # This means that a female fisher can move between 5-6 pixels per month or 30-36 pixels in each time step
  # And for simplicity, a male fisher can move 2*dist_mov (twice the distance of a female)
  # dist_mov relates to the number of cells (not quite right if fisher moving diagonally across a cell but works for our purposes)

  # fishers=t1
  # land=w1$land
  whoDFishers <- fishers[fishers$disperse=="D" & fishers$age>0,]$who
  disperseInd <- turtle(fishers, who = whoDFishers) # fishers who are dispersing (i.e., kits)

  # Have each fisher move 1 step in random heading
  # The landscape is not wrapped (torus = FALSE); or try with torus wrapped and OUT=TRUE
  # and the fishers can disperse outside of the landscape (out=TRUE)
  disperseInd <- right(disperseInd, angle = runif(n = NLcount(disperseInd), min = 0, max = 360))
  # patchHere(land, disperseInd)
  disperseInd <- fd(turtle(disperseInd, who=whoDFishers),dist=dist_mov, land, torus=torus, out=out) # all kits move 1 cell # this doesn't allow for dispersing fishers
  # patchHere(land, disperseInd)

  # if any dispersing fishers have exited the worlds extent, remove them from the simulation
  fisher.location <- as.data.frame(patchHere(land, disperseInd))
  fisher.location$who <- disperseInd$who
  out.of.bounds.fisher <- fisher.location[is.na(fisher.location$pxcor),]$who

  disperseInd <- die(disperseInd, who=out.of.bounds.fisher) # remove fishers who traveled outside worlds extent from dispersing object
  fishers <- die(fishers, who=out.of.bounds.fisher) # remove fishers who traveled outside worlds extent from main object

  ###--- for dispersing FEMALES
  # determine patch information for dispersing females
  # only run the loop if there are dispersing females
  if(NLcount(disperseInd[disperseInd$sex=="F"])!=0){
    disperseIndF <- turtle(disperseInd, who = disperseInd[disperseInd$sex=="F",]$who) # identify those dispersing (i.e., female kits)
    disperseHabitatF <- of(land, agents=patchHere(land, disperseIndF)) # identify habitat quality of current location
    disperseHabitatF[is.na(disperseHabitatF)] <- 0 # any NA habitat (i.e., outside of world is NOT suitable)
    dispersePatchF <- patchHere(land, disperseIndF) # the coordinates for the cells
    dispersePatchF[is.na(dispersePatchF)] <- 0

    # run loop to determine if females are dispersing
    # if the female kit finds a good quality cell (1) that is NOT occupied by another female (occupancy==1) can stay,
    # otherwise (if habitat = 0 OR occupancy>1) kit keeps moving
    # "D" = disperse; "E" = establish territory

    for(k in 1:NLcount(disperseIndF)){
      # check how many fishers are currently on the cell
      # k=3
      disperseInd.patch.occ <- turtlesOn(world = land, turtles = disperseIndF[k],
                                         agents = patch(land, dispersePatchF[k,1], dispersePatchF[k,2]))

      if(disperseHabitatF[k]==1 & NLcount(disperseInd.patch.occ)==1){ # if the habitat is good and there is only one turtle on the patch
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
    disperseHabitatM[is.na(disperseHabitatM)] <- 0 # any NA habitat (i.e., outside of world is NOT suitable)

    # if the male kit finds a good quality cell (1) and there aren't other established males within 2 cells (x and y direction) can stay,
    # otherwise (if habitat = 0 OR nearby male fishers>0) kit keeps moving
    for(k in 1:NLcount(disperseIndM)){

      # check if any established male fisher is nearby (within 2 cells east:west and 2 cells north:south to consider within male territory)
      # k=2
      dispserseInd.neighbour <- turtlesAt(land, turtles = fishers[fishers$sex=="M" & fishers$disperse=="E"], agents = turtle(disperseIndM, who = disperseIndM[k]$who),
                                          dx=c(-1:1), dy=c(-1:1), torus=torus)

      if(disperseHabitatM[k]==1 & NLcount(dispserseInd.neighbour)==0){ # if the habitat is good and there are no other established male territories nearby
        disperseIndM <- NLset(turtles = disperseIndM, agents = turtle(disperseIndM, who = disperseIndM[k]$who), var = "disperse", val = "E") # then establish
      } else {
        disperseIndM <- NLset(turtles = disperseIndM, agents = turtle(disperseIndM, who = disperseIndM[k]$who), var = "disperse", val = "D") # otherwise keep dispersing

        # move dispersing male one more cell (can move 2 cells for every 1 cell female moves); keep the male moving in the same direction, moving one more dist.mov distance
        disperseIndTMP <- fd(turtle(disperseIndM, who=disperseIndM[k]$who),dist=dist_mov, land, torus=torus, out=out)
        # patchHere(land, disperseIndTMP) # the coordinates for the cells

        # if any dispersing fishers have exited the worlds extent, remove them from the simulation
        fisher.location <- as.data.frame(patchHere(land, disperseIndTMP))
        fisher.location$who <- disperseIndTMP$who
        out.of.bounds.fisher <- fisher.location[is.na(fisher.location$pxcor),]$who

        disperseIndTMP <- die(disperseIndTMP, who=out.of.bounds.fisher) # remove fishers who traveled outside worlds extent from dispersing object
        fishers <- die(fishers, who=out.of.bounds.fisher) # remove fishers who traveled outside worlds extent from main object

        disperseHabitatTMP <- of(land, agents=patchHere(land, disperseIndTMP))
        disperseHabitatTMP[is.na(disperseHabitatTMP)] <- 0

        dispserseInd.neighbourTMP <- turtlesAt(land, fishers[fishers$sex=="M" & fishers$disperse=="E"], agents = turtle(disperseIndM, who = disperseIndTMP$who),
                                               dx=c(-1:1), dy=c(-1:1), torus=torus)
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


disperse_FEMALE <- function(land=land, fishers=fishers, dist_mov=1.0, out=TRUE, torus=TRUE) {
  # Only want fishers without established territories to move
  # Assume female fisher can move ~35 km in a month, and that each pixel is 5.5 km in length or 7.8 km in diameter
  # For ease of calculations, assume a dist_mov of 1.0 is one pixel
  # This means that a female fisher can move between 5-6 pixels per month or 30-36 pixels in each time step
  # dist_mov relates to the number of cells (not quite right if fisher moving diagonally across a cell but works for our purposes)

  # fishers=tmp$t0
  # land=tmp$land
  whoDFishers <- fishers[fishers$disperse=="D" & fishers$age>0,]$who
  disperseInd <- turtle(fishers, who = whoDFishers) # fishers who are dispersing (i.e., kits)

  # Have each fisher move 1 step in random heading
  # The landscape is not wrapped (torus = FALSE); or try with torus wrapped and OUT=TRUE
  # and the fishers can disperse outside of the landscape (out=TRUE)
  disperseInd <- right(disperseInd, angle = runif(n = NLcount(disperseInd), min = 0, max = 360))
  # patchHere(land, disperseInd)
  disperseInd <- fd(turtle(disperseInd, who=whoDFishers),dist=dist_mov, land, torus=torus, out=out) # all kits move 1 cell # this doesn't allow for dispersing fishers
  # patchHere(land, disperseInd)

  # if any dispersing fishers have exited the worlds extent, remove them from the simulation
  fisher.location <- as.data.frame(patchHere(land, disperseInd))
  fisher.location$who <- disperseInd$who
  out.of.bounds.fisher <- fisher.location[is.na(fisher.location$pxcor),]$who

  disperseInd <- die(disperseInd, who=out.of.bounds.fisher) # remove fishers who traveled outside worlds extent from dispersing object
  fishers <- die(fishers, who=out.of.bounds.fisher) # remove fishers who traveled outside worlds extent from main object

  ###--- for dispersing FEMALES
  # determine patch information for dispersing females
  # only run the loop if there are dispersing females
  if(NLcount(disperseInd)!=0){
    disperseIndF <- turtle(disperseInd, who = disperseInd$who) # identify those dispersing (i.e., female kits)
    disperseHabitatF <- of(land, agents=patchHere(land, disperseIndF)) # identify habitat quality of current location
    disperseHabitatF[is.na(disperseHabitatF)] <- 0 # any NA habitat (i.e., outside of world is NOT suitable)
    dispersePatchF <- patchHere(land, disperseIndF) # the coordinates for the cells
    dispersePatchF[is.na(dispersePatchF)] <- 0

    # run loop to determine if females are dispersing
    # if the female kit finds a good quality cell (1) that is NOT occupied by another female (occupancy==1) can stay,
    # otherwise (if habitat = 0 OR occupancy>1) kit keeps moving
    # "D" = disperse; "E" = establish territory

    for(k in 1:NLcount(disperseIndF)){
      # check how many fishers are currently on the cell
      # k=1
      disperseInd.patch.occ <- turtlesOn(world = land, turtles = disperseIndF[k],
                                         agents = patch(land, dispersePatchF[k,1], dispersePatchF[k,2]))

      if(disperseHabitatF[k]==1 & NLcount(disperseInd.patch.occ)==1){ # if the habitat is good and there is only one turtle on the patch
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


  return(fishers)

}

