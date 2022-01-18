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
set_up_world <- function(nfishers=nfishers, xlim=c(1,10), ylim=c(1,10), prophab=0.5){
  # nfishers=10 # must be divisible by 2 for equal portions of male and female to start

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

  fishers_start <- as.data.frame(sampleStratified(rtmp, size=nfishers, xy=TRUE, )) %>% filter(layer==1)
  fishers_start <- as.matrix(fishers_start[c("x","y")])

  # Start with a landscape of adult females and males, all with established territories
  # Place the fishers on "good" habitat

  t0 <- createTurtles(n = nfishers, coords=fishers_start, breed="adult")

  # assign 50:50 sex ratio, all females have an established territory, males do if >2 cells from another male
  t0 <- turtlesOwn(turtles = t0, tVar = c("sex"), tVal = rep(c("F","M"), each=nfishers/2))
  t0 <- turtlesOwn(turtles = t0, tVar = c("shape"), tVal = rep(c(16,15), each=nfishers/2)) # females are circles, males are squares
  t0 <- turtlesOwn(turtles = t0, tVar = c("disperse"), tVal = c(rep("E", each=nfishers)))
  t0 <- turtlesOwn(turtles = t0, tVar = c("mate_avail"), tVal = c(rep("NA", each=nfishers)))
  t0 <- turtlesOwn(turtles = t0, tVar = c("repro"), tVal = 0)

  # if a male is within 2 cells (i.e., same one or 8 adjacent, must be dispersing "D", otherwise can stay)
  male.fishers <- t0[t0$sex=="M" & t0$disperse=="E"]
  for(k in 1:NLcount(male.fishers)){
    #k=1; rm(k)
    patchAt(land, male.fishers[male.fishers$who==5,], dx=c(1:1), dy=c(1:1))
    nearby.male <- turtlesAt(land, turtles = t0,
                             agents = turtle(male.fishers, who = male.fishers[k]$who),
                             dx=c(-1:1), dy=c(-1:1), torus = FALSE)

    if(NLcount(nearby.male)==1){ # if there are no established male territories nearby then can stay
      t0 <- NLset(turtles = t0, agents = turtle(male.fishers, who = male.fishers[k]$who), var = "disperse", val = "E") # able to stay
    } else {
      t0 <- NLset(turtles = t0, agents = turtle(male.fishers, who = male.fishers[k]$who), var = "disperse", val = "D") # not able to stay

    }
  }

  # create a random age for the fishers
  # randomly assign females with ages from 2.5-8 and males with ages from 2.5-4
  # keep in mind that time steps are 6 months so have ages in 6 month increments
  # the oldest a female fisher can be is 8 or 16 time steps
  # the youngest time step for an adult is 2.5 or 5 time steps (juvenile = up to 2 years of 4 time steps)
  yrs.adult <- c(sample(5:16, nfishers/2, replace=TRUE), sample(5:8, nfishers/2, replace=TRUE))
  t0 <- turtlesOwn(turtles=t0, tVar = c("age"), tVal = yrs.adult/2)

  # Visualize the turtles on the landscape with their respective color
  plot(land)
  points(t0, pch = t0$shape, col = of(agents = t0, var = "color"))

  return(list("land"=land, "t0"=t0, "actual.prop.hab"=actual.prop.hab))

}

set_up_world(nfishers=10)
###--- REPRODUCE
find_mate <- function(fishers, dx=c(-4:4), dy=c(-4:4)){
  # fishers=t0; rm(k)
  whoMFishers <- fishers[fishers$sex=="F" & fishers$age>1 & fishers$disperse=="E"]$who

  if(length(whoMFishers)!=0){

    # if the female finds finds an established males within 2 cells (x and y direction) can mate,
    # otherwise (if nearby male fishers==0) no chance for mating
    for(k in 1:length(whoMFishers)){
      nearby.male <- turtlesAt(land, turtles = fishers[fishers$sex=="M" & fishers$age>0.5], agents = turtle(fishers, who = whoMFishers[k]),
                               dx=dx, dy=dy, torus = FALSE)
      nearby.male

      if(NLcount(nearby.male)>0){ # if there are established male territories nearby (i.e., potential mate(s))
        fishers <- NLset(turtles = fishers, agents = turtle(fishers, who = whoMFishers[k]), var = "mate_avail", val = "Y") # able to mate
      } else {
        fishers <- NLset(turtles = fishers, agents = turtle(fishers, who = whoMFishers[k]), var = "mate_avail", val = "N") # otherwise no potential to mate

      }
    }
  }
  return(fishers)
}


denning <- function(fishers, denLCI=denLCI, denUCI=denUCI) {

  # Random selection for which adult females reproduce, based on denning mean and SD (Central Interior)
  # fishers=t0; rm(fishers)
  whoFishers <- of(agents = fishers, var = c("who","breed","sex","mate_avail")) # "who" of the fishers before they reproduce
  whoAFFishers <- whoFishers %>% filter(breed=="adult" & sex=="F" & mate_avail=="Y") %>% dplyr::select(who)

  repro <- rbinom(n = nrow(whoAFFishers), size=1, prob=denLCI:denUCI) # prob can be a range - confidence intervals?
  whoAFFishers$repro <- repro

  fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=whoAFFishers$who), var = "repro", val = whoAFFishers$repro)

  return(fishers)
}


kits_produced <- function(fishers, ltrM=ltrM, ltrSD=ltrSD) {

  # Random selection for which adult females reproduce, based on denning mean and SD (Central Interior)
  # fishers=t1; rm(fishers)
  whoFishers <- of(agents = fishers, var = c("who","mate_avail","repro")) # "who" of the fishers before they reproduce
  reproWho <- whoFishers %>% filter(repro==1) %>% dplyr::select(who) # "who" of fishers which reproduce
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
    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=fishers[fishers$sex=="F",]$who), var = "shape", val = 16) # assign circle shape to females
    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=fishers[fishers$sex=="M",]$who), var = "shape", val = 15) # assign square shape to males
    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=offspring$who), var = "disperse", val = "D")
    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=offspring$who), var = "age", val = 0) # just born so time step 0
    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=offspring$who), var = "mate_avail", val = "NA") # just born so time step 0
    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=offspring$who), var = "repro", val = 0) # just born so time step 0

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

survive <- function(fishers, surv_estimates=lwdh_surv_estimates, Fpop="C") {

  # fishers=t1
  survFishers <- of(agents = fishers, var = c("who","breed","sex","disperse","age")) # "who" of the fishers before they reproduce
  survFishers$Cohort <- toupper(paste0(rep(Fpop,times=nrow(survFishers)),survFishers$sex,substr(survFishers$breed,1,1)))

  # survFishers <- left_join(survFishers,surv_estimates %>% dplyr::select(-Time_step, -age_6mnths, -Time),
  #                          by=c("Cohort", "age"))

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
  oldM <- survFishers %>% filter(sex=="M" & age>4) # "who" of male fishers who die of 'old age' (i.e., > 4 yrs)
  oldF <- survFishers %>% filter(sex=="F" & age>8) # "who" of female fishers who die of 'old age' (i.e., > 8 yrs)
  dispersing <- survFishers %>% filter(disperse=="D" & age>2) # "who" of dispersing fishers over 2

  fishers <- die(fishers, who=c(dieWho$who, oldM$who, oldF$who, dispersing$who))

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

disperse <- function(land=land, fishers=fishers, dist_mov=1.0) {
  # Only want fishers without established territories to move
  # Assume female fisher can move ~35 km in a month, and that each pixel is 5.5 km in length or 7.8 km in diameter
  # Assume a male fisher can move ~70 km in a month
  # For ease of calculations, assume a dist_mov of 1.0 is one pixel
  # This means that a female fisher can move between 5-6 pixels per month or 30-36 pixels in each time step
  # And for simplicity, a male fisher can move 2*dist_mov (twice the distance of a female)
  # dist_mov relates to the number of cells (not quite right if fisher moving diagonally across a cell but works for our purposes)

  # fishers=t2
  whoDFishers <- fishers[fishers$disperse=="D" & fishers$age>0,]$who
  disperseInd <- turtle(fishers, who = whoDFishers) # fishers who are dispersing (i.e., kits)

  # Have each fisher move 1 step in random heading
  # The landscape is not wrapped (torus = FALSE)
  # and the fishers can disperse outside of the landscape (out=TRUE)
  disperseInd <- right(disperseInd, angle = runif(n = NLcount(disperseInd), min = 0, max = 360))
  disperseInd <- fd(turtle(fishers, who=whoDFishers),dist=dist_mov, land, torus = FALSE, out = FALSE) # all kits move 1 cell # this doesn't allow for dispersing fishers
  dispersePatchF <- patchHere(land, disperseInd) # the coordinates for the cells

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
                                          dx=c(-1:1), dy=c(-1:1), torus = FALSE)

      if(disperseHabitatM[k]==1 & NLcount(dispserseInd.neighbour)==0){ # if the habitat is good and there are no other established male territories nearby
        disperseIndM <- NLset(turtles = disperseIndM, agents = turtle(disperseIndM, who = disperseIndM[k]$who), var = "disperse", val = "E") # then establish
      } else {
        disperseIndM <- NLset(turtles = disperseIndM, agents = turtle(disperseIndM, who = disperseIndM[k]$who), var = "disperse", val = "D") # otherwise keep dispersing

        # move dispersing male one more cell (can move 2 cells for every 1 cell female moves); keep the male moving in the same direction, moving one more dist.mov distance
        disperseIndTMP <- fd(turtle(disperseIndM, who=disperseIndM[k]$who),dist=dist_mov, land, torus = FALSE, out = TRUE)
        # patchHere(land, disperseIndTMP) # the coordinates for the cells

        disperseHabitatTMP <- of(land, agents=patchHere(land, disperseIndTMP))
        disperseHabitatTMP[is.na(disperseHabitatTMP)] <- 0

        dispserseInd.neighbourTMP <- turtlesAt(land, fishers[fishers$sex=="M" & fishers$disperse=="E"], agents = turtle(disperseIndM, who = disperseIndTMP$who),
                                               dx=c(-1:1), dy=c(-1:1), torus = FALSE)
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

