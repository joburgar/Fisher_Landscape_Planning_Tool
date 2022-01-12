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

# denLCI=repro.CI$drC[3], denUCI=repro.CI$drC[4], ltrM=repro.CI$lsC[1], ltrSD=repro.CI$lsC[2])

###--- REPRODUCE
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
  whoFishers <- of(agents = fishers, var = c("who","repro")) # "who" of the fishers before they reproduce
  whoRFishers <- whoFishers %>% filter(repro==1) %>% dplyr::select(who)

  repro <- rbinom(n = nrow(whoAFFishers), size=1, prob=denLCI:denUCI) # prob can be a range - confidence intervals?
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
    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=offspring$who), var = "age", val = 0) # just born so time step 0
    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=offspring$who), var = "mated", val = "N") # just born so time step 0

  }

  return(fishers)
}
###--- SURVIVE
# Have the fisher survive one time step depending on their age and cohort
# Use the survival function output from Eric's latest survival analysis
# Cohorts are broken down by population (Boreal / Central Interior), sex (M/F), and ageclass (A/J)
# create a function that runs each time step to determine the probability of a fisher surviving to the next time step

survive <- function(fishers, surv_estimates=km_surv_estimates, Fpop="C") {

  # fishers=t2
  survFishers <- of(agents = fishers, var = c("who","breed","sex","age")) # "who" of the fishers before they reproduce
  survFishers$Cohort <- toupper(paste0(rep(Fpop,times=nrow(survFishers)),survFishers$sex,substr(survFishers$breed,1,1)))

  survFishers <- left_join(survFishers,surv_estimates %>% dplyr::select(-Time_step, -age_6mnths, -Time),
                           by=c("Cohort", "age"))
  survFishers[is.na(survFishers)] <- 0
  survFishers$live <- NA

  for(i in 1:nrow(survFishers)){
    survFishers[i,10] <- rbinom(n=1, size=1, prob=survFishers[i,8]:survFishers[i,9])
  }

  dieWho <- survFishers %>% filter(live!=TRUE) # "who" of fishers which die

  fishers <- die(fishers, who=dieWho$who)
  return(fishers)
}

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

