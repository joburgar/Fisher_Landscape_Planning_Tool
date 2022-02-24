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
# 00b_Female_IBM_functions.R
# script for Individual Based Model functions - female only model
# reproduce, survive, and disperse
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 02-Dec-2021, revised 15-Feb-2022
#####################################################################################

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

  # create values and assign as adult (females) with established territories
  t0 <- turtlesOwn(turtles = t0, tVar = c("shape"), tVal =16) # females are circles, males are squares
  t0 <- turtlesOwn(turtles = t0, tVar = c("disperse"), tVal = c(rep("E", each=nfishers)))
  t0 <- turtlesOwn(turtles = t0, tVar = c("repro"), tVal = 0)

  # have fishers randomly assigned a year between 2.5 and 1 year less than max life span
  yrs.adult <- (sample(5:((maxAgeFemale-1)*2), nfishers, replace=TRUE))/2
  t0 <- turtlesOwn(turtles=t0, tVar = c("age"), tVal = yrs.adult)

  # Visualize the turtles on the landscape with their respective color
  plot(land)
  points(t0, pch = t0$shape, col = of(agents = t0, var = "color"))

  return(list("land"=land, "t0"=t0, "actual.prop.hab"=actual.prop.hab))

}


###--- SET-UP WORLD with actual aoi - for female only IBM
set_up_REAL_world_FEMALE <- function(nFemales=nFemales, maxAgeFemale=maxAgeFemale,raoi=raoi){
  # nFemales = 10
  # maxAgeFemale = 9

  nfishers = nFemales

  # Each cell is assumed to be the same size as one fisher territory
  # Cell values are either 0 or 1
  # Assume 0 = habitat unsuitable for fisher territory; 1 = suitable fisher habitat
  # Upload the raster (binary for habitat)
  cells.good.habitat <- sum(raoi@data@values)
  total.cells <- raoi@ncols * raoi@nrows
  actual.prop.hab <- cells.good.habitat / total.cells

  land <- raster2world(raoi)
  as.matrix(land@pCoords)

  # for some reason NetLogoR world matrices are set up differently from rasters
  # need to flip, change coordinates (-1), and keep in mind that NL worlds are col by row
  habM <- as.matrix(land@.Data)
  habMflipped <- habM[nrow(habM):1,]

  mHabitat <- which(habMflipped==1, arr.ind=TRUE)
  tmpMatrix <- matrix(1, nrow=nrow(mHabitat), ncol=ncol(mHabitat))
  NLmHabitat <- mHabitat - tmpMatrix

  fishers_start <- as.data.frame(NLmHabitat)
  colnames(fishers_start) <- c("pycor", "pxcor")
  fishers_start$rank <- rank(round(runif(cells.good.habitat, min=100, max=999)))
  fishers_start <- fishers_start %>% filter(rank <= nFemales) %>% dplyr::select(-rank)
  fishers_start <- fishers_start[c("pxcor","pycor")]

  fishers_start <- as.matrix(fishers_start)
  # Start with a landscape of adult females and males, all on "good" habitat
  t0 <- createTurtles(n = nfishers, coords=fishers_start, breed="adult")

  # create values and assign as adult (females) with established territories
  t0 <- turtlesOwn(turtles = t0, tVar = c("shape"), tVal =16) # females are circles, males are squares
  t0 <- turtlesOwn(turtles = t0, tVar = c("disperse"), tVal = c(rep("E", each=nfishers)))
  t0 <- turtlesOwn(turtles = t0, tVar = c("repro"), tVal = 0)

  # have fishers randomly assigned a year between 2.5 and 1 year less than max life span
  yrs.adult <- (sample(5:((maxAgeFemale-1)*2), nfishers, replace=TRUE))/2
  t0 <- turtlesOwn(turtles=t0, tVar = c("age"), tVal = yrs.adult)

  # Visualize the turtles on the landscape with their respective color
  plot(land)
  points(t0, pch = t0$shape, col = of(agents = t0, var = "color"))

  return(list("land"=land, "t0"=t0, "actual.prop.hab"=actual.prop.hab))

}

# tmp <- set_up_world_FEMALE(nFemales=10, maxAgeFemale=9,xlim=c(1,10), ylim=c(1,10), prophab=0.7)

###--- REPRODUCE
repro_FEMALE <- function(fishers=fishers, repro_estimates=repro.CI, Fpop="C") {

  # Random (binomial) selection for which adult females reproduce, based on denning rates confidence intervals
  # fishers=tApr; fishers=tmp$t0; rm(fishers)
  whoFishers <- of(agents = fishers, var = c("who","breed")) # "who" of the fishers before they reproduce
  whoAFFishers <- whoFishers[whoFishers$breed=="adult",]$who

  denLCI=repro.CI[repro.CI$Pop==Fpop & repro.CI$Param=="L95CI",]$dr
  denUCI=repro.CI[repro.CI$Pop==Fpop & repro.CI$Param=="U95CI",]$dr

  # repro <- as.integer(rbernoulli(n=length(whoAFFishers), p=c(denLCI:denUCI))) # prob can be a range - use confidence intervals
  repro <- rbinom(n = length(whoAFFishers), size=1, prob=denLCI:denUCI) # prob can be a range - use confidence intervals
  fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=whoAFFishers), var = "repro", val = repro)

  # Random selection for which adult females reproduce, based on denning mean and SD (Central Interior)
  whoFishers <- as.data.frame(of(agents = fishers, var = c("who","repro"))) # "who" of the fishers before they reproduce
  reproWho <- whoFishers[whoFishers$repro==1,]$who # "who" of fishers which reproduce

  ltrM=repro.CI[repro.CI$Pop==Fpop & repro.CI$Param=="mean",]$ls
  ltrSD=repro.CI[repro.CI$Pop==Fpop & repro.CI$Param=="sd",]$ls

  # if there is at least one fisher reproducing
  # have those fishers have offspring, based on the mean and sd of empirical data
  if (length(reproWho) != 0) {
    fishers <- hatch(turtles = fishers, who = reproWho, n=round(rnorm(n=1, mean=ltrM, sd=ltrSD)/2),breed="juvenile") # litter size based on empirical data (divided by 2 for female only model)

    # assign all of the offsprig as dispersing, change repro and age values to reflect newborn kits rather than their moms
    allFishers <- of(agents=fishers, var="who")
    offspring <- allFishers[!(allFishers %in% whoFishers$who)]

    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=offspring), var = "disperse", val = "D")
    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=offspring), var = "age", val = 0) # just born so time step 0
    fishers <- NLset(turtles = fishers, agents = turtle(fishers, who=offspring), var = "repro", val = 0) # just born not yet reproductive
  }

  return(fishers)
}

# tmp2 <- repro_FEMALE(fishers=tmp$t0)
# age.val <- of(agents=tmp2, var=c("age"))+0.5
# tmp2 <- NLset(turtles = tmp2, agents=turtle(tmp2, who=tmp2$who),var="age", val=age.val)

###--- DISPERSE
disperse_FEMALE <- function(land=land, fishers=fishers, dist_mov=1.0, out=TRUE, torus=TRUE) {
  # Only want fishers without established territories to move
  # Assume female fisher can move ~35 km in a month, and that each pixel is 5.5 km in length or 7.8 km in diameter
  # For ease of calculations, assume a dist_mov of 1.0 is one pixel
  # This means that a female fisher can move between 5-6 pixels per month or 30-36 pixels in each time step
  # dist_mov relates to the number of cells (not quite right if fisher moving diagonally across a cell but works for our purposes)

  # fishers=tmp2
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

# tmp3 <- disperse_FEMALE(land=land, fishers=tmp2)
# tmp3
# plot(land)
# points(tmp3, pch = tmp3$shape, col = of(agents = tmp3, var = "color"))

###--- SURVIVE
# Have the fisher survive one time step depending on their age and cohort
# Use the survival function output from Eric's latest survival analysis
# Cohorts are broken down by population (Boreal / Central Interior), sex (M/F), and ageclass (A/J)
# create a function that runs each time step to determine the probability of a fisher surviving to the next time step
# also need to kill off any fishers that are over 8 years (female) and 4 years (male)
# *** UPDATE - not enough fishers were surviving when using age survival probabilities, changed to cohort level probabilities

survive_FEMALE <- function(fishers=fishers, surv_estimates=rf_surv_estimates, Fpop="C", maxAgeFemale=9) {

  # fishers=t1; fishers=tmp3
  survFishers <- of(agents = fishers, var = c("who","breed","disperse","age")) # "who" of the fishers at start of survival
  survFishers$Cohort <- toupper(paste0(rep(Fpop,times=nrow(survFishers)),rep("F",times=nrow(survFishers)),survFishers$sex,substr(survFishers$breed,1,1)))

  survFishers <- as.data.frame(left_join(survFishers,surv_estimates,by=c("Cohort")))

  survFishers[is.na(survFishers)] <- 0
  survFishers$live <- NA

  for(i in 1:nrow(survFishers)){
    # i=1; rm(i)
    if(survFishers[i,]$age!=0){ # can't kill off juveniles that haven't reached 6 months
      # survFishers[i,]$live <- as.integer(rbernoulli(n=1, p=c(survFishers[i,]$L95CL:survFishers[i,]$U95CL)))
      survFishers[i,]$live <- rbinom(n=1, size=1, prob=survFishers[i,]$L95CL:survFishers[i,]$U95CL)
    }
  }

  dieWho <- survFishers %>% filter(live==0) # "who" of fishers which die, based on probability
  oldF <- survFishers %>% filter(age>maxAgeFemale) # "who" of female fishers who die of 'old age' (i.e., > 8 yrs)
  dispersing <- survFishers %>% filter(disperse=="D" & age>2) # "who" of dispersing fishers over 2

  fishers <- die(fishers, who=c(dieWho$who, oldF$who, dispersing$who))
  return(fishers)
}

# tmp4 <- survive_FEMALE(fishers=tmp3)
# tmp4

