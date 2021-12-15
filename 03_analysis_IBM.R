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

# Create a square landscape of 9 by 9 cells (81 cells total)
# Each cell is assumed to be the same size as one fisher territory
# Cell values are randomly chosen either 1 or 2
# Assume 1 = habitat unsuitable for fisher territory; 2 = suitable fisher habitat
# Create the patches
w1 <- createWorld(minPxcor = 0, maxPxcor = 4, minPycor = 0, maxPycor = 4,
                  data = runif(25))
plot(w1)




land <- createWorld(minPxcor = 1, maxPxcor = 9,
                    minPycor = 1, maxPycor = 9,
                    sample(c(1, 2), 81, replace = TRUE))
plot(land) # visualize the landscape

# Create the turtles
# 10 moving individuals (start with all female)
# Place the turtles at different locations within the world
nturtle = 10 # number of turtles
t1 <- createTurtles(n = nturtle, coords = randomXYcor(land, n=nturtle), breed = c(rep("fisher", nturtle)))

# For 2-sex model, add a "sex" variable in t1
# t1 <- turtlesOwn (turtles = t1, tVar = "sex",
#                   tVal = c("M", "M", "M", "M", "M", "F", "F", "F", "F", "F"))

# Visualize the turtles on the landscape with their respective color
points(t1, pch = 16, col = of(agents = t1, var = "color"))

# Define a variable
distRate <- 0.5

# MODEL
plot(land) # visualize the landscape

for(i in 1:10){ # run the model 10 times
  # Identify the cells the turtles are on
  cellTurtle <- patchHere(world = land, turtles = t1)
  ?patchHere
  # And the values of these cells
  distMove <- of(world = land, agents = cellTurtle)
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
           world = land, torus = FALSE, out = FALSE)
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
  points(t1, pch = 16, col = of(agents = t1, var = "color"))
}
