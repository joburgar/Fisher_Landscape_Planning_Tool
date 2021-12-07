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
# 00_BBN_tutorial.R
# script to run through Bayesian Belief Network tutorials
# collated by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 09-Nov-2021
#####################################################################################
rversion <- R.Version()
.libPaths(paste0("C:/Program Files/R/R-",rversion$major,".",rversion$minor,"/library")) # to ensure reading/writing libraries from C drive

# Bauduin, S., McIntire, E.J.B., Chubaty, A.M.
# NetLogoR: A package to build and run spatially explicit agentbased models in R
# Ecography
# R script example of a spatially explicit agent-based model using NetLogoR

library(NetLogoR)
#install.packages("nnls")
#install.packages("lcmix", repos="http://R-Forge.R-project.org")
#install.packages("MASS")
library(lcmix)
library(MASS)
# AGENTS
# Create a square landscape of 9 by 9 cells (81 cells total)
# Cell values are randomly chosen either 1 or 2
land <- createWorld(minPxcor = 1, maxPxcor = 9,
                    minPycor = 1, maxPycor = 9,
                    sample(c(1, 2), 81, replace = TRUE))
plot(land) # visualize the landscape
# Create three moving individuals (three turtles)
# Place the turtles in the middle of the landscape just created
t1 <- createTurtles(n = 3, world = land)
# Visualize the turtles on the landscape with their respective color
points(t1, pch = 16, col = of(agents = t1, var = "color"))
# Define a variable
distRate <- 0.5
# MODEL
for(i in 1:10){ # run the model 10 times
  # Identify the cells the turtles are on
  cellTurtle <- patchHere(world = land, turtles = t1)
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
  # The turtles t1 move with a step length of distMoveRan (one value
  each)
# The landscape is not a torus (torus = FALSE)
# and the turtles cannot move outside of the landscape (out =
FALSE)
t1 <- fd(turtles = t1, dist = distMoveRan,
         world = land, torus = FALSE, out = FALSE)
# Then the turtles rotate with a multivariate normal turn angle,
# based on the mean of the group, correlated at 0.8
meanHeading <- mean(of(agents = t1, var = "heading"))
Sigma <- matrix(rep(0.8 * meanHeading, length = nrow(t1) *
                      nrow(t1)),
                ncol = nrow(t1))
diag(Sigma) <- meanHeading
angleInd = mvrnorm(n = 1, mu = rep(meanHeading, nrow(t1)), Sigma =
                     Sigma)
# Turtles rotate to the right if angleInd > 0
# or to the left if angleInd < 0
t1 <- right(turtles = t1, angle = angleInd)
# Visualize the turtles' new position
points(t1, pch = 16, col = of(agents = t1, var = "color"))
}
