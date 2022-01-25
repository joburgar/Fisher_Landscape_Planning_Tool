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
# 04_output.R
# script to produce outputs of Individual Based Models (IBMs) for fisher
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 25-Jan-2022
#####################################################################################
version$major
version$minor
R_version <- paste0("R-",version$major,".",version$minor)

.libPaths(paste0("C:/Program Files/R/",R_version,"/library")) # to ensure reading/writing libraries from C drive
tz = Sys.timezone() # specify timezone in BC

# Load Packages
list.of.packages <- c("tidyverse", "NetLogoR","nnls","lcmix","MASS","Cairo","PNWColors")
# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

source("00_IBM_functions.R")
#####################################################################################
# Create 3 sets of 100 simulations - vary the proportion of habitat and survival
# Low, medium and high habitat = 0.5, 0.6, and 0.7 (same world set up, get actual values)
# Low, medium and high survival = 0.7, 0.8, 0.9

load("out/IBM_noescape.RData")
w1 <- IBM_noescape[[1]]; w1$actual.prop.hab # 0.45
w2 <- IBM_noescape[[2]]; w2$actual.prop.hab # 0.56
w3 <- IBM_noescape[[3]]; w3$actual.prop.hab # 0.71

###--- plot the simulated landbases
Cairo(file="out/IBM_noescape_w1.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
plot(w1$land, main="Simulated Landbase\n45% Suitable Habitat")
points(w1$t0, pch = w1$t0$shape, col = of(agents = w1$t0, var = "color"))
dev.off()

Cairo(file="out/IBM_noescape_w2.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
plot(w2$land, main="Simulated Landbase\n56% Suitable Habitat")
points(w2$t0, pch = w2$t0$shape, col = of(agents = w2$t0, var = "color"))
dev.off()

Cairo(file="out/IBM_noescape_w3.PNG",type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
plot(w3$land, main="Simulated Landbase\n71% Suitable Habitat")
points(w3$t0, pch = w3$t0$shape, col = of(agents = w3$t0, var = "color"))
dev.off()

# Run 100 simulations for each, save as objects
# Calculate mean # of animals per cell at 10 years for each simulation to produce a heat map
# Create a figure with mean number of animals (+/- SE) for each time step and graph for each simulation

# IBM.w1.surv7.sim100 <- IBM_escape[[4]]
# IBM.w1.surv8.sim100 <- IBM_escape[[5]]
# IBM.w1.surv9.sim100 <- IBM_escape[[6]]
# IBM.w2.surv7.sim100 <- IBM_escape[[7]]
# IBM.w2.surv8.sim100 <- IBM_escape[[8]]
# IBM.w2.surv9.sim100 <- IBM_escape[[9]]
# IBM.w3.surv7.sim100 <- IBM_escape[[10]]
# IBM.w3.surv8.sim100 <- IBM_escape[[11]]
# IBM.w3.surv9.sim100 <- IBM_escape[[12]][[100]][[23]]

IBM_NLcount <- vector('list',9*100*23)
for(c in 1:20700){
  for(i in 4:12){
  for(j in 1:100){
      for(k in 1:23){
        if(length(IBM_noescape[[i]][[j]][[k]])==1){
          IBM_NLcount[c] <- 0
        } else {
          IBM_NLcount[c] <- NLcount(IBM_noescape[[i]][[j]][[k]])
        }
      }
    }
  }
}

