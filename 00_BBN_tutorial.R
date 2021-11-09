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

.libPaths("C:/Program Files/R/R-4.1.1/library") # to ensure reading/writing libraries from C drive

# Load Packages
list.of.packages <- c("tidyverse","bnstruct","qgraph","DAAG","visNetwork","tictoc")

# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
rm(list.of.packages, new.packages) # for housekeeping

# graph package no longer available from CRAN, need this workaround
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("graph")
# BiocManager::install("Rgraphviz")

#####################################################################################
# citation("bnstruct")

?BNDataset
# There are two ways to build a BNDataset:
# using two files containing respectivelyheader informations and data, and manually providing the data table and the
# related header informations (variable names, cardinality and discreteness).
# The key information needed are:
# 1.  the data;
# 2.  the state of variables (discrete or continuous);
# 3.  the names of the variables;
# 4.  the cardinalities of the variables (if discrete), or the number of levels they have to be quantized into (if continuous)

child.dataset <- child() # already in class BNDataset
asia.dataset <- asia()

options(max.print=200, width=60)
raw.data(child.dataset)

dataset <- bootstrap(child.dataset, num.boots=100)

tic() # to compute model run time
dataset.with.imputed.samples <- bootstrap(child.dataset, num.boots=100, imputation=TRUE)
toc()

show(dataset.with.imputed.samples)
complete.subset <- complete(dataset.with.imputed.samples, c(1,4,7))
show(complete.subset)

# for (i in 1:num.boots(dataset))+   print( boot(dataset, i) )

net <- learn.network(dataset.with.imputed.samples, algo="sem", scoring.func="AIC")
plot(net)

# simple example
dataset <- child()

# learning with available cases analysis, MMHC, BDeu

net <- learn.network(dataset)
plot(net)

# learning with imputed data, MMHC, BDeu
imp.dataset <- impute(dataset)
net <- learn.network(imp.dataset, use.imputed.data = TRUE)
plot(net)

# SEM, BDeu using previous network as starting point
tic()
net <- learn.network(dataset, algo = "sem",
                     scoring.func = "BDeu",
                     initial.network = net,struct.threshold = 10,
                     param.threshold = 0.001)
toc()
plot(net)

# we update the probabilities with EM from the raw dataset
# starting from the first network
engine  <- InferenceEngine(net)
results <- em(engine, dataset)
updated.engine  <- results$InferenceEngine
updated.dataset <- results$BNDataset
