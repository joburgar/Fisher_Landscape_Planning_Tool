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
# script to produce outputs (through functions) of Individual Based Models (IBMs) for fisher
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 25-Jan-2022
#####################################################################################
version$major
version$minor
R_version <- paste0("R-",version$major,".",version$minor)

.libPaths(paste0("C:/Program Files/R/",R_version,"/library")) # to ensure reading/writing libraries from C drive
tz = Sys.timezone() # specify timezone in BC

# Load Packages
list.of.packages <- c("tidyverse", "NetLogoR","nnls","lcmix","MASS","Cairo","PNWColors", "ggplot2",
                      "sf","raster","rgdal","data.table")
# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

################################################################################
###--- output FUNCTIONS

# measures of uncertainty
se <- function(x) sqrt(var(x)/length(x))
LCL <- function(x) quantile(x, probs=0.05)
UCL <- function(x) quantile(x, probs=0.95)

# initial world plot
setup_plot <- function(sim_out=sim_out, name_out=name_out){
  # sim_out=scenario1 name_out="canBex1"
  Cairo(file=paste0("out/",name_out,"_",round(sim_out[[1]]$actual.prop.hab*100),"hab_setup.PNG"),type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
  plot(sim_out[[1]]$land, main=c(paste0("Simulated Landbase"),paste0(round(sim_out[[1]]$actual.prop.hab*100),"% Suitable Habitat")))
  points(sim_out[[1]]$t0, pch = sim_out[[1]]$t0$shape, col = of(agents = sim_out[[1]]$t0, var = "color"))
  dev.off()

  Cairo(file=paste0("out/",name_out,"_",round(sim_out[[1]]$actual.prop.hab*100),"hab_setup_nofishers.PNG"),type="png",width=2200,height=2000,pointsize=12,bg="white",dpi=300)
  plot(sim_out[[1]]$land, main=c(paste0("Simulated Landbase"),paste0(round(sim_out[[1]]$actual.prop.hab*100),"% Suitable Habitat")))
  dev.off()
}


# grab output from one set of 100 simulations
sim_output <- function(sim_out=sim_out, sim=sim, numsims=numsims, yrs_sim=yrs_sim){
  # sim_out=scenario1; sim=2; numsims=100; yrs_sim=10
  num.runs <- yrs_sim + 1

  ABM.df <- as.data.frame(array(NA,c(numsims,num.runs)))
  colnames(ABM.df) <- paste0("TimeStep_",str_pad(seq_len(num.runs),2,pad="0"))

  # Suggesting using data.table and lapply
  # Rows = replicates (n = 100)
  # Columns = time steps (n = 12)

  Reps <- 1:numsims
  timeSteps <- 1:num.runs

  ABM.df <- rbindlist(lapply(Reps, function(rps){
    ABM.df_ts <- rbindlist(lapply(timeSteps, function(ts){
      DT <- as.array(sim_out[[sim]][[rps]][[ts]])

      if (length(DT) != 0){
        nAdults = as.numeric(table(DT$breed)["adult"])
        nJuvenile = as.numeric(table(DT$breed)["juvenile"])
      } else {
        nAdults = 0
        nJuvenile = 0
      }
      tb <- data.table(Sim = paste0("Sim",str_pad(sim,2,pad="0")),
                       Run = rps,
                       TimeStep = paste0("TimeStep_",str_pad(ts,
                                                             2, pad="0")),
                       Count = nAdults)
      return(tb)
    }))
  }))

}


# grab output from all sets of simulations - use if saved more than one scenario in a list
pop_output <- function(sim_out=sim_out, sim_order=c(4:6), numsims=100, yrs_sim=10){
  # sim_out=Boreal_escape_FEMALE_binom
  # sim_order=c(4:6)
  num.runs <- yrs_sim + 1

  ABM.df <- as.data.frame(array(NA,c(num.runs*numsims*length(sim_order),4)))
  colnames(ABM.df) <- c("Run","TimeStep","Count","Sim")

  a=1
  b=num.runs*numsims
  for(i in 4:6){
    ABM.df[a:b,] <- sim_output(sim_out=sim_out, sim=i, numsims=numsims, yrs_sim=yrs_sim)
    a=a+(num.runs*numsims)
    b=b+(num.runs*numsims)
  }
return(ABM.df)
}


ABM_fig_1sim <- function(sim_out=sim_out, numsims=100, yrs_sim=10, Fpop=Fpop){
  # sim_out=scenario1; numsims=100; yrs_sim=10; Fpop="B"
  ABM.df <- sim_output(sim_out=sim_out, sim=2, numsims=numsims, yrs_sim=yrs_sim)
  ABM.df$Pop <- rep(Fpop, each=dim(ABM.df)[1])

  ABM.df$Pcnthab <- sim_out[[1]]$actual.prop.hab
  ABM.TS.mean <- ABM.df %>% dplyr::select(-Run) %>% pivot_wider(names_from=TimeStep, values_from=Count, values_fn=mean)
  ABM.TS.mean$Param <- "Mean"

  ABM.TS.se <- ABM.df %>% dplyr::select(-Run) %>% pivot_wider(names_from=TimeStep, values_from=Count, values_fn=se)
  ABM.TS.se$Param <- "SE"
  ABM.TS.LCL <- ABM.df %>% dplyr::select(-Run) %>% pivot_wider(names_from=TimeStep, values_from=Count, values_fn=LCL)
  ABM.TS.LCL$Param <- "LCL"

  ABM.TS.UCL <- ABM.df %>% dplyr::select(-Run) %>% pivot_wider(names_from=TimeStep, values_from=Count, values_fn=UCL)
  ABM.TS.UCL$Param <- "UCL"

  ABM.TS <- rbind(ABM.TS.mean, ABM.TS.se, ABM.TS.LCL, ABM.TS.UCL)

  ABM.TS.df <- ABM.TS %>% pivot_longer(cols = 4:(3+yrs_sim+1),names_to = "TimeStep",values_to = "Value" )
  ABM.TS.df <- ABM.TS.df %>% pivot_wider(names_from = Param, values_from = Value)

  ABM.TS.use <- ABM.TS.df %>% filter(!TimeStep %in% c("TimeStep_01"))

  ABM.TS.use$TimeStepNum <- as.numeric(substr(ABM.TS.use$TimeStep,10,11))

  pal_col <- pnw_palette(name="Starfish",n=7,type="discrete")

  Fpop_name <- ifelse(Fpop=="C","Columbian","Boreal")
  fishers_to_start <- nrow(sim_out[[1]]$t0)

  sim.TS.plot <- ggplot(data = ABM.TS.use) +
    theme_bw() + theme(strip.background = element_rect(fill = "white", colour = "white")) +
    theme(panel.grid = element_blank())+
    geom_ribbon(aes(x = TimeStepNum, ymin = LCL, ymax = UCL), fill = "#2c6184") +
    geom_vline(xintercept =7, col="darkgrey", lty=4) +
    geom_hline(yintercept = 0, col="grey", lty=4) +
    geom_line(aes(x = TimeStepNum, y = Mean)) +
    theme(axis.text.x = element_blank()) +
    xlab(expression(paste("Annual Predictions Starting at T"[0]))) +
    ylab("Number of Adult Female Fishers (Mean + 95% Confidence Intervals)")+
    ggtitle(paste0("Simulating ",yrs_sim," Years of ",Fpop_name," Fisher Populations,\nStarting with ",
                   round(unique(ABM.TS.df$Pcnthab*100)),"% Suitable Habitat and ",
                   fishers_to_start," Adult Female Fishers"))

  sim.TS.plot_se <- ggplot(data = ABM.TS.use) +
    theme_bw() + theme(strip.background = element_rect(fill = "white", colour = "white")) +
    theme(panel.grid = element_blank())+
    geom_vline(xintercept = "TimeStep_06", col="grey", lty=4) +
    geom_hline(yintercept = 0, col="grey", lty=4) +
    geom_point(aes(x = TimeStep, y = Mean), size=2) +
    geom_errorbar(aes(x = TimeStep, y = Mean, ymin=Mean-SE, ymax= Mean+SE),
                  width=.2, position=position_dodge(0.05)) +
    theme(axis.text.x = element_blank()) +
    xlab(expression(paste("Annual Predictions Starting at T"[0]))) +
    ylab("Number of Adult Female Fishers (Mean \u00B1 1 SE)")+ # \u00B1 is ± in unicode
    ggtitle(paste0("Simulating ",yrs_sim," Years of ",Fpop_name," Fisher Populations,\nStarting with ",
                   round(min(ABM.TS.df$Pcnthab*100)),"% Suitable Habitat and ",
                   fishers_to_start," Adult Female Fishers"))


  return(list(ABM.TS.df=ABM.TS.df, sim.TS.plot=sim.TS.plot, sim.TS.plot_se=sim.TS.plot_se))

}


################################################################################

### Create heatmaps for the outputs
# keep in mind that WGS84 lat/long espg = 4326; BC Albers espg = 3005; NAD83 / UTM zone 10N espg = 26910

# find coordinates for each fisher at 11 year mark (initial year + yrs+sim)
# create a heat map based on number of times fisher territory is selected
# uses presence/absence for female fishers on pixel (either 1 or 0) for each run
# then uses sum to count how many times each pixel is selected (out of numsims run, i.e., 100)
# the mean and se relate to the number of pixels (i.e., territories) selected per simulation

heatmap_output <- function(sim_out=sim_out, sim_order=sim_order, numsims=100, yrs_sim=10, TS=TS, name_out=name_out, dir_name=dir_name){
  # sim_out=C.w1_real.FEMALE; sim_order=2; numsims=100; yrs_sim=10; TS=12; name_out="QTSA_ex2"

  TS_full=paste0("TimeStep_",TS)

  # find out how many runs had at least one female adult fisher alive at end)
  tmp <- sim_output(sim_out=sim_out, sim=sim_order, numsims=numsims, yrs_sim=yrs_sim)

  fishers_to_start <- tmp %>% filter(TimeStep=="TimeStep_01") %>% summarise(numAF=mean(Count))

  Nozero.runs <- tmp %>% filter(TimeStep==TS_full) %>%
    group_by(Sim) %>%
    filter(Count!=0)

  tmp2 <- Nozero.runs %>% dplyr::select(Run)
  nozerosims <- tmp2$Run

  # set extent of raster the same as extent of initial world
  rw <- world2raster(sim_out[[1]]$land)

  r <- raster()
  r <- setExtent(r, rw, keepres=TRUE)

  r_list=list()

  # for simulations where at least one fisher survived
  for(i in 1:length(nozerosims)){
    ftmp <- as.data.frame(patchHere(sim_out[[1]]$land, sim_out[[sim_order]][[nozerosims[i]]][[TS]]))
    ftmp$Fisher <- 1
    ftmp.sf <- st_as_sf(ftmp, coords = c("pxcor", "pycor"))
    ftmp.sfp <- st_buffer(ftmp.sf, dist=.1)

    # r_list[[i]] <- rasterize(ftmp.sfp, r, field="Fisher", fun=rFun, background=0) # interim work around until terra and new raster package uploaded
    r_list[[i]] <- rasterize(ftmp.sfp, r, field="Fisher", background=0) # interim work around until terra and new raster package uploaded
  }

  r_zeroes <- raster()
  r_zeroes <- setExtent(r_zeroes, rw, keepres=TRUE)
  values(r_zeroes) <- 0

  r_zeroes_list=list()

  if(length(nozerosims)!=100){
    for(i in 1:(100-length(nozerosims))){
      r_zeroes_list[[i]] <- r_zeroes
    }
  }

  r_stack = stack(r_list, r_zeroes_list)
  r_stackApply <- stackApply(r_stack, indices=1, fun=sum)

  writeRaster(r_stackApply, file=paste0("out/",dir_name,"/rSim_",name_out,"_",round(sim_out[[sim_order-3]]$actual.prop.hab*100),"hab.tif"), bylayer=TRUE, overwrite=TRUE)

  Fisher_Nmean <- mean(r_stackApply@data@values)
  Fisher_Nse <- se(r_stackApply@data@values)

  suitable_habitat <- sum(sim_out[[1]]$land)
  total_habitat <- dim(sim_out[[1]]$land)[1]*dim(sim_out[[1]]$land)[2]

  Cairo(file=paste0("out/",dir_name,"/rHeatmap_",name_out,"_hab.PNG"), type="png", width=2200, height=2000,pointsize=15,bg="white",dpi=300)
  plot(r_stackApply, oma=c(2, 3, 5, 2))
  mytitle = paste0("Estimated Fisher Territories over ",numsims," Simulations")
  mysubtitle1 = paste0("Starting with ",fishers_to_start$numAF," fishers and ",round(suitable_habitat/total_habitat*100),"% habitat")
  mysubtitle2 = paste0("predicted ",round(Fisher_Nmean)," \u00B1 ",round(Fisher_Nse)," (mean \u00B1 1 SE) established fisher territories after ",yrs_sim," years.")
  mtext(side=3, line=3, at=-0.07, adj=0, cex=1, mytitle)
  mtext(side=3, line=2, at=-0.07, adj=0, cex=0.8, mysubtitle1)
  mtext(side=3, line=1, at=-0.07, adj=0, cex=0.8, mysubtitle2)
  dev.off()

  return(list(raster=r_stackApply, Fisher_Nmean=Fisher_Nmean, Fisher_Nse=Fisher_Nse, nozerosims=nozerosims))

}


################################################################################
# Run 100 simulations for each, save as objects
# Calculate mean # of adult females per cell at 10 years for each simulation to produce a heat map
# Create a figure with mean number of adult females (+/- SE or 95% CIs) for each time step and graph for each simulation

###--- load real world simulations (list of 1 run for 1 population from analysis script)
load("out/canBex.FEMALE.RData")

scenario1 <- canBex.FEMALE[[1]]
scenario2 <- canBex.FEMALE[[2]]
scenario3 <- canBex.FEMALE[[3]]
scenario4 <- canBex.FEMALE[[4]]

length(canBex.FEMALE)

### create output data from scenario output
scenario_output_function <- function(sim_out=sim_out, name_out=name_out,
                                     numsims=100, yrs_sim=10, Fpop="B",
                                     sim_order=2, TS=11){
  dir_name <- substr(name_out,1,nchar(name_out)-1)

  scenario_setup <- setup_plot(sim_out=sim_out, name_out=name_out)

  scenario_out <- ABM_fig_1sim(sim_out=sim_out, numsims=numsims, yrs_sim=yrs_sim, Fpop=Fpop)

  Cairo(file=paste0("out/",dir_name,"/IBM_",name_out,"_MeanSE.PNG"),type="png",width=3000,height=2200,pointsize=15,bg="white",dpi=300)
  print(scenario_out$sim.TS.plot_se)
  dev.off()

  Cairo(file=paste0("out/",dir_name,"/IBM_",name_out,"_MeanCL.PNG"),type="png",width=3000,height=2200,pointsize=15,bg="white",dpi=300)
  print(scenario_out$sim.TS.plot)
  dev.off()

  Cairo(file=paste0("out/",dir_name,"/IBM_",name_out,"_Saoi.PNG"),type="png",width=3000,height=2200,pointsize=15,bg="white",dpi=300)
  plot(sim_out[[1]]$land, legend=FALSE, main="Simulated Fisher Established Territories within Area of Interest")
  points(sim_out[[1]]$t0, pch = sim_out[[1]]$t0$shape, col = of(agents = sim_out[[1]]$t0, var = "color"))
  dev.off()

  scenario_heatmap <- heatmap_output(sim_out=sim_out, sim_order=sim_order, numsims=numsims, yrs_sim=yrs_sim, TS=TS,
                                     name_out=name_out, dir_name=dir_name)

  return(list(scenario_setup=scenario_setup, scenario_out=scenario_out, scenario_heatmap=scenario_heatmap))
}


dir.create(file.path(paste(getwd(),"canBex",sep="/out/")), recursive = TRUE)

for(i in 1:length(canBex.FEMALE)){
  scenario_output_function(sim_out=canBex.FEMALE[[i]], name_out=paste0("canBex",i))
}
