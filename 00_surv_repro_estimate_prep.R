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
# 00_surv_repro_estimate_prep.R
# script to prepare probabilities, taken from Rich and Eric's survival paper analyses
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 09-Dec-2021
#####################################################################################
version$major
version$minor
R_version <- paste0("R-",version$major,".",version$minor)

.libPaths(paste0("C:/Program Files/R/",R_version,"/library")) # to ensure reading/writing libraries from C drive
tz = Sys.timezone() # specify timezone in BC

# Load Packages
list.of.packages <- c("tidyverse", "survival","ggfortify")
# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
#####################################################################################

###--- REPRODUCTIVE PARAMETERS
# Code to determine confidence intervals from Mean, Standard deviation, and N
# https://bookdown.org/logan_kelly/r_practice/p09.html

CI_from_meanSDn <- function(mean=mean, sd=sd, n=n, alpha=0.5){
  sample.mean <- mean
  # print(sample.mean)

  sample.n <- n
  sample.sd <- sd
  sample.se <- sample.sd/sqrt(sample.n)
  # print(sample.se)

  alpha <- alpha
  degrees.freedom = sample.n - 1
  t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
  # print(t.score)

  margin.error <- t.score * sample.se
  lower.bound <- sample.mean - margin.error
  upper.bound <- sample.mean + margin.error
  # print(c(lower.bound,upper.bound))

  return(c(lower.bound, upper.bound))
}

# use on reproductive parameters from Rich and Eric's survival paper
# for Central Interior - denning rate mean = 0.54; sd = 0.41; n = 37
# for Central Interior - litter size mean = 1.7; sd = 0.73; n = 14
# for Boreal - denning rate mean = 0.75; sd = 0.39; n = 22
# for Boreal - litter size mean = 2.6; sd = 0.70; n = 18
drCImean <- 0.54; drCIsd <- 0.41; drCIn <- 37
lsCImean <- 1.7; lsCIsd <- 0.73; lsCIn = 14
drBmean <- 0.75; drBsd <- 0.39; drBn <- 22
lsBmean <- 2.6; lsBsd <- 0.70; lsBn <- 18


drCI <- CI_from_meanSDn(mean=drCImean, sd=drCIsd, n=drCIn)
lsCI <- CI_from_meanSDn(mean=lsCImean, sd=lsCIsd, n=lsCIn)
drB <- CI_from_meanSDn(mean=drBmean, sd=drBsd, n=drBn)
lsB <- CI_from_meanSDn(mean=lsBmean, sd=lsBsd, n=lsBn)


repro.CI <- cbind(drCI, lsCI, drB, lsB)
rownames(repro.CI) <- c("L95CI", "U95CI")
write.csv(as.data.frame(repro.CI), "data/repro.CI.csv")

###--- SURVIVE
# Have the fisher survive one time step depending on their age and cohort
# Use the survival function output from Eric's latest survival analysis
# Cohorts are broken down by population (Boreal / Central Interior), sex (M/F), and ageclass (A/J)
# Write the function with defaults to Central Interior but ability to update as necessary

###--- SURVIVE
# load Eric Lofroth's survival data (recieved Dec 2021)
load("./data/fisher_survival.RData")
library(survival)
?survfit
summary(fishersurvival)
# Call: survfit(formula = fishersurv ~ Population + Sex + Ageclass, data = fisher1)
library(ggfortify)
cox_full <- coxph(fishersurv ~ Population + Sex + Ageclass, data=fisher1)
cox_nosex <- coxph(fishersurv ~ Population + Ageclass, data=fisher1)
cox_nopop <- coxph(fishersurv ~ Sex + Ageclass, data=fisher1)
cox_noac <- coxph(fishersurv ~ Population + Sex, data=fisher1)
cox_pop <- coxph(fishersurv ~ Population, data=fisher1)
cox_sex <- coxph(fishersurv ~ Sex, data=fisher1)
cox_ac <- coxph(fishersurv ~ Ageclass, data=fisher1)
anova(cox_full, cox_nosex, cox_noac, cox_ac, cox_pop, cox_nopop, cox_sex)

summary(cox_full)
summary(cox_nosex)
cox_fit <- survfit(cox_full)
autoplot(cox_fit)

summary(fishersurv)
fisher1 %>% group_by(Population, Sex, Ageclass) %>% summarise(mean(DaysMonitored))
fisher1$rownum <- as.numeric(rownames(fisher1))
fisher_adult <- fisher1 %>% filter(Ageclass=="Adult")
fisher_adult$rownum
glimpse(fisher_adult)
fishersurv_adult <- fishersurv[fisher_adult$rownum,]

km_full <- survfit(formula = fishersurv ~ Population + Sex + Ageclass, data=fisher1)
km_adult <- survfit(formula = fishersurv_adult ~ Population + Sex, data=fisher_adult)

# no difference if running survfit with all data or just adult, so run with all data
# ignore subadult output after 2 years as not biologically meaningful
# extract survival (mean and SE) for each cohort (Pop + Sex + Ageclass) in 6 month increments (365/2=182.5)
# use these values as the probabilty of a fisher surviving in the IBM survival function
surv_times <- c(182, 365, 548, 730, 912, 1095, 1278, 1460, 1642, 1825, 2008, 2190, 2372, 2555, 2738, 2920)

km_full_summary <- summary(km_full, times=surv_times, extend=TRUE)

quantile(km_full)
class(km_full)

str(km_full_summary)
length(surv_times) # 16 times (every 6 months for 8 years)
length(km_full_summary$strata) # 8 strata, repeating 16 values

# [1] "Population=Boreal, Sex=Female, Ageclass=Adult   "           "Population=Boreal, Sex=Female, Ageclass=Subadult"           "Population=Boreal, Sex=Male  , Ageclass=Adult   "
# [4] "Population=Boreal, Sex=Male  , Ageclass=Subadult"           "Population=Central Interior, Sex=Female, Ageclass=Adult   " "Population=Central Interior, Sex=Female, Ageclass=Subadult"
# [7] "Population=Central Interior, Sex=Male  , Ageclass=Adult   " "Population=Central Interior, Sex=Male  , Ageclass=Subadult"

km_surv_estimates <- as.data.frame(cbind(km_full_summary$surv, km_full_summary$std.err, km_full_summary$lower, km_full_summary$upper))
colnames(km_surv_estimates) <- c("Surv","SE","L95CL","U95CL")
km_surv_estimates$Cohort <- rep(c("BFA","BFJ","BMA","BMJ","CFA","CFJ","CMA","CMJ"),each=16) # add in the group, e.g., BFA = Boreal Population, Female, Adult)
km_surv_estimates$Time <- rep(surv_times, times=8) # add in the time interval for the estimates
km_surv_estimates$Time_step <- rep(seq_len(16),times=8) # add in the time step (as per the IBM - every six months = 1 time step)
km_surv_estimates$Use <- 1 # default to 1 (use) to set up case_when function
km_surv_estimates <- km_surv_estimates %>% mutate(Use = case_when(grepl("J", Cohort) & Time_step > 4 ~ 0, TRUE ~ Use)) # only consider juvenile time steps for 2 years

###--- NEED TO CONSIDER THAT Time_step should make intuitive sense for adults and juveniles. Created "age" and "age_6mnths" categories where adults become adult at time_step 5 so 5-16 = adult and 1-4 = juvenile
km_surv_estimates <- km_surv_estimates %>% mutate(age_6mnths = case_when(grepl("A", Cohort) ~ as.numeric(Time_step) + 4, TRUE ~ as.numeric(Time_step))) # only consider juvenile time steps for 2 years
km_surv_estimates$age <- km_surv_estimates$age_6mnths / 2
glimpse(km_surv_estimates)

write.csv(km_surv_estimates, "data/km_surv_estimates.csv", row.names = FALSE)

