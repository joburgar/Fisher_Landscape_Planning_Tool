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
list.of.packages <- c("tidyverse", "survival","ggfortify","Cairo")
# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
#####################################################################################

###--- REPRODUCTIVE PARAMETERS
# Code to determine confidence intervals from Mean, Standard deviation, and N
# https://bookdown.org/logan_kelly/r_practice/p09.html

CI_from_meanSDn <- function(mean=mean, sd=sd, n=n, alpha=0.05){
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

  return(c(sample.mean, sample.sd, lower.bound, upper.bound))
}

# use on reproductive parameters from Rich and Eric's survival paper
# for Central Interior - denning rate mean = 0.54; sd = 0.41; n = 37
# for Central Interior - litter size mean = 1.7; sd = 0.73; n = 14
# for Boreal - denning rate mean = 0.75; sd = 0.39; n = 22
# for Boreal - litter size mean = 2.6; sd = 0.70; n = 18
drCmean <- 0.54; drCsd <- 0.41; drCn <- 37
lsCmean <- 1.7; lsCsd <- 0.73; lsCn = 14
drBmean <- 0.75; drBsd <- 0.39; drBn <- 22
lsBmean <- 2.6; lsBsd <- 0.70; lsBn <- 18


drC <- CI_from_meanSDn(mean=drCmean, sd=drCsd, n=drCn)
lsC <- CI_from_meanSDn(mean=lsCmean, sd=lsCsd, n=lsCn)
drB <- CI_from_meanSDn(mean=drBmean, sd=drBsd, n=drBn)
lsB <- CI_from_meanSDn(mean=lsBmean, sd=lsBsd, n=lsBn)


repro.CI <- cbind(drC, lsC, drB, lsB)
rownames(repro.CI) <- c("mean", "sd", "L50CI", "U50CI")
write.csv(as.data.frame(repro.CI), "data/repro_CI.csv")

SD_from_meanCIn <- function(mean=mean, CI=c(LCI, UCI), n=n, alpha=0.05){
  # SD = sqrt(n) * (UCI-LCI)/3.92 # assumes 95% CI (for 90% CI, should be 3.29 rather than 3.92)
  degrees.freedom = n - 1
  t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
  # print(t.score)

  sample.SD <- sqrt(n) * (CI[2]-CI[1])/t.score
  sample.SE <- sample.SD/(sqrt(n))

  return(list(mean=mean, sample.SD=sample.SD, sample.SE=sample.SE))
}

CFA_surv <- SD_from_meanCIn(mean=0.79, CI=c(0.61, 0.97), n=30)
CFJ_surv <- SD_from_meanCIn(mean=0.41, CI=c(0.06, 0.76), n=11)

BFA_surv <- SD_from_meanCIn(mean=0.86, CI=c(0.69, 1.00), n=18)
BFJ_surv <- SD_from_meanCIn(mean=0.54, CI=c(0.26, 1.00), n=11)


###--- SURVIVE
# Have the fisher survive one time step depending on their age and cohort
# Use the survival function output from Eric's latest survival analysis
# Cohorts are broken down by population (Boreal / Central Interior), sex (M/F), and ageclass (A/J)
# Write the function with defaults to Central Interior but ability to update as necessary

###--- SURVIVE
# load Eric Lofroth's survival data (recieved Dec 2021)
load("./data/fisher_survival.RData")
library(survival)

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
summary(km_full)


# no difference if running survfit with all data or just adult, so run with all data
# ignore subadult output after 2 years as not biologically meaningful
# extract survival (mean and SE) for each cohort (Pop + Sex + Ageclass) in 6 month increments (365/2=182.5)
# use these values as the probabilty of a fisher surviving in the IBM survival function
surv_times <- c(182, 365, 548, 730, 912, 1095, 1278, 1460, 1642, 1825, 2008, 2190, 2372, 2555, 2738, 2920)
surv_times <- c(365, 730, 1095, 1460, 1825, 2190, 2555, 2920)

km_full_summary <- summary(km_full, times=surv_times, extend=TRUE)

quantile(km_full)
class(km_full)

str(km_full_summary)
length(surv_times) # 16 times (every 6 months for 8 years); or 8 times (annually for 8 years)
length(km_full_summary$strata) # 8 strata, repeating 16 values; or 64 (8 strata repeating 8 values)

# [1] "Population=Boreal, Sex=Female, Ageclass=Adult   "           "Population=Boreal, Sex=Female, Ageclass=Subadult"           "Population=Boreal, Sex=Male  , Ageclass=Adult   "
# [4] "Population=Boreal, Sex=Male  , Ageclass=Subadult"           "Population=Central Interior, Sex=Female, Ageclass=Adult   " "Population=Central Interior, Sex=Female, Ageclass=Subadult"
# [7] "Population=Central Interior, Sex=Male  , Ageclass=Adult   " "Population=Central Interior, Sex=Male  , Ageclass=Subadult"

km_surv_estimates <- as.data.frame(cbind(km_full_summary$surv, km_full_summary$std.err, km_full_summary$lower, km_full_summary$upper))
colnames(km_surv_estimates) <- c("Surv","SE","L95CL","U95CL")
# km_surv_estimates$Cohort <- rep(c("BFA","BFJ","BMA","BMJ","CFA","CFJ","CMA","CMJ"),each=16) # add in the group, e.g., BFA = Boreal Population, Female, Adult)
km_surv_estimates$Cohort <- rep(c("BFA","BFJ","BMA","BMJ","CFA","CFJ","CMA","CMJ"),each=8) # add in the group, e.g., BFA = Boreal Population, Female, Adult)

km_surv_estimates$Time <- rep(surv_times, times=8) # add in the time interval for the estimates
# km_surv_estimates$Time_step <- rep(seq_len(16),times=8) # add in the time step (as per the IBM - every six months = 1 time step)
km_surv_estimates$Time_step <- rep(seq_len(8),times=8) # add in the time step (as per the IBM - every six months = 1 time step)

km_surv_estimates$Use <- 1 # default to 1 (use) to set up case_when function
# km_surv_estimates <- km_surv_estimates %>% mutate(Use = case_when(grepl("J", Cohort) & Time_step > 4 ~ 0, TRUE ~ Use)) # only consider juvenile time steps for 2 years
km_surv_estimates <- km_surv_estimates %>% mutate(Use = case_when(grepl("J", Cohort) & Time_step > 2 ~ 0, TRUE ~ Use)) # only consider juvenile time steps for 2 years
km_surv_estimates <- km_surv_estimates %>% mutate(Use = case_when(grepl("A", Cohort) & Time_step < 3 ~ 0, TRUE ~ Use)) # only consider juvenile time steps for 2 years

###--- NEED TO CONSIDER THAT Time_step should make intuitive sense for adults and juveniles. Created "age" and "age_6mnths" categories where adults become adult at time_step 5 so 5-16 = adult and 1-4 = juvenile
# km_surv_estimates <- km_surv_estimates %>% mutate(age_6mnths = case_when(grepl("A", Cohort) ~ as.numeric(Time_step) + 4, TRUE ~ as.numeric(Time_step))) # only consider juvenile time steps for 2 years
km_surv_estimates$age <- km_surv_estimates$Time_step
glimpse(km_surv_estimates)

# subset to estimates needed for survival function
km_surv_estimates <- km_surv_estimates %>% filter(Use==1 & age<8.5) %>% dplyr::select(-Use)
glimpse(km_surv_estimates)

# data check - to make sure it makes sense for each age class
km_surv_estimates %>% group_by(Cohort) %>% summarise(max(age))
km_surv_estimates %>% filter(grepl("J", Cohort))
km_surv_estimates %>% filter(grepl("A", Cohort))

km_surv_estimates$Fpop <- substr(km_surv_estimates$Cohort,1,1)
km_surv_estimates_F <- km_surv_estimates %>% arrange(Fpop, age) %>% filter(!grepl("M", Cohort))

write.csv(km_surv_estimates_F, "data/km_surv_estimates_F.csv", row.names = FALSE)

############################################
se <- function(x) sqrt(var(x)/length(x))

age_df <- read.csv("data/RF_age_distribution.csv")

Cairo(file=paste0("out/age_dist_plot.PNG"), type="png", width=4000, height=3000,pointsize=15,bg="white",dpi=300)
ggplot(age_df, aes(x=as.factor(Age), y=F_Age_Distribution))+
  geom_bar(stat="identity") +
  facet_wrap(~Total_Pop, scales = "free_y")+
  xlab("Age") + ylab("Female Age Distribution") +
  theme(legend.position = "bottom", legend.title =element_blank())
dev.off()

propF <- age_df %>% group_by(Total_Pop) %>% summarise(sum(F_Age_Distribution))
colnames(propF)[2] <- "Total_F"
propF$Prop_F <- propF$Total_F / propF$Total_Pop # 71% female in pop (overall)

age_df$F_Age_Distribution_prop <- age_df$F_Age_Distribution / (age_df$Total_Pop*0.71)

Cairo(file=paste0("out/age_dist_prop_plot.PNG"), type="png", width=4000, height=3000,pointsize=15,bg="white",dpi=300)
ggplot(age_df, aes(x=as.factor(Age), y=F_Age_Distribution_prop))+
  geom_bar(stat="identity") +
  facet_wrap(~Total_Pop, scales = "free_y")+
  xlab("Age") + ylab("Female Age Distribution")+
  theme(legend.position = "bottom", legend.title =element_blank())
dev.off()

age_prop <- age_df %>% group_by(Age) %>% summarise(mean=mean(F_Age_Distribution_prop), se=se(F_Age_Distribution_prop))

Cairo(file=paste0("out/age_dist_meanprop_plot.PNG"), type="png", width=2000, height=1800,pointsize=15,bg="white",dpi=300)
ggplot(age_prop, aes(x=as.factor(Age), y=mean, label=round(mean, digits=2)))+
  geom_bar(stat="identity") +
  geom_text(hjust=0, vjust=-0.5)+
  # geom_linerange(aes(ymin=mean-se, ymax=mean+se))+
  xlab("Age") + ylab("Proportion of Females (mean)")+
  theme(legend.position = "bottom", legend.title =element_blank())
dev.off()


Cairo(file=paste0("out/age_dist_RFprop_plot.PNG"), type="png", width=2000, height=1800,pointsize=15,bg="white",dpi=300)
ggplot(age_df %>% filter(Total_Pop==571), aes(x=as.factor(Age), y=F_Age_Distribution, label=round(F_Age_Distribution_prop, digits=2)))+
  geom_bar(stat="identity") +
  geom_text(hjust=0, vjust=-0.5)+
  xlab("Age") + ylab("Proportion of Females (RF final data)")+
  theme(legend.position = "bottom", legend.title =element_blank())
dev.off()


##################
# survival options
surv <- read.csv("data/updated_surv_estimates_9May2022.csv")
survF <- surv %>% filter(grepl("F", Cohort))
# survF$Age_class <- ifelse(grepl("A",survF$Cohort),"A","J")
# survF$Fpop <- ifelse(grepl("B",survF$Cohort),"Boreal","Columbian")
survFLEX <- as.data.frame(matrix(NA, 20, 2))
colnames(survFLEX) <- c("Age","Fpop")
survFLEX$Age <- rep(0:9,2)
survFLEX$Fpop <- rep(c("Boreal","Columbian"), each=10)

survFLEX <- survFLEX %>% mutate(Surv=case_when(Age==0 & Fpop=="Boreal" ~ 1,
                                               Age==1 & Fpop=="Boreal" ~ survF[survF$Cohort=="BFJ",]$Surv,
                                               Age>1 & Fpop=="Boreal" ~ survF[survF$Cohort=="BFA",]$Surv,
                                               Age==0 & Fpop=="Columbian" ~ 1,
                                               Age==1 & Fpop=="Columbian" ~ survF[survF$Cohort=="CFJ",]$Surv,
                                               Age>1 & Fpop=="Columbian" ~ survF[survF$Cohort=="CFA",]$Surv))


survFLEX <- survFLEX %>% mutate(SE=case_when(Age==0 & Fpop=="Boreal" ~ 0,
                                             Age==1 & Fpop=="Boreal" ~ survF[survF$Cohort=="BFJ",]$SE,
                                             Age>1 & Fpop=="Boreal" ~ survF[survF$Cohort=="BFA",]$SE,
                                             Age==0 & Fpop=="Columbian" ~ 0,
                                             Age==1 & Fpop=="Columbian" ~ survF[survF$Cohort=="CFJ",]$SE,
                                             Age>1 & Fpop=="Columbian" ~ survF[survF$Cohort=="CFA",]$SE))

survFLEX$Adj <- ifelse(survFLEX$Age<4,1, ifelse(survFLEX$Age<7,2, ifelse(survFLEX$Age<10,3)))
survFLEX$SE_adj <- survFLEX$SE * survFLEX$Adj
survFLEX$SurvLSE <- survFLEX$Surv - survFLEX$SE_adj
survFLEX$SurvHSE <- survFLEX$Surv + survFLEX$SE_adj
survFLEX$SurvHSE <- case_when(survFLEX$SurvHSE >1 ~ 1,TRUE ~ survFLEX$SurvHSE)


Cairo(file=paste0("out/surv_SE_adj_plot.PNG"), type="png", width=2000, height=1800,pointsize=15,bg="white",dpi=300)
ggplot(survFLEX, aes(x=as.factor(Age), y=Surv))+
  geom_point() +
  geom_linerange(aes(ymin=SurvLSE, ymax=SurvHSE))+
  xlab("Age") + ylab("Survival")+
  facet_wrap(~Fpop)+
  theme(legend.position = "bottom", legend.title =element_blank())
dev.off()
