library(lme4)
library(lmerTest)
library(interactions)
library(car)
library(emmeans)
library(tidyverse)
library(viridis)
library(extrafont)

# ------------------------------------------------------------------------------
############
# Get data #
############
# ------------------------------------------------------------------------------

# Read in data
d1 <- read.csv("D:/bCFS_EEG_Reanalysis/results/erponsets_alltrials.csv")
d1$OnsetType <- rep("Mode",nrow(d1))

d2 <- read.csv("D:/bCFS_EEG_Reanalysis/results/erponsets_control_alltrials.csv")
d2$OnsetType <- rep("Last",nrow(d2))

d <- rbind(d1,d2)

# Convert to factors
d$Subject <- as.factor(d$Subject)
d$Emotion <- as.factor(d$Emotion)
d$Expectation <- as.factor(d$Expectation)
d$Condition <- as.factor(d$Condition)

# ------------------------------------------------------------------------------

# Mean-centre variables
d$mcAnxiety  <- d$Anxiety  - mean(d$Anxiety)

# Z-score subject data
d <- d %>%
  group_by(Subject,OnsetType) %>%
  mutate(zPOonsetSL = scale(POonsetSL),
         zCPonsetSL = scale(CPonsetSL),
         zPOonsetRL = scale(POonsetRL),
         zCPonsetRL = scale(CPonsetRL),
         zPOdrift = scale(POdrift),
         zCPdrift = scale(CPdrift),
         znondecision = scale(nondecision),
         zdrift = scale(drift),
         zboundary = scale(boundary)) %>%
  as.data.frame()

# ------------------------------------------------------------------------------
##############
# Run Models #
##############
# ------------------------------------------------------------------------------

# Parietal-occipital, stimulus-locked
m1 <- lmer(POonsetSL ~ Emotion*Expectation + (1|Subject),d,REML=F)
summary(m1)
cat_plot(m1,pred=Emotion,modx=Expectation)

m2 <- lmer(CPonsetSL ~ Emotion*Expectation + (1|Subject),d,REML=F)
summary(m2)
cat_plot(m2,pred=Emotion,modx=Expectation)
  
# Parietal-occipital, response-locked
  