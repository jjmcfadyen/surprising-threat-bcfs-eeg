library(lmerTest)
library(lsmeans)
library(ggplot2)
library(viridis)
library(dplyr)

setwd("D:/Scratch/bCFS_EEG_Reanalysis/results")

exp1 <- read.csv("trial_data_exp1.csv", header = TRUE)
exp2 <- read.csv("trial_data_exp2.csv", header = TRUE)
exp1 <- exp1[,c(1:15,17)]

exp1_subjects <- unique(exp1$Subject) # how many subjects in exp 1?
exp2$Subject = exp2$Subject + tail(exp1_subjects,n=1) # add this number to exp. 2 so it goes up to 63 (no repeats of subject numbers)

CSV <- rbind(exp1,exp2)
CSV$Experiment <- c(rep("1",dim(exp1)[1]),
                    rep("2",dim(exp2)[1]))

# change category numbers to words
CSV$Gender[which(CSV$Gender == 1)] = "male"
CSV$Gender[which(CSV$Gender == 2)] = "female"
CSV$Order[which(CSV$Order == 1)] = "neutralfirst"
CSV$Order[which(CSV$Order == 2)] = "fearfulfirst"
CSV$Emotion[which(CSV$Emotion == 1)] = "neutral"
CSV$Emotion[which(CSV$Emotion == 2)] = "fearful"
CSV$Expectation[which(CSV$Expectation == 1)] = "expected"
CSV$Expectation[which(CSV$Expectation == 2)] = "unexpected"
CSV$Subject = as.character(CSV$Subject)

# mean-centre predictors
CSV$STAI = CSV$STAI - mean(CSV$STAI)
CSV$ExpTrial = CSV$ExpTrial - mean(CSV$ExpTrial)
CSV$Block = CSV$Block - mean(CSV$Block)

# compare models
m0 <- lmer(RT ~ (1+Block|Subject) + Block + Experiment + STAI, data=CSV, REML=FALSE)
m1 <- lmer(RT ~ Emotion + (1+Emotion+Block|Subject) + Block + Experiment + STAI, data=CSV, REML=FALSE)
m2 <- lmer(RT ~ Expectation + (1+Expectation+Block|Subject) + Block + Experiment + STAI, data=CSV, REML=FALSE)
m3 <- lmer(RT ~ Emotion + Expectation + (1+Emotion+Expectation+Block|Subject) + Block + Experiment + STAI, data=CSV, REML=FALSE)
m4 <- lmer(RT ~ Emotion*Expectation + (1+Emotion*Expectation+Block|Subject) + Block + Experiment + STAI, data=CSV, REML=FALSE)

anova(m0,m1,m2,m3,m4)

# WINNING MODEL = M4 (interaction model)
# WINNING MODEL = M4 (interaction model)
LSM <- as.data.frame(lsmeans(m4, pairwise~Expectation|Emotion))
actual_means <- CSV %>% group_by(Emotion, Expectation) %>% summarise(mean = mean(RT), sd = sd(RT)) # get actual data means
round(actual_means$mean,3)
round(actual_means$sd,3)

# Plot model fit

CSV$fit <- predict(m4)

library(ggalt)

ggplot(CSV, aes(y=fit, x=RT)) + 
  geom_point(size=0.5, alpha=.5, show.legend = FALSE) +
  stat_bkde2d(aes(fill=..level..),geom="polygon") +
  geom_abline(intercept=0,slope=1) + 
  scale_fill_viridis() + 
  theme_classic() + 
  coord_cartesian(xlim=c(.5,10),ylim=c(1,7))

library(MuMIn)
r.squaredGLMM(m4)

# Plot fixed effects
M <- summary(m4)
M <- as.data.frame(M$coefficients)
coeff <- data.frame(Coefficient = rownames(M), Estimate = M$Estimate, SE = M$`Std. Error`)
coeff <- coeff[2:dim(coeff)[1],]
coeff$upper <- coeff$Estimate + coeff$SE
coeff$lower <- coeff$Estimate - coeff$SE

ggplot(data=coeff,aes(x=Coefficient,y=Estimate)) + 
  geom_point(size=4) + 
  geom_errorbar(data=coeff,aes(x=Coefficient,ymin=lower,ymax=upper),width=.25, size=1) + 
  theme_classic() + 
  scale_x_discrete(labels=c("Block","Emotion","Interaction","Expectation","Anxiety")) + 
  geom_abline(slope=0,intercept=0)


## Plot LSM
CSV$Condition = rep(0,dim(CSV)[1])
CSV$Condition[CSV$Emotion == "neutral" & CSV$Expectation == "expected"] = 1
CSV$Condition[CSV$Emotion == "neutral" & CSV$Expectation == "unexpected"] = 2
CSV$Condition[CSV$Emotion == "fearful" & CSV$Expectation == "expected"] = 3
CSV$Condition[CSV$Emotion == "fearful" & CSV$Expectation == "unexpected"] = 4

CSV$Condition = as.character(CSV$Condition)

LSM <- LSM[c(which(LSM$lsmeans.Emotion == "neutral" & LSM$lsmeans.Expectation == "expected"), # put in right order
             which(LSM$lsmeans.Emotion == "neutral" & LSM$lsmeans.Expectation == "unexpected"),
             which(LSM$lsmeans.Emotion == "fearful" & LSM$lsmeans.Expectation == "expected"),
             which(LSM$lsmeans.Emotion == "fearful" & LSM$lsmeans.Expectation == "unexpected")),]
LSM$Condition <- 1:4
LSM$Upper <- LSM$lsmeans.lsmean + LSM$lsmeans.SE
LSM$Lower <- LSM$lsmeans.lsmean - LSM$lsmeans.SE

CSV$Subject <- sprintf("%02d", as.numeric(CSV$Subject))

mean_colour = "#FFC300"
ggplot() + 
  geom_point(data=CSV, aes(x = Condition, y = fit, color = Subject), position="jitter", alpha=.1) + 
  stat_summary(data=CSV, aes(x = Condition, y = fit, fill = as.character(Subject)), fun.y=mean, 
               size=4, pch=24, color="black", alpha=.7, position_dodge(width=.5), geom = "point") + 
  geom_errorbar(data=LSM, aes(x=Condition,ymin=Lower,ymax=Upper),width=0,size=1,color="black") + 
  geom_line(data=LSM, aes(x=Condition,y=lsmeans.lsmean),color="black", size=1.3) + 
  geom_point(data=LSM, aes(x=Condition,y=lsmeans.lsmean), size=5, pch=21, fill=mean_colour, color="black", stroke=2) +
  scale_color_viridis(discrete=TRUE) + scale_fill_viridis(discrete=TRUE) + 
  theme_classic() + 
  coord_cartesian(ylim=c(1,6.5))