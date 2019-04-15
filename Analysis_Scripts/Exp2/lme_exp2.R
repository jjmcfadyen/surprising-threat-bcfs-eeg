library(lme4)
library(lsmeans)
library(ggplot2)
library(reshape2)
library(merTools)
#library(sjPlot)

setwd("G:/Jessica/Private/Experiments/2_Experiment2_CFS/Final/Exp2_EEG/behavioural_data_refined")

CSV <- read.csv("trial_data.csv", header = TRUE)

# mean-centre trial number
CSV$ExpTrial = CSV$ExpTrial - mean(CSV$ExpTrial)

# change category numbers to words
CSV$Gender[which(CSV$Gender == 1)] = "male"
CSV$Gender[which(CSV$Gender == 2)] = "female"
CSV$Order[which(CSV$Order == 1)] = "neutralfirst"
CSV$Order[which(CSV$Order == 2)] = "fearfulfirst"
CSV$Emotion[which(CSV$Emotion == 1)] = "neutral"
CSV$Emotion[which(CSV$Emotion == 2)] = "fearful"
CSV$Expectation[which(CSV$Expectation == 1)] = "expected"
CSV$Expectation[which(CSV$Expectation == 2)] = "unexpected"

# make models (2 by 2 repeated-measures ANOVA)
m0 <- lmer(RT ~ Gender + Order + (1|Subject) + mcSTAI, data=CSV, REML=FALSE)

m1 <- lmer(RT ~ Emotion + Gender + Order + (1|Subject) + (1|mcSTAI) + (1|ExpTrial), data=CSV, REML=FALSE)
m2 <- lmer(RT ~ Expectation + Gender + Order + (1|Subject) + (1|mcSTAI) + (1|ExpTrial), data=CSV, REML=FALSE)
m3 <- lmer(RT ~ Emotion + Expectation + Gender + Order + (1|Subject) + (1|mcSTAI) + (1|ExpTrial), data=CSV, REML=FALSE)
m4 <- lmer(RT ~ Emotion*Expectation + Gender + Order + (1|Subject) + (1|mcSTAI) + (1|ExpTrial), data=CSV, REML=FALSE)

anova(m0,m1,m2,m3,m4)

lsmeans(m4,pairwise~Expectation|Emotion)

CSV$fit <- predict(m4)
ggplot(CSV, aes(y=fit, x=RT, colour=interaction(Emotion,Expectation))) + geom_point()

write.csv(cbind(predictInterval(m0),predictInterval(m1),predictInterval(m2),predictInterval(m3),predictInterval(m4)),"model_fits.csv", row.names=FALSE)

rbeta <- ranef(m4)
fbeta <- fixef(m4)

write.csv(rbeta,"random_m4.csv",row.names=FALSE)
write.csv(fbeta,"fixed_m4.csv",row.names=FALSE)

ggplot(CSV, aes(y=RT, x=interaction(Emotion,Expectation), colour=Subject, group=Subject)) + geom_line() + scale_color_gradientn(colours=rainbow(5))

library(lattice)
dotplot(ranef(m4,condVar=TRUE))
dotplot(fixef(m4,condVar=TRUE))

R <- as.numeric(as.character(unlist(rbeta$Subject[[1]])))
print(paste(mean(R),sd(R)/sqrt(length(R)),sep=" "))

R <- as.numeric(as.character(unlist(rbeta$mcSTAI[[1]])))
print(paste(mean(R),sd(R)/sqrt(length(R)),sep=" "))

R <- as.numeric(as.character(unlist(rbeta$ExpTrial[[1]])))
print(paste(mean(R),sd(R)/sqrt(length(R)),sep=" "))
