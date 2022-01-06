library(lme4)
library(lmerTest)
library(interactions)
library(car)
library(emmeans)
library(tidyverse)
library(viridis)
library(ggpubr)

# ------------------------------------------------------------------------------
############
# Get data #
############
# ------------------------------------------------------------------------------

# Read in data
d <- read.csv("D:/bCFS_EEG_Reanalysis/results/erponsets_alltrials.csv")

# Convert to factors
d$Subject <- as.factor(d$Subject)
d$Emotion <- as.factor(d$Emotion)
d$Expectation <- as.factor(d$Expectation)
d$Condition <- as.factor(d$Condition)
d$Electrode <- as.factor(d$Electrode)
d$onsetType <- as.factor(d$onsetType)

# Create new variables
d$onsetPcnt <- d$onsetSL / d$RT

# ------------------------------------------------------------------------------

# Mean-centre variables
d$mcAnxiety  <- d$Anxiety  - mean(d$Anxiety)

# Z-score subject data
d <- d %>%
  group_by(Subject,Electrode,onsetType) %>%
  mutate(zonsetSL = scale(onsetSL),
         zonsetRL = scale(onsetRL),
         zonsetPcnt = scale(onsetPcnt),
         zSlope = scale(slope),
         zdrift = scale(drift),
         znondecision = scale(nondecision),
         zboundary = scale(boundary),
         zRT = scale(RT)) %>%
  as.data.frame()
 
# ------------------------------------------------------------------------------
################
# Descriptives #
################
# ------------------------------------------------------------------------------

pd <- d %>%
  group_by(Subject,onsetType) %>%
  summarise(onsetSL=mean(onsetSL)) %>%
  group_by(onsetType) %>%
  summarise(mean=round(mean(onsetSL),3),
            std=round(sd(onsetSL),3),
            min=round(min(onsetSL),3),
            max=round(max(onsetSL),3)) %>%
  as.data.frame()
pd

pd <- d %>%
  group_by(Subject,onsetType) %>%
  summarise(onsetPcnt=mean(onsetPcnt)) %>%
  group_by(onsetType) %>%
  summarise(mean=round(mean(onsetPcnt*100),2),
            std=round(sd(onsetPcnt*100),2),
            min=round(min(onsetPcnt*100),2),
            max=round(max(onsetPcnt*100),2)) %>%
  as.data.frame()
pd

pd <- d %>%
  group_by(Subject,onsetType) %>%
  summarise(onsetRL=mean(onsetRL)) %>%
  group_by(onsetType) %>%
  summarise(mean=round(mean(onsetRL),3),
            std=round(sd(onsetRL),3),
            min=round(min(onsetRL),3),
            max=round(max(onsetRL),3)) %>%
  as.data.frame()
pd

# ------------------------------------------------------------------------------
################
# Correlations #
################
# ------------------------------------------------------------------------------

pd <- d %>%
  group_by(Subject,Condition,onsetType) %>%
  summarise(Emotion=unique(Emotion),
            Expectation=unique(Expectation),
            zonsetSL=mean(zonsetSL),
            zonsetRL=mean(zonsetRL),
            zSlope=mean(zSlope),
            zRT=mean(zRT),
            znondecision=mean(znondecision),
            zdrift=mean(zdrift),
            zboundary=mean(zboundary),
            Anxiety=mean(Anxiety)) %>%
  as.data.frame()
pd <- pd %>% filter(!is.na(zonsetSL))
  
# # calculate Neutral - Fearful
# pd_emotion <- pd %>%
#   group_by(Subject,onsetType,Emotion) %>%
#   summarise_at(c("zonsetSL","zonsetRL","zSlope","zRT","znondecision","zdrift","zboundary"),mean) %>%
#   mutate(zonsetSL=zonsetSL-lag(zonsetSL),
#          zonsetRL=zonsetRL-lag(zonsetRL),
#          zSlope=zSlope-lag(zSlope),
#          zRT=zRT-lag(zRT),
#          znondecision=znondecision-lag(znondecision),
#          zdrift=zdrift-lag(zdrift),
#          zboundary=zboundary-lag(zboundary)) %>% 
#   select(-Emotion) %>%
#   filter(!is.na(zonsetSL)) %>%
#   as.data.frame
# 
# # calculate Expected - Unexpected
# pd_expectation <- pd %>%
#   group_by(Subject,onsetType,Expectation) %>%
#   summarise_at(c("zonsetSL","zonsetRL","zSlope","zRT","znondecision","zdrift","zboundary"),mean) %>%
#   mutate(zonsetSL=zonsetSL-lag(zonsetSL),
#          zonsetRL=zonsetRL-lag(zonsetRL),
#          zSlope=zSlope-lag(zSlope),
#          zRT=zRT-lag(zRT),
#          znondecision=znondecision-lag(znondecision),
#          zdrift=zdrift-lag(zdrift),
#          zboundary=zboundary-lag(zboundary)) %>% 
#   select(-Expectation) %>%
#   filter(!is.na(zonsetSL)) %>%
#   as.data.frame
# 
# # calculate Neutral pred - Fearful pred
# pd_interaction <- pd %>%
#   group_by(Subject,onsetType,Emotion) %>%
#   mutate(zonsetSL=zonsetSL-lag(zonsetSL),
#          zonsetRL=zonsetRL-lag(zonsetRL),
#          zSlope=zSlope-lag(zSlope),
#          zRT=zRT-lag(zRT),
#          znondecision=znondecision-lag(znondecision),
#          zdrift=zdrift-lag(zdrift),
#          zboundary=zboundary-lag(zboundary)) %>% 
#   select(-Expectation) %>%
#   filter(!is.na(zonsetSL)) %>%
#   group_by(Subject,onsetType) %>%
#   mutate(zonsetSL=zonsetSL-lag(zonsetSL),
#          zonsetRL=zonsetRL-lag(zonsetRL),
#          zSlope=zSlope-lag(zSlope),
#          zRT=zRT-lag(zRT),
#          znondecision=znondecision-lag(znondecision),
#          zdrift=zdrift-lag(zdrift),
#          zboundary=zboundary-lag(zboundary)) %>% 
#   select(-Emotion) %>%
#   filter(!is.na(zonsetSL)) %>%
#   as.data.frame

# make function to do correlations
summarisecorr <- function(testname,thisonset,thisdata) {

  thisd <- filter(thisdata,onsetType==thisonset)
  
  if (onsetType!="late") {
    p1 <- cor.test(thisd$zonsetSL,thisd$znondecision)
    p2 <- cor.test(thisd$zonsetSL,thisd$zdrift)
    p3 <- cor.test(thisd$zonsetSL,thisd$zboundary)
  } else {
    p1 <- cor.test(thisd$zonsetRL*(-1),thisd$znondecision)
    p2 <- cor.test(thisd$zonsetRL*(-1),thisd$zdrift)
    p3 <- cor.test(thisd$zonsetRL*(-1),thisd$zboundary)
  }
  
  stats <- list(p1,p2,p3)
  ntests <- length(stats)
  
  df <- data.frame(Name=rep(testname,ntests),
                   parameter=c("nondecision","drift","boundary"),
                   p=rep(NA,ntests),
                   r=rep(NA,ntests),
                   df=rep(NA,ntests))
  for (i in 1:ntests) {
    df$p[i] <- stats[[i]]$p.value
    df$r[i] <- stats[[i]]$estimate
    df$df[i] <- stats[[i]]$parameter
  }
  return(df)
}


# ALL CONDITIONS POOLED
# 0 to early change point
corrstats <- summarisecorr("initial","early",pd)

# early to late change point
corrstats <- rbind(corrstats,
                   summarisecorr("middle","diff",pd))

# late change point to response onset
corrstats <- rbind(corrstats,
                   summarisecorr("late","late",pd))

# correct for multiple comparisons
corrstats$mcp <- corrstats$p * nrow(corrstats)
corrstats$sig <- corrstats$mcp < .05


# plot
thispd <- pd

thispd$onsetType <- gsub("early","A_early",thispd$onsetType)
thispd$onsetType <- gsub("diff","B_diff",thispd$onsetType)
thispd$onsetType <- gsub("late","C_late",thispd$onsetType)

idx <- thispd$onsetType=="C_late"
thispd$zonsetSL[idx] <- thispd$zonsetRL[idx] * (-1)

p1 <- ggplot(data=thispd,aes(x=zonsetSL,y=znondecision,group=onsetType)) + 
  geom_point() +
  geom_smooth(method=lm,color="black") + 
  facet_wrap(~onsetType) + 
  theme_void() + 
  coord_cartesian(xlim=c(-0.4,0.4),ylim=c(-2,2)) + 
  geom_segment(aes(x=-.4,xend=.4,y=-2,yend=-2)) +
  geom_segment(aes(x=-.4,xend=-.4,y=-2,yend=2))

p2 <- ggplot(data=thispd,aes(x=zonsetSL,y=zdrift,group=onsetType)) + 
  geom_point() +
  geom_smooth(method=lm,color="black") + 
  facet_wrap(~onsetType) + 
  theme_void() + 
  coord_cartesian(xlim=c(-0.4,0.4),ylim=c(-2,2)) + 
  geom_segment(aes(x=-.4,xend=.4,y=-2,yend=-2)) +
  geom_segment(aes(x=-.4,xend=-.4,y=-2,yend=2))

p3 <- ggplot(data=thispd,aes(x=zonsetSL,y=zboundary,group=onsetType)) + 
  geom_point() +
  geom_smooth(method=lm,color="black") + 
  facet_wrap(~onsetType) + 
  theme_void() + 
  coord_cartesian(xlim=c(-0.4,0.4),ylim=c(-2,2)) + 
  geom_segment(aes(x=-.4,xend=.4,y=-2,yend=-2)) +
  geom_segment(aes(x=-.4,xend=-.4,y=-2,yend=2))

ggarrange(p1,p2,p3,ncol=1,nrow=3)

# # EMOTION
# # 0 to early change point
# corrstats <- summarisecorr("initial","early",pd_emotion)
# 
# # early to late change point
# corrstats <- rbind(corrstats,
#                    summarisecorr("middle","diff",pd_emotion))
# 
# # late change point to response onset
# corrstats <- rbind(corrstats,
#                    summarisecorr("late","late",pd_emotion))
# 
# # correct for multiple comparisons
# corrstats$mcp <- corrstats$p * nrow(corrstats)
# corrstats$sig <- corrstats$mcp < .05
# 
# 
# # EXPECTATION
# # 0 to early change point
# corrstats <- summarisecorr("initial","early",pd_expectation)
# 
# # early to late change point
# corrstats <- rbind(corrstats,
#                    summarisecorr("middle","diff",pd_expectation))
# 
# # late change point to response onset
# corrstats <- rbind(corrstats,
#                    summarisecorr("late","late",pd_expectation))
# 
# # correct for multiple comparisons
# corrstats$mcp <- corrstats$p * nrow(corrstats)
# corrstats$sig <- corrstats$mcp < .05
# 
# 
# # INTERACTION
# # 0 to early change point
# corrstats <- summarisecorr("initial","early",pd_interaction)
# 
# # early to late change point
# corrstats <- rbind(corrstats,
#                    summarisecorr("middle","diff",pd_interaction))
# 
# # late change point to response onset
# corrstats <- rbind(corrstats,
#                    summarisecorr("late","late",pd_interaction))
# 
# # correct for multiple comparisons
# corrstats$mcp <- corrstats$p * nrow(corrstats)
# corrstats$sig <- corrstats$mcp < .05

# ------------------------------------------------------------------------------
##############
# Run Models #
##############
# ------------------------------------------------------------------------------

# convert long to wide
md <- d %>%
  group_by(Subject,Trial,Electrode) %>%
  mutate(
    early_onsetSL = rep(onsetSL[onsetType=="early"],3),
    diff_onsetSL = rep(onsetSL[onsetType=="diff"],3),
    late_onsetSL = rep(onsetSL[onsetType=="late"],3),
    early_zonsetSL = rep(zonsetSL[onsetType=="early"],3),
    diff_zonsetSL = rep(zonsetSL[onsetType=="diff"],3),
    late_zonsetSL = rep(zonsetSL[onsetType=="late"],3),
    early_onsetRL = rep(onsetRL[onsetType=="early"],3),
    diff_onsetRL = rep(onsetRL[onsetType=="diff"],3),
    late_onsetRL = rep(onsetRL[onsetType=="late"],3),
    early_zonsetRL = rep(zonsetRL[onsetType=="early"],3),
    diff_zonsetRL = rep(zonsetRL[onsetType=="diff"],3),
    late_zonsetRL = rep(zonsetRL[onsetType=="late"],3),
    early_onsetPcnt = rep(onsetPcnt[onsetType=="early"],3),
    diff_onsetPcnt = rep(onsetPcnt[onsetType=="diff"],3),
    late_onsetPcnt = rep(onsetPcnt[onsetType=="late"],3)) %>%
  filter(onsetType=="early") %>%
  select(-onsetType) %>%
  ungroup() %>%
  as.data.frame()
md[,grepl("_onset",names(md))] <- lapply(md[,grepl("_onset",names(md))],function(x) scale(x,center=TRUE,scale=FALSE)) %>% as.data.frame()



m1 <- lmer(early_onsetSL ~ Emotion*Expectation + Electrode + early_onsetPcnt + (1|Subject),
                md,REML=F)
summary(m1)
vif(m1)
cat_plot(m1,pred=Emotion,modx=Expectation)
em <- emmeans(m1,pairwise~Emotion*Expectation)

m2 <- lmer(late_onsetSL ~ Emotion*Expectation + early_onsetSL + Electrode + early_onsetPcnt + late_onsetPcnt + (1|Subject),
           md,REML=F)
summary(m2)
vif(m2)
cat_plot(m2,pred=Emotion,modx=Expectation)
em <- emmeans(m2,pairwise~Emotion*Expectation)

