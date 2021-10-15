library(lme4)
library(lmerTest)
library(interactions)
library(car)
library(emmeans)
library(tidyverse)
library(viridis)
library(extrafont)

# ------------------------------------------------------------------------------

experiment <- 3 # 1 or 2, or 3 to combine both
st <- "last" # "all", "first", or "last"

# ------------------------------------------------------------------------------
############
# Get data #
############
# ------------------------------------------------------------------------------

# Read in data
if (experiment==3){
  
  d1 <- read.csv("D:/bCFS_EEG_Reanalysis/data/exp1_behavioural.csv")
  d1$Experiment <- rep(1,nrow(d1))
  
  d2 <- read.csv("D:/bCFS_EEG_Reanalysis/data/exp2_behavioural.csv")
  d2$Subject <- d2$Subject + max(d1$Subject)
  d2$Experiment <- rep(2,nrow(d2))
  
  d <- rbind(d1,d2)
  d$Experiment <- as.factor(d$Experiment)
  
} else {
  d <- read.csv(paste("D:/bCFS_EEG_Reanalysis/data/exp",experiment,"_behavioural.csv",sep=""))
}

# ------------------------------------------------------------------------------

# Add info
d$Response[d$Response=="NaN"] = NA
d$Response <- as.factor(d$Response)

d$StartBlock <- rep(NA,nrow(d))
for (subject in unique(d$Subject)){
  tmp <- filter(d,Subject==subject)
  d$StartBlock[d$Subject==subject] = rep(tmp$Emotion[1],nrow(tmp))
}

d$BlockType <- rep(NA,nrow(d))
d$BlockType[d$Condition==1 | d$Condition==4] = "neutral"
d$BlockType[d$Condition==2 | d$Condition==3] = "fearful"

# ------------------------------------------------------------------------------

# Select data
if (st!="all") {
  d <- d %>% filter(Type==st | Type=="deviant")
}

d <- d %>% filter(Acc==1,Outlier==0,RT>.5)

if (experiment==2){
  d <- d %>% filter(Subject!=3)
} else if (experiment==3) {
  d <- d %>% filter(Subject!=(3+max(d1$Subject)),StateAnxiety!="NaN")
}

# Convert to factors
d$Subject <- as.factor(d$Subject)
d$Emotion <- as.factor(d$Emotion)
d$Expectation <- as.factor(d$Expectation)
d$Condition <- as.factor(d$Condition)

# ------------------------------------------------------------------------------

# Mean-centre variables
anxietyMeans <- c(mean(d$STAI),mean(d$StateAnxiety),mean(d$TraitAnxiety)) # in case of back-transforming later
d$STAI <- d$STAI - mean(d$STAI)
d$StateAnxiety  <- d$StateAnxiety  - mean(d$StateAnxiety)
d$TraitAnxiety  <- d$TraitAnxiety  - mean(d$TraitAnxiety)

if (experiment==3){
  d <- d %>% group_by(Experiment) %>% mutate(zRT = scale(RT))
  d$RT <- d$zRT
}

# Add "mismatch" response time for the deviants
d$Mismatch <- rep(NA,nrow(d))
for (s in unique(d$Subject)){
  idx <- d$Subject==s
  d$Mismatch[idx & d$Condition==2] <- d$RT[idx & d$Condition==2] - mean(d$RT[idx & d$Condition==3],na.rm=TRUE)
  d$Mismatch[idx & d$Condition==4] <- d$RT[idx & d$Condition==4] - mean(d$RT[idx & d$Condition==1],na.rm=TRUE)
}

# ------------------------------------------------------------------------------
##############
# Run Models #
##############
# ------------------------------------------------------------------------------
 
if (experiment < 3) {
  
  m0 <- lmer(RT ~ (1|Subject),d,REML=F)
  m <- lmer(RT ~ Emotion*Expectation + (1|Subject),
            d,REML=F)
  anova(m0,m)
  
} else {
  
  # d$AnxietyGroup <- rep(NA,nrow(d))
  # d$AnxietyGroup[d$TraitAnxiety < median(d$TraitAnxiety)] <- "low"
  # d$AnxietyGroup[d$TraitAnxiety >= median(d$TraitAnxiety)] <- "high"
  # d$AnxietyGroup <- as.factor(d$AnxietyGroup)
  
  m0 <- lmer(RT ~ Experiment + (1|Subject),d,REML=F)
  m1 <- lmer(RT ~ Emotion*Expectation + Experiment + (1|Subject),
             d,REML=F)
  m <- lmer(RT ~ Emotion*Expectation*TraitAnxiety + Experiment + (1|Subject),
            d,REML=F)
  anova(m0,m,m1)
  
}

summary(m)
vif(m)

# Create table for manuscript
S <- summary(m)
S <- as.data.frame(S$coefficients)

T <- data.frame(Col1 = rownames(S),
                Estimate = round(S$Estimate,3),
                SEM = round(S$`Std. Error`,3),
                t_upper = round(S$`t value`,3),
                p_lower = round(S$`Pr(>|t|)`,3))

if (experiment < 3) {

  e <- emmeans(m,pairwise~Emotion*Expectation)
  e <- as.data.frame(e$emmeans)
  T <- rbind(T,data.frame(Col1 = paste(e$Expectation,e$Emotion),
                          Estimate = round(e$emmean,3),
                          SEM = round(e$SE,3),
                          t_upper = round(e$`asymp.LCL`,3),
                          p_lower = round(e$`asymp.UCL`,3)))
  
  T <- T[c(1,2,3,4,6,8,5,7),]
  print(T, row.names = FALSE)
  
  emmeans(m,pairwise~Emotion)
  emmeans(m,pairwise~Expectation)
  emmeans(m,pairwise~Emotion*Expectation)
  
  cat_plot(m,pred=Emotion,modx=Expectation)
} else {
  
  e <- emmeans(m,pairwise~Emotion*Expectation*AnxietyGroup)
  e <- as.data.frame(e$emmeans)
  T <- rbind(T,data.frame(Col1 = paste(e$Expectation,e$Emotion,e$AnxietyGroup),
                          Estimate = round(e$emmean,3),
                          SEM = round(e$SE,3),
                          t_upper = round(e$`asymp.LCL`,3),
                          p_lower = round(e$`asymp.UCL`,3)))
  
  T <- T[c(1,2,3,4,5,6,7,8,14,16,13,15,10,12,9,11),]
  print(T, row.names = FALSE)
  
  emmeans(m,pairwise~Emotion*Expectation*AnxietyGroup)
  
  interact_plot(m,pred=TraitAnxiety,modx=Emotion,mod2=Expectation,
                interval=TRUE,vary.lty=FALSE,colors=c("#FF5D00","#00C9FF")) + 
    theme_classic()
  
  pred <- jtools::make_predictions(m,pred="TraitAnxiety",
                                   at=list(Emotion=unique(d$Emotion),Expectation=unique(d$Expectation)))
  write.csv
  
}

# ------------------------------------------------------------------------------

# Plot

# font_import()
loadfonts(device = "win")
setwd("C:/Users/jmcfadyen/OneDrive/Manuscripts/bCFS_EEG_Study/2021")

if (experiment < 3) {
  
  cmap <- c("#00E8FF","#336EFF","#FFBA30","#FF0000")
  
  md <- as.data.frame(emmeans(m,pairwise~Emotion*Expectation)$emmeans)
  md$Condition <- as.factor(c(3,1,4,2))
  
  # EMM
  g <- ggplot(data=md, aes(x=Condition,y=emmean)) + 
    geom_errorbar(aes(ymin=asymp.LCL,ymax=asymp.UCL), width=.1,size=1) + 
    geom_point(aes(color=Condition),size=5) + 
    theme(
      text=element_text(size=14, family="Arial"),
      panel.background  = element_blank(),
      axis.line = element_line(color="black",size=1),
      axis.text.x = element_text(color="black",size=12),
      axis.text.y = element_text(color="black",size=12),
      axis.title.x = element_text(margin=margin(t=20)),
      axis.title.y = element_text(margin=margin(r=20)),
      axis.ticks.length = unit(7, "points"),
      legend.position="none"
      ) + 
    scale_color_manual(values=cmap) + 
    ylab("Estimated Marginal Mean\nResponse Time (seconds)") + xlab("Condition") + 
    scale_x_discrete(labels=c("1"="Expected\nNeutral",
                              "2"="Unexpected\nNeutral",
                              "3"="Expected\nFearful",
                              "4"="Unexpected\nFearful"))
  
  if (experiment==1){
    g <- g + coord_cartesian(ylim=c(2.8,3.8)) + scale_y_continuous(breaks=seq(2.6,3.8,0.2))
  } else if (experiment==2){
    g <- g + coord_cartesian(ylim=c(1.7,2)) + scale_y_continuous(breaks=seq(1.7,2,.1))
  }
  
  ggsave(paste("exp",experiment,"_emmeans.svg",sep=""),g,width=150,height=122,units="mm",dpi=300)
  
  # Densities
  pd <- d %>% group_by(Subject,Condition) %>% summarise(RT=median(RT))
  
  g <- ggplot(data=pd, aes(x=Condition,y=RT,color=Condition,fill=Condition)) + 
    geom_violin(alpha=.15,scale="count",size=1) + 
    geom_jitter(width=.075,size=2,alpha=.4) +
    geom_boxplot(alpha=.33, width=.25,outlier.size=3,outlier.stroke=0,size=1) + 
    theme(
      text=element_text(size=14, family="Arial"),
      panel.background  = element_blank(),
      axis.line = element_line(color="black",size=1),
      axis.text.x = element_text(color="black",size=12),
      axis.text.y = element_text(color="black",size=12),
      axis.title.x = element_text(margin=margin(t=20)),
      axis.title.y = element_text(margin=margin(r=20)),
      axis.ticks.length = unit(7, "points"),
      legend.position="none"
    ) + 
    scale_color_manual(values=cmap) + scale_fill_manual(values=cmap) +
    ylab("Median Response Time (seconds)") + xlab("Condition") + 
    scale_x_discrete(labels=c("1"="Expected\nNeutral",
                              "2"="Unexpected\nNeutral",
                              "3"="Expected\nFearful",
                              "4"="Unexpected\nFearful"))
  
  if (experiment==1){
    g <- g + coord_cartesian(ylim=c(1.5,6)) + scale_y_continuous(breaks=seq(1.5,6,0.5))
  } else if (experiment==2) {
    g <- g + coord_cartesian(ylim=c(1,2.5)) + scale_y_continuous(breaks=seq(1,2.5,0.25))
  }
  
  ggsave(paste("exp",experiment,"_rtmedians.svg",sep=""),g,width=150,height=122,units="mm",dpi=300)
  
} else if (experiment==3) {
  
  medians <- data.frame(Median=c(median(d$TraitAnxiety),
                                 median(d$TraitAnxiety[d$AnxietyGroup=="low"]),
                                 median(d$TraitAnxiety[d$AnxietyGroup=="high"])),
                        Group=c("all","low","high"))
  
  # HISTOGRAM
  g <- ggplot(data=d,aes(x=TraitAnxiety,fill=AnxietyGroup)) + 
    geom_histogram(aes(y=..density..),alpha=0.5,position="identity",color="black") + 
    geom_vline(data=medians,aes(xintercept=Median,color=Group)) + 
    theme(
      text=element_text(size=14, family="Arial"),
      panel.background  = element_blank(),
      axis.line = element_line(color="black",size=1),
      axis.text.x = element_text(color="black",size=12),
      axis.text.y = element_text(color="black",size=12),
      axis.title.x = element_text(margin=margin(t=20)),
      axis.title.y = element_text(margin=margin(r=20)),
      axis.ticks.length = unit(7, "points"),
      legend.position="none"
    ) + 
    coord_cartesian(xlim=c(-20,25),ylim=c(0,0.12)) + 
    scale_x_continuous(breaks=seq(-20,25,5)) + scale_y_continuous(breaks=seq(0,0.12,0.02))
    
  ggsave("anxiety_histogram.svg",g,width=150,height=122,units="mm",dpi=300)
  
  # EMMEANS
  md <- as.data.frame(emmeans(m,pairwise~Emotion*Expectation*AnxietyGroup)$emmeans)
  md$Condition <- as.factor(rep(c(3,1,4,2),2))
  
  cmap <- c("#00E8FF","#336EFF","#FFBA30","#FF0000")
  
  # EMM
  g <- ggplot(data=md, aes(x=Condition,y=emmean)) + 
    geom_errorbar(aes(ymin=asymp.LCL,ymax=asymp.UCL), width=.1,size=1) + 
    geom_point(aes(color=Condition),size=5) + 
    facet_wrap(~AnxietyGroup) +
    theme(
      text=element_text(size=14, family="Arial"),
      panel.background  = element_blank(),
      axis.line = element_line(color="black",size=1),
      axis.text.x = element_text(color="black",size=12),
      axis.text.y = element_text(color="black",size=12),
      axis.title.x = element_text(margin=margin(t=20)),
      axis.title.y = element_text(margin=margin(r=20)),
      axis.ticks.length = unit(7, "points"),
      legend.position="none"
    ) + 
    scale_color_manual(values=cmap) + 
    ylab("Estimated Marginal Mean\nResponse Time (seconds)") + xlab("Condition") + 
    scale_x_discrete(labels=c("1"="Expected\nNeutral",
                              "2"="Unexpected\nNeutral",
                              "3"="Expected\nFearful",
                              "4"="Unexpected\nFearful")) +
    coord_cartesian(ylim=c(-0.3,0.3)) + 
    scale_y_continuous(breaks=seq(-0.3,0.3,0.1))
  
  ggsave("allexp_emmeans.svg",g,width=150,height=122,units="mm",dpi=300)
}


