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
  d$Experiment <- rep(experiment,nrow(d))
}

# Ignore subject who had too many missed trials
rmsub <- c()
for (e in unique(d$Experiment)) {
  
  thisd <- d %>% 
    filter(Experiment==e) %>% 
    group_by(Subject) %>% 
    summarise(n = n(),
              missed = sum(is.na(RT)),
              proportion = round(missed/n,2)) %>% 
    as.data.frame()
  
  z <- as.vector(abs(scale(thisd$proportion)))
  
  if (!all(is.na(z))){
    rmsub <- c(rmsub,as.vector(thisd$Subject[z>4]))
  }
}

d <- filter(d,Subject!=rmsub)

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
  m1 <- lmer(RT ~ Emotion + (1|Subject),d,REML=F)
  m2 <- lmer(RT ~ Expectation + (1|Subject),d,REML=F)
  m3 <- lmer(RT ~ Emotion + Expectation + (1|Subject),d,REML=F)
  m4 <- lmer(RT ~ Emotion*Expectation + (1|Subject),d,REML=F)
  anova(m0,m1,m2,m3,m4)
  
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
  
  emtrends(m, pairwise~Emotion*Expectation, var="TraitAnxiety")
  
  # save subject means for Bayes analysis in JASP
  sm <- d %>% 
    group_by(Subject,Emotion,Expectation) %>%
    summarise(RT=mean(RT),
              TraitAnxiety=unique(TraitAnxiety)) %>%
    as.data.frame()
  sm <- data.frame(EN=sm$RT[sm$Expectation=="expected" & sm$Emotion=="neutral"],
                   UN=sm$RT[sm$Expectation=="unexpected" & sm$Emotion=="neutral"],
                   EF=sm$RT[sm$Expectation=="expected" & sm$Emotion=="fearful"],
                   UF=sm$RT[sm$Expectation=="unexpected" & sm$Emotion=="fearful"],
                   Anxiety=sm$TraitAnxiety[sm$Expectation=="expected" & sm$Emotion=="neutral"])
  sm$NeutPred <- sm$UN-sm$EN
  sm$FearPred <- sm$UF-sm$EF
  write.csv(sm,"exp3.csv")
  
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
  
  ss <- sim_slopes(m,pred=TraitAnxiety,modx=Emotion,mod2=Expectation)
  
  interact_plot(m,pred=TraitAnxiety,modx=Emotion,mod2=Expectation,
                interval=TRUE,vary.lty=FALSE,colors=c("#FF5D00","#00C9FF")) + 
    theme_classic()
  
  pred <- jtools::make_predictions(m,pred="TraitAnxiety",
                                   at=list(Emotion=unique(d$Emotion),Expectation=unique(d$Expectation)))
  
  T <- data.frame()
  for (emotion in c("neutral","fearful")) {
    for (expectation in c("expected","unexpected")) {
      tmp <- filter(pred,Emotion==emotion,Expectation==expectation)
      tmp <- tmp[order(tmp$TraitAnxiety),]
      thisname <- paste(substr(str_to_title(expectation),1,1),substr(str_to_title(emotion),1,1),sep="")
      thisT <- data.frame(Anxiety=tmp$TraitAnxiety,
                          upper=tmp$ymax,
                          mean=tmp$RT,
                          lower=tmp$ymin)
      colnames(thisT) <- c("Anxiety",
                           paste("upper",thisname,sep=""),
                           thisname,
                           paste("lower",thisname,sep=""))
      if (nrow(T)==0){
        T <- thisT
      } else {
        T <- cbind(T,thisT)
      }
    }
  }
  write.csv(T,"anxslopes.csv")
  
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

# --------------------------------------------------------------------------------------
# DDM

experiment <- 3

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
  d$Experiment <- rep(experiment,nrow(d))
}

# Ignore subject who had too many missed trials
rmsub <- c()
for (e in unique(d$Experiment)) {
  
  thisd <- d %>% 
    filter(Experiment==e) %>% 
    group_by(Subject) %>% 
    summarise(n = n(),
              missed = sum(is.na(RT)),
              proportion = round(missed/n,2)) %>% 
    as.data.frame()
  
  z <- as.vector(abs(scale(thisd$proportion)))
  
  if (!all(is.na(z))){
    rmsub <- c(rmsub,as.vector(thisd$Subject[z>4]))
  }
}

if (!all(is.null(rmsub))){
  d <- filter(d,Subject!=rmsub)
}

if (experiment==3){
  d <- d %>% group_by(Experiment) %>% mutate(zRT = scale(RT))
  d$RT <- d$zRT
}


# RUN EZ-DDM
logit <- function(x){
  y <- log(x/(1-x))
  return(y)
}

sign <- function(x){
  if (x>0){
    return(1)
  } else if (x==0){
    return(0)
  } else {
    return(-1)
  }
}

ezddm <- function(p,MRT,VRT,s){
  if (p==0 | p==1){
    error("Accuracy cannot be 0% or 100%")
  }
  s2     = s*s;
  L      = logit(p);
  x      = L*(L*p*p - L*p + p - 0.5) / VRT;
  v      = sign(p-0.5)*s*x^(1/4);
  a      = s2*logit(p)/v;
  y      = -v*a/s2;
  MDT    = (a/(2*v))*(1-exp(y))/(1+exp(y));
  Ter    = MRT - MDT;
  return (c(v,a,Ter));
}



T <- data.frame()
for (subject in unique(d$Subject)){
  
  thissum <- d %>% filter(Subject==subject,!is.na(RT)) %>% group_by(Emotion,Expectation) %>% summarise(acc=mean(Acc))
  if (any(thissum$acc==0) | any(thissum$acc==1)){
    print(paste("Skipping subject ",subject,sep=""))
  } else {
  
    for (condition in 1:4){
      
      thisd <- filter(d,Subject==subject)
      if (condition==1){
        thisd <- filter(thisd,Emotion=="neutral",Expectation=="expected")
      } else if (condition==2){
        thisd <- filter(thisd,Emotion=="neutral",Expectation=="unexpected")
      } else if (condition==3){
        thisd <- filter(thisd,Emotion=="fearful",Expectation=="expected")
      } else if (condition==4){
        thisd <- filter(thisd,Emotion=="fearful",Expectation=="unexpected")
      }
      thisd <- filter(thisd,!is.na(thisd$RT))
      
      params <- ezddm(mean(thisd$Acc),mean(thisd$RT),var(thisd$RT*1000)/(1000*1000),0.1)
      
      thisT <- data.frame(Subject=subject,
                          Experiment=thisd$Experiment[1],
                          TraitAnxiety=thisd$TraitAnxiety[1],
                          Emotion=thisd$Emotion[1],
                          Expectation=thisd$Expectation[1],
                          RT=mean(thisd$RT),
                          Drift=params[1],
                          Boundary=params[2],
                          Nondecision=params[3])
      
      if (nrow(T)==0){
        T <- thisT
      } else {
        T <- rbind(T,thisT)
      }
    }
  }
}

T$Condition <- rep(NA,nrow(T))
T$Condition[T$Emotion=="neutral" & T$Expectation=="expected"] = 1
T$Condition[T$Emotion=="neutral" & T$Expectation=="unexpected"] = 2
T$Condition[T$Emotion=="fearful" & T$Expectation=="expected"] = 3
T$Condition[T$Emotion=="fearful" & T$Expectation=="unexpected"] = 4

# z-score
T <- T %>%
  group_by(Subject) %>%
  mutate(Drift=scale(Drift),
         Boundary=scale(Boundary),
         Nondecision=scale(Nondecision))

Tsum <- T %>%
  group_by(Condition) %>%
  summarise(
    n = n(),
    sRT=sd(RT,na.rm=TRUE)/sqrt(n),
    RT=mean(RT,na.rm=TRUE),
    sDrift=sd(Drift,na.rm=TRUE)/sqrt(n),
    Drift=mean(Drift,na.rm=TRUE),
    sBoundary=sd(Boundary,na.rm=TRUE)/sqrt(n),
    Boundary=mean(Boundary,na.rm=TRUE),
    sNondecision=sd(Nondecision,na.rm=TRUE)/sqrt(n),
    Nondecision=mean(Nondecision,na.rm=TRUE)
            )

ggplot() + 
  geom_point(data=T,aes(x=Condition,y=Drift)) + 
  geom_errorbar(data=Tsum,aes(x=Condition,ymin=Drift-sDrift,ymax=Drift+sDrift),width=0.15) + 
  geom_point(data=Tsum,aes(x=Condition,y=Drift),size=6) 

ggplot() + 
  geom_point(data=T,aes(x=Condition,y=Boundary)) + 
  geom_errorbar(data=Tsum,aes(x=Condition,ymin=Boundary-sBoundary,ymax=Boundary+sBoundary),width=0.15) + 
  geom_point(data=Tsum,aes(x=Condition,y=Boundary),size=6) 

ggplot() + 
  geom_point(data=T,aes(x=Condition,y=Nondecision)) + 
  geom_errorbar(data=Tsum,aes(x=Condition,ymin=Nondecision-sNondecision,ymax=Nondecision+sNondecision),width=0.15) + 
  geom_point(data=Tsum,aes(x=Condition,y=Nondecision),size=6) 

ggplot() + 
  geom_point(data=T,aes(x=Condition,y=RT)) + 
  geom_point(data=Tsum,aes(x=Condition,y=RT),size=6) 
  

csvtable <- data.frame()
for (condition in 1:4) {
  tmp <- filter(T,Condition==condition)
  tmp <- tmp[order(tmp$Subject),]
  thisname <- paste("C",condition,sep="")
  thisT <- data.frame(Subject=tmp$Subject,
                      Experiment=tmp$Experiment,
                      Anxiety=tmp$TraitAnxiety,
                      drift=tmp$Drift,
                      boundary=tmp$Boundary,
                      nondecision=tmp$Nondecision)
  colnames(thisT) <- c("Subject","Experiment","Anxiety",
                       paste("drift_",thisname,sep=""),
                       paste("boundary_",thisname,sep=""),
                       paste("nondecision_",thisname,sep=""))
  if (nrow(csvtable)==0){
    csvtable <- thisT
  } else {
    csvtable <- cbind(csvtable,thisT[,seq(4,ncol(thisT))])
  }
}

write.csv(csvtable,"ddmparams.csv",row.names=FALSE)




# T$Drift <- scale(T$Drift,center=TRUE,scale=FALSE)
# T$Boundary <- scale(T$Boundary,center=TRUE,scale=FALSE)
# T$Nondecision <- scale(T$Nondecision,center=TRUE,scale=FALSE)
# T$TraitAnxiety <- scale(T$TraitAnxiety,center=TRUE,scale=FALSE)

m_drift <- lmer(Drift ~ Emotion*Expectation + (1|Subject),T,REML=F)
summary(m_drift)

m_boundary <- lmer(Boundary ~ Emotion*Expectation + (1|Subject),T,REML=F)
summary(m_boundary)

m_nondecision <- lmer(Nondecision ~ Emotion*Expectation + (1|Subject),T,REML=F)
summary(m_nondecision)

m_rt <- lmer(RT ~ Emotion*Expectation + (1|Subject),T,REML=F)
summary(m_rt)




