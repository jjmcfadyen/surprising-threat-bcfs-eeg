library(tidyverse)

#############################################
### LOAD DATA ###############################
#############################################

# Read in long data
setwd("D:/Scratch/bCFS_EEG_Reanalysis/results")
CSV <- read.csv("trial_data_exp1.csv")

# Transform & clean data for LBA modelling
CSV$Condition = rep(1,dim(CSV)[1])
CSV$Condition[CSV$Emotion == 1 & CSV$Expectation == 1] = 1
CSV$Condition[CSV$Emotion == 1 & CSV$Expectation == 2] = 2
CSV$Condition[CSV$Emotion == 2 & CSV$Expectation == 1] = 3
CSV$Condition[CSV$Emotion == 2 & CSV$Expectation == 2] = 4

D <- CSV %>% mutate(blocknumber=Block,
                       trialnumber=ExpTrial,
                       s=factor(as.numeric(as.factor(Subject))),
                       St=factor(Orientation,levels=1:2,labels=c("SL","SR")),
                       Cond=factor(Condition,levels=1:4,labels=c("EN","UN","EF","UF")),
                       R=factor(Response,levels=1:2,labels=c("RL","RR")),
                       RT=RT) %>%
  arrange(s,blocknumber,trialnumber) %>%
  select(s,St,Cond,R,RT,blocknumber,trialnumber)

D <- D %>% filter(RT > .5, RT < 10)

# check if any subjects don't have stimulus orientation information
keep_subjects <- unique(D$s[which(!is.na(D$St))])
D <- D %>% mutate(s=factor(s,levels=keep_subjects,labels=as.cha
D <- D[which(D$s == keep_subjects),]racter(keep_subjects)))

#############################################
### MODEL ###################################
#############################################

# Load in Dynamic Models of Choice
setwd("D:/DMC/DMC-MBN18/dmc")
source("dmc.R")
load_model("lba","lba_B.R")

#--------------------------------
#set up model to fit
factors=list(St=c("SL","SR"),Cond=c("EN","UN","EF","UF"))
responses=c("RL","RR")
match.map=list(M=list(SL="RL",SR="RR"))
consts <-c(sd_v=1,st0=1)
p.map=list(A="Cond",B="Cond",mean_v="Cond",sd_v="1",t0="1",st0="1")
model <- model.dmc(type="norm",constants=consts,p.map=p.map,
                   match.map=match.map,factors=factors,responses=responses)

data.model <- data.model.dmc(D, model)

#--------------------------------
#set priors
pop.mean <- c(rep(1,length(attr(model,"p.vector"))-1),0.2)
names(pop.mean) <- names(attr(model,"p.vector"))

pop.prior <- prior.p.dmc(
  dists = rep("tnorm",length(pop.mean)),
  p1=pop.mean,
  p2=c(rep(.1,length(grep("A",names(pop.mean)))+length(grep("B",names(pop.mean)))),
       rep(.2,length(grep("mean_v",names(pop.mean)))),   
       rep(.1,length(grep("sd_v",names(pop.mean)))),
       rep(.05,length(grep("t0",names(pop.mean))))),
  lower=c(rep(0,length(grep("A",names(pop.mean)))+length(grep("B",names(pop.mean)))),
          rep(NA,length(grep("mean_v",names(pop.mean)))),   
          rep(0,length(grep("sd_v",names(pop.mean)))),
          rep(.1,length(grep("t0",names(pop.mean))))),
  upper=c(rep(NA,length(grep("A",names(pop.mean)))+length(grep("B",names(pop.mean)))),
          rep(NA,length(grep("mean_v",names(pop.mean)))),   
          rep(NA,length(grep("sd_v",names(pop.mean)))),
          rep(1,length(grep("t0",names(pop.mean)))))
)

mean.prior <- prior.p.dmc(
  dists = rep("tnorm",length(pop.mean)),
  p1=pop.mean,                           
  p2=c(rep(1,length(grep("A",names(pop.mean)))+length(grep("B",names(pop.mean)))),
       rep(2,length(grep("mean_v",names(pop.mean)))),   
       rep(1,length(grep("sd_v",names(pop.mean)))),
       rep(1,length(grep("t0",names(pop.mean))))),
  lower=c(rep(0,length(grep("A",names(pop.mean)))+length(grep("B",names(pop.mean)))),
          rep(NA,length(grep("mean_v",names(pop.mean)))),   
          rep(0,length(grep("sd_v",names(pop.mean)))),
          rep(.1,length(grep("t0",names(pop.mean))))),
  upper=c(rep(NA,length(grep("A",names(pop.mean)))+length(grep("B",names(pop.mean)))),
          rep(NA,length(grep("mean_v",names(pop.mean)))),   
          rep(NA,length(grep("sd_v",names(pop.mean)))),
          rep(1,length(grep("t0",names(pop.mean)))))
)


scale.prior <- prior.p.dmc(
  dists = rep("beta", length(pop.mean)),
  p1=(pop.mean>-1)*1, #hack to get p1 to be a vector of 1's with variable names
  p2=rep(1,length(pop.mean))
)

pp.prior <- list(mean.prior, scale.prior) 

#-----------------------------
# Generate starting values

starting_samples <- h.samples.dmc(nmc=100,pop.prior,data.model,thin=10,pp.prior=pp.prior)

#load LBA
load_model("lba","lba_B.R")

#Gets rid of bad chains
system.time({unstuck_samples  <- h.run.unstuck.dmc(starting_samples, p.migrate = .05, cores = 2)})

#Runs until gelman diag is below 1.1 for each chain
system.time({converged_samples <- h.run.converge.dmc(h.samples.dmc(nmc=100, samples=unstuck_samples), 
                                                     nmc=100,cores=2)})
