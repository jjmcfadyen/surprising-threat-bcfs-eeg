library(tidyverse)
library(ggplot2)

#############################################
### LOAD DATA ###############################
#############################################

min_rt = .5
max_rt = 10

# Read in long data
setwd("D:/Scratch/bCFS_EEG_Reanalysis/results")
CSV <- read.csv("trial_data_exp1.csv")

# Transform & clean data for LBA modelling
CSV$Condition = rep(1,dim(CSV)[1])
CSV$Condition[CSV$Emotion == 1 & CSV$Expectation == 1] = 1
CSV$Condition[CSV$Emotion == 1 & CSV$Expectation == 2] = 2
CSV$Condition[CSV$Emotion == 2 & CSV$Expectation == 1] = 3
CSV$Condition[CSV$Emotion == 2 & CSV$Expectation == 2] = 4

# check if any subjects don't have stimulus orientation information
na_subj <- unique(CSV$Subject[which(is.na(CSV$Orientation))])
CSV <- subset(CSV,Subject != na_subj)

# arrange data
D <- CSV %>% mutate(blocknumber=Block,
                       trialnumber=ExpTrial,
                       s=factor(as.numeric(as.factor(Subject))),
                       St=factor(Orientation,levels=1:2,labels=c("left","right")),
                       # Condition=factor(Condition,levels=1:4,labels=c("EN","UN","EF","UF")),
                       Emotion=factor(Emotion,levels=1:2,labels=c("Neutral","Fearful")),
                       Expectation=factor(Expectation,levels=1:2,labels=c("Expected","Unexpected")),
                       R=factor(Response,levels=1:2,labels=c("LEFT","RIGHT")),
                       RT=RT) %>%
  arrange(s,blocknumber,trialnumber) %>%
  select(s,St,Emotion,Expectation,R,RT,blocknumber,trialnumber)

# remove outlier RTs
D <- D %>% filter(RT > min_rt, RT < max_rt)

# Add accuracy info
D$Acc <- rep(FALSE,dim(D)[1])
D$Acc[D$St == "left" & D$R == "LEFT"] = TRUE
D$Acc[D$St == "right" & D$R == "RIGHT"] = TRUE

#############################################
### MODEL ###################################
#############################################

# Load in Dynamic Models of Choice
setwd("D:/DMC/DMC-MBN18/dmc")
source("dmc.R")
load_model("lba","lba_B.R")

factors=list(St=c("left","right"),Emotion=c("Neutral","Fearful"),Expectation=c("Expected","Unexpected"))
responses=c("LEFT","RIGHT")
match.map=list(M=list(left="LEFT",right="RIGHT"))

consts=c(sd_v.false=0,
         mean_v.Expected.Neutral.false=0,
         mean_v.Unexpected.Neutral.false=0,
         mean_v.Expected.Fearful.false=0,
         mean_v.Unexpected.Fearful.false=0)

p.map <- list(A="1",B=c("Emotion","Expectation"),mean_v=c("Expectation","Emotion","M"),sd_v="M",t0="1",st0="1")
model <- model.dmc(type="norm",
                   constants=consts,
                   p.map=p.map,
                   match.map=match.map,
                   factors=factors,
                   responses=responses)

data.model <- data.model.dmc(D, model)


# #--------------------------------
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
          rep(0,length(grep("mean_v",names(pop.mean)))),   
          rep(0,length(grep("sd_v",names(pop.mean)))),
          rep(.1,length(grep("t0",names(pop.mean))))),
  upper=c(rep(NA,length(grep("A",names(pop.mean)))+length(grep("B",names(pop.mean)))),
          rep(NA,length(grep("mean_v",names(pop.mean)))),   
          rep(NA,length(grep("sd_v",names(pop.mean)))),
          rep(max_rt/2,length(grep("t0",names(pop.mean)))))
)
par(mfcol=c(2,length(names(pop.prior))/2)); for (i in names(pop.prior)) plot.prior(i,pop.prior)

mean.prior <- prior.p.dmc(
  dists = rep("tnorm",length(pop.mean)),
  p1=pop.mean,                           
  p2=c(rep(1,length(grep("A",names(pop.mean)))+length(grep("B",names(pop.mean)))),
       rep(2,length(grep("mean_v",names(pop.mean)))),   
       rep(1,length(grep("sd_v",names(pop.mean)))),
       rep(1,length(grep("t0",names(pop.mean))))),
  lower=c(rep(0,length(grep("A",names(pop.mean)))+length(grep("B",names(pop.mean)))),
          rep(0,length(grep("mean_v",names(pop.mean)))),   
          rep(0,length(grep("sd_v",names(pop.mean)))),
          rep(.1,length(grep("t0",names(pop.mean))))),
  upper=c(rep(NA,length(grep("A",names(pop.mean)))+length(grep("B",names(pop.mean)))),
          rep(NA,length(grep("mean_v",names(pop.mean)))),   
          rep(NA,length(grep("sd_v",names(pop.mean)))),
          rep(max_rt/2,length(grep("t0",names(pop.mean)))))
)
par(mfcol=c(2,length(names(mean.prior))/2)); for (i in names(mean.prior)) plot.prior(i,mean.prior)

scale.prior <- prior.p.dmc(
  dists = rep("beta", length(pop.mean)),
  p1=(pop.mean>-1)*1, #hack to get p1 to be a vector of 1's with variable names
  p2=rep(1,length(pop.mean))
)
par(mfcol=c(2,length(names(scale.prior))/2)); for (i in names(scale.prior)) plot.prior(i,scale.prior)

pp.prior <- list(mean.prior, scale.prior) 

#-----------------------------
# Generate starting values
setwd("D:/Scratch/bCFS_EEG_Reanalysis/results")
starting_samples <- h.samples.dmc(nmc=100,pop.prior,data.model,thin=10,pp.prior=pp.prior)
save(starting_samples,file="starting_samples.RData")

#Gets rid of bad chains
system.time({unstuck_samples  <- h.run.unstuck.dmc(starting_samples, p.migrate = .05, cores = 4)})
save(unstuck_samples,file="unstuck_samples.RData")

#Runs until gelman diag is below 1.1 for each chain
system.time({converged_samples <- h.run.converge.dmc(h.samples.dmc(nmc=100, samples=unstuck_samples),
                                                     nmc=100,cores=4)})
save(converged_samples,file="converged_samples.RData")

gelman.diag.dmc(converged_samples)

plot.dmc(converged_samples[[1]],pll.chain=TRUE)
for (i in 1:1) {
  plot.dmc(converged_samples,density=FALSE,smooth=FALSE,subject=i,layout=c(2,5))
  Sys.sleep(0.1) 
}

#############################################
### CONVERGENCE #############################
#############################################

#get multivariate potential scale reduction factor (rhat)
gelman.diag.dmc(converged_samples)

# get effective sample size
effectiveSize.dmc(converged_samples)

#generate traceplots for each subject
Nsubj=length(converged_samples)
for(s in 1:Nsubj){
  plot.dmc(converged_samples[[s]])
  Sys.sleep(0.5)
}

#############################################
### POSTERIOR PREDICTIVES ###################
#############################################

#---------------------------------------------------------------------------------
#Here, we generate posterior predictives and plot in two ways. The first is via the method
#provided by dmc. 


#generate postior predictives dmc object
posterior_predictives_dmc <- h.post.predict.dmc(converged_samples)
save(posterior_predictives_dmc,file="posterior_predictives_dmc.RData")

#plot posterior cdfs
plot.pp.dmc(posterior_predictives_dmc)

#This summary function calculates the mean choice probabilities and quantiles collapsed across subjects.
#It removes a cell from consideration in the quantile estimate if it has fewer than 5 observations.
detach(package:plyr)
collapse<-function(data,design){
  #Remove quantiles from cells with fewer than five observations
  quant = data %>%
    group_by(St,Emotion,Expectation) %>%
    mutate(ncell=length(s)) %>%
    group_by(St,Emotion,Expectation,R) %>%
    mutate(nobs=length(s),prop=nobs/mean(ncell))
    
    
    summarise(nobs=length(s),
              prop=nobs/mean(ncell),
              q1=quantile(RT,0.1),
              q3=quantile(RT,0.3),
              q5=quantile(RT,0.5),
              q7=quantile(RT,0.7),
              q9=quantile(RT,0.9))
  quant[quant$nobs<5,c('q1','q3','q5','q7','q9')]<-NA

  quant2=left_join(design,quant, by = c("s", "St", "Emotion", "Expectation", "R"))
  quant2$prop[is.na(quant2$prop)]<-0

  collapsed=quant2 %>%
    group_by(St,Emotion,Expectation,R) %>%
    summarise(prop.m=mean(prop),
              prop.upper=prop.m+sd(prop)/sqrt(length(R)),
              prop.lower=prop.m-sd(prop)/sqrt(length(R)),
              q1.m=mean(q1,na.rm=T),
              q1.upper=q1.m+sd(q1,na.rm=T)/sqrt(sum(!is.na(R))),
              q1.lower=q1.m-sd(q1,na.rm=T)/sqrt(sum(!is.na(R))),
              q3.m=mean(q3,na.rm=T),
              q3.upper=q3.m+sd(q3,na.rm=T)/sqrt(sum(!is.na(R))),
              q3.lower=q3.m-sd(q3,na.rm=T)/sqrt(sum(!is.na(R))),
              q5.m=mean(q5,na.rm=T),
              q5.upper=q5.m+sd(q5,na.rm=T)/sqrt(sum(!is.na(R))),
              q5.lower=q5.m-sd(q5,na.rm=T)/sqrt(sum(!is.na(R))),
              q7.m=mean(q7,na.rm=T),
              q7.upper=q7.m+sd(q7,na.rm=T)/sqrt(sum(!is.na(R))),
              q7.lower=q7.m-sd(q7,na.rm=T)/sqrt(sum(!is.na(R))),
              q9.m=mean(q9,na.rm=T),
              q9.upper=q9.m+sd(q9,na.rm=T)/sqrt(sum(!is.na(R))),
              q9.lower=q9.m-sd(q9,na.rm=T)/sqrt(sum(!is.na(R))))
  return(collapsed)
}

Nsubj = length(converged_samples)
full_design = expand.grid(s=factor(1:Nsubj),
                          St=unique(D$St),
                          Emotion=unique(D$Emotion),
                          Expectation=unique(D$Expectation),
                          R=c('LEFT','RIGHT'))

observed_data=collapse(D,full_design)

names(observed_data)[c(5,8,11,14,17,20)]<-c('prop.mid','q1.mid','q3.mid','q5.mid','q7.mid','q9.mid')
observed_data$source='Data'

#Generate Posterior Predictives 
model = attributes(converged_samples[[1]]$data)$model

#Get sample parameter values for each subject
mcmc.list=list()
#use = matrix(NA,Nsubj,Nsamp)
for(s in 1:Nsubj){
  mcmc.list[[s]]=as.matrix(theta.as.mcmc.list(converged_samples[[s]]))
}
Nsamp <- 100
Nreps <- 60
use <- sample(1:dim(mcmc.list[[s]])[1],Nsamp)

#for each iteration generate data based on the sampled value
iter.list=list()
ctr1=0
# create progress bar
pb <- txtProgressBar(min = 0, max = Nsamp, style = 3)
for(i in 1:Nsamp){
  p.mat = matrix(NA,Nsubj,dim(converged_samples[[1]]$theta)[2])
  colnames(p.mat) = dimnames(converged_samples[[1]]$theta)[[2]]
  rownames(p.mat) = c("p.vector",rep("",dim(p.mat)[1]-1))
  
  #create nsubject x nparameter matrix of sampled values for that iteration 
  for(s in 1:Nsubj){
    p.mat[s,] = mcmc.list[[s]][use[i],]
  }
  
  #Simulate data based on parameter matrix for that iteration
  sim=h.simulate.dmc(model,ps=p.mat,ns=Nsubj,n=Nreps)
  
  #get RT quantile and choice proportion for each subject in each condition
  ctr1=ctr1+1
  iter.list[[ctr1]] <- collapse(sim,full_design) %>%
    mutate(source="Model",
           iter=ctr1)
  setTxtProgressBar(pb, i) #increment progress bar
}

#create and save alternate posterior predictives object 
posterior_predictives_alt = bind_rows(iter.list)
save(posterior_predictives_alt,file="posterior_predictives_alt.RData")

#calculate CIs by taking the quantiles of the different samples
predicted_data = posterior_predictives_alt %>% group_by(St,Emotion,Expectation,R) %>%
  summarise(prop.mid = quantile(prop.m,0.5),
            prop.lower = quantile(prop.m,0.025),
            prop.upper = quantile(prop.m,0.975),
            q1.mid = quantile(q1.m,0.5,na.rm=T),
            q1.lower = quantile(q1.m,0.025,na.rm=T),
            q1.upper = quantile(q1.m,0.975,na.rm=T),
            q3.mid = quantile(q3.m,0.5,na.rm=T),
            q3.lower = quantile(q3.m,0.025,na.rm=T),
            q3.upper = quantile(q3.m,0.975,na.rm=T),
            q5.mid = quantile(q5.m,0.5,na.rm=T),
            q5.lower = quantile(q5.m,0.025,na.rm=T),
            q5.upper = quantile(q5.m,0.975,na.rm=T),
            q7.mid = quantile(q7.m,0.5,na.rm=T),
            q7.lower = quantile(q7.m,0.025,na.rm=T),
            q7.upper = quantile(q7.m,0.975,na.rm=T),
            q9.mid = quantile(q9.m,0.5,na.rm=T),
            q9.lower = quantile(q9.m,0.025,na.rm=T),
            q9.upper = quantile(q9.m,0.975,na.rm=T),
            source='Model')

#get data in the right format to plot
plot_data = rbind(observed_data,predicted_data) %>%
  ungroup() %>%
  mutate(Emotion = factor(Emotion,levels=c("Neutral","Fearful"),labels=c("Neutral","Fearful")),
         Expectation = factor(Expectation,levels=c("Expected","Unexpected"),labels=c("Expected","Unexpected")),
         St = factor(St,levels=c("left","right"),labels=c("Left","Right")))
plot_data$Condition = rep(0,dim(plot_data)[1])
plot_data$Condition[plot_data$Emotion == "Neutral" & plot_data$Expectation == "Expected"] = "EN"
plot_data$Condition[plot_data$Emotion == "Neutral" & plot_data$Expectation == "Unexpected"] = "UN"
plot_data$Condition[plot_data$Emotion == "Fearful" & plot_data$Expectation == "Expected"] = "EF"
plot_data$Condition[plot_data$Emotion == "Fearful" & plot_data$Expectation == "Unexpected"] = "UF"

plot_data$Condition = as.factor(plot_data$Condition)
levels(plot_data$Condition) = c("EN","UN","EF","UF")

#plot RT data and save
plot_data_left <- plot_data[plot_data$St == "Left" & plot_data$R == "LEFT",]
plot_data_right <- plot_data[plot_data$St == "Right" & plot_data$R == "RIGHT",] 
qidx <- grep("q",names(plot_data))
plot_data_avg <- plot_data_left
plot_data_avg[,qidx] = (as.matrix(plot_data_left[,qidx]) + as.matrix(plot_data_right[,qidx])) / 2

par(mfrow=c(1,3))

plot_fit <- function(data,colour1,colour2,alpha){ # colour1 = data, colour2 = model
  ggplot(data=data,aes(x=Condition,group=source,colour=source)) +
    geom_line(aes(y=q1.mid),size=0.75,alpha=alpha) +
    geom_point(aes(y=q1.mid),size=2,alpha=alpha) +
    geom_errorbar(aes(ymin=q1.lower,ymax=q1.upper),width=0.1,size=0.75) +
    geom_line(aes(y=q3.mid),size=0.75,alpha=alpha) +
    geom_point(aes(y=q3.mid),size=2,alpha=alpha) +
    geom_errorbar(aes(ymin=q3.lower,ymax=q3.upper),width=0.1,size=0.75) +
    geom_line(aes(y=q5.mid),size=0.75,alpha=alpha) +
    geom_point(aes(y=q5.mid),size=2,alpha=alpha) +
    geom_errorbar(aes(ymin=q5.lower,ymax=q5.upper),width=0.1,size=0.75) +
    geom_line(aes(y=q7.mid),size=0.75,alpha=alpha) +
    geom_point(aes(y=q7.mid),size=2,alpha=alpha) +
    geom_errorbar(aes(ymin=q7.lower,ymax=q7.upper),width=0.1,size=0.75) +
    geom_line(aes(y=q9.mid),size=0.75,alpha=alpha) +
    geom_point(aes(y=q9.mid),size=2,alpha=alpha) +
    geom_errorbar(aes(ymin=q9.lower,ymax=q9.upper),width=0.1,size=0.75) +
    theme_minimal() +ylim(c(1.5,5)) + 
    # scale_color_viridis(discrete=TRUE) + scale_fill_viridis(discrete=TRUE) + 
    scale_color_manual(values = c(colour1,colour2)) + 
    labs(x="Condition",y="Response Time (seconds)",colour="Data Type")
}

plot_fit(plot_data_avg,"#900C3F","#FFC300",.8)


#############################################
### PARAMETERS ##############################
#############################################

#extract parameters from samples object into mcmc.list object
Nsubj = length(converged_samples)
Nsamp = dim(converged_samples[[1]]$theta)[1]*dim(converged_samples[[1]]$theta)[3]
mcmc.list=list()
use = matrix(NA,Nsubj,Nsamp)
for(s in 1:Nsubj){
  mcmc.list[[s]]=data.frame(as.matrix(theta.as.mcmc.list(converged_samples[[s]])))
  mcmc.list[[s]]$s = s
  mcmc.list[[s]]$iter = 1:Nsamp
}

parameters <- bind_rows(mcmc.list)

top_row <- parameters %>% 
  select(s,iter,B.Neutral.Expected:mean_v.Unexpected.Fearful.true) %>%
  gather(key,value,B.Neutral.Expected:mean_v.Unexpected.Fearful.true) %>%
  group_by(key,iter) %>%
  summarise(value.mean = mean(value))
top_row_B <- extract(top_row[grep("B.",top_row$key),],
          col=key,into=c('Parameter','Emotion','Expectation'),regex="(.+)\\.(.+)\\.(.+)")
top_row_mean_v <- extract(top_row[grep("mean_v.",top_row$key),],
          col=key,into=c('Parameter','Expectation','Emotion','R'),regex="(.+)\\.(.+)\\.(.+)\\.(.+)")

#get CIs for each comparison in top row
B_summary <- top_row_B %>%
  group_by(Parameter,Emotion,Expectation) %>%
  summarise(mean = mean(value.mean),
            lower = quantile(value.mean,0.025),
            upper = quantile(value.mean,0.975))

V_summary <- top_row_mean_v %>%
  group_by(Parameter,Emotion,Expectation) %>%
  summarise(mean = mean(value.mean),
            lower = quantile(value.mean,0.025),
            upper = quantile(value.mean,0.975))

# colours <- c("#2191F1","#650CA2","#FF0049","#FF9A00")
colours <- c("#FFC300","#DA125F")

plot_violin <- function(data,name,x_summary){
  ggplot() + 
  geom_violin(data=as.data.frame(data),
              aes(x=Emotion,y=value.mean,fill=Expectation),
              position=position_dodge(1)) + 
  geom_pointrange(data=x_summary,
               aes(x=Emotion,y=mean,ymin=lower,ymax=upper,group=Expectation),
               position=position_dodge(1),color="black") + 
  theme_minimal() + 
  scale_fill_manual(values = colours) + 
  labs(x="Emotion",y=name) + 
  scale_x_discrete(limits=c("Neutral","Fearful"))
}

plot_violin(top_row_B,"Threshold",B_summary)
plot_violin(top_row_mean_v,"Drift Rate",V_summary)


#Create plot of drift rates. Here we normalise on each subject's 
#mean drift across conditions and response alternatives. This isolates
#the within subject variance in drift and therefore the pure effect of the
#manipulations
plot_data <- parameters %>% 
  select(s,iter,B.Neutral.Expected:mean_v.Unexpected.Fearful.true) %>%
  gather(key,value,B.Neutral.Expected:mean_v.Unexpected.Fearful.true) %>%
  group_by(iter,s) %>%
  mutate(value.mc = value-mean(value))

plot_B <- extract(plot_data[grep("B.",plot_data$key),],
                  col=key,into=c('Parameter','Emotion','Expectation'),regex="(.+)\\.(.+)\\.(.+)")
plot_V <- extract(plot_data[grep("mean_v.",plot_data$key),],
                  col=key,into=c('Parameter','Expectation','Emotion','R'),regex="(.+)\\.(.+)\\.(.+)\\.(.+)")
plot_V <- select(plot_V,s:Emotion,value:value.mc)

plot_data <- rbind(plot_B,plot_V)

sum_data <- plot_data %>%
  group_by(Parameter,Emotion,Expectation) %>%
  summarise(lower = quantile(value.mc,0.025),
            median = quantile(value.mc,0.5),
            upper = quantile(value.mc,0.975))

B_sum <- sum_data[sum_data$Parameter == "B",]
V_sum <- sum_data[sum_data$Parameter == "mean_v",]

plot_data_bysubj <- plot_data %>% 
  group_by(s,Parameter,Emotion,Expectation) %>%
  mutate(subj_value.mc = mean(value.mc)) %>%
  select(s,Parameter,Emotion,Expectation,subj_value.mc)

plot_mc <- function(sumdata,name){
ggplot() + 
  geom_violin(data=plot_data,aes(x=Emotion,fill=Expectation,y=value.mc),
              position=position_dodge(1)) +
  geom_pointrange(data=sumdata,aes(x=Emotion,y=median,group=Expectation,ymin=lower,ymax=upper),
                  position=position_dodge(1),color="black",size=1) + 
  theme_classic() + 
  scale_x_discrete(limits=c("Neutral","Fearful")) + 
  labs(x="Emotion",y=name)
}
  
plot_mc(B_sum,"Threshold")
plot_mc(V_sum,"Drift Rate")

subject_plot_data <- plot_data %>% 
  group_by(s,Parameter,Emotion,Expectation) %>% 
  summarise(subject_value.mc = mean(value.mc))

plot_subject <- function(param,name){
  ggplot(data=subject_plot_data[subject_plot_data$Parameter == param,],
       aes(x=Emotion,fill=Expectation,y=subject_value.mc)) + 
  geom_dotplot(binaxis='y',position=position_dodge(1),dotsize=.5,stackdir='center') +
  stat_summary(fun.data=mean_sdl,position=position_dodge(1)) +
  theme_classic() + 
  scale_x_discrete(limits=c("Neutral","Fearful")) + 
  labs(x="Emotion",y=name)
}

plot_subject("B","Threshold")
plot_subject("mean_v","Drift Rate")

#############################################
### BAYES FACTORS ###########################
#############################################

#number of samples from the prior
nsamples = 1E5 

#for each subject, sample from prior for threshold and drift parameters
prior_list = list()
for(i in 1:Nsubj){
  priors_tmp = data.frame(
    B.Neutral.Expected = rtnorm(n = nsamples,
                         mean=rtnorm(n=nsamples,mean=1,sd=1,lower=0,upper=Inf),
                         sd=runif(n=nsamples,min=0,max=1),
                         lower = 0,
                         upper = Inf),
    B.Neutral.Unexpected = rtnorm(n = nsamples,
                         mean=rtnorm(n=nsamples,mean=1,sd=1,lower=0,upper=Inf),
                         sd=runif(n=nsamples,min=0,max=1),
                         lower = 0,
                         upper = Inf),
    B.Fearful.Expected = rtnorm(n = nsamples,
                         mean=rtnorm(n=nsamples,mean=1,sd=1,lower=0,upper=Inf),
                         sd=runif(n=nsamples,min=0,max=1),
                         lower = 0,
                         upper = Inf),
    B.Fearful.Unexpected = rtnorm(n = nsamples,
                          mean=rtnorm(n=nsamples,mean=1,sd=1,lower=0,upper=Inf),
                          sd=runif(n=nsamples,min=0,max=1),
                          lower = 0,
                          upper = Inf),
    mean_v.Neutral.Expected = rtnorm(n = nsamples,
                                     mean=rtnorm(n=nsamples,mean=1,sd=2,lower=-Inf,upper=Inf),
                                     sd=runif(n=nsamples,min=0,max=1),
                                     lower = -Inf,
                                     upper = Inf),
    mean_v.Neutral.Unexpected = rtnorm(n = nsamples,
                                       mean=rtnorm(n=nsamples,mean=1,sd=2,lower=-Inf,upper=Inf),
                                       sd=runif(n=nsamples,min=0,max=1),
                                       lower = -Inf,
                                       upper = Inf),
    mean_v.Fearful.Expected = rtnorm(n = nsamples,
                                     mean=rtnorm(n=nsamples,mean=1,sd=2,lower=-Inf,upper=Inf),
                                     sd=runif(n=nsamples,min=0,max=1),
                                     lower = -Inf,
                                     upper = Inf),
    mean_v.Fearful.Unexpected = rtnorm(n = nsamples,
                                       mean=rtnorm(n=nsamples,mean=1,sd=2,lower=-Inf,upper=Inf),
                                       sd=runif(n=nsamples,min=0,max=1),
                                       lower = -Inf,
                                       upper = Inf)
    
  ) 
  priors_tmp$s = i
  priors_tmp$iter = 1:nsamples
  prior_list[[i]] = priors_tmp
  print(i)
}

priors = bind_rows(prior_list) %>%
  gather(key,value,B.Neutral.Expected:mean_v.Fearful.Unexpected) %>%
  extract(col=key,into=c('Parameter','Emotion','Expectation'),regex="(.+)\\.(.+)\\.(.+)")

levels(priors$Emotion) = c("Neutral","Fearful")
levels(priors$Expectation) = c("Expected","Unexpected")
levels(priors$Parameter) = c("B","mean_v")

#Get prior density for mean quality, quantity, and threshold at 0
prior_density_0 = c(B=NA,mean_v=NA)

for(parm in c('B','mean_v')){
  parm_priors = filter(priors,Parameter==parm) %>%
    group_by(iter,Emotion,Expectation) %>%
    summarise(value = mean(value))
  d_prior = approxfun(density(parm_priors$value),rule=2)
  prior_density_0[names(prior_density_0)==parm] = d_prior(0)
}

difference_parameters <- parameters %>%
  mutate(B.Neutral_vs_Fearful = (B.Neutral.Expected+B.Neutral.Unexpected)/2 - 
                                (B.Fearful.Expected+B.Fearful.Unexpected)/2,
         B.Expected_vs_Unexpected = (B.Neutral.Expected+B.Fearful.Expected)/2 - 
                                    (B.Neutral.Unexpected+B.Fearful.Unexpected)/2,
         B.EN_vs_UN = B.Neutral.Expected - B.Neutral.Unexpected,
         B.EF_vs_UF = B.Fearful.Expected - B.Fearful.Unexpected,
         B.EN_vs_EF = B.Neutral.Expected - B.Fearful.Expected,
         B.UN_vs_UF = B.Neutral.Unexpected - B.Fearful.Unexpected,
         mean_v.Neutral_vs_Fearful = (mean_v.Expected.Neutral.true+mean_v.Unexpected.Neutral.true)/2 -
                                     (mean_v.Expected.Fearful.true+mean_v.Unexpected.Fearful.true)/2,
         mean_v.Expected_vs_Unexpected = (mean_v.Expected.Neutral.true+mean_v.Expected.Fearful.true) -
                                        (mean_v.Unexpected.Neutral.true+mean_v.Unexpected.Fearful.true),
         mean_v.EN_vs_UN = mean_v.Expected.Neutral.true - mean_v.Unexpected.Neutral.true,
         mean_v.EF_vs_UF = mean_v.Expected.Fearful.true - mean_v.Unexpected.Fearful.true,
         mean_v.EN_vs_EF = mean_v.Expected.Neutral.true - mean_v.Expected.Fearful.true,
         mean_v.UN_vs_UF = mean_v.Unexpected.Neutral.true - mean_v.Unexpected.Fearful.true) %>%
  select(s,iter,B.Neutral_vs_Fearful:mean_v.UN_vs_UF) %>%
  gather(key,value,B.Neutral_vs_Fearful:mean_v.UN_vs_UF) %>%
  extract(col=key,into=c('Parameter','Inc'),regex="(.+)\\.(.+)") %>%
  mutate(Inc = factor(Inc,levels=c('Neutral_vs_Fearful','Expected_vs_Unexpected','EN_vs_UN','EF_vs_UF','EN_vs_EF','UN_vs_UF'),
                      labels=c("Neutral vs Fearful","Expected vs Unexpected","EN vs UN","EF vs UF","EN vs EF","UN vs UF")),
         Parm = factor(Parameter,levels=c('B','mean_v'),labels=c('Threshold','Drift Rate')))

bfs = difference_parameters %>% 
  group_by(Parameter,Inc) %>%
  summarise(bf=NA)

for(parm in c('B','mean_v')){
  for(inc in c("Neutral vs Fearful","Expected vs Unexpected","EN vs UN","EF vs UF","EN vs EF","UN vs UF")){
      
      #Get mean bfs
      tmp = filter(difference_parameters,Parameter==parm,Inc==inc) %>%
        group_by(iter) %>%
        summarise(value = mean(value))
      
      d_post = approxfun(density(tmp$value),rule=2)
      bf = prior_density_0[names(prior_density_0)==parm] / d_post(0)
      bfs[bfs$Parameter==parm & bfs$Inc == inc,'bf'] = bf
  }
}

bfs


