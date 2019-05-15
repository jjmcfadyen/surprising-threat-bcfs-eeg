require(dplyr)
require(rtdists)
require(ggplot2)

############################################################
### Load in data 
############################################################
setwd("D:/Scratch/bCFS_EEG_Reanalysis/results")

CSV <- read.csv("trial_data_exp2.csv", header = TRUE)

# change category numbers to words
CSV$Gender[which(CSV$Gender == 1)] = "male"
CSV$Gender[which(CSV$Gender == 2)] = "female"
CSV$Order[which(CSV$Order == 1)] = "neutralfirst"
CSV$Order[which(CSV$Order == 2)] = "fearfulfirst"
CSV$Emotion[which(CSV$Emotion == 1)] = "neutral"
CSV$Emotion[which(CSV$Emotion == 2)] = "fearful"
CSV$Expectation[which(CSV$Expectation == 1)] = "expected"
CSV$Expectation[which(CSV$Expectation == 2)] = "unexpected"
CSV$Response[which(CSV$Response == 1)] = "left"
CSV$Response[which(CSV$Response == 2)] = "right"
CSV$Orientation[which(CSV$Orientation == 1)] = "left"
CSV$Orientation[which(CSV$Orientation == 2)] = "right"

CSV$condition <- rep(0,dim(CSV)[1])
CSV$condition[which(CSV$Emotion=="neutral" & CSV$Expectation=="expected")] = "EN"
CSV$condition[which(CSV$Emotion=="neutral" & CSV$Expectation=="unexpected")] = "UN"
CSV$condition[which(CSV$Emotion=="fearful" & CSV$Expectation=="expected")] = "EF"
CSV$condition[which(CSV$Emotion=="fearful" & CSV$Expectation=="unexpected")] = "UF"


############################################################
### LBA
############################################################

### Load in DMC
setwd("D:/Scratch/bCFS_EEG_Reanalysis/scripts/surprising-threat-bcfs-eeg/Analysis_Scripts/DMC-MBN18")
source ("dmc/dmc.R")

load_model (dir_name="LBA",model_name="lba_B.R")

factors <- list(St=c("left","right"),Emotion="neutral","fearful",Expectation="expected","unexpected")
responses <- c("left","right")
match.map = list(M=list(left="left",right="right"))
consts <-c(sd_v=1,st0=0)
p.map=list(A="1",B=c("Inc","R"),mean_v=c("Diff","Inc","M"),sd_v="1",t0="1",st0="1")




