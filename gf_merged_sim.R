library(tidyverse)
library(stringr)
library(tictoc)

# setwd("/work/04734/dhbrand/stampede2/GitHub/drought_metagen/")
source("TwoStage_Package_Code.R")


GF <- read_tsv("MERGED_genefamilies.tsv") %>% as.data.frame
GF <- GF[-1,]
factors <- colnames(GF)[-1]
labels <- tibble(factors,label = c(rep("a",10),rep("b",10))) %>% as.data.frame

TwoStage_Package(GF,labels,"sigtest_gf_sim_tmm1.csv",1)
TwoStage_Package(GF,labels,"sigtest_gf_sim_tmm2.csv",2) 
# [1] "No significantly differentially abundant features!!!"
# Error in Twostage(X, Y, fileoutput) : object 'pvalue.mat' not found
X <- GF
Y <- labels
count <-  X.elast
category <-  Y

p.val = rep(NA,ncol(count))  #Declare initial pvalue

for (i in 1:ncol(count)) {
  dat = data.frame(testing=count[,1], groups=as.factor(category))
  if (sum(dat$testing) != 0) { #in case that all counts are zero
    fit = glm.nb(testing ~ groups, data=dat)
    p.val[i] = summary(fit)$coefficients[2,4] 
  }
summary(fit)
fit = glmnet(X.new, Y.new, family = "binomial", alpha = alpha.opt, lambda = cv$lambda.1se)
summary(fit)
enet.val <- colnames(X.elast)
GF.no.binary <- GF[!(GF[,1]) %in% enet.val,]
TwoStage_Package(GF.no.binary,labels,"sigtest_gf_sim_tmm1.csv",1)
