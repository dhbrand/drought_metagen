library(tidyverse)
library(stringr)
library(tictoc)

setwd("/work/04734/dhbrand/stampede2/GitHub/drought_metagen/")
source("TwoStage_Package_Code.R")

# MASS select masks dplyr select
select <- dplyr::select

GF <- read_tsv("DELIVER/GeneFamilies.merged.tsv", col_types = ("cddddddddddddddddddddd"))
factors <- read_csv("phyllo_factors.csv")
factors <- factors[c(1:8,13:21),]
# The imputed classes needed to run cv.glmnet
factors <- rbind(factors, c("PHYLLO30", "HF", "drought"),c("PHYLLO31", "HF", "drought"), c("PHYLLO32", "HF", "watered"), c("PHYLLO33", "HF", "drought"), c("PHYLLO34", "HF", "watered"),c("PHYLLO35", "CA", "drought"), c("PHYLLO36", "CA", "watered"), c("PHYLLO37", "CA", "drought"), c("PHYLLO38", "CA", "watered"),c("PHYLLO39", "DE", "drought"), c("PHYLLO40", "DE", "watered"), c("PHYLLO41", "DE", "drought"), c("PHYLLO42", "DE", "watered"))



subjects <- str_split(names(GF[,-1]), "_", simplify = TRUE)[,1]
names(GF) <- c("ID", subjects)

# Sites without maize removed
GF <- GF[-1,c(1:9,14:22)]

# impute 3rd drought location for HF
GF <- GF %>% mutate(PHYLLO30 = map2_dbl(.$PHYLLO28,.$PHYLLO29, ~ (.x + .y)/2))

# impute 4th and 5th watered and drought for HF
GF <- GF %>% mutate(PHYLLO31 = pmap_dbl(list(.$PHYLLO28, .$PHYLLO29, .$PHYLLO30), ~ (..1 + ..2 + ..3)/3), PHYLLO32 = pmap_dbl(list(.$PHYLLO09, .$PHYLLO26, .$PHYLLO27), ~ (..1 + ..2 + ..3)/3)) %>% 
  mutate(PHYLLO33 = pmap_dbl(list(.$PHYLLO28, .$PHYLLO29, .$PHYLLO30, .$PHYLLO31), ~ (..1 + ..2 + ..3 + ..4)/4), PHYLLO34 = pmap_dbl(list(.$PHYLLO09, .$PHYLLO26, .$PHYLLO27, .$PHYLLO31), ~ (..1 + ..2 + ..3 + ..4)/4))

# impute 4th and 5th watered and drought for CA
GF <- GF %>% mutate(PHYLLO35 = pmap_dbl(list(.$PHYLLO12, .$PHYLLO14, .$PHYLLO22), ~ (..1 + ..2 + ..3)/3), PHYLLO36 = pmap_dbl(list(.$PHYLLO11, .$PHYLLO13, .$PHYLLO21), ~ (..1 + ..2 + ..3)/3)) %>% 
  mutate(PHYLLO37 = pmap_dbl(list(.$PHYLLO12, .$PHYLLO14, .$PHYLLO22, .$PHYLLO35), ~ (..1 + ..2 + ..3 + ..4)/4), PHYLLO38 = pmap_dbl(list(.$PHYLLO11, .$PHYLLO13, .$PHYLLO21, .$PHYLLO36), ~ (..1 + ..2 + ..3 + ..4)/4))

# impute 4th and 5th watered and drought for DE
GF <- GF %>% mutate(PHYLLO39 = pmap_dbl(list(.$PHYLLO16, .$PHYLLO24, .$PHYLLO25), ~ (..1 + ..2 + ..3)/3), PHYLLO40 = pmap_dbl(list(.$PHYLLO10, .$PHYLLO15, .$PHYLLO23), ~ (..1 + ..2 + ..3)/3)) %>% 
  mutate(PHYLLO41 = pmap_dbl(list(.$PHYLLO16, .$PHYLLO24, .$PHYLLO25, .$PHYLLO39), ~ (..1 + ..2 + ..3 + ..4)/4), PHYLLO42 = pmap_dbl(list(.$PHYLLO10, .$PHYLLO15, .$PHYLLO23, .$PHYLLO40), ~ (..1 + ..2 + ..3 + ..4)/4))

setwd("output/")

city <- c("HF", "CA", "DE")
tic()
for (i in city) {
  Y <- factors %>% filter(city == i) %>% select(Sample_ID, treatment) %>% as.data.frame
  X <- GF %>% select(ID,pull(Y, Sample_ID)) %>% as.data.frame
  TwoStage_Package(X,Y,paste("sigtest_", i, "_tmm1.csv", sep = ""), 1)
  TwoStage_Package(X,Y,paste("sigtest_", i, "_tmm2.csv", sep = ""), 2)
}
print(toc())
# Y.hf <- factors %>% filter(city == "HF") %>% select(Sample_ID, treatment) %>% as.data.frame
# X.hf <- GF %>% select(ID,pull(Y.hf, Sample_ID)) %>% as.data.frame
# 
# TwoStage_Package(X.hf,Y.hf,"sigtest_hf_tmm1.csv",1)
# TwoStage_Package(X.hf,Y.hf,"sigtest_hf_tmm2.csv",2)
# 
# Y.ca <- factors %>% filter(city == "CA") %>% select(Sample_ID, treatment) %>% as.data.frame
# X.ca <- GF %>% select(ID,pull(Y.ca, Sample_ID)) %>% as.data.frame
# 
# TwoStage_Package(X.ca,Y.ca,"sigtest_ca_tmm1.csv",1)
# TwoStage_Package(X.ca,Y.ca,"sigtest_ca_tmm2.csv",2)
# 
# Y.de <- factors %>% filter(city == "DE") %>% select(Sample_ID, treatment) %>% as.data.frame
# X.de <- GF %>% select(ID,pull(Y.de, Sample_ID)) %>% as.data.frame
# 
# TwoStage_Package(X.de,Y.de,"sigtest_de_tmm1.csv",1)
# TwoStage_Package(X.de,Y.de,"sigtest_de_tmm2.csv",2)



# Find rows with 1 or 2 observations; do not need as ENNB transposes X
# GF1 <- GF %>% as_tibble %>% mutate_all(funs(.==0)) %>% reduce(`+`) %>% cbind(GF, Count = .)
# GF1 <- GF1 %>% filter(!(Count %in% 18:19)) %>% select(1:19)
# 
# GF2 <- GF %>% 
#   as_tibble %>% 
#   mutate_all(funs(.==0)) %>% 
#   reduce(`+`) %>% 
#   {!(. %in% c(18, 19))} %>%
#   magrittr::extract(GF, .)



 
# Repexp for stackoverflow help
# set.seed(123)
# library(Matrix)
# M1 <- as.matrix(rsparsematrix(100, 20, .1, rand.x = runif))
# x <- vector("integer")
# for(i in 1:dim(M1)[1]){
#   l <- (length(which(M1[i,] == 0)))
#   x <- c(x,l)
# }
# ind <- which(x == 19 | x == 20)
# M1 <- M1[-ind,]
# 
# M1 %>% mutate(zero_count = length(map(which(. == 0))))
# 
# M2 <- M1 %>% as_tibble %>% mutate_all(funs(.==0)) %>% reduce(`+`) %>% cbind(M1, Count = .)

# Example data from developers website
x <- read.csv("example_data/abundance.csv")
z <- read.csv("phenotype.csv")
TwoStage_Package(as.data.frame(GF),as.data.frame(Y),"sigtest_example.csv",2)

# Checking sparsity in rows
g <- vector("integer")
for(i in 1:20){
  l <- (length(which(x[i,2:21] == 0)))
  g <- c(g,l)
}
which(g==20)

# create test dataset for ENNB
x <- matrix(rnorm(100), nrow = 10)
feat <- NULL
for (i in 1:10) { feat <- c(feat,paste("f",i,sep = "")) }
subj <- c("id")
for (i in 1:10) { subj <- c(subj,paste("s",i,sep = "")) }
x <- data.frame(cbind(feat, x))
set_names(x,subj)
x <- data.frame(map_at(x,2:11,as.numeric))
y <- data.frame(cbind(subj[-1],sample(1:2,10, replace = T)))
TwoStage_Package(x,y,"sigtest_example.csv",1)

