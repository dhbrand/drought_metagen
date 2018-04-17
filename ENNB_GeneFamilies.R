library(tidyverse)
library(stringr)
library(tictoc)
require(glmnet)

source("TwoStage_Package_Code.R")
select <- dplyr::select
#x <- spec_tsv("DELIVER/GeneFamilies.merged.tsv");cols_condense(x)
GF <- read_tsv("DELIVER/GeneFamilies.merged.tsv", col_types = ("cddddddddddddddddddddd"))
factors <- read_csv("phyllo_factors.csv")
factors <- rbind(factors, c("PHYLLO30", "HF", "drought"),c("PHYLLO31", "HF", "drought"), c("PHYLLO32", "HF", "watered"), c("PHYLLO33", "HF", "drought"), c("PHYLLO34", "HF", "watered"))
factors <- factors[c(1:8,13:22),]
# str_length(names(GF[,-1]))
# str_length(names(factors[,1]))
# str_view(names(GF[,-1]), "_")

subjects <- str_split(names(GF[,-1]), "_", simplify = TRUE)[,1]
names(GF) <- c("ID", subjects)
GF <- GF[-1,c(1:9,14:22)]

# impute 3rd drought location for HF
GF <- GF %>% mutate(PHYLLO30 = map2_dbl(.$PHYLLO28,.$PHYLLO29, ~ (.x + .y)/2))

# impute 4th and 5th watered and drought for HF
GFm <- GF %>% mutate(PHYLLO31 = pmap_dbl(list(.$PHYLLO28, .$PHYLLO29, .$PHYLLO30), ~ (..1 + ..2 + ..3)/3), PHYLLO32 = pmap_dbl(list(.$PHYLLO09, .$PHYLLO26, .$PHYLLO27), ~ (..1 + ..2 + ..3)/3)) %>% 
  mutate(PHYLLO33 = pmap_dbl(list(.$PHYLLO28, .$PHYLLO29, .$PHYLLO30, .$PHYLLO31), ~ (..1 + ..2 + ..3 + ..4)/4), PHYLLO34 = pmap_dbl(list(.$PHYLLO09, .$PHYLLO26, .$PHYLLO27, .$PHYLLO31), ~ (..1 + ..2 + ..3 + ..4)/4))


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


Y.hf <- factors %>% filter(city == "HF") %>% select(Sample_ID, treatment) %>% as.data.frame
X.hf <- GFm %>% select(ID,pull(Y, Sample_ID)) %>% as.data.frame

TwoStage_Package(X.hf,Y.hf,"sigtest_example.csv",1)
# 
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


x <- read.csv("abundance.csv.csv")
z <- read.csv("phenotype.csv")
TwoStage_Package(as.data.frame(GF),as.data.frame(Y),"sigtest_example.csv",2)


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

# Find a way to determine lambda without CV
x <- matrix(rnorm(100), nrow = 9)
y <- sample(1:2,9, replace = T)
fit <- glmnet(x,y,family = "binomial")
plot(fit)
# confirm 6 observations 2 small for cv.glmnet
cv = cv.glmnet(x, y, family = "binomial", nfolds = 3)
tLL <- fit$nulldev - deviance(fit)
k <- fit$df
n <- fit$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
which(AICc == min(AICc))
BIC<-log(n)*k - tLL
which(BIC == min(BIC))
lambda <- c(AICc[which(AICc == min(AICc))], BIC[which(BIC == min(BIC))])

require(HDeconometrics)
fit1 <- ic.glmnet(x, y)
plot(fit1)
lambda.ic <- fit1$lambda
