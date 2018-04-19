library(tidyverse)
library(stringr)
library(tictoc)
library(glmnet)
library(edgeR)

setwd("/work/04734/dhbrand/stampede2/GitHub/drought_metagen/"


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
#GF <- GF[-1,]

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

Y.ca <- factors %>% filter(city == "CA") %>% select(Sample_ID, treatment) %>% as.data.frame
X.ca <- GF %>% select(ID,pull(Y.ca, Sample_ID)) %>% as.data.frame
X <- X.ca;Y <- Y.ca


TMMNorm = function(X, Y, TMM.option){
  factors = calcNormFactors(X,method="TMM") #Calculate normaization factors
  if(TMM.option==1){
    eff.lib.size = colSums(X)*factors;
    ref.lib.size = mean(eff.lib.size); #Use the mean of the effective library sizes as a reference library size
    X.output = sweep(X,MARGIN=2,eff.lib.size,"/")*ref.lib.size; #Normalized read counts
  } else if(TMM.option==2){ #Regenerate counts with a common dispersion.
    group = as.factor(Y);
    X.tmp = DGEList(counts=X, group=group);
    X.tmp = calcNormFactors(X.tmp);
    X.tmp = estimateCommonDisp(X.tmp);
    X.output = as.matrix(X.tmp$pseudo.counts);
  } else {
    stop("TMM.option must be either 1 or 2!");
  }
  return(X.output);
}


rownames(X) = X[,1]
X = X[,-1]
X = as.matrix(X)
Y = factor(Y[,2])
X.tmm = TMMNorm(X,Y,1)  # TMM normalization (edgeR)
X = t(X.tmm)
ngroups = nlevels(Y)


# find lambda using CLT
alpha = seq(0.01,0.1,by=0.01)
cv.mat = matrix(NA, nrow = length(alpha), ncol = 3)
colnames(cv.mat) = c("Alpha", "Lambda", "CVM")
cv.mat[,1] = alpha
npop1 = table(Y)[1]
npop2 = table(Y)[2]
N = npop1 + npop2
if (N <= 15) nfold = 3
if (N > 15  & N < 50) nfold = 5
if (N >= 50) nfold = 10

alpha.step <- NULL
lambda.step <- NULL

for (i in 1:100) {
    for (j in 1:length(alpha)) {
      try({
      cv = cv.glmnet(X, Y, family = "binomial", alpha = alpha[j], nfolds = nfold)
      ind = match(cv$lambda.min, cv$lambda)
      cv.mat[j,2] = cv$lambda[ind]
      cv.mat[j,3] = cv$cvm[ind]
      })
    }
  alpha.opt = cv.mat[cv.mat[,"CVM"] == min(cv.mat[,"CVM"]),1]
  lambda.opt <- cv.mat[cv.mat[,"CVM"] == min(cv.mat[,"CVM"]),2]
  alpha.step <- c(alpha.step, alpha.opt)
  lambda.step <- c(lambda.step, lambda.opt)
}

write_rds(alpha.step, "alphaCA.rds")
write_rds(lambda.step, "lambdaCA.rds")
