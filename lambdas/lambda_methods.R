library(glmnet)
library(edgeR)
library(HDeconometrics)
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




# Example data from developers website
X <- read.csv("example_data/abundance.csv")
Y <- read.csv("example_data/phenotype.csv")


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
      cv = try(cv.glmnet(X, Y, family = "binomial", alpha = alpha[j], nfolds = nfold))
      ind = match(cv$lambda.min, cv$lambda)
      cv.mat[j,2] = cv$lambda[ind]
      cv.mat[j,3] = cv$cvm[ind]
    }
  alpha.opt = cv.mat[cv.mat[,"CVM"] == min(cv.mat[,"CVM"]),1]
  lambda.opt <- cv.mat[cv.mat[,"CVM"] == min(cv.mat[,"CVM"]),2]
  alpha.step <- c(alpha.step, alpha.opt)
  lambda.step <- c(lambda.step, lambda.opt)
}
summary(alpha.step); summary(lambda.step)
