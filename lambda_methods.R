library(glmnet)
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
