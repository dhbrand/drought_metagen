library(tidyverse)

source("ennb.R")
source("ennb_no_glm_cv.R")
source("ennb_no_glm_cv_no_threshold_all_features.R")
source("ennb_no_threshold_all_features.R")

GF <- read_tsv("DELIVER/GeneFamilies.merged.tsv", col_types = ("cddddddddddddddddddddd"))
factors <- read_csv("phyllo_factors.csv")
factors <- factors[c(1:8, 13:21), ]
# The imputed classes needed to run cv.glmnet
factors <- rbind(factors, c("PHYLLO30", "HF", "drought"))
                 
                 , c("PHYLLO31", "HF", "drought"), c("PHYLLO32", "HF", "watered"), c("PHYLLO33", "HF", "drought"), c("PHYLLO34", "HF", "watered"), c("PHYLLO35", "CA", "drought"), c("PHYLLO36", "CA", "watered"), c("PHYLLO37", "CA", "drought"), c("PHYLLO38", "CA", "watered"), c("PHYLLO39", "DE", "drought"), c("PHYLLO40", "DE", "watered"), c("PHYLLO41", "DE", "drought"), c("PHYLLO42", "DE", "watered"))

setwd("test_output")

subjects <- str_split(names(GF[, -1]), "_", simplify = TRUE)[, 1]
names(GF) <- c("ID", subjects)

# Sites without maize removed
GF <- GF[-c(1:2), c(1:9, 14:22)]
# GF <- GF[-1,]

# impute 3rd drought location for HF
GF <- GF %>% mutate(PHYLLO30 = map2_dbl(.$PHYLLO28, .$PHYLLO29, ~ (.x + .y) / 2))

# impute 4th and 5th watered and drought for HF
GF <- GF %>%
  mutate(PHYLLO31 = pmap_dbl(list(.$PHYLLO28, .$PHYLLO29, .$PHYLLO30), ~ (..1 + ..2 + ..3) / 3), PHYLLO32 = pmap_dbl(list(.$PHYLLO09, .$PHYLLO26, .$PHYLLO27), ~ (..1 + ..2 + ..3) / 3)) %>%
  mutate(PHYLLO33 = pmap_dbl(list(.$PHYLLO28, .$PHYLLO29, .$PHYLLO30, .$PHYLLO31), ~ (..1 + ..2 + ..3 + ..4) / 4), PHYLLO34 = pmap_dbl(list(.$PHYLLO09, .$PHYLLO26, .$PHYLLO27, .$PHYLLO31), ~ (..1 + ..2 + ..3 + ..4) / 4))

# impute 4th and 5th watered and drought for CA
GF <- GF %>%
  mutate(PHYLLO35 = pmap_dbl(list(.$PHYLLO12, .$PHYLLO14, .$PHYLLO22), ~ (..1 + ..2 + ..3) / 3), PHYLLO36 = pmap_dbl(list(.$PHYLLO11, .$PHYLLO13, .$PHYLLO21), ~ (..1 + ..2 + ..3) / 3)) %>%
  mutate(PHYLLO37 = pmap_dbl(list(.$PHYLLO12, .$PHYLLO14, .$PHYLLO22, .$PHYLLO35), ~ (..1 + ..2 + ..3 + ..4) / 4), PHYLLO38 = pmap_dbl(list(.$PHYLLO11, .$PHYLLO13, .$PHYLLO21, .$PHYLLO36), ~ (..1 + ..2 + ..3 + ..4) / 4))

# impute 4th and 5th watered and drought for DE
GF <- GF %>%
  mutate(PHYLLO39 = pmap_dbl(list(.$PHYLLO16, .$PHYLLO24, .$PHYLLO25), ~ (..1 + ..2 + ..3) / 3), PHYLLO40 = pmap_dbl(list(.$PHYLLO10, .$PHYLLO15, .$PHYLLO23), ~ (..1 + ..2 + ..3) / 3)) %>%
  mutate(PHYLLO41 = pmap_dbl(list(.$PHYLLO16, .$PHYLLO24, .$PHYLLO25, .$PHYLLO39), ~ (..1 + ..2 + ..3 + ..4) / 4), PHYLLO42 = pmap_dbl(list(.$PHYLLO10, .$PHYLLO15, .$PHYLLO23, .$PHYLLO40), ~ (..1 + ..2 + ..3 + ..4) / 4))



Y.hf <- factors %>% filter(city == "HF") %>% dplyr::select(Sample_ID, treatment) %>% as.data.frame()
X.hf <- GF %>% dplyr::select(ID, pull(Y.hf, Sample_ID)) %>% as.data.frame()
X.hf <- X.hf[which(rowSums(X.hf[, 2:7]) != 0), ]

lambda_hf <- mean(read_rds("../lambdas/lambdaHF.rds"), na.rm = T)
alpha_hf <- mean(read_rds("../lambdas/alphaHF.rds"), na.rm = T)

ennb_no_glm_cv(X = X.hf,Y = Y.hf, fileoutput =  "sigtest_hf_tmm1_alpha05.csv", TMM_option = 1,threshold =  .05, lambda.opt = lambda_hf, alpha.opt = alpha_hf)
ennb_no_glm_cv(X = X.hf,Y = Y.hf, fileoutput =  "sigtest_hf_tmm2_alpha05.csv", TMM_option = 2,threshold =  .05, lambda.opt = lambda_hf, alpha.opt = alpha_hf)
ennb_no_glm_cv_no_threshold(X = X.hf, Y = Y.hf, fileoutput =  "fulltest_hf_tmm1.csv", TMM_option = 1, lambda.opt = lambda_hf, alpha.opt = alpha_hf)
ennb_no_glm_cv_no_threshold(X = X.hf, Y = Y.hf, fileoutput =  "fulltest_hf_tmm2.csv", TMM_option = 2, lambda.opt = lambda_hf, alpha.opt = alpha_hf)

Y.ca <- factors %>% filter(city == "CA") %>% dplyr::select(Sample_ID, treatment) %>% as.data.frame() 
X.ca <- GF %>% dplyr::select(ID, pull(Y.ca, Sample_ID)) %>% as.data.frame() %>%
  dplyr::select(1:2,4,6,3,5,7)
  dplyr::select(PHYLLO11, PHYLLO13, PHYLLO21, PHYLLO12, PHYLLO14, PHYLLO22)
X.ca <- X.ca[which(rowSums(X.ca[, 2:7]) != 0), ]

lambda_ca <- mean(read_rds("../lambdas/lambdaCA.rds"), na.rm = T)
alpha_ca <- mean(read_rds("../lambdas/alphaCA.rds"), na.rm = T)

ennb_no_glm_cv(X = X.ca,Y = Y.ca, fileoutput =  "sigtest_ca_tmm1_alpha05.csv", TMM_option = 1,threshold =  .05,lambda.opt =lambda_ca,alpha.opt = alpha_ca)
ennb_no_glm_cv(X = X.ca,Y = Y.ca, fileoutput =  "sigtest_ca_tmm2_alpha05.csv", TMM_option = 2,threshold =  .05,lambda.opt =lambda_ca,alpha.opt = alpha_ca)
ennb_no_glm_cv_no_threshold(X = X.ca,Y = Y.ca, fileoutput =  "fulltest_ca_tmm1.csv", TMM_option = 1, lambda.opt =lambda_ca,alpha.opt = alpha_ca)
ennb_no_glm_cv_no_threshold(X = X.ca,Y = Y.ca, fileoutput =  "fulltest_ca_tmm2.csv", TMM_option = 2, lambda.opt =lambda_ca,alpha.opt = alpha_ca)

Y.de <- factors %>% filter(city == "DE") %>% dplyr::select(Sample_ID, treatment) %>% as.data.frame()
X.de <- GF %>% dplyr::select(ID, pull(Y.de, Sample_ID)) %>% as.data.frame() %>%
  dplyr::select(1:3,5,4,6:7)
X.de <- X.de[which(rowSums(X.de[, 2:7]) != 0), ]

lambda_de <- mean(read_rds("../lambdas/lambdaDE.rds"), na.rm = T)
alpha_de <- mean(read_rds("../lambdas/alphaDE.rds"), na.rm = T)

ennb_no_glm_cv(X = X.de, Y = Y.de, fileoutput =  "sigtest_de_tmm1_alpha001.csv", TMM_option = 1,threshold =  .001,lambda.opt = lambda_de, alpha.opt = alpha_de)
ennb_no_glm_cv(X = X.de, Y = Y.de, fileoutput =  "sigtest_de_tmm2_alpha001.csv", TMM_option = 2,threshold =  .001,lambda.opt = lambda_de, alpha.opt = alpha_de)
ennb_no_glm_cv_no_threshold(X = X.de,Y = Y.de, fileoutput =  "fulltest_de_tmm1.csv", TMM_option = 1, lambda.opt =lambda_de, alpha.opt = alpha_de)
ennb_no_glm_cv_no_threshold(X = X.de,Y = Y.de, fileoutput =  "fulltest_de_tmm2.csv", TMM_option = 2, lambda.opt =lambda_de, alpha.opt = alpha_de)


ca <- colSums(X.ca == 0) / dim(X.ca)[1]
ca
hf <- colSums(X.hf == 0) / dim(X.hf)[1]
hf
de <- colSums(X.de == 0) / dim(X.de)[1]
de
mean(ca[-1])
mean(hf[-1])
mean(de[-1])

fit0 <- glmnet(t(X.de[, -1]), Y.de[, 2], family = "binomial", alpha = 0)
fit.5 <- glmnet(t(X.de[, -1]), Y.de[, 2], family = "binomial", alpha = 0.5)
fit1 <- glmnet(t(X.de[, -1]), Y.de[, 2], family = "binomial", alpha = 1)

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
X <- read.csv("../example_data/abundance.csv")
Y <- read.csv("../example_data/phenotype.csv")
TwoStage_Package(as.data.frame(GF), as.data.frame(Y), "sigtest_example.csv", 2)

# Checking sparsity in rows
g <- vector("integer")
for (i in 1:20) {
  l <- (length(which(x[i, 2:21] == 0)))
  g <- c(g, l)
}
which(g == 20)

# create test dataset for ENNB
x <- matrix(rnorm(100), nrow = 10)
feat <- NULL
for (i in 1:10) {
  feat <- c(feat, paste("f", i, sep = ""))
}
subj <- c("id")
for (i in 1:10) {
  subj <- c(subj, paste("s", i, sep = ""))
}
x <- data.frame(cbind(feat, x))
set_names(x, subj)
x <- data.frame(map_at(x, 2:11, as.numeric))
y <- data.frame(cbind(subj[-1], sample(1:2, 10, replace = T)))
TwoStage_Package(x, y, "sigtest_example.csv", 1)

ggplot(sig.df, aes(p.val)) +
  geom_histogram(bins = 10)
