library(tidyverse)
library(magrittr)
source("TwoStage_Package_Code.R")

gof <- read_tsv("GOFuncIDs.cts_75.Mapping.padded.main.cumulative.summary_table.tsv")[,-2]
GOF <- gof %>% t() #%>% as.data.frame()
genes <- rownames(GOF)
genes <- genes[-1]
GOF <- GOF %>% as_data_frame()
GOF <- GOF %>% set_colnames(as.character(GOF[1,])) #%>% GOF[-1,]
GOF <- GOF[-1,]
GOF[["Sample_ID"]] <- genes
GOF <- GOF %>% dplyr::select(Sample_ID, everything())

factors <- read_csv("phyllo_factors.csv")#, col_types = cols("treatment" = col_factor(levels = NULL)))

Y.hf <-  (factors %>% filter(city == "CA")) %>% dplyr:::select(id = Sample_ID, treatment) 
Y.hf[,2] <- factor(as.matrix(Y.hf[,2]))
Y.hf[,2] <- factor(as.matrix(Y.hf[,2]))
Y.hf <- as.data.frame(Y.hf)
#Y.hf %>% map_if(is.factor, as.character) %>% as_data_frame -> Y.hf

Y.hf[,1] <- factor(as.matrix(Y.hf[,1]))
X.hf <- GOF %>% dplyr:::select(Sample_ID, which(names(.) %in% pull(Y.hf, id)))
X.hf <- X.hf %>% map_at(2:10, as.numeric) %>% as_data_frame() 
X.hf$Sample_ID <- as.factor(X.hf$Sample_ID)
# X.hf <- X.hf %>% mutate(total = rowSums(X.hf[,2:10]))
# X.hf <- X.hf %>% filter(!(total == 0)) %>% dplyr::select(1:10)

TwoStage_Package(as.data.frame(X.hf), Y.hf, fileoutput = "sigtest.csv", TMM.option = 1)

Y <- factor(Y.hf[,2])
ngroups <- nlevels(Y)
attr(Y.hf, "spec") <- NULL
attributes(Y.hf)
