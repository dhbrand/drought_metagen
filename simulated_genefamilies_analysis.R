library(tidyverse)
library(magrittr)
library(rlang)

source("TwoStage_Package_Code.R")

dat_calls <- read_tsv("DNApolymAIII/Calls_genefamilies-cpm.tsv")[-1,] %>% dplyr::select(1, Annotation = "EstID", label = "EstEC", 4:34) %>% replace_na(list(label = "unknown")) %>% as.data.frame 

factors <- cbind(names(dat_calls)[4:33], c(rep("A", 10), rep("B", 10), rep("C", 10))) %>% as.data.frame

TwoStage_Package(dat_id, factors, "gf_sim_id_tmm1.csv", 1)

# Sig EC labels
# 2.7.7.7, 2.3.1.16, 2.4.2.14

# Group A vs B
gf_a_b <- dat_calls %>% filter(Call %in% c("A", "B")) %>% dplyr::select(1:23)
factors_a_b <- factors %>% filter(V2 %in% c("A", "B"))
TwoStage_Package(gf_a_b[,-c(1,3)], factors_a_b, "gf_a_b_tmm1.csv", 1, .001)

# sig_a_b <- read_csv("gf_a_b_tmm1.csv")
# a_b <- left_join(sig_a_b, dat_calls, by = "Annotation") %>% dplyr::select(1:9)
# sig_ec <- c("2.7.7.7", "2.3.1.16", "2.4.2.14")
# sig_ec_found <- intersect(sig_ec, unique(a_b$label))
# ec_else <- setdiff(unique(a_b$label), sig_ec)
# tp <- sum(table(a_b$label)[sig_ec_found])
# fp <- sum(table(a_b$label)[ec_else])
# fn <- sum(table(gf_a_b$label)[sig_ec_found]) - tp
# tn <- nrow(gf_a_b) - sum(tp, fp, fn)
# confusion_matrix <- cbind("positive" = c(tp, fp), "negative" = c(tn, fn)) 
# rownames(confusion_matrix) <- (c("true", "false"))
# fdr <-  tn / (tn + fn)
# confusion_matrix; paste("FDR = ",fdr)


# Group A vs C
gf_a_c <- dat_calls %>% filter(Call %in% c("A", "C")) %>% dplyr::select(1:12, 23:33)
factors_a_c <- factors %>% filter(V2 %in% c("A", "C"))
TwoStage_Package(gf_a_c[,-c(1,3)], factors_a_c, "gf_a_c_tmm1.csv", 1, .001)

# sig_a_c <- read_csv("gf_a_c_tmm1.csv")
# a_c <- left_join(sig_a_c, dat_calls, by = "Annotation") %>% dplyr::select(1:9)
# sig_ec <- c("2.7.7.7", "2.3.1.16", "2.4.2.14")
# sig_ec_found <- intersect(sig_ec, unique(a_c$label))
# ec_else <- setdiff(unique(a_c$label), sig_ec)
# tp <- sum(table(a_c$label)[sig_ec_found])
# fp <- sum(table(a_c$label)[ec_else])
# fn <- sum(table(gf_a_c$label)[sig_ec_found]) - tp
# tn <- nrow(gf_a_c) - sum(tp, fp, fn)
# confusion_matrix <- cbind("positive" = c(tp, fp), "negative" = c(tn, fn)) 
# rownames(confusion_matrix) <- (c("true", "false"))
# fdr <-  tn / (tn + fn)
# confusion_matrix; paste("FDR = ",fdr)


# Group B vs C
gf_b_c <- dat_calls %>% filter(Call %in% c("B", "C")) %>% dplyr::select(1:3, 13:32)
factors_b_c <- factors %>% filter(V2 %in% c("B", "C"))
TwoStage_Package(gf_b_c[,-c(1,3)], factors_b_c, "gf_b_c_tmm1.csv", 1, .001)

# sig_b_c <- read_csv("gf_b_c_tmm1.csv")
# b_c <- left_join(sig_b_c, dat_calls, by = "Annotation") %>% dplyr::select(1:9)
# sig_ec <- c("2.7.7.7", "2.3.1.16", "2.4.2.14")
# sig_ec_found <- intersect(sig_ec, unique(b_c$label))
# ec_else <- setdiff(unique(b_c$label), sig_ec)
# tp <- sum(table(b_c$label)[sig_ec_found])
# fp <- sum(table(b_c$label)[ec_else])
# fn <- sum(table(gf_b_c$label)[sig_ec_found]) - tp
# tn <- nrow(gf_b_c) - sum(tp, fp, fn)
# confusion_matrix <- cbind("positive" = c(tp, fp), "negative" = c(tn, fn)) 
# rownames(confusion_matrix) <- (c("true", "false"))
# fdr <-  tn / (tn + fn)
# confusion_matrix; paste("FDR = ",fdr)

calc_confusion_mat <- function(df){
  new_df <- enquo(df)
  groups <- unlist(str_split(as.character(get_expr(new_df)), "_", n = 3))
  sig <- read_csv(paste0(get_expr(new_df),"_tmm1.csv"))
  sig_label <- left_join(sig, dat_calls, by = "Annotation") %>% dplyr::select(1:9)
  sig_ec <- c("2.7.7.7", "2.3.1.16", "2.4.2.14")
  sig_ec_found <- intersect(sig_ec, unique(sig_label$label))
  ec_else <- setdiff(unique(sig_label$label), sig_ec)
  tp <- sum(table(sig_label$label)[sig_ec_found])
  fp <- sum(table(sig_label$label)[ec_else])
  fn <- sum(table(df$label)[sig_ec_found]) - tp
  tn <- nrow(df) - sum(tp, fp, fn)
  confusion_matrix <- cbind("positive" = c(tp, fp), "negative" = c(tn, fn)) 
  rownames(confusion_matrix) <- (c("true", "false"))
  fdr <-  fp / (tp + fp)
  cat(paste("output for", groups[2],"&", groups[3], sep = " "), "\n")
  print(confusion_matrix)
  cat(paste("FDR = ",fdr), "\n\n")
}
calc_confusion_mat(gf_b_c)
calc_confusion_mat(gf_a_c)
calc_confusion_mat(gf_a_b)

