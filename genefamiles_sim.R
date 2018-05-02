setwd("~/GitHub/drought_metagen/SIM_CONFUSION/")
gf_files <- list.files(path = "./", pattern = "genefamilies.tsv", recursive = TRUE, include.dirs = TRUE)
data <- map(gf_files, read_tsv)
all_data <- data.frame()
for (i in seq_along(data)){
  all_data <- rbind(all_data,data[[i]][,1])
}
row_track <- 1
for (i in seq_along(data)) {
  col <- i + 1
  if (i == 1) {
    all_data[row_track:(row_track+nrow(data[[i]])-1),col] <- data[[i]][,2]
    row_track <- row_track + nrow(data[[i]]) - 1
  } else {
    all_data[(row_track+1):(row_track+nrow(data[[i]])),col] <- data[[i]][,2]
    row_track <- row_track + nrow(data[[i]])
  }
  
  
}
all_data1 <- all_data %>% mutate_all(funs(replace(., is.na(.), 0)))

all_data2 <- all_data1[(which(all_data1[,1] != "UNMAPPED")),]
