list.of.packages <- c("dplyr", "ggplot2", "gridExtra", "parallel", "data.table", "jsonlite", "purrr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
