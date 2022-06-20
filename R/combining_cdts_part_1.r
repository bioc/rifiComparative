#extract output from rifi_statistics (first condition)
setwd("path_to_the_first_condition")
inp_s <- get(load("stats_se.rda"))
inp_s <- as.data.frame(rowRanges(inp_s))

#extract output from rifi_statistics (second condition)
setwd("path_to_the_second_condition")
inp_f <- get(load("stats_se.rda"))
inp_f <- as.data.frame(rowRanges(inp_f))

setwd("path_new_directory_data_comparison")
#add condition to both dataframes
inp_s$cdt <- "cdt1"
inp_f$cdt <- "cdt2"

# Important: intensity of both conditions should be normalized together
# prior to comparison

##############################################################################
#'* Run the differential expression at probe or bin level of*
#'* intensity at time 0. The result is log2FC of both condition* 
#'* Pick-up the columns:  log2FC, p_value adjusted, position and strand:*
#'* call the first 2 columns: logFC_int and P.Value and save the object as dc*'

inp_s <- left_join(inp_s[,-c(1:4)], dc, by = c("position","strand"))
inp_f <- left_join(inp_f[,-c(1:4)], dc, by = c("position", "strand"))
##############################################################################

#save objects
save(inp_s, file="inp_s.rda")
save(inp_f, file="inp_f.rda")

