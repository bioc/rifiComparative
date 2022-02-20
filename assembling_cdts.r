#extract output from rifi_statistics
#
setwd("path_to_the_first_condition")
load("probe_sta.rdata")
inp_s <- probe_sta

#treated data
setwd("path_to_the_second_condition")
load("probe_sta.rdata")
inp_f <- probe_sta

setwd("path_new_directory_data_comparison/")
#add condition to probe_sta dataframe
inp_s$cdt <- "condition1"
inp_f$cdt <- "condition2"
#Important: intensity of both conditions should be normalized together
# prior to comparison
data_combined <- rbind(inp_s, inp_f)


