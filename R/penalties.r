# looking into the data from 2 cdts, the segments comparative approach is a complicated process since the segments are different for one cdt to another.
# We developed a workflow to make the comparative easier. We segmented the genome upon HL difference between 2 cdts and Intensity log2FC on the same way. 
# The approach used to segment the genome is dynamic programming.
# HL is segmented using difference between 2 cdts. Intensity was segmented using log2FC.

source("score_fun_ave.r")
source("make_pen.r")
source("fragment_HL_pen.r")
source("fragment_inty_pen.r")

setwd("path_new_directory_data_comparison/")
load("data_combined_se.rda")
load("df_comb_se.rda")

#' Title
#'
#' @param data 
#'
#' @return
#' @export
#'
#' @examples

penalties <- function(data){
# calculate the difference of half-life from both conditions. 
# difference is referred to distance 
data[,"distance_HL"] <-
    data[, "half_life.cdt1"] - data[, "half_life.cdt2"]

# find the best penalties set for half-life fragmentation using dynamic programming 
# on half-life distance 
pen_HL <- make_pen(
    probe = data,
    FUN = fragment_HL_pen,
    cores = 60,
    logs = as.numeric(rep(NA, 8)),
    dpt = 1,
    smpl_min = 10,
    smpl_max = 50,
    sta_pen = 0.5,
    end_pen = 4.5,
    rez_pen = 9,
    sta_pen_out = 0.5,
    end_pen_out = 3.5,
    rez_pen_out = 7
)

##########################DP for log2FC(intensity) penalties###################
#add log2FC(intensity) to the data frame
data[,"distance_int"] <- data[,"logFC_int"]

# dynamic programming on log2FC(intensity) to find the best penalties
pen_int <- make_pen(
    probe = data,
    FUN = fragment_inty_pen,
    cores = 60,
    logs = as.numeric(rep(NA, 8)),
    dpt = 1,
    smpl_min = 10,
    smpl_max = 50,
    sta_pen = 0.5,
    end_pen = 4.5,
    rez_pen = 9,
    sta_pen_out = 0.5,
    end_pen_out = 3.5,
    rez_pen_out = 7
)
return(list(pen_HL, pen_int))
}

pen_HL <- penalties(data)[[1]]
pen_int <- penalties(data)[[2]]
save(pen_HL, file="pen_HL.rda")
save(pen_int, file="pen_int.rda")

