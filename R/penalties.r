# =========================================================================
#  penalties           
# -------------------------------------------------------------------------
#' penalties wraps conveniently all penalty steps
#' 
#' penalties finds the best set of penalties for half-life and intensity
#' fragmentation using dynamic programming. The segmentation of the HL uses the
#' difference between 2 conditions and the segmentation of the intensity uses 
#' the log2FC. 
#' 
#' The function uses 4 functions: 
#' 
#' score_fun_ave.r
#' 
#' make_pen.r
#' 
#' fragment_HL_pen.r
#' 
#' fragment_inty_pen.r
#'
#' @param data data frame with the joined columns from differential 
#' expression and output of rifi_stats.
#' @param cores integer: the number of assigned cores for the task. It needs to 
#' be increased in case of big data.
#' 
#' @return A list of data frame and penalties, The first element is data frame
#' with 2 more variables, second and third are HL and intensity penalties
#' respectively.
#'
#' @examples
#' data(df_comb_minimal) 
#' penalties_df <- penalties(df_comb_minimal)[[1]]
#' pen_HL <- penalties(df_comb_minimal)[[2]]
#' pen_int <- penalties(df_comb_minimal)[[3]]
#' @export


penalties <- function(data, cores = 2){
# calculate the difference of half-life from both conditions. 
# difference is referred to distance 

    data[,"distance_HL"] <-
        data[, "half_life.cdt1"] - data[, "half_life.cdt2"]
    
# find the best penalties set for half-life fragmentation using dynamic
# programming on half-life distance 

    pen_HL <- make_pen(
        probe = data,
        FUN = fragment_HL_pen,
        cores = cores,
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
    
    ##DP for log2FC(intensity) penalties
    #add log2FC(intensity) to the data frame
    data[,"distance_int"] <- data[,"logFC_int"]
    
    # dynamic programming on log2FC(intensity) to find the best penalties
    pen_int <- make_pen(
        probe = data,
        FUN = fragment_inty_pen,
        cores = cores,
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
    return(list(data, pen_HL, pen_int))
}


