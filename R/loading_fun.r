# =========================================================================
# loading_fun           
# -------------------------------------------------------------------------
#' loading_fun loads the data of all conditions
#'
#' loading_fun extract outputs from rifi_stats of all conditions and join
#' each data to the differential expression table. The differential gene
#' expression at time 0 needs to be run separately.
#' The columns log2FC, p_value adjusted, position and strand are extracted and
#' saved to a data frame. loading_fun joins the differential gene expression
#' table and the output from rifi statistics into one data frame.
#'
#' @param data1 data frame from rifi_stats of one condition
#' @param data2 data frame from rifi_stats of other condition
#' @param data3 data frame from differential gene expression at time 0
#'
#' @return A list of two data frames with joined columns from differential
#' expression and output of rifi_stats with the corresponding columns: ID with
#' position, strand, intensity, probe_TI, flag, position_segment, delay,
#' half_life, TI_termination_factor, delay_fragment, velocity_fragment,
#' intercept, slope, HL_fragment, HL_mean_fragment, intensity_fragment,
#' intensity_mean_fragment, TU, TI_termination_fragment,
#' TI_mean_termination_factor, seg_ID, pausing_site, iTSS_I, ps_ts_fragment,
#' event_ps_itss_p_value_Ttest, p_value_slope, delay_frg_slope, velocity_ratio,
#' event_duration, event_position, FC_HL, FC_fragment_HL, p_value_HL,
#' FC_intensity, FC_fragment_intensity, p_value_intensity, FC_HL_intensity,
#' FC_HL_intensity_fragment, FC_HL_adapted, synthesis_ratio,
#' synthesis_ratio_event, p_value_Manova, p_value_TI, cdt (condition),
#' logFC_int (log2FC(intensity)), P.Value of log2FC(intensity).
#'
#' @examples
#' data(stats_se_cdt1)
#' data(stats_se_cdt2)
#' data(differential_expression)
#' inp_s <-
#' loading_fun(stats_se_cdt1, stats_se_cdt2, differential_expression)[[1]]
#' inp_f <-
#' loading_fun(stats_se_cdt1, stats_se_cdt2, differential_expression)[[2]]
#' 
#' @export

loading_fun <- function(data1, data2, data3) {
    inp_s <- as.data.frame(rowRanges(data1))
    
    #extract output from rifi_statistics (second condition)
    #set the correct path: setwd("path to the second data")
    inp_f <- as.data.frame(rowRanges(data2))
    
    #set the correct path for the analysis
    #add condition to both dataframes
    inp_s$cdt <- "cdt1"
    inp_f$cdt <- "cdt2"
    
    # Important: intensity of both conditions should be normalized together
    # prior to comparison
    # Run the differential expression at probe or bin level of
    # intensity at time 0. The result is log2FC of both condition
    # Pick-up the columns:  log2FC, p_value adjusted, position and strand
    # call the first 2 columns: logFC_int and P.Value and save the object as dc
    
    # get the differential expression object
    inp_s <-
        left_join(inp_s[, -c(1:4)], data3, by = c("position", "strand"))
    inp_f <-
        left_join(inp_f[, -c(1:4)], data3, by = c("position", "strand"))
    
    return(list(inp_s, inp_f))
}
