# =========================================================================
#  joining_data_row           Joining two data frames by row
# -------------------------------------------------------------------------
#' 
#' 'joining_data_row': joins two data frames from different conditions by row
#'
#' @param input1 data frame from SummarizedExperiment output of rifi_stats 
#' from rifi package of the first condition. 
#' @param input2 data frame from SummarizedExperiment output of rifi_stats 
#' from rifi package of the second condition.
#'
#' @return the data frame with joined rows from both conditions with the 
#' corresponding columns: ID with position, strand, intensity, probe_TI, flag, 
#' position_segment, delay, half_life, TI_termination_factor, delay_fragment, 
#' velocity_fragment, intercept, slope, HL_fragment, HL_mean_fragment, 
#' intensity_fragment, intensity_mean_fragment, TU, TI_termination_fragment, 
#' TI_mean_termination_factor, seg_ID, pausing_site, iTSS_I, ps_ts_fragment, 
#' event_ps_itss_p_value_Ttest, p_value_slope, delay_frg_slope, velocity_ratio, 
#' event_duration, event_position, FC_HL, FC_fragment_HL, p_value_HL, 
#' FC_intensity, FC_fragment_intensity, p_value_intensity, FC_HL_intensity, 
#' FC_HL_intensity_fragment, FC_HL_adapted, synthesis_ratio, 
#' synthesis_ratio_event, p_value_Manova, p_value_TI, cdt (condition), 
#' logFC_int (log2FC(intensity)), P.Value of log2FC(intensity)
#'
#' @examples
#' data(inp_s)
#' data(inp_f)
#' data_combined_minimal <- 
#' joining_data_row(input1 = inp_s, input2 = inp_f)
#' @export
 

joining_data_row <- function(input1, input2){
    data <- rbind(input1, input2)
    return(data)
}
