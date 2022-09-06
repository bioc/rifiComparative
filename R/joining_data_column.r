# =========================================================================
# joining_data_column           Joining two data frames by column
# -------------------------------------------------------------------------
#' 
#''joining_data_column': joins two data frames from different conditions by
#' column.
#'
#' @param data data frame with joined columns from both conditions
#'
#' @return the data frame with joined columns from both conditions with the 
#' corresponding columns: strand, position, ID, intensity.cdt1, position_segment,
#' half_life.cdt1, TI_termination_factor.cdt1", HL_fragment.cdt1, 
#' intensity_fragment.cdt1, TI_termination_fragment.cdt1, logFC_int, P.Value, 
#' intensity.cdt2, half_life.cdt2, TI_termination_factor.cdt2, HL_fragment.cdt2,
#' intensity_fragment.cdt2, TI_termination_fragment.cdt2
#' cdt1: first condition, cdt2: second condition.
#'
#' @examples
#' data(data_combined_minimal)
#' df_comb_minimal <- joining_data_column(data = data_combined_minimal)
#' 
#' @export

joining_data_column <- function(data){
    sta_cdt1 <- data %>%
        filter(cdt == "cdt1")
    sta_cdt1 <- sta_cdt1[, -c(4, 6, 7, 9, 12:15, 17, 19:20, 22:47)]
    sta_cdt2 <- data %>%
        filter(cdt == "cdt2")
    sta_cdt2 <- sta_cdt2[, -c(4, 6, 7, 9, 12:15, 17, 19:20, 22:49)]
    df_comb <- full_join(
        sta_cdt1,
        sta_cdt2,
        by = c("position", "strand",
               "ID", "position_segment"),
        suffix = c(".cdt1", ".cdt2")
    )
    return(df_comb)
} 


