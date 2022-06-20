setwd("path_new_directory_data_comparison/")

###joining both data by column
#' Title
#'
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
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

df_comb <- joining_data_column(data = data_combined)
save(df_comb, file="df_comb_se.rda")
