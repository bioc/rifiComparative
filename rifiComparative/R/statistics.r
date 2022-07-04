# ==============================================================================
# statistics                Statistical prediction between 2 segments
# ------------------------------------------------------------------------------
# 
# 'statistics': check segment significance using t-test 
#' 'statistics' uses t-test to check HL and intensity segments significance. The 
#' function return the data frame with p_value p_value adjusted. 
#' The function used is t_test_function.
#'
#' @param data data frame of the joined column
#' 
#' @return a list of two data frames, the first one contains all segments with 
#' p_value and p_value adjusted. The second one removes the duplicated segments 
#' from intensity and could be saved as excel file. 
#'
#' @examples
#' data(df_comb_minimal)
#' df_comb_minimal <- statistics(data= fragment_int)[[1]]
#' df_comb_uniq_minimal <- statistics(data= fragment_int)[[2]]
#' write_xlsx(df_comb_uniq_minimal, "df_comb_uniq_minimal.xlsx")
#' @export

statistics <- function(data){
    
    #eliminate outlier from HL fragment
    frag_HL <-
        unique(data[,"HL_comb_fragment"][
            grep(paste0("Dc_\\d+$"), data[,"HL_comb_fragment"])])
    
    #eliminate outlier from intensity fragment
    frag_int <-
        unique(data[,"intensity_comb_fragment"][
            grep(paste0("I_\\d+$"), data[,"intensity_comb_fragment"])])
    
    data <-
        t_test_function(
            data = data,
            par = "HL",
            par1 = "half_life"
        )
    p_adjusted <-
        as.numeric(p.adjust(as.numeric(data[ ,"p_value_distance_HL"]),
                            method = "fdr"))
    data <-
        tibble::add_column(data, formatC(p_adjusted, format = "e", digits = 2),
                           .after = ncol(data)-1)
    
    colnames(data)[ncol(data)] <- "p_adjusted_HL"
    
    data <-
        t_test_function(
            data = data,
            par = "intensity",
            par1 = "intensity"
        )

        p_adjusted <-
            as.numeric(p.adjust(as.numeric(
                data[ ,"p_value_distance_intensity"]), method = "fdr"))
        data <-
            tibble::add_column(data, formatC(p_adjusted, format = "e", 
                                             digits = 2), .after = ncol(data)-1)
        
        colnames(data)[ncol(data)] <- "p_adjusted_intensity"
        
        df_comb_uniq <-
            data[!duplicated(data[,"intensity_comb_fragment"]), ]
        
        return(c(data, df_comb_uniq))
        }
