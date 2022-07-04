# =========================================================================
#  fragmentation       Conveniently wraps all fragmentation steps
# -------------------------------------------------------------------------
#' 'fragmentation': fragments the half-life and intensity into segments
#' 'fragmentation': fragments the half-life and intensity into segments using 
#' the penalties output.
#'
#' @param data data frame: data frame combined data by column
#' @param pen_HL list: list of the penalties set optimal for the fragmentation
#' for half-life
#' @param pen_int list: list of the penalties set optimal for the fragmentation
#' for intensity
#' 
#' @return two data frames with half-life and intensity fragments and the mean 
#' of the coefficient fragment based. 
#'
#' @examples
#' data(df_comb_minimal)
#' data(pen_HL)
#' data(pen_int)
#  df_comb_minimal <- fragmentation(data=df_comb_minimal, pen_HL, pen_int)
#' @export

 
fragmentation <- function(data, pen_HL, pen_int){
#I. Dynamic Programming for HL
    data <- fragment_HL(
    probe = data,
    cores = 60,
    pen = pen_HL[[1]][[9]],
    pen_out = pen_HL[[1]][[10]]
    )

#II. Dynamic Programming for intensity
    data <- fragment_inty(
    probe = data,
    cores = 60,
    pen = pen_int[[1]][[9]],
    pen_out = pen_int[[1]][[10]]
    )
    return(data)
}


