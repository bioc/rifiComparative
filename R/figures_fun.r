# =========================================================================
# figures_fun                    
# -------------------------------------------------------------------------
#' 'figures_fun': generates several plots 
#' 
#' 'figures_fun' plots at one the density of HL, the HL category as histogram, 
#' log2FC of decay and synthesis rate and their heatscatter. Scatter plot of HL
#' and volcano plot. The function uses the four output generated previously.
#' 
#' The functions used are:
#' 
#'  plot_decay_synt: plots the changes in RNA decay rates versus the changes 
#'  in RNA synthesis rates 
#'    
#'  plot_heatscatter: plots the changes in RNA decay rates versus the changes 
#'  in RNA synthesis rates with density.
#'  
#'  plot_volcano: plots statistical significance  versus magnitude of change .
#'  
#'  plot_histogram: plot a histogram of probe/bin half-life categories from 
#'  2 to 20 minutes in both conditions.
#'  
#'  plot_density: plots the probe/bin half-life density in both conditions.
#'  
#'  plot_scatter: plots of the bin/probe half-life in one condition1 vs.
#'  condition2.
#' 
#' extract the object of rifi_statistics from both conditions. 
#' The differential gene expression at time 0 needs to be run separately. 
#' The columns log2FC, p_value adjusted, position and strand are extracted and 
#' saved to a data frame. loading_fun_fig joins the differential gene expression 
#' table and the output from rifi statistics into one data frame.
#'
#' @param data.1 data frame output of statistic
#' @param data.2 data frame joining two outputs from rifi_stats by row
#' @param input.1 data frame joining two outputs from rifi_stats by column
#' @param input.2 data frame of differential gene expression at time 0
#' @param cdt1 string for the first condition
#' @param cdt2 string for the second condition
#' @param y integer to break the scaling in scatter plot for y_axis
#' @param x integer to break the scaling in scatter plot for x_axis
#' @param limits vector to limit the scaling in scatter plot for both axis
#'
#' @return several plots
#'
#' @examples
#' data(data_combined_minimal)
#' data(df_comb_minimal)
#' data(differential_expression) 
#' data(df_mean_minimal)
#' figures_fun(data.1 = df_mean_minimal, data.2 = data_combined_minimal, 
#' input.1 = df_comb_minimal, input.2 = differential_expression, cdt1 = "sc", 
#' cdt2 = "fe") 
#' @export


figures_fun <- function(data.1, data.2, input.1, input.2, cdt1, cdt2, y = 30, 
                        x = 30, limits = c(0, 20)) {
    
    # log2FC(decay.rate) vs. log2FC(synthesis.rate)
    # data.1 = df_mean_minimal
    plot_decay_synt(data.1)
    
    # plot heatscatter log2FC(decay.rate) vs. log2FC(synthesis.rate)
    # data.1 = df_mean_minimal
    plot_heatscatter(data.1)
    
    # plot density of HL
    # data.2 = data_combined_minimal
    plot_density(data.2, cdt1 = cdt1, cdt2 = cdt2)
    
    # plot histogram of HL classification
    # input.1 = df_comb_minimal
    plot_histogram(input.1, cdt1 = cdt1, cdt2 = cdt2)
    
    # plot scatter of HL
    # input.1 = df_comb_minimal
    plot_scatter(input.1)
    
    # plot volcano plot
    # input.2 = differential_expression
    plot_volcano(input.2)
    
}


