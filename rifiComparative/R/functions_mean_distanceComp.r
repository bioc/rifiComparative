p_value_function <- function(data){
    p_value_hl <- na.omit(unique(data[,"p_value_distance_HL"]))
    p_value_int <- na.omit(unique(data[,"p_value_distance_intensity"]))
return(c(p_value_hl, p_value_int))
}

eliminate_outlier_hl <- function(data){
    output <- data[grep(paste0("Dc_\\d+", "$"), data[,"HL_comb_fragment"]),]
    return(output)
}

eliminate_outlier_int <- function(data){
    output <- data[grep(paste0("I_\\d+", "$"),
                        data[,"intensity_comb_fragment"]),]
    return(output)
}

#extract the mean of logged values intensity and the length of the fragments implicated
mean_length_int <- function(data){
    output <- data %>%
        summarise(mean(distance_int), length(unique(intensity_comb_fragment)))
    colnames(output) <- c("mean", "length")
    return(output)
}

#extract the mean of logged values HL and the length of the fragments implicated
mean_length_hl <- function(data){
    output <- data %>%
        summarise(mean(distance_HL_log), length(unique(HL_comb_fragment)))
    colnames(output) <- c("mean", "length")
    return(output)
}

calculating_rate <- function(data){
    #log2FC of decay_rate
    dr <- log2((log(2)/mean(data$half_life.cdt1))/
             (log(2)/mean(data$half_life.cdt2)))
    #log2 FC of steady-state
    int <- log2(mean(data$intensity.cdt1)/
                   mean(data$intensity.cdt2))
    return(list(dr, int))
}
