
##joining the data by row
joining_data_row <- function(input.1, input.2){
    data <- rbind(input1, input2)
    return(data)
}
data_combined <- joining_data_row(input.1 = inp_s, input.2 = inp_f)
save(data_combined, file = "data_combined_se.rda")