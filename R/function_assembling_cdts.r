##########################preparation of intensity DP##########################
#get the fit1 from normalized data (sc and fe) from Limma at probe level
#load fit from limma, contrast of 2 conditions, intensity is replaced

fit2_df <- get(load("fit2_sc_fe_notAveraged.rda"))

#' Title
#'
#' @param data 
#' @param input.1 
#' @param input.2 
#'
#' @return list containing outputs
#' @export
#'
#' @examples
#' 
intensity_normalization_position <- function(data, input.1, input.2){ 

    data <- data[, c(8, 9, 10, 58, 61)]
    colnames(data)[c(3,4)] <- c("strand", "logFC_int")

#replace orientation (+ with - and vice versa as in limma orientation is opposite)
plus <- which(data$strand == "+")
minus <- which(data$strand == "-")
data$strand[plus] <- "-"
data$strand[minus] <- "+"
positive <- data[ ,c(2:5)] %>%
    filter(strand == "+") %>%
    arrange(end)
colnames(positive)[1] <- "position"
negative <- data[ ,c(1,3:5)] %>%
    filter(strand == "-") %>%
    arrange(start)
colnames(negative)[1] <- "position"
data <- rbind(positive, negative)

#save(data, file="Data_SC_iron_fit1_positionSE.rdata")
dc <- get(load("Data_SC_iron_fit1_positionSE.rdata"))
#join the data together
input.2[,"ID"] <- as.numeric(input.2[,"ID"])

input.2 <- left_join(input.2[,-c(1:4)], dc, by = c("position", "strand"))

input.1[,"ID"] <- as.numeric(input.1[,"ID"])
input.1 <- left_join(input.1[,-c(1:4)], dc, by = c("position","strand"))

 return(list(input.1, input.2))
}