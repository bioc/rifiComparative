HL_category <- function(data, indice) {
    l1 <- length(which(data[, paste0("half_life.", indice)] <= 2))
    l2 <- length(which(data[, paste0("half_life.", indice)] > 2 &
                           data[, paste0("half_life.", indice)] <= 4))
    l3 <- length(which(data[, paste0("half_life.", indice)] > 4 &
                           data[, paste0("half_life.", indice)] <= 8))
    l4 <- length(which(data[, paste0("half_life.", indice)] > 8 &
                           data[, paste0("half_life.", indice)] <= 10))
    l5 <- length(which(data[, paste0("half_life.", indice)] > 10 &
                           data[, paste0("half_life.", indice)] <= 20))
    l6 <- length(which(data[, paste0("half_life.", indice)] > 20))
    return(c(l1, l2, l3, l4, l5, l6))
}