setwd("path_to_your_dataframe")
load("summary_cdts.rdata")

time_v <- c(2,4,8,10,20)
indice <- c("cdt1", "cdt2")

HL_count <- function(data, indice, time_v) {
    l1 <- length(which(data[, paste0("HL_", indice)] <= time_v[1]))
    l <- c()
    for(i in 1:(length(time_v)-1)){
        l <- c(l, length(which(data[, paste0("HL_", indice)] > time_v[i] &
                               data[, paste0("HL_", indice)] <= time_v[i+1]))) 
        }
    l2 <- length(which(data[, paste0("HL_", indice)] > time_v[length(time_v)]))
    return(c(l1, l, l2))
}

hist_function <- function(data, indice){
hist1 <- HL_count(data, indice = condition[1], time_v = time_v)
hist2 <- HL_count(data, indice = condition[2], time_v = time_v)

category <- c("<=2", "2-4", "4-8", "8-10", "10-20", ">20")

#add category to the dataframe
hist.bothC <- cbind.data.frame(hist1, hist2, category)

hist.bothC[, 1:2] <- sapply(hist.bothC[, 1:2], function(x)
    x / sum(x))

colnames(hist.bothC)[1:2] <- c(indice[1], indice[2])
hist.bothC <- melt(hist.bothC)

pdf("histogram_count")
hs <- ggplot(hist.bothC, aes(x = category, y = value,
                          fill = variable)) +
    geom_bar(
        color = "black",
        stat = "identity",
        position = position_dodge2(),
        orientation = "x",
        show.legend = TRUE
    ) +
    scale_x_discrete(limits = category) +
    labs(x = "\nHalf-life [min]", y = "Relative frequency % \n",
         title = paste("Half-life classification in", indice[1], "and", indice[2])) +
    scale_fill_brewer(palette = 'Purples') +
    theme(
        legend.position = c(.85, .9),
        legend.title = element_blank(),
        axis.text = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(colour = "black"),
        axis.line = element_line(colour = "dimgrey"),
        text = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = .5)
    ) +
    theme_bw()
print(hs)
dev.off()
}
