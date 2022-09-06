#plot HL categories
plot_histogram <- function(data, cdt1, cdt2) {
    
    png("histogram_count.png", 900, 900, res = 200)
    
    hist1 <- HL_category(data, "cdt1")
    hist2 <- HL_category(data, "cdt2")
    
    category <- c("<=2", "2-4", "4-8", "8-10", "10-20", ">20")
    
    #add category to the dataframe
    hist_sf <- cbind.data.frame(hist1, hist2, category)
    
    #Value Standardization
    hist_sf[, 1:2] <- sapply(hist_sf[, 1:2], function(x)
        x / sum(x))
    
    colnames(hist_sf)[1:2] <- c(cdt1, cdt2)
    hist_sf <- reshape2::melt(hist_sf, id.vars = "category")
    
    h <- ggplot(hist_sf, aes(x = category, y = value,
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
             title = "Half-life classification") +
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
    print(h)
    dev.off()
}
