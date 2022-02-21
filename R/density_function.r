#density plot plots density of half-life in both conditions.
#the data is the combined output from rifi_statistics
#value is the maximum value of half-life to be plotted
#conditions are your conditions
density_plot <- function(data, value, condition){
pdf("density_HL.pdf")
ggplot(filter(data, HL_mean_fragment < value),
       aes(HL_mean_fragment, color = cdt)) +
    geom_density(alpha = 0.2) +
    geom_vline(data = data, 
               mapping = aes(xintercept = 
                                 mean(na.omit(HL_mean_fragment)), color = cdt)) +
    scale_color_discrete(labels = c(condition[1], condition[2])) +
    scale_x_continuous(breaks = seq(0, value, by = 4), limits = c(0, value)) +
    labs(x = "Half-life [min]\n") +
    ggtitle(paste("Decay density in", condition[1], "vs", condition[2])) +
    theme(
        panel.background = element_blank(),
        axis.line = element_line(colour = "dimgrey"),
        legend.title = element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = .5)
    )
dev.off()
}