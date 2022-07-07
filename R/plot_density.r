# density plot on HL
plot_density <- function(data, cdt1, cdt2) {
    
    png("density_HL.png", 900, 900, res = 200)
    
    d <- ggplot(filter(data, HL_mean_fragment < 30),
           aes(HL_mean_fragment, color = cdt)) +
        geom_density(alpha = 0.2) +
        geom_vline(data = data,
                   mapping = aes(xintercept =
                                     mean(na.omit(
                                         HL_mean_fragment
                                     )), color = cdt)) +
        scale_color_discrete(labels = c(cdt1, cdt2)) +
        scale_x_continuous(breaks = seq(0, 30, by = 4), limits = c(0, 30)) +
        labs(x = "Half-life [min]\n") +
        ggtitle("Decay density") +
        theme(
            panel.background = element_blank(),
            axis.line = element_line(colour = "dimgrey"),
            legend.title = element_blank(),
            plot.title = element_text(
                size = 12,
                face = "bold",
                hjust = .5
            )
        )
    print(d)
    dev.off()
    
}