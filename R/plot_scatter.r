#plot scatter of HL from both conditions at probe/bin level
plot_scatter <- function(data,
                         y = 30,
                         x = 30,
                         limits = c(0, 20)) {
    png("scatter_plot_HL.png", 900, 900, res = 200)
    
    s <- ggplot(data, aes(x = half_life.cdt1, y = half_life.cdt2)) +
        geom_point(alpha = .5) +
        scale_y_continuous(breaks = seq(0, y, by = 2), limits = limits) +
        scale_x_continuous(breaks = seq(0, x, by = 2), limits = limits) +
        geom_smooth(
            method = 'lm',
            se = TRUE,
            span = .01,
            colour = "yellowgreen",
            show.legend = TRUE
        ) +
        labs(x = "\nHalf-life in standard conditions [min]",
             y = "Half-life in iron depletion [min]\n") +
        theme(
            panel.background = element_blank(),
            legend.position = "none",
            axis.text = element_text(size = 13, face = "bold"),
            axis.text.x = element_text(colour = "black"),
            axis.line = element_line(colour = "dimgrey"),
            text = element_text(size = 12)
        )
    print(s)
    dev.off()
}
