plot_volcano <- function(input) {
    ##load data with log2FC intensity by segment
    png("volcano_plot.png", 900, 900, res = 200)
    v <- input %>%
        ggplot(aes(x = logFC_int,
                   y = -log10(P.Value))) +
        geom_point() +
        geom_point(
            data = . %>% filter(-log10(P.Value) > 3 & logFC_int < 0),
            color = "cornflowerblue",
            size = .5
        ) +
        geom_point(
            data = . %>% filter(-log10(P.Value) > 3 & logFC_int > 0),
            color = "firebrick",
            size = .5
        ) +
        labs(x = "log2FC(mRNA)", y = "-log10(P.Value)",
             title = "cdt2 vs. cdt1") +
        theme(title = element_text(face = "bold"))
    print(v)
    dev.off()
}