plot_decay_synt <- function(data) {
    
    colnames(data)[11:13] <- c("fc_dr", "ratio_dr_ss", "fc_sr")
    
    #plot log2FC decay rate vs log2FC synthesis rate highlight some genes of
    #interest.
    png("Decay_rate_vs_Synthesis_rate.png", 900, 900, res = 200)
    
    ds <- data %>%
        ggplot(aes(x = fc_dr, y = fc_sr)) +
        geom_point(col = "magenta", size = .5) +
        geom_point(
            data = . %>% filter(fc_dr >= .5 & fc_sr <= -.5),
            color = "yellow",
            size = .5
        ) +
        
        #synthesis rate
        geom_hline(yintercept = median(na.omit(data$fc_sr)), show.legend = T) +
        geom_hline(
            colour = "grey60",
            linetype = "dashed",
            size = 1,
            yintercept = .5,
            show.legend = T
        ) +
        geom_hline(
            colour = "grey60",
            linetype = "dashed",
            size = 1,
            yintercept = -.5,
            show.legend = T
        ) +
        geom_hline(colour = "grey60",
                   yintercept = 0,
                   show.legend = T) +
        
        #decay rate
        geom_vline(xintercept = median(na.omit(data$fc_dr)), show.legend = T) +
        geom_vline(
            colour = "grey80",
            linetype = "dashed",
            size = 1,
            xintercept = -.5,
            show.legend = T
        ) +
        geom_vline(
            colour = "grey80",
            linetype = "dashed",
            size = 1,
            xintercept = +.5,
            show.legend = T
        ) +
        geom_vline(colour = "grey80",
                   xintercept = 0,
                   show.legend = T) +
        
        #mRNA
        geom_abline(
            intercept = 0.7,
            slope = median(na.omit(data$intensity_FC)),
            show.legend = T
        ) +
        geom_abline(
            colour = "grey70",
            linetype = "dashed",
            size = 1,
            slope = 1,
            intercept = .7,
            show.legend = T
        ) +
        geom_abline(
            colour = "grey70",
            linetype = "dashed",
            size = 1,
            slope = 1,
            intercept = -.7,
            show.legend = T
        ) +
        geom_abline(colour = "grey70", slope = 1) +
        
        scale_y_continuous(breaks = seq(-5, 5, by = 1), limits = c(-5, 5)) +
        scale_x_continuous(breaks = seq(-5, 5, by = 1), limits = c(-5, 5)) +
        labs(x = "log2FC(Decay rate)", y = "log2FC(Synthesis rate)") +
        # geom_text_repel(data = . %>% filter(fc_dr >= .5 & fc_sr <= -.5),
        #                 aes(label=gene), size=4, cex=3,
        #                 arrow = arrow(length = unit(0.01, "npc"),
        #                               type = "closed", ends = "first"),
        #                 colour="blue", segment.colour="blue")+
        theme(panel.background = element_blank())
    print(ds)
    dev.off()
}
