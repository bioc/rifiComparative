plot_heatscatter <- function(data) {
    
    colnames(data)[11:13] <- c("fc_dr", "ratio_dr_ss", "fc_sr")
    
    png("mRNA_Synthesis-Decay_Compensation.png", 1200, 900, res = 200)
    
    h <- heatscatter(
        data$fc_dr,
        data$fc_sr,
        ylab = "log2 (Synthesis rate fold Fe-/Fe+)",
        xlab = "log2 (Decay rate fold Fe-/Fe+)",
        xlim = c(-3, 3),
        ylim = c(-3, 3),
        add.contour = TRUE,
        cor = TRUE,
        daltonize = T,
        cvd = "d",
        size = TRUE,
        separate = T,
        add.quartiles = TRUE,
        simulate = F,
        colpal = "spectral",
        color.contour = "black",
        main = "Fe- versus Fe+"
    )
    print(h)
    dev.off()
}