#strand_function plots HL and intensity segments of both conditions.

strand_function <-
    function(data,
             data_p,
             data_n,
             data_p_c,
             data_n_c,
             Strand,
             condition,
             frag,
             i,
             fontface,
             axis_text_y_size,
             axis_title_y_size
             ) {
        #first plot for intensity segments
        for (j in seq_along(Strand)) {
            df <- data.frame()
            if (Strand[j] == "+") {
                df <- data_p
                df_c <- data_p_c
            } else{
                df <- data_n
                df_c <- data_n_c
            }
            if (nrow(df) != 0 & nrow(df_c) != 0) {
                p1 <-
                    ggplot(df, aes(x = get('position'))) +
                    scale_x_continuous(limits = c(frag[i], frag[i + 1])) +
                    labs(y = "Intensity [log2FC]") +
                    theme_bw() +
                    background_grid(major = "xy", minor = "none") +
                    theme(
                        legend.title = element_blank(),
                        axis.title.y = element_text(colour = 5,
                                                    size = axis_title_y_size),
                        axis.text.y = element_text(
                            angle = 90,
                            hjust = 1,
                            size = axis_text_y_size
                        ),
                        axis.text.x = element_blank(),
                        axis.title.x = element_blank(),
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor.x = element_blank(),
                        legend.position = "none",
                        plot.margin = margin(.1, .2, .1, .2, "cm"),
                        panel.border = element_blank()
                    )
                
                #######################segment plot###################
                #first plot for half-life segments
                #select segments without outliers
                df <- indice_function(df, "HL_fragment")
                # increase the limit to 20 in case 3 or more probes/bins have a HL
                # above 10
                Limit_h_df1 <-
                    limit_function(df, "half_life", ind = 1)
                if (Limit_h_df1 == 20) {
                    Breaks_h <- seq(0, Limit_h_df1, by = 4)
                } else{
                    Breaks_h <- seq(0, Limit_h_df1, by = 2)
                }
                df1.h <- secondaryAxis(df, "half_life", ind = 1)
                #in case only one bin is available and the HL is above 20
                if (all(df$half_life > 20)) {
                    df$half_life <- 20
                }
                p2 <-
                    ggplot(df, aes(x = get('position'))) +
                    scale_x_continuous(limits = c(frag[i], frag[i + 1])) +
                    labs(y = "Half-life [min]") +
                    theme_bw() +
                    background_grid(major = "xy", minor = "none") +
                    theme(
                        legend.title = element_blank(),
                        legend.position = "none",
                        axis.title.y = element_text(colour = 6,
                                                    size = axis_title_y_size),
                        axis.text.x = element_text(size = 6),
                        axis.title.x = element_blank(),
                        axis.text.y = element_text(
                            angle = 90,
                            hjust = 1,
                            size = axis_text_y_size
                        ),
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor.x = element_blank(),
                        plot.margin = margin(.1, .2, .1, .2, "cm"),
                        panel.border = element_blank()
                    )
                #add the second axis for half-life segments plot
                if (length(unique(df$half_life)) == 1) {
                    if (is.na(unique(df$half_life))) {
                        p2 <- p2 +
                            geom_blank() +
                            scale_y_continuous(
                                limits = c(0, Limit_h_df1),
                                breaks = Breaks_h,
                                sec.axis = sec_axis( ~ . * 1,
                                                     name = "Half-life [min]",
                                                     breaks = Breaks_h)
                            )
                    }
                }
               ######################intensity plot################
                    #add a reference to outliers probes or bins
                    #plot distance and new fragments results of DP of both
                    #conditions
                    df_c <-
                        indice_function(df_c, "intensity_comb_fragment")
                    df1_wo <-
                        meanPosition(df_c %>%
                                         filter(indice==1),
                                     "intensity_comb_fragment")
                    if (nrow(df_c %>% filter(indice == 1)) != 0) {
                        p1 <- p1 +
                            geom_point(
                                data = filter(
                                    df_c,
                                    indice == 1,
                                    distance_int < 0,
                                    distance_int > -5
                                ),
                                aes(x = position, y = distance_int),
                                col = 6,
                                size = .5
                            ) +
                            geom_point(
                                data = filter(
                                    df_c,
                                    indice == 1,
                                    distance_int > 0,
                                    distance_int < 5
                                ),
                                aes(x = position, y = distance_int),
                                col = 4,
                                size = .5
                            ) +
                            geom_line(data =  filter(df_c, indice == 1),
                                      aes(
                                          x = get('position'),
                                          y = get('intensity_mean_comb_fragment'),
                                          col = get('intensity_comb_fragment')
                                      )
                            ) +
                            geom_text(
                                data = df1_wo,
                                aes(
                                    x = get('meanPosi'),
                                    y = get('intensity_mean_comb_fragment'),
                                    label = get('intensity_comb_fragment')
                                ),
                                size = 1.3,
                                check_overlap = TRUE
                            )
                        if (nrow(df_c %>% filter(
                            indice == 1,
                            p_value_distance_intensity < 0.05
                        )) != 0) {
                            if (unique(df_c$strand) == "+") {
                                df_c_mean <- arrange_byGroup(
                                    df_c %>% filter(
                                        indice == 1,
                                        p_value_distance_intensity < 0.05
                                    ),
                                    "intensity_comb_fragment"
                                )
                            } else{
                                df_c_mean <- df_c %>% filter(indice == 1,
                                                             p_value_distance_intensity < 0.05)
                                df_c_mean <-
                                    df_c_mean[!duplicated(df_c_mean$p_value_distance_intensity), ]
                            }
                            p1 <- p1 +
                                geom_text(
                                    data = df_c_mean,
                                    aes(
                                        x = get('position'),
                                        y = get('intensity_mean_comb_fragment')
                                    ),
                                    label = "**",
                                    fontface = fontface,
                                    size = 2,
                                    check_overlap = TRUE
                                )
                            
                        }
                    }
                    if (nrow(df_c %>% filter(get('indice') == 2)) != 0) {
                        p1 <- p1 +
                            geom_point(
                                data = df_c %>%
                                    filter(
                                        get('indice') == 2 & distance_int > 0 &
                                            distance_int < 10
                                    ),
                                aes(x = position, y = distance_int),
                                col = 2,
                                shape = 17,
                                size = .5
                            ) +
                            geom_point(
                                data = df_c %>%
                                    filter(
                                        get('indice') == 2 & distance_int < 0 &
                                            distance_int > -10
                                    ),
                                aes(x = position, y = distance_int),
                                col = 3,
                                shape = 12,
                                size = .5
                            )
                    }
                    #######################HL plot###################
                    #plot distance and new fragments results of DP of both
                    #conditions
                    #dismiss outliers
                    df_c <-
                        indice_function(df_c, "HL_comb_fragment")
                    #add fragment label 
                    df1_wo <-
                        meanPosition(df_c %>%
                                         filter(indice==1),
                                     "HL_comb_fragment")
                    if (nrow(df_c %>% filter(indice == 1)) != 0) {
                        p2 <- p2 +
                            geom_point(
                                data = filter(
                                    df_c,
                                    indice == 1,
                                    distance_HL < 0,
                                    distance_HL < 5
                                ),
                                aes(x = position, y = distance_HL),
                                col = 6,
                                size = .5
                            ) +
                            geom_point(
                                data = filter(
                                    df_c,
                                    indice == 1,
                                    distance_HL > 0,
                                    distance_HL > -5
                                ),
                                aes(x = position, y = distance_HL),
                                col = 4,
                                size = .5
                            ) +
                            geom_line(data =  filter(df_c, indice == 1),
                                      aes(
                                          x = get('position'),
                                          y = get('HL_mean_comb_fragment'),
                                          col = get('HL_comb_fragment')
                                      )
                            )+
                            geom_text(
                                data = df1_wo,
                                aes(
                                    x = get('meanPosi'),
                                    y = get('HL_mean_comb_fragment'),
                                    label = get('HL_comb_fragment')
                                ),
                                size = 1.3,
                                check_overlap = TRUE
                            )
                            
                        if (nrow(df_c %>% filter(p_value_distance_HL < 0.05)) != 0) {
                            if (unique(df_c$strand) == "+") {
                                df_c_mean <- arrange_byGroup(
                                    df_c %>% filter(
                                        indice == 1,
                                        p_value_distance_HL < 0.05
                                    ),
                                    "HL_comb_fragment"
                                )
                            } else{
                                df_c_mean <- df_c %>% filter(indice == 1,
                                                             p_value_distance_HL < 0.05)
                                df_c_mean <-
                                    df_c_mean[!duplicated(df_c_mean$p_value_distance_HL), ]
                            }
                            p2 <- p2 +
                                geom_text(
                                    data = df_c_mean,
                                    aes(
                                        x = get('position'),
                                        y = get('HL_mean_comb_fragment')
                                    ),
                                    label = "**",
                                    fontface = fontface,
                                    size = 2,
                                    check_overlap = TRUE
                                )
                            
                        }
                    }
                    if (nrow(df_c %>% filter(get('indice') == 2)) != 0) {
                        p2 <- p2 +
                            geom_point(
                                data = df_c %>%
                                    filter(
                                        get('indice') == 2 & distance_HL > 0 &
                                            distance_HL < 10
                                    ),
                                aes(x = position, y = distance_HL),
                                col = 2,
                                shape = 17,
                                size = .5
                            ) +
                            geom_point(
                                data = df_c %>%
                                    filter(
                                        get('indice') == 2 & distance_HL < 0 &
                                            distance_HL > -10
                                    ),
                                aes(x = position, y = distance_HL),
                                col = 3,
                                shape = 12,
                                size = .5
                            )
                    }
                
            }
            if (Strand[j] == "+") {
                p <- list(p1, p2)
            } else {
                p4 <- p1 + coord_trans(y = "reverse")
                p5 <- p2 + coord_trans(y = "reverse")
                p.1 <- list(p4, p5)
            }
        }
        p <- c(p, p.1)
        return(p)
    }
