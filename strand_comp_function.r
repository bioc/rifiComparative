#strand_function plots delay, HL and intensity segments of both conditions.
#strand_function plots events as pausing sites and iTSS_I.
strand_function <-
    function(data,
             data_p,
             data_n,
             Strand,
             condition,
             coverage,
             frag,
             i,
             fontface,
             HL_threshold_1,
             HL_threshold_2,
             HL_threshold_color_1,
             axis_text_y_size,
             axis_title_y_size,
             p_value_int,
             p_value_event,
             p_value_hl,
             event_duration_ps,
             event_duration_itss
           ) {
        #first plot for intensity segments
        for (j in seq_along(Strand)) {
            df <- data.frame()
            if (Strand[j] == "+") {
                df <- data_p
            }else{
                df <- data_n
            }
        if (nrow(df) != 0) {
        p1 <-
            ggplot(df, aes(
                x = get('position'),
                y = get('intensity'),
                col = cdt
            )) +
            scale_x_continuous(limits = c(frag[i], frag[i + 1])) +
            scale_y_continuous(
                trans = 'log2',
                labels = label_log2_function,
                limits = c(NA, NA),
                sec.axis = sec_axis(~ . * 1, name = "Coverage",
                                    labels = label_square_function)
            ) +
            labs(y = "Intensity [A.U]") +
            theme_bw() +
            background_grid(major = "xy", minor = "none") +
            theme(
                legend.title = element_blank(),
                axis.title.y = element_text(colour = 5, size = axis_title_y_size),
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
        
        #######################segment plot positive strand###################
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
            ggplot(df, aes(
                x = get('position'),
                y = get('half_life'),
                col = cdt
            )) +
            scale_x_continuous(limits = c(frag[i], frag[i + 1])) +
            scale_y_continuous(
                limits = c(0, Limit_h_df1),
                breaks = Breaks_h,
                sec.axis = sec_axis(~ . * 1, name = "Half-life [min]", 
                                    breaks = Breaks_h)
            ) +
            labs(y = "Half-life [min]") +
            theme_bw() +
            background_grid(major = "xy", minor = "none") +
            theme(
                legend.title = element_blank(),
                legend.position = "none",
                axis.title.y = element_text(colour = 6, 
                                            size = axis_title_y_size),
                axis.text.x = element_blank(),
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
                        sec.axis = sec_axis(~ . * 1, 
                                            name = "Half-life [min]",
                                            breaks = Breaks_h)
                    )
            }
        }
        #select segments without outliers
        df <- indice_function(df, "delay_fragment")
        #increase the limit to 20 in case 3 or more probes/bins have a delay
        #above 10
        Limit_df1 <- limit_function(df, "delay", ind = 1)
        if (Limit_df1 == 20) {
            Breaks_d <- seq(0, Limit_df1, by = 4)
        } else{
            Breaks_d <- seq(0, Limit_df1, by = 2)
        }
        #first plot for delay segments
        df1.d <- secondaryAxis(df, "delay", ind = 1)
        #in case only one bin is available and the delay is above 20
        if (all(df$delay > 20)) {
            df$delay <- 20
        }
        p3 <-
            ggplot(df, aes(
                x = get('position'),
                y = get('delay'),
                col = cdt
            )) +
            scale_x_continuous(limits = c(frag[i], frag[i + 1])) +
            scale_y_continuous(
                limits = c(0, Limit_df1),
                breaks = Breaks_d,
                sec.axis = sec_axis(~ . * 1, name = "Delay [min]", breaks =
                                        Breaks_d)
            ) +
            labs(y = "Delay [min]") +
            theme_bw() +
            background_grid(major = "xy", minor = "none") +
            theme(
                legend.title = element_blank(),
                axis.title.x = element_blank(),
                legend.position = "none",
                axis.title.y = element_text(colour = 4, 
                                            size = axis_title_y_size),
                axis.text.y = element_text(
                    angle = 90,
                    hjust = 1,
                    size = axis_text_y_size
                ),
                axis.text.x = element_text(size = 6),
                panel.grid.major.x = element_blank(),
                plot.title = element_blank(),
                plot.margin = margin(.1, .2, .2, .2, "cm"),
                panel.border = element_blank()
            )
        #add the second axis for delay segments plot
        if (length(unique(df$delay)) == 1) {
            if (is.na(unique(df$delay))) {
                p3 <- p3 +
                    scale_y_continuous(
                        limits = c(0, Limit_df1),
                        breaks = Breaks_d,
                        sec.axis = sec_axis(~ . * 1, name = "Delay [min]",
                                            breaks = Breaks_d)
                    )
            }
        }
        #add the segments to the plot base and check
        #intensity is plotted independently if delay/HL data are present or not.
        if (length(na.omit(df$delay)) == 0) {
            p1 <- p1 +
                geom_point(size = .5)
            p2 <- p2
            p3 <- p3
        } else if (length(na.omit(df$delay)) < 3) {
            p1 <- p1 +
                geom_point(size = .5)
            p2 <- p2 +
                geom_point(size = .5)
            p3 <- p3 +
                geom_point(size = .5)
        } else{
            ######################intensity plot positive strand################
            #add a reference to outliers probes or bins
            df <-
                indice_function(df, "intensity_fragment")
            #intensities/HL and delay with NA are not plotted, high intensities
            #are plotted with a different shape,
            #outliers probes/bins are plotted with a different color, each
            #segment has a color, a line is added for each segment
            #indicating the mean. "FC*" indicate significant t-test and each
            #segment is labeled.
            for (o in seq_along(condition)) {
                if (nrow(df %>% filter(indice == 1 &
                                           cdt == condition[o])) != 0) {
                    p1 <- p1 +
                        geom_point(
                            data = df %>% 
                                filter(indice == 1 & cdt == condition[o]),
                            aes(fill = condition[o]),
                            size = .5
                        )
                    #eliminate outliers probes or bins
                    df1_wo <- df %>% 
                        filter(indice == 1)
                    #assign a mean position for the plot
                    df1_wo <-
                        meanPosition(df1_wo, "intensity_fragment")
                    if (length(df1_wo$intensity_fragment) != 0) {
                        p1 <- p1 +
                            geom_line(
                                data = df %>%
                                    filter(indice == 1 &
                                               cdt == condition[o]),
                                aes(
                                    x = get('position'),
                                    y = get('intensity_mean_fragment'),
                                    fill = get('intensity_fragment')
                                )
                            ) 
                    }
                }
            }
            #in case of coverage is available.
            if (coverage == 1) {
                p1 <- p1 +
                    geom_line(data = tmp.c1,
                              aes(x = get('position'), y = coverage),
                              col = col_coverage)
            }
            #######################HL plot positive strand###################
            #plot bins and mean fragments
            df <- indice_function(df, "HL_fragment")
            for (o in seq_along(condition)) {
                if (nrow(df %>% filter(indice == 1 &
                                           cdt == condition[o])) != 0) {
                    p2 <- p2 +
                        geom_point(
                            data = df %>%
                                filter(indice == 1 & cdt == condition[o]),
                            aes(fill = get('HL_fragment')),
                            size = .5
                        )
                    #eliminate outliers probes or bins
                    df1_wo <-
                        df %>% filter(indice == 1 & cdt == condition[o])
                    #assign a mean position for the plot
                    df1_wo <-
                        meanPosition(df1_wo, "HL_fragment")
                    if (length(df1_wo$HL_fragment) != 0) {
                        p2 <- p2 +
                            geom_line(
                                data = df %>%
                                    filter(indice == 1 & cdt == condition[o]),
                                aes(
                                    x = get('position'),
                                    y = get('HL_mean_fragment'),
                                    fill = get('HL_fragment')
                                )
                            ) 
                    }
                }
            }
            ########################delay plot positive strand#################
            #assign the indice to df1
            df <- indice_function(df, "delay_fragment")
            #plot delay bins
            for (o in seq_along(condition)) {
                if (nrow(df %>% 
                         filter(indice == 1 & cdt == condition[o])) != 0) {
                    p3 <- p3 +
                        geom_point(
                            data = df %>% 
                                filter(indice == 1 & cdt == condition[o]),
                            aes(y = delay,
                                fill = get('delay_fragment')),
                            size = .5
                        )
                }
                
                #plot the regression line from delay
                #plot of the delay or delay predicted is not always a line close to
                #the real values as delay data is very noisy,
                #the regression line could be far from the real data. In this cases,
                #an internal lm fit from ggplot is applied.
                #in case of slope of delay fragments is 0, the intercept is plotted
                #otherwise the delay (delay.p) is predicted from the slope and
                #plotted.
                if (length(na.omit(df$delay)) > 2) {
                    df.c <-
                        regr(
                            df %>% filter(cdt == condition[o]),
                            ind = 1,
                            data = data %>% filter(cdt == condition[o])
                        )
                    #delay_mean column is added to draw a line in case of velocity
                    #is NA or equal to 60.
                    if (nrow(df.c) != 0) {
                        if(Strand[j] == "+"){
                        p3 <- p3 +
                            geom_line(data = df.c,
                                      aes(
                                          y = get('predicted_delay'),
                                          fill = get('delay_fragment')
                                      ),
                                      size = .4)
                        }else{
                            p3 <- p3 +
                                geom_smooth(
                                    data = df.c %>%
                                        filter(get('indice') == 1),
                                    formula = y ~ x,
                                    aes(
                                        y = get('delay.p_rev'),
                                        fill = get('delay_fragment')
                                    ),
                                    se = FALSE,
                                    size = .4,
                                    method = "lm"
                                )
                        }
                      
                    }
                    
                }
            }
        }
        
        ####################pausing site positive strand##############
        #plot pausing sites events
        for(o in seq_along(condition)){
            if (nrow(df %>%
                 filter(get('pausing_site') == "+" &
                        cdt == condition[o])) != 0) {
            #select pausing site
            df1_ps <- df %>%
                filter(get('pausing_site') == "+" & cdt == condition[o])
            #select pausing site duration
            df1_ps <-
                df1_ps %>%
                filter(na.omit(get('event_duration')) >= event_duration_ps)
            if (nrow(df1_ps %>%
                     filter(get('event_ps_itss_p_value_Ttest')
                            < p_value_event)) != 0) {
                df1_ps_s <-
                    df1_ps %>%
                    filter(get('event_ps_itss_p_value_Ttest')
                           < p_value_event)
                p3 <-
                    my_segment_T(
                        p3,
                        data = df1_ps_s,
                        "PS*",
                        y = 0,
                        yend = 3,
                        dis = 50,
                        ytext = 3.8,
                        color = "orange",
                        linetype = "dashed",
                        df = "pausing",
                        fontface = fontface
                    )
            }
            if (nrow(df1_ps %>%
                     filter(get('event_ps_itss_p_value_Ttest')
                            > p_value_event)) != 0) {
                df1_ps_b <-
                    df1_ps %>%
                    filter(get('event_ps_itss_p_value_Ttest')
                           > p_value_event)
                p3 <-
                    my_segment_T(
                        p3,
                        data = df1_ps_b,
                        "PS",
                        y = 0,
                        yend = 3,
                        dis = 50,
                        ytext = 3.8,
                        color = "orange",
                        linetype = "dashed",
                        df = "pausing",
                        fontface = fontface
                    )
            }
            if (nrow(df1_ps %>%
                     filter(is.na(
                         get('event_ps_itss_p_value_Ttest')
                     ))) != 0) {
                df1_ps_s <- df1_ps %>%
                    filter(is.na('event_ps_itss_p_value_Ttest'))
                p3 <-
                    my_segment_T(
                        p3,
                        data = df1_ps_s,
                        "PS",
                        y = 0,
                        yend = 3,
                        dis = 50,
                        ytext = 3.8,
                        color = "orange",
                        linetype = "dashed",
                        df = "pausing",
                        fontface = fontface
                    )
                }
            }
        }
        ####################iTSS_I positive strand###################
        #plot internal starting sites events
        if (nrow(df %>%
                 filter(get('iTSS_I') == "+")) != 0) {
            #select iTSS
            df1_itss <- df %>%
                filter(get('iTSS_I') == "+")
            #select iTSS duration
            df1_itss <-
                df1_itss %>%
                filter(na.omit(get('event_duration'))
                       <= event_duration_itss)
            if (nrow(df1_itss %>%
                     filter(get('event_ps_itss_p_value_Ttest')
                            < p_value_event)) != 0) {
                df1_itss_s <-
                    df1_itss %>%
                    filter(get('event_ps_itss_p_value_Ttest')
                           < p_value_event)
                p3 <-
                    my_segment_NS(
                        p3,
                        data = df1_itss_s,
                        "iTSS*",
                        y = 0,
                        yend = 3.2,
                        dis = 10,
                        ytext = 3.6,
                        color = 4,
                        linetype = "dotted",
                        fontface = fontface
                    )
            }
            if (nrow(df1_itss %>%
                     filter(get('event_ps_itss_p_value_Ttest')
                            > p_value_event)) != 0) {
                df1_itss_b <-
                    df1_itss %>%
                    filter(get('event_ps_itss_p_value_Ttest')
                           > p_value_event)
                p3 <-
                    my_segment_NS(
                        p3,
                        data = df1_itss_b,
                        "iTSS",
                        y = 0,
                        yend = 3.2,
                        dis = 10,
                        ytext = 3.6,
                        color = 4,
                        linetype = "dotted",
                        fontface = fontface
                    )
            }
            if (nrow(df1_itss %>%
                     filter(is.na('event_ps_itss_p_value_Ttest'))) != 0) {
                df1_itss_b <-
                    df1_itss %>%
                    filter(is.na('event_ps_itss_p_value_Ttest'))
                p3 <-
                    my_segment_NS(
                        p3,
                        data = df1_itss_b,
                        "iTSS",
                        y = 0,
                        yend = 3.2,
                        dis = 10,
                        ytext = 3.6,
                        color = 4,
                        linetype = "dotted",
                        fontface = fontface
                    )
            }
        }
        ####################FC positive strand###################
        #add FC for intensity ratio test if p_value is significant
        #select the last row for each segment and add 40 nucleotides in case
        #of negative strand for a nice plot
        df <- indice_function(df, "intensity_fragment")
        df1_wo <- df[which(df$indice == 1), ]
        #select the last row for each segment and add 40 nucleotides in case
        #of negative strand for a nice plot
        df1_wo_pvalue <-
            arrange_byGroup(df1_wo, "p_value_intensity")
        if (nrow(df1_wo_pvalue) != 0) {
            df1_p_val_int <-
                df1_wo_pvalue[which(df1_wo_pvalue$p_value_intensity
                                    < p_value_int), ]
            if (nrow(df1_p_val_int) != 0) {
                p1 <- p1 +
                    geom_text(
                        data = df1_p_val_int,
                        aes(
                            x = get('position'),
                            y = get('intensity_mean_fragment')
                        ),
                        label = "FC*",
                        fontface = fontface,
                        size = 1,
                        check_overlap = TRUE
                    )
            }
            df1_p_val_int <-
                df1_wo_pvalue[which(df1_wo_pvalue$p_value_intensity
                                    > p_value_int), ]
            if (nrow(df1_p_val_int) != 0) {
                p1 <- p1 +
                    geom_text(
                        data = df1_p_val_int,
                        aes(
                            x = get('position'),
                            y = get('intensity_mean_fragment')
                        ),
                        label = "FC",
                        fontface = fontface,
                        size = 1,
                        check_overlap = TRUE
                    )
            }
        }
        #add FC for HL ratio lower than HL_threshold upon p_value significance
        df <- indice_function(df, "HL_fragment")
        df1_wo <- df[which(df$indice == 1), ]
        df1_hl <- arrange_byGroup(df, "FC_HL")
        df1_p_val_hl <-
            df1_hl[which(df1_hl$FC_HL >= HL_threshold_1), ]
        if (nrow(df1_p_val_hl) != 0) {
            df1_p_val_hl_p1 <-
                df1_p_val_hl[which(df1_p_val_hl$p_value_HL <= p_value_hl), ]
            if (nrow(df1_p_val_hl_p1) != 0) {
                p2 <- my_segment_NS(
                    p2,
                    data = df1_p_val_hl_p1,
                    "HL*",
                    y = 0,
                    yend = 3,
                    dis = 10,
                    ytext = 3.4,
                    color = HL_threshold_color_1,
                    linetype = "dashed",
                    fontface = fontface
                )
            }
            df1_p_val_hl_p2 <-
                df1_p_val_hl[which(df1_p_val_hl$p_value_HL > p_value_hl), ]
            if (nrow(df1_p_val_hl_p2) != 0) {
                p2 <- my_segment_NS(
                    p2,
                    data = df1_p_val_hl_p2,
                    "HL",
                    y = 0,
                    yend = 3,
                    dis = 10,
                    ytext = 3.4,
                    color = HL_threshold_color_1,
                    linetype = "dashed",
                    fontface = fontface
                )
            }
        }
        #add FC for HL ratio higher than HL_threshold upon p_value significance
        df1_p_val_hl <-
            df1_hl[which(df1_hl$FC_HL <= HL_threshold_2), ]
        if (nrow(df1_p_val_hl) != 0) {
            df1_p_val_hl_p1 <-
                df1_p_val_hl[which(df1_p_val_hl$p_value_HL <= p_value_hl), ]
            if (nrow(df1_p_val_hl_p1) != 0) {
                p2 <- my_segment_NS(
                    p2,
                    data = df1_p_val_hl_p1,
                    "HL*",
                    y = 0,
                    yend = 3,
                    dis = 10,
                    ytext = 3.4,
                    color = HL_threshold_color_1,
                    linetype = "dashed",
                    fontface = fontface
                )
            }
            df1_p_val_hl_p2 <-
                df1_p_val_hl[which(df1_p_val_hl$p_value_HL > p_value_hl), ]
            if (nrow(df1_p_val_hl_p2) != 0) {
                p2 <- my_segment_NS(
                    p2,
                    data = df1_p_val_hl_p2,
                    "HL",
                    y = 0,
                    yend = 3,
                    dis = 10,
                    ytext = 3.4,
                    color = HL_threshold_color_1,
                    linetype = "dashed",
                    fontface = fontface
                )
            }
        }
        df1_p_val_hl <-
            df1_hl[which(df1_hl$FC_HL > HL_threshold_2 &
                             df1_hl$FC_HL < HL_threshold_1), ]
        if (nrow(df1_p_val_hl) != 0) {
            df1_p_val_hl_p <-
                df1_p_val_hl[which(df1_p_val_hl$p_value_HL <= p_value_hl), ]
            if (nrow(df1_p_val_hl_p) != 0) {
                p2 <- my_segment_NS(
                    p2,
                    data = df1_p_val_hl_p,
                    "HL*",
                    y = 0,
                    yend = 3,
                    dis = 10,
                    ytext = 3.4,
                    color = HL_threshold_color_1,
                    linetype = "dashed",
                    fontface = fontface
                        )
                    }
               }
        }
            if(Strand[j] == "+"){
                p <- list(p1, p2, p3)
            }else{
                p4 <- p1 + coord_trans(y = "reverse")
                p5 <- p2 + coord_trans(y = "reverse")
                p6 <- p3 + coord_trans(y = "reverse")
                p.1 <- list(p4, p5, p6)
            }
        }
        p <- c(p, p.1)
        return(p)
    }
