#annotation_plot function plots genome annotation from a gff preprocessed file
#and the transcription units annotated from rifi framework. 
#Since the genome annotation from gff file is same for both conditions, 
#only one box is plotted for locus_tag and gene respectively but 2 boxes for TUs 
#according to the condition.
#annotation_plot plots as well events as termination and iTSS_II.
annotation_plot <-
  function(data_p,
           data_n,
           condition = condition,
           tmp.1,
           tmp.2,
           frag,
           i,
           an = an,
           region = region,
           color_region = color_region,
           fontface = fontface,
           color_text.1 = color_text.1,
           color_TU = color_TU,
           scaling_TU = scaling_TU,
           Alpha = Alpha,
           size_tu = size_tu,
           termination_threshold =
             termination_threshold,
           iTSS_threshold =
             iTSS_threshold,
           p_value_manova =
             p_value_manova,
           pos.1 = pos.1,
           pos.2 = pos.2) {
    #create an annotation data frame upon strand orientation
    an.p <- an[which(an$strand == "+"),]
    an.p <- as.character(unique(an.p$region))
    an.m <- an[which(an$strand == "-"),]
    an.m <- as.character(unique(an.m$region))
    breaks <- seq(frag[i], frag[i + 1], 1000)
    
    df.annt <- data.frame()
    segment_data_t.1 <- data.frame()
    segment_data_g.1 <- data.frame()
    segment_data_g.2 <- data.frame()
    segment_data_g.p1 <- data.frame()
    segment_data_g.p2 <- data.frame()
    #TU_annotation and gene_annot_function functions are used to annotate
    #TUs and genes on the positive strand.
    
    tryCatch({
      #only in case the dataframe from the positive strand is not empty,
      #the TU are annotated
      for (o in seq_along(condition)) {
        if (length(which(data_p$cdt == condition[o])) != 0) {
          segment_data_t.1 <-
            TU_annotation(data_p[which(data_p$cdt == condition[o]),],
                          "TU",
                          2 + scaling_TU[o],
                          5 + scaling_TU[o],
                          color_TU[o],
                          cdt = condition[o])
          segment_data_t.1 <-
            segment_data_t.1[grep("_NA", segment_data_t.1$annotation,
                                  invert = TRUE), ]
        }
        if (nrow(segment_data_t.1) == 0) {
            segment_data_t.1 <- empty_boxes(
              ystart = 2 + scaling_TU[o],
              yend = 5 + scaling_TU[o],
              frag = frag,
              i = i,
              strand = "+"
            )
          } 
        df.annt <- rbind(df.annt, segment_data_t.1)
            # Function to find TUs split into 2 pages
            tu_border <-
              function_TU_arrow(data_p[which(data_p$cdt == condition[o]),],
                                tmp.1,
                                frag = frag,
                                i = i)
            # Eliminate the TU split on the border to avoid arrow plot
            if (length(tu_border) != 0) {
              segment_data_t.1.tuLeft <-
                segment_data_t.1[-which(tu_border ==
                                          segment_data_t.1$annotation), ]
            } else {
              segment_data_t.1.tuLeft <- segment_data_t.1
        }
      }
      if (length(an.p) != 0) {
        for (k in seq_along(an.p)) {
          segment_data_g.1 <-
            gene_annot_function(
              pos.1 = pos.1,
              pos.2 = pos.2,
              yint = 13,
              yend = 16,
              annot = an,
              Strand = "+",
              reg = an.p[k],
              cdt = condition[o],
              col = color_region[which(an.p[k] == region)]
            )
          segment_data_g.2 <-
            gene_annot_function(
              pos.1 = pos.1,
              pos.2 = pos.2,
              yint = 9.2,
              yend = 12.2,
              annot = an,
              Strand = "+",
              reg = an.p[k],
              cdt = condition[o],
              col = color_region[which(an.p[k] == region)]
            )
          segment_data_g.p1 <-
            rbind(segment_data_g.p1, segment_data_g.1)
          segment_data_g.p2 <-
            rbind(segment_data_g.p2, segment_data_g.2)
        }
      }
    }, warning = function(war) {
      
    },
    error = function(err) {
    })
    
    segment_data_t.2 <- data.frame()
    segment_data_g.1 <- data.frame()
    segment_data_g.2 <- data.frame()
    segment_data_g.m1 <- data.frame()
    segment_data_g.m2 <- data.frame()
    #TU and genes annotation in the negative strand
    tryCatch({
      for (o in seq_along(condition)) {
        if (length(which(data_n$cdt == condition[o])) != 0) {
           segment_data_t.2 <-
            TU_annotation(data_n[which(data_n$cdt == condition[o]),],
                          "TU",
                          -2 - scaling_TU[o],
                          -5 - scaling_TU[o],
                          color_TU[o],
                          cdt = condition[o])
          segment_data_t.2 <-
            segment_data_t.2[grep("_NA", segment_data_t.2$annotation,
                                  invert = TRUE), ]
        }
        if (nrow(segment_data_t.2) == 0) {
            segment_data_t.2 <- empty_boxes(
                ystart = -2 + scaling_TU[o],
                yend = -5 + scaling_TU[o],
                frag = frag,
                i = i,
                strand = "-"
            )
           # segment_data_t.2 <- empty_boxes(-3, -6, frag = frag, i = i)
        } 
        df.annt <- rbind(df.annt, segment_data_t.2)
            # Function to find TUs split into 2 pages
            tu_border <-
              function_TU_arrow(data_n[which(data_n$cdt == condition[o]),],
                                tmp.2,
                                frag = frag,
                                i = i)
            # Eliminate the TU split on the border to avoid arrow plot
            if (length(tu_border) != 0) {
              segment_data_t.2.tuLeft <-
                segment_data_t.2[-which(tu_border ==
                                          segment_data_t.2$annotation), ]
            } else {
              segment_data_t.2.tuLeft <- segment_data_t.2
            }
      }
      if (length(an.m) != 0) {
        for (k in seq_along(an.m)) {
          segment_data_g.1 <-
            gene_annot_function(
              pos.1 = pos.1,
              pos.2 = pos.2,
              yint = -13,
              yend = -16,
              annot = an,
              Strand = "-",
              reg = an.m[k],
              cdt = condition[o],
              col = color_region[which(an.m[k] == region)]
            )
          segment_data_g.2 <-
            gene_annot_function(
              pos.1 = pos.1,
              pos.2 = pos.2,
              yint = -9.2,
              yend = -12.2,
              annot = an,
              Strand = "-",
              reg = an.m[k],
              cdt = condition[o],
              col = color_region[which(an.m[k] == region)]
            )
          segment_data_g.m1 <-
            rbind(segment_data_g.m1, segment_data_g.1)
          segment_data_g.m2 <-
            rbind(segment_data_g.m2, segment_data_g.2)
        }
      }
    }, warning = function(war) {
      
    },
    error = function(err) {
      
    })
    
    df.annt <-
      rbind(
        df.annt,
        segment_data_g.p1,
        segment_data_g.p2,
        segment_data_g.m1,
        segment_data_g.m2
      )
    
    ###################################plot#############################
    if (nrow(df.annt) != 0) {
      p7 <- ggplot(data = df.annt) +
        scale_x_continuous(limits = c(frag[i], frag[i + 1]),
                           breaks = breaks) +
        scale_y_continuous(limits = c(-16, 16), expand = c(0, 0)) +
        geom_hline(yintercept = 0) +
        theme_minimal(base_size = 11) +
        theme(
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none",
          axis.ticks.y =  element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_blank(),
          panel.border = element_blank()
        )
      
      #################add arrows to TUs, both strands separately##########
      if (length(which(df.annt$region == "r")) != 0) {
      p7 <- p7 +
      annotate(
        "segment",
        x = frag[i + 1],
        xend = frag[i],
        y = unique((
          df.annt[which(df.annt$region == "r" &
                          df.annt$strand == "+"), "ystart"] +
            df.annt[which(df.annt$region == "r" &
                            df.annt$strand == "+"), "yend"]
        ) / 2),
        yend = unique((
          df.annt[which(df.annt$region == "r" &
                          df.annt$strand == "+"), "ystart"] +
            df.annt[which(df.annt$region == "r" &
                            df.annt$strand == "+"), "yend"]
        ) / 2),
        size = .1,
        arrow = my_arrow(1, "open")
      )+
      annotate(
        "segment",
        xend = frag[i + 1],
        x = frag[i],
        y = unique((
          df.annt[which(df.annt$region == "r" &
                          df.annt$strand == "-"), "ystart"] +
            df.annt[which(df.annt$region == "r" &
                            df.annt$strand == "-"), "yend"]
        ) / 2),
        yend = unique((
          df.annt[which(df.annt$region == "r" &
                          df.annt$strand == "-"), "ystart"] +
            df.annt[which(df.annt$region == "r" &
                            df.annt$strand == "-"), "yend"]
        ) / 2),
        size = .1,
        arrow = my_arrow(1, "open")
      )
      }
      ###############################################################
       #add rectangle for each box
      if (nrow(df.annt) != 0) {
        p7 <- p7 +
          geom_rect(
            data = df.annt,
            aes(
              xmin = get('xstart'),
              ymin = get('ystart'),
              xmax = get('xend'),
              ymax = get('yend')
            ),
            alpha = Alpha,
            fill = df.annt$color
          ) +
          geom_text(
            data = df.annt[grep("UTR", df.annt$region, invert = TRUE), ],
            aes((get('xstart') + get('xend')) / 2, (get('ystart') +
                                                      get('yend')) / 2,
                label = get('annotation')),
            size = size_tu,
            fontface = fontface,
            color = color_text.1,
            check_overlap = TRUE
          )
      }
   
      if (nrow(df.annt %>% filter(region == "r")) != 0) {
        ####################### add small ticks to TUs positive strand#########
     if(nrow(df.annt[which(df.annt$region == "r" &
                      df.annt$strand == "+"),]) != 0){
      for(o in seq_along(condition)){
        df.annt_p <- df.annt[which(df.annt$region == "r" &
                        df.annt$strand == "+"),]
        p7 <- p7 +
          annotate(
            "segment",
            x = unique(df.annt_p[which(df.annt_p$cdt == condition[o]),
                               "xstart"]) - 10,
            xend = unique(df.annt_p[which(df.annt_p$cdt == condition[o]),
                                  "xstart"]) - 10,
            y = unique(df.annt_p[which(df.annt_p$cdt == condition[o]),
                               "ystart"]),
            yend = unique(df.annt_p[which(df.annt_p$cdt == condition[o]),
                                  "yend"]) + .3,
            size = .2
          ) +
          annotate(
            "segment",
            xend = unique(df.annt_p[which(df.annt_p$cdt == condition[o]),
                                  "xstart"]) - 5,
            x = unique(df.annt_p[which(df.annt_p$cdt == condition[o]),
                               "xstart"]) + 60,
            y = unique(df.annt_p[which(df.annt_p$cdt == condition[o]),
                               "yend"]) + .3,
            yend = unique(df.annt_p[which(df.annt_p$cdt == condition[o]),
                                  "yend"]) + .3,
            size = .2,
            arrow = my_arrow(.4, "open")
          )
      }
     }
        ####################### add small ticks to TUs negative strand#########
     if(nrow(df.annt[which(df.annt$region == "r" &
                              df.annt$strand == "-"),]) != 0){ 
      for(o in seq_along(condition)){
          df.annt_n <- df.annt[which(df.annt$region == "r" &
                                       df.annt$strand == "-"),]
          p7 <- p7 +
            annotate(
              "segment",
              x = unique(df.annt_n[which(df.annt_n$cdt == condition[o]),
                                   "xend"]) + 10,
              xend = unique(df.annt_n[which(df.annt_n$cdt == condition[o]),
                                      "xend"]) + 10,
              y = unique(df.annt_n[which(df.annt_n$cdt == condition[o]),
                                   "ystart"]),
              yend = unique(df.annt_n[which(df.annt_n$cdt == condition[o]),
                                      "yend"]) - .3,
              size = .2
            ) + 
            annotate(
              "segment",
              xend = unique(df.annt_n[which(df.annt_n$cdt == condition[o]),
                                      "xend"]) + 5,
              x = unique(df.annt_n[which(df.annt_n$cdt == condition[o]),
                                   "xend"]) - 60,
              y = unique(df.annt_n[which(df.annt_n$cdt == condition[o]),
                                   "yend"]) - .3,
              yend = unique(df.annt_n[which(df.annt_n$cdt == condition[o]),
                                      "yend"]) - .3,
              size = .2,
              arrow = my_arrow(.4, "open")
            )
        }  
      }
    }
        
      ############################events plot##############################
        df1_syR_T <- NA
      for(o in seq_along(condition)){
        if (length(which(!is.na(data_p[which(data_p$cdt == condition[o]),
                                       "synthesis_ratio"]))) != 0) {
          df1_syR <- data_p %>%
            filter(cdt == condition[o] & !is.na(synthesis_ratio))
          #in case last position matches with an event which needs to be
          #on the next page of the plot.
          if (last(df1_syR$position) == frag[c(i + 1)]) {
            fc_seg <-
              df1_syR[which(df1_syR$position ==
                              frag[c(i + 1)]), "FC_HL_intensity_fragment"]
            df1_syR[which(df1_syR$FC_HL_intensity_fragment == fc_seg),
                    c("synthesis_ratio_event",
                      "p_value_Manova")] <- NA
          }
          if (length(which(
            df1_syR$synthesis_ratio < termination_threshold &
            !is.na(df1_syR$FC_HL_intensity_fragment)
          )) != 0) {
            df1_syR_T <-
              df1_syR[which(
                df1_syR$synthesis_ratio < termination_threshold &
                  !is.na(df1_syR$FC_HL_intensity_fragment)
              ), ]
            df1_syR_T <-
              arrange_byGroup(df1_syR_T, "FC_HL_intensity_fragment")
            if (length(which(df1_syR_T$p_value_Manova < p_value_manova)) != 0) {
              df1_syR_T.m <-
                df1_syR_T %>%
                filter(get('p_value_Manova') < p_value_manova)
              p7 <-
                my_segment_T(
                  p7,
                  data = df1_syR_T.m,
                  "Ter*",
                  y = 2 + scaling_TU[o],
                  yend = 5 +  scaling_TU[o],
                  dis = 50,
                  ytext = 2.1,
                  color = 2,
                  linetype = "solid",
                  df = "termination",
                  fontface = fontface
                )
            }
            if (length(which(df1_syR_T$p_value_Manova > p_value_manova)) != 0) {
              df1_syR_T.t <-
                df1_syR_T %>%
                filter(get('p_value_Manova') > p_value_manova)
              p7 <-
                my_segment_T(
                  p7,
                  data = df1_syR_T.t,
                  "Ter",
                  y = 2 + scaling_TU[o],
                  yend = 5 +  scaling_TU[o],
                  dis = 50,
                  ytext = 2.1,
                  color = 2,
                  linetype = "solid",
                  df = "termination",
                  fontface = fontface
                )
            }
            if (length(which(is.na(df1_syR_T$p_value_Manova))) != 0) {
              df1_syR_T.t <- df1_syR_T %>%
                filter(is.na(get('p_value_Manova')))
              p7 <-
                my_segment_T(
                  p7,
                  data = df1_syR_T.t,
                  "Ter",
                  y = 2 + scaling_TU[o],
                  yend = 5 +  scaling_TU[o],
                  dis = 50,
                  ytext = 2.1,
                  color = 2,
                  linetype = "solid",
                  df = "termination",
                  fontface = fontface
                )
            }
          }
          df1_syR_T <- NA
          #plot New_start event from synthesis_ratio_event
          if (length(which(
            df1_syR$synthesis_ratio > iTSS_threshold &
            !is.na(df1_syR$FC_HL_intensity_fragment)
          )) != 0) {
            df1_syR_T <-
              df1_syR[which(
                df1_syR$synthesis_ratio > iTSS_threshold &
                  !is.na(df1_syR$FC_HL_intensity_fragment)
              ), ]
            df1_syR_T <-
              arrange_byGroup(df1_syR_T, "FC_HL_intensity_fragment")
            if (length(which(df1_syR_T$p_value_Manova < p_value_manova)) != 0) {
              df1_syR_T.m <-
                df1_syR_T %>%
                filter(get('p_value_Manova') < p_value_manova)
              p7 <-
                my_segment_NS(
                  p7,
                  data = df1_syR_T.m,
                  "NS*",
                  y = 2 + scaling_TU[o],
                  yend = 5 +  scaling_TU[o],
                  dis = 10,
                  ytext = 2.1,
                  color = "coral4",
                  linetype = "solid",
                  fontface = fontface
                )
            }
            if (length(which(df1_syR_T$p_value_Manova > p_value_manova)) != 0) {
              df1_syR_T.t <-
                df1_syR_T %>%
                filter(get('p_value_Manova') > p_value_manova)
              p7 <-
                my_segment_NS(
                  p7,
                  data = df1_syR_T.t,
                  "NS",
                  y = 2 + scaling_TU[o],
                  yend = 5 + scaling_TU[o],
                  dis = 10,
                  ytext = 2.1,
                  color = "coral4",
                  linetype = "solid",
                  fontface = fontface
                )
            }
            if (length(which(is.na(df1_syR_T$p_value_Manova))) != 0) {
              df1_syR_T.t <- df1_syR_T %>%
                filter(is.na(get('p_value_Manova')))
              p7 <-
                my_segment_NS(
                  p7,
                  data = df1_syR_T.t,
                  "NS",
                  y = 2 + scaling_TU[o],
                  yend = 5 + scaling_TU[o],
                  dis = 10,
                  ytext = 2.1,
                  color = "coral4",
                  linetype = "solid",
                  fontface = fontface
                )
              }
            }
         }
      }
    }
    
    if (nrow(df.annt %>% filter(region == "r" & strand == "-")) != 0) {
         df2_syR_T <- NA
         for(o in seq_along(condition)){
          if (length(which(!is.na(data_n[which(data_n$cdt == condition[o]),
                                         "synthesis_ratio"]))) != 0) {
            df2_syR <- data_n %>%
              filter(cdt == condition[o] & !is.na(synthesis_ratio))
          if (length(which(
            df2_syR$synthesis_ratio < termination_threshold &
            !is.na(df2_syR$FC_HL_intensity_fragment)
          )) != 0) {
            df2_syR_T <-
              df2_syR[which(
                df2_syR$synthesis_ratio < termination_threshold &
                  !is.na(df2_syR$FC_HL_intensity_fragment)
              ),]
            df2_syR_T <-
              df2_syR_T[!duplicated(df2_syR_T$FC_fragment_intensity),]
            if (length(which(df2_syR_T$p_value_Manova < p_value_manova)) != 0) {
              df2_syR_T.m <-
                df2_syR_T %>%
                filter(get('p_value_Manova') < p_value_manova)
              p7 <-
                my_segment_T(
                  p7,
                  data = df2_syR_T.m,
                  "Ter*",
                  y = -2 - scaling_TU[o] ,
                  yend = -5 - scaling_TU[o],
                  dis = 50,
                  ytext = -2.1,
                  color = 2,
                  linetype = "solid",
                  df = "termination",
                  fontface = fontface
                )
            }
            if (length(which(df2_syR_T$p_value_Manova > p_value_manova)) != 0) {
              df2_syR_T.t <-
                df2_syR_T %>%
                filter(get('p_value_Manova') > p_value_manova)
              p7 <-
                my_segment_T(
                  p7,
                  data = df2_syR_T.t,
                  "Ter",
                  y = -2 - scaling_TU[o] ,
                  yend = -5 - scaling_TU[o],
                  dis = 50,
                  ytext = -2.1,
                  color = 2,
                  linetype = "solid",
                  df = "termination",
                  fontface = fontface
                )
            }
            if (length(which(is.na(df2_syR_T$p_value_Manova))) != 0) {
              df2_syR_T.t <- df2_syR_T %>%
                filter(is.na(get('p_value_Manova')))
              p7 <-
                my_segment_T(
                  p7,
                  data = df2_syR_T.t,
                  "Ter",
                  y = -2 - scaling_TU[o] ,
                  yend = -5 - scaling_TU[o],
                  dis = 50,
                  ytext = -2.1,
                  color = 2,
                  linetype = "solid",
                  df = "termination",
                  fontface = fontface
                )
            }
          }
          df2_syR_T <- NA
          #plot New_start event from synthesis_ratio_event column
          if (length(which(
            df2_syR$synthesis_ratio > iTSS_threshold &
            !is.na(df2_syR$FC_HL_intensity_fragment)
          )) != 0) {
            df2_syR_T <-
              df2_syR[which(
                df2_syR$synthesis_ratio > iTSS_threshold &
                  !is.na(df2_syR$FC_HL_intensity_fragment)
              ),]
            df2_syR_T <-
              df2_syR_T[!duplicated(df2_syR_T$FC_fragment_intensity),]
            if (length(which(df2_syR_T$p_value_Manova < p_value_manova)) != 0) {
              df2_syR_T.m <-
                df2_syR_T %>%
                filter(get('p_value_Manova') < p_value_manova)
              p7 <-
                my_segment_NS(
                  p7,
                  data = df2_syR_T.m,
                  "NS*",
                  y = -2 - scaling_TU[o] ,
                  yend = -5 - scaling_TU[o],
                  dis = 10,
                  ytext = -2.1,
                  color = "coral4",
                  linetype = "solid",
                  fontface = fontface
                )
            }
            if (length(which(df2_syR_T$p_value_Manova > p_value_manova)) != 0) {
              df2_syR_T.t <-
                df2_syR_T %>%
                filter(get('p_value_Manova') > p_value_manova)
              p7 <-
                my_segment_NS(
                  p7,
                  data = df2_syR_T.t,
                  "NS",
                  y = -2 - scaling_TU[o] ,
                  yend = -5 - scaling_TU[o],
                  dis = 10,
                  ytext = -2.1,
                  color = "coral4",
                  linetype = "solid",
                  fontface = fontface
                )
            }
            if (length(which(is.na(df2_syR_T$p_value_Manova))) != 0) {
              df2_syR_T.t <- df2_syR_T %>%
                filter(is.na(get('p_value_Manova')))
               p7 <-
                my_segment_NS(
                  p7,
                  data = df2_syR_T.t,
                  "NS",
                  y = -2 - scaling_TU[o] ,
                  yend = -5 - scaling_TU[o],
                  dis = 10,
                  ytext = -2.1,
                  color = "coral4",
                  linetype = "solid",
                  fontface = fontface
                )
              }
            }
          }
        }
      }

    return(p7)
  }
