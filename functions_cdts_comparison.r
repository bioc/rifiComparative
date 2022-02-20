function_df <- function(indice,
                        dt3,
                        input,
                        an,
                        i,
                        data,
                        locus_t,
                        conditions) {
  for (o in seq_along(indice)) {
    ind = indice[[o]]
    data3 = dt3[which(dt3$cdt == conditions[[o]]), ]
    inp = input[which(input$cdt == conditions[[o]]), ]
    
    #################################delay##############################
    # grep rows with position flanking the gene positions from
    # data combined df
    inp.1 <- inp[which(inp$strand %in% an$strand &
                         between(inp$position, an$start, an$end)),]
    #grep delay fragments
    d <-
      unique(inp.1$delay_fragment[grep(paste0("\\D_", "\\d+", "$"),
                                       inp.1$delay_fragment)])
    
    #grep delay fragments for those fragments with only delay outliers
    if (length(d) == 0) {
      d <- unique(inp.1$delay_fragment[grep(paste0("\\D_", "\\d+", "_0"),
                                            inp.1$delay_fragment)])
    }
    data[i, paste0("delay_frg_", ind)] <- as.numeric(length(d))
    ##################################HL################################
    #grep HL fragments except outliers
    dc <-
      unique(inp.1$HL_fragment[grep(paste0("\\Dc_", "\\d+", "$"),
                                    inp.1$HL_fragment)])
    # Outliers are catched for those cases such sRNA
    if (length(dc) == 0) {
      dc <- unique(inp.1$HL_fragment[grep(paste0("\\Dc_", "\\d+", "_0"),
                                          inp.1$HL_fragment)])
    }
    if (length(dc) != 0) {
      data[i, paste0("HL_frg_", ind)] <-
        length(dc)
      #mean of half-life gene based without outliers
      data[i, paste0("HL_mean_", ind)] <-
        round(as.numeric(mean(inp.1[grep(paste(paste0("\\Dc_", "\\d+", "$"),
                                               collapse = "|"),
                                         inp.1$HL_fragment),
                                    "half_life"])), 2)
      
      HL_t_test <-
        round(as.numeric(inp.1[grep(paste(paste0("\\Dc_", "\\d+", "$"),
                                          collapse = "|"),
                                    inp.1$HL_fragment),
                               "half_life"]), 2)
      #assemble HL data for t-test
      data[i, paste0("HL_P_value_", ind)] <-
        paste0(HL_t_test, collapse = ",")
      
      #assemble the mean of half_life of the fragments
      data[i, paste0("HL_", ind)] <-
        paste0(round(unique(inp.1[grep(paste(paste0("\\Dc_", "\\d+", "$"),
                                             collapse = "|"),
                                       inp.1$HL_fragment),
                                  "HL_mean_fragment"]), 2),
               collapse = "|")
    }
    ############################intensity################################
    #number of intensity fragment in sc
    I <-
      unique(inp.1$intensity_fragment[grep(paste0("\\I_", "\\d+", "$"),
                                           inp.1$intensity_fragment)])
    if (length(I) == 0) {
      I <-
        unique(inp.1$intensity_fragment[grep(paste0("\\I_", "\\d+", "_0"),
                                             inp.1$intensity_fragment)])
    }
    
    if (length(I) != 0) {
      data[i, paste0("int_frg_", ind)] <-
        length(I)
      data[i, paste0("int_mean_", ind)] <-
        round(as.numeric(mean(inp.1[grep(paste(paste0("\\I_", "\\d+", "$"),
                                               collapse = "|"),
                                         inp.1$intensity_fragment),
                                    "intensity"])), 0)
      
      int_t_test <-
        round(as.numeric(inp.1[grep(paste(paste0("\\I_", "\\d+", "$"),
                                          collapse = "|"),
                                    inp.1$intensity_fragment),
                               "intensity"]), 2)
      #assemble HL data for t-test
      data[i, paste0("int_P_value_", ind)] <-
        paste0(int_t_test, collapse = ",")
      
      data[i, paste0("int_", ind)] <-
        paste0(round(unique(inp.1[grep(paste(paste0("\\I_", "\\d+", "$"),
                                             collapse = "|"),
                                       inp.1$intensity_fragment),
                                  "intensity_mean_fragment"]), 0),
               collapse = "|")
      
    }
    ##############################events#################################
    tryCatch({
      ev <-
        data3[grep(paste0("^", locus_t, "$"), data3$locus_tag),]
      data[i, paste0("paus_", ind)] <-
        ifelse(nrow(ev[grep("ps", ev$event),]) != 0,
               paste0(floor(ev[grep("ps", ev$event),
                               "event_position"]), collapse = "|"), 0)
      data[i, paste0("ter_", ind)] <-
        ifelse(nrow(ev[grep("Termination", ev$event),]) != 0,
               paste0(floor(ev[grep("Termination", ev$event),
                               "event_position"]), collapse = "|"), 0)
      data[i, paste0("iTSS_I_", ind)] <-
        ifelse(nrow(ev[grep("iTSS_I", ev$event), ]) != 0,
               paste0(floor(ev[grep("iTSS_I", ev$event),
                               "event_position"]), collapse = "|"), 0)
      data[i, paste0("iTSS_II_", ind)] <-
        ifelse(nrow(ev[grep("iTSS_II", ev$event), ]) != 0,
               paste0(floor(ev[grep("iTSS_II", ev$event),
                               "event_position"]), collapse = "|"), 0)
    },
    error = function(e) {
      
    })
  }
  return(data)
}

function_tTest <- function(data, parameter, indice)
{
  data[, paste0("P.Value_", parameter)] <- NA
  data[, paste0("adj.P.Val_", parameter)] <- NA
  for (p in seq_len(nrow(data))) {
    cd.1 <- as.numeric(unlist(str_split(data[p, paste0(parameter, "_P_value_", indice[1])], ",")))
    cd.2 <- as.numeric(unlist(str_split(data[p, paste0(parameter, "_P_value_", indice[2])], ",")))
    if (length(cd.1) < 3 || length(cd.2) < 3) {
      data[p, paste0("P.Value_", parameter)] <- NA
    } else{
      data[p, paste0("P.Value_", parameter)] <- t.test(cd.1, cd.2)[[3]]
    }
  }
  data[, paste0("adj.P.Val_", parameter)]  <- p.adjust(data[, paste0("P.Value_", parameter)], method = "fdr")
  return(data)
}


function_manova <- function(data, parameter, indice)
{
  data[, paste0("P.Value_", parameter)] <- NA
  data[, paste0("adj.P.Val_", parameter)] <- NA
  for (p in seq_len(nrow(data))) {
    cd.1 <- as.numeric(unlist(str_split(data[p, paste0(parameter, "_P_value_", indice[1])], ",")))
    cd.2 <- as.numeric(unlist(str_split(data[p, paste0(parameter, "_P_value_", indice[2])], ",")))
    if (length(cd.1) < 3 || length(cd.2) < 3) {
      data[p, paste0("P.Value_", parameter)] <- NA
    } else{
      data[p, paste0("P.Value_", parameter)] <- t.test(cd.1, cd.2)[[3]]
    }
  }
  data[, paste0("adj.P.Val_", parameter)]  <- p.adjust(data[, paste0("P.Value_", parameter)], method = "fdr")
  return(data)
}