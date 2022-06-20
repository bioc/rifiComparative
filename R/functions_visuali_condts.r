#strand_selection used to select the strand
strand_selection <- function(data, Strand) {
  data %>%
    filter(get('strand') == Strand) %>%
    arrange(get('position'))
}

#splitGenome_function used to split the genome into same fragments
splitGenome_function <- function(x, gLength) {
  fl_v <- c()
  for (i in seq_along(unique(x))) {
    fl_v <- c(fl_v, x[which(x == i)][1])
  }
  fl_v <- na.omit(fl_v)
  fl_v <- fl_v[c(1, fl_v)]
  names(fl_v)[1] <- 1
  fl_v[1] <- 0
  frag <- as.numeric(as.character(names(fl_v)))
  frag <- c(frag, last(gLength))
  frag <- frag[!duplicated(frag)]
  return(frag)
}

#indice_function assign to distinguish between outlier, terminal or bin/probe belonging to the segment
indice_function <- function(data, col) {
  data[, "indice"] <- NA
  pos <- grep(paste0("\\D_\\d+", "$"), data[, col])
  pos_0 <- grep(paste0("\\D_\\d+\\_O", "$"), data[, col])
  pos_1 <- grep(paste0("\\D_\\d+\\_T", "$"), data[, col])
  data[pos, "indice"] <- 1
  data[pos_0, "indice"] <- 2
  data[pos_1, "indice"] <- 3
  return(data)
}

#TU_annotation plots TUs
TU_annotation <- function(data, annot, ystart, yend, color_tu, cdt) {
  pos.1_v <- c()
  pos.2_v <- c()
  annot_tmp <- (unique(data[, annot]))
  for (j in seq_along(annot_tmp)) {
    if (length(annot_tmp) == 0) {
      break ()
    }
    pos <- which(data[, annot] == annot_tmp[j])
    pos.1 <- data$position[pos[1]]
    pos.1_v <- c(pos.1_v, pos.1)
    pos.2 <- data$position[last(pos)]
    pos.2_v <- c(pos.2_v, pos.2)
  }
  segment_data <- data.frame(
    xstart = pos.1_v,
    xend = pos.2_v,
    ystart = rep(ystart, times = length(pos.1_v)),
    yend = yend,
    region = rep("r", times = length(pos.1_v)),
    color = color_tu,
    annotation = annot_tmp,
    gene = rep("g", times = length(pos.1_v)),
    cdt = cdt,
    stringsAsFactors = FALSE
  )
  if (unique(data$strand) == "-") {
    segment_data$strand <- "-"
  } else{
    segment_data$strand <- "+"
  }
  return(na.omit(segment_data))
}

#gene_annot_function plots genes or locus_tag
gene_annot_function <-
  function(pos.1,
           pos.2,
           annot,
           Strand,
           yint,
           yend,
           reg,
           col,
           cdt) {
    annot_tmp <- annot[which(annot[, "strand"] %in% Strand), ]
    annot_tmp <-
      annot_tmp[order(annot_tmp[, "start"], decreasing = FALSE), ]
    if (Strand == "-") {
      annot_tmp <- annot_tmp[between(annot_tmp[, "start"], pos.1, pos.2), ]
    } else{
      annot_tmp <- annot_tmp[between(annot_tmp[, "end"], pos.1, pos.2), ]
    }
    if (!is_empty(reg)) {
      annot_tmp <- annot_tmp[which(annot_tmp[, "region"] %in% reg), ]
      if (!is_empty(annot_tmp[, "region"])) {
        pos_v <- c()
        color_v <- c()
        for (j in seq_along(annot_tmp[, "region"])) {
          pos <- which(annot_tmp[, "region"][j] == reg)
          region <- reg[pos]
          color <- col[pos]
          pos_v <- c(pos_v, pos)
          color_v <- c(color_v, color)
        }
        segment_data <- data.frame(
          xstart = annot_tmp$start,
          xend = annot_tmp$end,
          ystart = rep(yint, times = length(color_v)),
          yend = rep(yend, times = length(color_v)),
          region = annot_tmp$region,
          color = color_v,
          annotation = annot_tmp$locus_tag,
          gene = annot_tmp$gene,
          cdt = cdt,
          strand = Strand,
          stringsAsFactors = FALSE
        )
      }
    } else{
      segment_data <- data.frame(
        xstart = annot_tmp$start,
        xend = annot_tmp$end,
        ystart = yint,
        yend = yend,
        region = rep(NA, times = length(annot_tmp$start)),
        color = rep(NA, times = length(annot_tmp$start)),
        annotation = annot_tmp$locus_tag,
        gene = annot_tmp$gene,
        strand = Strand,
        stringsAsFactors = FALSE
      )
    }
    segment_data[which(is.na(segment_data$gene)), "gene"] <-
      segment_data[which(is.na(segment_data$gene)), "annotation"]
    segment_data$xend <- ifelse(segment_data$xend > pos.2 - 2000, 
                                pos.2 - 2000, segment_data$xend)
    segment_data$xstart <- ifelse(segment_data$xstart < pos.1 + 1999, 
                                pos.1 + 2000, segment_data$xstart)
    return(segment_data)
  }

#label_log2_function converts intensity scale to exponential 2
label_log2_function <-
  function(x)
    parse(text = paste0('2^', log(x, 2)))

label_square_function <- function(x)
  round(sqrt(x), 0)

#secondaryAxis designs the secondray axis
secondaryAxis <-  function(data, parameter, ind) {
  data <- data %>%
    filter(get('indice') == ind)
  if (nrow(data) > 1 & length(data[, parameter] > 20) >= 1) {
    data[which(data[, parameter] > 20), parameter] <- 20
  }
  return(data)
}

#arrange_byGroup used to slice the rows by parameter
arrange_byGroup <- function(input, parameter) {
  slice <- slice
  if (unique(input$strand) == "+") {
    input <- as.data.frame(input %>%
                          group_by(input[, parameter]) %>%
                          arrange(get('position')) %>%
                            dplyr::slice(n()))
  }else{
    input <- as.data.frame(input %>%
                            group_by(input[, parameter]) %>%
                            arrange(get('position')) %>%
                             dplyr::slice(n()))
    input[which(input$pausing_site == "+"), "position"] <-
      input[which(input$pausing_site == "+"), "position"] + 40
    input[which(input$iTSS_I == "+"), "position"] <-
      input[which(input$iTSS_I == "+"), "position"] + 40
  }
  return(input)
}

#regr fits the delay in case of negative strand. regr uses the coefficients extracted from slop and intercept in case of positive strand.

regr <- function(input, ind, data) {
  predicted_delay <- NA
  df_bind <- data.frame()
  input <- indice_function(input, "delay_fragment")
  data <- indice_function(data, "delay_fragment")
  #outliers and slope equal to 0 are not selected from the short data and
  #the original data
  df <-
    input %>%
    filter(get('indice') == ind) %>%
    filter(get('slope') != 0)
  df.1 <-
    data %>%
    filter(get('indice') == ind) %>%
    filter(get('slope') != 0)
  #slope equal to 0 are not further taken in consideration
  if (nrow(df) != 0) {
    if (unique(input$strand) == "+") {
      df$predicted_delay <-
        df[, "position"] * df[, "slope"] + df[, "intercept"]
    } else{
      #in case df has less than 3 rows, no regression is estimated.
      #df is returned empty
      tryCatch({
      if(nrow(df) <= 2){
        df <- data.frame()
      }else{
      #positions are adjusted to 0 or subtracted from genome length, the fit
      #is the same.
      df$position_adjusted <- abs(max(df$position) - df$position) + 1
      #loop into unique delay fragments, extract and store the coefficients
      delayF <- unique(df$delay_fragment)
      for (i in seq_along(delayF)) {
        df.f <- df[which(df$delay_fragment %in% delayF[i]), ]
        if (nrow(df.f) <= 2) {
          next ()
        }
        df.o <- df.1[which(df.1$delay_fragment %in% delayF[i]), ]
        coef.lm <- lm(rev(delay) ~ position_adjusted, df.f)
        df.f$interp_adjusted <- coef.lm[[1]][[1]]
        slope_adjusted <- coef.lm[[1]][[2]]
        df.f$slope_adjusted <- slope_adjusted
        df.f$delay.p <- df.f$slope_adjusted *
          df.f$position_adjusted + df.f$interp_adjusted
        if (nrow(df.o) != nrow(df.f) & slope_adjusted > 0) {
          #check the fragment on the original data
          df.o <- df.1[which(df.1$delay_fragment %in% delayF[i]), ]
          df.o$position_adjusted <-
            abs(max(df.o$position) - df.o$position) + 1
          coef.lm <- lm(rev(delay) ~ position_adjusted, df.o)
          df.o$interp_adjusted <- coef.lm[[1]][[1]]
          df.o$slope_adjusted <- coef.lm[[1]][[2]]
          df.o$delay.p <- df.o$slope_adjusted *
            df.o$position_adjusted + df.o$interp_adjusted
          df.o <- df.o[which(df.o$position %in% df$position), ]
          df.f <- df.o
        }
        df_bind <- rbind(df_bind, df.f)
          }
        }
      }, error=function(e){}
      )
      #revert the delay predicted to adjust it to the positions of the dataframe
      if(nrow(df_bind) != 0){
      df_bind <- as.data.frame(df_bind %>%
                                 group_by(get('delay_fragment')) %>%
                                 mutate(delay.p_rev = rev(get('delay.p'))))
      df <- df_bind
      }else{
        df <- data.frame()
      }
    }
  }
  return(df)
}

#add an arrow to TU 
my_arrow <- function(Unit, type) {
  return(arrow(
    type = type,
    ends = "first",
    length = unit(Unit, "mm")
  ))
}

#delay_mean adds a mean of all delay bins in case no fit is possible
delay_mean <- function(data, parameter, ind) {
  df <- data %>%
    filter(get('indice') == ind)
  df <-
    as.data.frame(tapply(data[, parameter], data$delay_fragment, mean))
  df <- cbind.data.frame(rownames(df), df[, 1])
  colnames(df) <- c("delay_fragment", "mean_delay")
  newDF <- left_join(data, df, by = "delay_fragment")
  return(newDF$mean_delay)
}

#my_segment_T plots terminal events
my_segment_T <-
  function(p,
           data,
           label,
           y,
           yend,
           dis,
           ytext,
           color,
           linetype,
           df,
           fontface) {
    p <- p +
      annotate(
        "segment",
        x = unique(data$position),
        xend = unique(data$position),
        y = y,
        yend = yend,
        size = .4,
        color = color,
        linetype = linetype
      ) +
      annotate(
        "segment",
        x = unique(data$position) + dis,
        xend = unique(data$position) - dis,
        y = yend,
        yend = yend,
        size = .2,
        color = color,
        lineend = "butt"
      ) +
      geom_text(
        data = data,
        aes(x = get('position'), y = ytext),
        label = label,
        fontface = fontface,
        color = color,
        size = 1.2
      )
    if (df == "pausing") {
      p <- p +
        geom_segment(
          data = data,
          xend = unique(data$position) - 15,
          x = unique(data$position) - 15,
          y = yend,
          yend = yend + .5,
          size = .1,
          color = color,
          lineend = "butt",
          inherit.aes = TRUE
        ) +
        geom_segment(
          data = data,
          xend = unique(data$position) + 15,
          x = unique(data$position) + 15,
          y = yend,
          yend = yend + .5,
          size = .1,
          color = color,
          lineend = "butt",
          inherit.aes = TRUE
        )
    }
    return(p)
  }

#my_segment_NS plots iTSS_II events
my_segment_NS <-
  function(p,
           data,
           label,
           y,
           yend,
           dis,
           ytext,
           color,
           linetype,
           fontface) {
    if (unique(data$strand) == "+") {
      p <- p +
        annotate(
          "segment",
          x = unique(data$position) - dis,
          xend =
            unique(data$position) - dis,
          y = y,
          yend = yend,
          size = .4,
          color = color,
          linetype = linetype
        ) +
        annotate(
          "segment",
          xend = unique(data$position),
          x = unique(data$position) + 60,
          y = yend,
          yend = yend,
          size = .2,
          arrow = my_arrow(.4, "open")
        ) +
        geom_text(
          data = data,
          aes(x = get('position'), y = ytext),
          label = label,
          fontface = fontface,
          color = color,
          size = 1
        )
    } else if (unique(data$strand) == "-") {
      p <- p +
        annotate(
          "segment",
          x = unique(data$position) + dis,
          xend = unique(data$position) + dis,
          y = y,
          yend = yend,
          size = .2,
          color = color,
          linetype = linetype
        ) +
        annotate(
          "segment",
          xend = unique(data$position) + 5,
          x = unique(data$position) - 60,
          y = yend,
          yend = yend,
          size = .2,
          arrow = my_arrow(.4, "open")
        ) +
        geom_text(
          data = data,
          aes(x = get('position'), y = ytext),
          label = label,
          fontface = fontface,
          color = color,
          size = 1
        )
    }
    return(p)
  }

#limit_function limits the scaling in delay and HL
limit_function <- function(data, parameter, ind) {
  data <- data %>%
    filter(get('indice') == ind)
  if (length(which(data[, parameter] > 10)) >= 3 &
      length(which(data[, parameter] <= 20)) >= 3) {
    Limit <- 20
  } else if (nrow(data) < 3 &
             length(which(data[, parameter] > 10)) >= 1) {
    Limit <- 20
  } else{
    Limit <- 10
  }
  return(Limit)
}

#empty_boxes plots empty boxes in case no genes or TUs are available
empty_boxes <- function(ystart,
                        yend,
                        frag = frag,
                        i = i) {
  data <- data.frame(
    xstart = frag[i],
    xend = frag[i + 1],
    ystart = ystart,
    yend = yend,
    region = NA,
    color = adjustcolor("white", alpha.f = 0.2),
    annotation = "",
    gene = NA,
    strand = NA
  )
  return(data)
}

# function_TU_arrow finds TUs split into 2 pages
function_TU_arrow <- function(data,
                              big_data,
                              frag = frag,
                              i = i) {
  tu <- unique(data$TU)
  tu <- tu[grep("_NA", tu, invert = TRUE)]
  tu_v <- c()
  for (j in seq_along(tu)) {
    pos <- big_data[which(big_data$TU == tu[j]), "position"]
    if (unique(data$strand == "+")) {
      if (j == 1 & pos[1] < frag[i]) {
        tu_v <- c(tu_v, tu[j])
      }
    } else{
      if (last(pos) > frag[i + 1]) {
        tu_v <- c(tu_v, tu[j])
      }
    }
  }
  return(tu_v)
}

#slope_function limits the slope 
slope_function <- function(data) {
  data[which(data$slope < 0.0001), "slope"] <- 0
  return(data)
}

#velo_function limits the velocity in case of infinite
velo_function <- function(data) {
  data$velocity_fragment[is.infinite(data$velocity_fragment)] <- NA
  return(data)
}


TI_frag_threshold <- function(data, TI_threshold) {
  frg <- data[, "TI_mean_termination_factor"]
  names(frg) <- data[, "TI_termination_fragment"]
  TIs  <- c()
  for (m in seq_len(length(frg) - 1)) {
    val <- frg[m] / frg[m + 1]
    if (val < TI_threshold) {
      TIs <- c(TIs, names(frg[m + 1]))
    }
  }
  return(TIs)
}

meanPosition <- function(input, parameter) {
  input <- as.data.frame(input %>%
                           group_by(input[, parameter]) %>%
                           mutate(meanPosi = mean(get('position'))))
  input <- input[!duplicated(input[, parameter]), ]
  return(input)
}

