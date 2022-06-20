fragment_HL <- function(probe, cores = 1, pen=pen, pen_out=pen_out) {
    num_args <- list(pen, pen_out)
    names(num_args) <- c("pen", "pen_out")    
    registerDoMC(cores)
    
    probe <- probe[with(probe, order(-xtfrm(probe$strand), probe$position)), ]
    
    probe[probe$strand == "-", ] <-
        probe[probe$strand == "-", ][order(probe[probe$strand == "-", ]$position,
                                           decreasing = TRUE), ]
    
    probe[, "HL_comb_fragment"] <- NA
    probe[, "HL_mean_comb_fragment"] <- NA
    
    # the position is not relevant in this case
    tmp_df <-
        data.frame(
            ID = probe$ID,
            val = probe$distance_HL,
            seg = probe$position_segment
        )
    
    # the fragmentation is performed on the delay_fragments, independent on
    # if they are (terminal) outliers or NAs.
    tmp_df[, "seg"] <- gsub("_O|_NA", "", tmp_df$seg)
    
    # although _NA will most probably will be dismissed here, because there
    # should be no half-life, if there is no delay
    tmp_df <- na.omit(tmp_df)
    
    unique_seg <- unlist(unique(tmp_df$seg))
    
    count <- 1
    
    # II. Dynamic Programming: the scoring function is interpreted
    
    frags <- foreach(k = seq_along(unique_seg)) %dopar% {
        section <- tmp_df[which(tmp_df$seg == unique_seg[k]), ]
        
        best_frags <- c()
        best_names <- c()
        
        if (nrow(section) > 1) {
            # only segments with more than one value are grouped...*
            
            for (i in 2:nrow(section)) {
                tmp_score <-
                    score_fun_ave(section[seq_len(i), "val"],
                                  section[seq_len(i), "ID"], pen)
                tmp_name <- names(tmp_score)
                if (i > 3) {
                    # fragments of the size at least 4 are allowed for half-life
                    for (j in (i - 1):3) {
                        tmp_val <- section[j:i, "val"]
                        tmp_ID <- section[j:i, "ID"]
                        tmp <-
                            score_fun_ave(tmp_val, tmp_ID, pen_out) + pen + best_frags[j - 2]
                        tmp_score <- c(tmp_score, tmp)
                        tmp_n <- paste0(best_names[j - 2], "|", names(tmp))
                        tmp_name <- c(tmp_name, tmp_n)
                    }
                }
                pos <- which(tmp_score == min(tmp_score))[1]
                tmp_score <- tmp_score[pos]
                tmp_name <- tmp_name[pos]
                best_frags <- c(best_frags, tmp_score)
                best_names <- c(best_names, tmp_name)
            }
        } else {
            #* ...all segments with less than two values are grouped automatically
            tmp_score <-
                score_fun_ave(section[, "val"], section[, "ID"], pen_out)
            tmp_name <- names(tmp_score)
            best_names <- c(best_names, tmp_name)
        }
        
        best_names[length(best_names)]
    }
    
    # III. Fill the dataframe
    
    for (k in seq_along(frags)) {
        na <- strsplit(frags[[k]], "\\|")[[1]]
        
        for (i in seq_along(na)) {
            tmp_trgt <- strsplit(na[i], "_")[[1]][1]
            trgt <- strsplit(tmp_trgt, ",")[[1]]
            if (length(strsplit(na[i], "_")[[1]]) == 3) {
                tmp_outl <- strsplit(na[i], "_")[[1]][3]
                outl <- strsplit(tmp_outl, ",")[[1]]
                trgt <- trgt[-which(trgt %in% outl)]
                rows <- match(outl, probe[, "ID"])
                nam <- paste0("Dc_", count, "_O")
                probe[rows, "HL_comb_fragment"] <- nam
                # the value assigned here is the mean (compare velocity and intercept
                # for delay)
                probe[rows, "HL_mean_comb_fragment"] <-
                    as.numeric(strsplit(na[i], "_")[[1]][2])
            }
            rows <- match(trgt, probe[, "ID"])
            nam <- paste0("Dc_", count)
            probe[rows, "HL_comb_fragment"] <- nam
            probe[rows, "HL_mean_comb_fragment"] <-
                as.numeric(strsplit(na[i], "_")[[1]][2])
            count <- count + 1
        }
    }
    row_NA <- which(is.na(probe$HL_comb_fragment))
    if (length(row_NA) > 0) {
        group <- c(row_NA[1])
        for (i in seq_along(row_NA)) {
            if (is.na(probe$HL_comb_fragment[row_NA[i] + 1])) {
                group <- c(group, row_NA[i] + 1)
            } else if (gsub("_O", "", probe$HL_comb_fragment[group[1] - 1]) ==
                       gsub("_O", "", probe$HL_comb_fragment[group[length(group)] + 1])) {
                probe$HL_comb_fragment[group] <-
                    paste0(gsub("_O", "", probe$HL_comb_fragment[group[1] - 1]), "_NA")
                probe$HL_mean_comb_fragment[group] <-
                    probe$HL_mean_comb_fragment[group[1] - 1]
                group <- row_NA[i + 1]
            } else if (probe$HL_comb_fragment[group[1] - 1] != probe$HL_comb_fragment[
                group[length(group)] + 1]) {
                probe$HL_comb_fragment[group] <-
                    paste0("Dc_", as.numeric(
                        gsub("Dc_|_O", "", probe$HL_comb_fragment[group[1] - 1])) + 0.5, "_NA")
                group <- row_NA[i + 1]
            }
        }
    }
    probe$HL_comb_fragment[is.na(probe$HL_comb_fragment)] <- paste0("Dc_", as.numeric(
        gsub("Dc_|_O", "", probe$HL_comb_fragment[!is.na(probe$HL_comb_fragment)]
             [length(probe$HL_comb_fragment[!is.na(probe$HL_comb_fragment)])])) + .5, "_NA")
    probe
}
