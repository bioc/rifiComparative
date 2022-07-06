score_fun_ave <-
    function(y, z, pen, n_out = min(10, max(1, 0.2 * length(z)))) {
        n_out <-
            min(n_out, length(z) - 3) # the number of allowed outliers is decided
        tmp_ave <- mean(y) # the mean is calculated
        tmp_dif <- abs(y - tmp_ave) # the difference to the mean is calculated
        names(tmp_dif) <-
            z # tmp_dif and y are named by z to decide for outliers
        names(y) <- z
        out <-
            sort(tmp_dif, decreasing = TRUE)[1] # out is is the sorted vector of
        #differences to the mean
        tmp_score <-
            sum(tmp_dif) # the sum of the differences to the mean is used as score
        tmp_mean <- tmp_ave # the mean is cached
        if (n_out >= 1) {
            # checks if more than 0 outliers are allowed
            for (i in seq_len(n_out)) {
                # the loop iterates over the number of allowed outliers
                tmp_out <- names(out) # all but the i worst are selected
                tmp_y <- y[!names(y) %in% tmp_out]
                tmp_ave_o <- mean(tmp_y) # new mean
                tmp_dif_o <- abs(tmp_y - tmp_ave_o) # new difference to mean
                out <- c(out, sort(tmp_dif_o, decreasing = TRUE)[1])
                tmp_score_o <-
                    sum(tmp_dif_o) + pen * i # new score with outliers penalty
                tmp_score <- c(tmp_score, tmp_score_o) # the score(s) are cached
                tmp_mean <- c(tmp_mean, tmp_ave_o) # the mean(s) are cached
            }
        }
        mi <-
            which(tmp_score == min(tmp_score))[1] # the lowest score is selected
        nam <- paste0(z, collapse = ",") # all names (z) are collected
        nam <- paste0(nam, "_", tmp_mean[mi]) # the mean is pasted to the name
        if (mi > 1) {
            # checks if outliers were found
            outlier <-
                paste(names(out)[seq_len(mi - 1)], collapse = ",") # pastes all
            #outliers together
            nam <- paste0(nam, "_", outlier) # pastes the outliers to the mean
        }
        res <- tmp_score[mi] # the lowest score
        names(res) <- nam # the name
        res
    }
