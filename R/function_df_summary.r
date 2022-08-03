#the function tries to cover all annotations. Positions 1 and 2 cover those segments starting before the beginning of the gene/TU annotated. 
annotation_function_df <-
  function(feature, pos.1, pos.2, strand, data_annotation) {
    data_annotation[, "region"]  <-
      as.character(data_annotation[, "region"])
    data_annotation <-
      data_annotation[which(data_annotation[, "strand"] %in% strand),]
    positions.1 <- data_annotation[between(data_annotation$end, pos.1, pos.2), "end"]
    positions.2 <- data_annotation[between(data_annotation$start, pos.1, pos.2),
                                   "start"]
    positions <- unique(c(positions.1, positions.2))
    positions <- positions[order(positions, decreasing = F)]
      features  <-
            unique(c(data_annotation[between(data_annotation$start, positions[1],
                                    last(positions)), feature],
             data_annotation[between(data_annotation$end, positions[1],
                              last(positions)), feature]))
 
    return(features)
  }


