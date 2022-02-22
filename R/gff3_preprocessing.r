source("functions_gff.r")
gff3_preprocess <- function(path) {
    ge_size <- genomeSize(path)
    #select columns region, start, end, strand and annotation
    inp <- readGFF(path)
    inp <- inp[-1,]
    tmp <- inp[, c("type", "start", "end", "strand", "gene", "locus_tag")]
    tmp <-
        tmp[grep("^gene$",
                 tmp$type,
                 invert = TRUE), ]
    colnames(tmp)[1] <- "region"
    #replace gene with NA with locus_tag in case NA is not recognized as NA
    tmp[which(tmp$gene == "NA"), "gene"] <-
        tmp[which(tmp$gene == "NA"), "locus_tag"]
    #replace gene with NA with locus_tag
    tmp[is.na(tmp$gene), "gene"] <- tmp[is.na(tmp$gene), "locus_tag"]
    if (length(which(tmp$region %in% "antisense_RNA")) != 0) {
        tmp[which(tmp$region == "antisense_RNA"), "region"] <- "asRNA"
    }
    return(list(tmp, ge_size))
}
