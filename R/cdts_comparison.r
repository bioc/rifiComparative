# cdts_comparison compares two data outputs from Rifi framework.
# 
#' @param specie string: organism name.
#' @param dt1 dataframe: the first element from rifi_summary list.
#' @param annot data.table: data.table generated from gff3_preprocessing function
#' @param indice vector: vector of conditions as cdt1 and cdt2.
#' @param dt3 dataframe: the third element from rifi_summary list.
#' @param input dataframe: joined dataframes generated from rifi_statistics 
#' of both conditions. A column is added indicating the condition.
#
# To use the function just add the corresponding parameters
# df_cdt <- cdts_comparison(
# specie = "Synechocystis_sp_PCC6803",
# dt1 = dt1,
# annot = annot_g_c,
# indice = c("cdt1", "cdt2"),
# dt3 = dt3,
# input = data_combined,
# conditions= c("condition 1 name", "condition 2 name")
# )

cdts_comparison <- function(specie,
                            dt1,
                            annot,
                            indice,
                            dt3,
                            input,
                            conditions) {
    df_cdt <- data.frame()
    #select all annotation features e.g. "CDS", "ncRNA"
    l <- unique(annot$region)
    for (k in seq_along(l)) {
        #create an empty dataframe where all segments features are gathered
        df <- data.frame(
            locus_tag = character(),
            gene = character(),
            feature = character(),
            strand = character(),
            delay_frg_cdt1 = numeric(0),
            delay_frg_cdt2 = numeric(0),
            HL_frg_cdt1 = numeric(0),
            HL_frg_cdt2 = numeric(0),
            HL_cdt1 = numeric(0),
            HL_cdt2 = numeric(0),
            HL_mean_cdt1 = numeric(0),
            HL_mean_cdt2 = numeric(0),
            log2FC_HL = numeric(0),
            int_frg_cdt1 = numeric(0),
            int_frg_cdt2 = numeric(0),
            int_cdt1 = numeric(0),
            int_cdt2 = numeric(0),
            int_mean_cdt1 = numeric(0),
            int_mean_cdt2 = numeric(0),
            log2FC_int = numeric(0),
            log2FC_HL_int = numeric(0),
            paus_cdt1 = numeric(0),
            paus_cdt2 = numeric(0),
            iTSS_I_cdt1 = numeric(0),
            iTSS_I_cdt2 = numeric(0),
            start = numeric(0),
            end = numeric(0),
            iTSS_II_cdt1 = numeric(0),
            iTSS_II_cdt2 = numeric(0),
            ter_cdt1 = numeric(0),
            ter_cdt2 = numeric(0),
            HL_P_value_cdt1 = numeric(0),
            HL_P_value_cdt2 = numeric(0),
            int_P_value_cdt1 = numeric(0),
            int_P_value_cdt2 = numeric(0),
            stringsAsFactors = F
        )
        #select annotation region by category
        annotation <- annot[which(annot$region %in% l[[k]]),]
        # remove duplicated locus_tag with the same features (region).
        # This happen when a gene has many transcripts starting from
        # different 5'UTR.
        # select the locus tag in case of duplication
        lt_dup <-
            annotation[duplicated(annotation$locus_tag), "locus_tag"]
        
        if (length(lt_dup) != 0) {
            lt_ndup <- annotation[-which(annotation$locus_tag %in% lt_dup),]
            
            ann <- data.frame()
            for (j in seq_along(lt_dup)) {
                an <- annotation[which(annotation$locus_tag == lt_dup[j]), ]
                an$start <- an$start[1]
                an$end <- last(an$end)
                ann <- rbind(ann, an)
                ann <- ann[!duplicated(ann),]
            }
            annotation <- rbind(lt_ndup, ann)
            annotation <-
                annotation[order(annotation$start, decreasing = F),]
        }
        #gather locus tag  for next loop
        genes <- unique(na.omit(annotation$locus_tag))
        #loop into all locus tag from same region
        for (i in seq_along(genes)) {
            print(paste0("i: ", i))
            #assign the locus tag
            locus_t <- as.character(genes[i])
            #extract the IRanges of the locus_tag from the annotation df
            an <-
                annot[grep(paste0("^", locus_t, "$"), annot$locus_tag), ]
            an <- an[which(an$region == l[[k]]), ]
            #in case locus_tag is collapsed with others
            #add locus_tag to df
            df[i, "locus_tag"] <- locus_t
            #add gene to df
            df[i, "gene"] <- paste(unique(an$gene), collapse = "|")
            df[i, "feature"] <-
                paste(unique(an$region), collapse = "|")
            
            #only in case of synechosystis, the annotation is customized and the first
            #gene has 2 annotations as its present on the start and at the end of
            #the chromosome
            if (specie == specie) {
                an <- as.data.frame(an) %>%
                    group_by(locus_tag) %>%
                    distinct(region, .keep_all = T)
                if (i == 1 & an$region == "CDS") {
                    an <- data.frame(
                        region = "CDS",
                        start = 1,
                        #(continuation 3573271)
                        end = 772,
                        #(continuation 3573470)
                        strand = "+",
                        locus_tag = locus_t,
                        stringsAsFactors = F
                    )
                }
            }
            #fill in the dataframe
            df[i, "start"] <- an[, "start"]
            df[i, "end"] <- an[, "end"]
            df[i, "strand"] <- an[, "strand"]
            #grep locus tag from df1
            lg <-
                dt1[grep(paste0("^", locus_t, "$"), dt1$locus_tag), ]
            lg <- lg[which(lg$feature_type == l[[k]]),]
            #function_df assemble half-life, intensity fragment corresponding 
            #to the flanking region of the locus tag
            if (nrow(lg) != 0) {
                df <- function_df (
                    indice = indice,
                    dt3 = dt3,
                    input = data_combined,
                    an = an,
                    i = i,
                    data = df,
                    locus_t = locus_t,
                    conditions = conditions
                )
            }
        }
        #the dataframes generated are assembled into df_cdt
        df_cdt <- rbind(df_cdt, df)
    }
    #assign log2FC for half-life, intensity and ratio of both
    df_cdt[, "log2FC_HL"] <-
        log2(df_cdt[, paste0("HL_mean_", indice[1])] /
                 df_cdt[, paste0("HL_mean_", indice[2])])
    df_cdt[, "log2FC_int"] <-
        log2(df_cdt[, paste0("int_mean_", indice[1])] /
                 df_cdt[, paste0("int_mean_", indice[2])])
    df_cdt[, "log2FC_HL_int"] <-
        df_cdt[, "log2FC_HL"] / df_cdt[, "log2FC_int"]
    #apply t-test to half-life and intensity fold change
    df_cdt <-
        function_tTest(data = df_cdt,
                       parameter = "HL",
                       indice = indice)
    df_cdt <-
        function_tTest(data = df_cdt,
                       parameter = "int",
                       indice = indice)
    #apply manova test to ratio of two fold changes
    ##eliminate the columns not needed anymore
    df_cdt[, c(
        paste0("HL_P_value_", indice[1]),
        paste0("HL_P_value_", indice[2]),
        paste0("int_P_value_", indice[1]),
        paste0("int_P_value_", indice[2])
    )] <- list(NULL)
    #add an ID column for an easy reference
    df_cdt <- df_cdt %>%
        add_column(ID = 1:nrow(df_cdt), .before = "locus_tag")
    #substitute indice conditions with the conditions indicated
    colnames(df_cdt) <-
        mgsub(colnames(df_cdt), c("cdt1", "cdt2"), conditions)
    return(df_cdt)
}
write_xlsx(df_cdt, "summary_cdts.xlsx")
save(df_cdt, file = "summary_cdts.rdata")
