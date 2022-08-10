# =========================================================================
# adjusting_HLToInt  
# -------------------------------------------------------------------------
#' adjusting_HLToInt Creates one table merging HL and intensity fragments
#' with genome annotation
#'
#'
#' 'adjusting_HLToInt' merges HL and intensity segments adapting the positions 
#' to each other and combining to the genome annotation.
#' To make HL and intensity segments comparable, log2FC(HL) is used to generate 
#' the data frame instead of distance.
#' The fragments should have a significant p_value from t-test at least from one
#' segmentation, either HL or intensity.
#'
#' The functions used are:
#'
#' 1. p_value_function extracts and return the p_values of HL and intensity
#' segments respectively.
#'
#' 2. eliminate_outlier_hl eliminates outliers from HL fragments.
#'
#' 3. eliminate_outlier_int eliminates outliers from intensity fragments.
#'
#' 4. mean_length_int calculates the mean of the log2FC(intensity) fragments
#' adapted to HL_fragments and their lengths
#'
#' 5. mean_length_hl calculates the mean of log2FC(HL) fragments adapted to the
#' intensity fragments and their lengths.
#'
#' 6. calculating_rate calculates decay rate and log2FC(intensity). Both are
#' used to calculate synthesis rate.

#' @param data data frame: data frame combined data by column
#' @param Strand string: either "+" or "-"
#' @param annotation data frame: data frame from processed gff3 file.
#'
#' @return The data frame with the corresponding columns:
#'    \describe{
#'      \item{position:}{Integer, position of the first fragment}
#'      \item{region:}{String, region annotation covering the fragments}
#'      \item{gene:}{String, gene annotation covering the fragments}
#'      \item{locus_tag:}{String, locus_tag annotation covering the fragments}
#'      \item{strand:}{Boolean. The bin/probe specific strand (+/-)}
#'      \item{fragment_HL:}{String, HL fragments}
#'      \item{fragment_int:}{String, intensity fragments}
#'      \item{position_frg_int:}{Integer, position of the first fragment and
#'      the last position of the last fragment}
#'      \item{mean_HL_fragment:}{Integer, mean of the HL of the fragments
#'      involved}
#'      \item{mean_int_fragment:}{Integer, mean of the intensity of the
#'      fragments involved}
#'      \item{log2FC(decay_rate):}{Integer, log2FC(decay(condition1)/
#'      decay(condition2))}
#'      \item{log2FC(synthesis_rate):}{Integer, sum of log2FC(decay_rate) and
#'      log2FC(intensity)}
#'      \item{intensity_FC:}{Integer, log2FC(mean(intensity(condition1))/mean(
#'      intensity(condition2)))}
#'      \item{Log2FC(HL)+Log2FC(int):}{Integer, sum of log2FC(decay_rate) and
#'      log2FC(intensity)}
#'      \item{p_value:}{String, indicated by "*" means at least one fragment
#'      either HL fragment or intensity fragment has a significant p_value}
#'     }
#'
#' @examples
#' data(stats_df_comb_minimal)
#' data(annot_g)
#' df_mean_minimal <- adjusting_HLToInt(data = stats_df_comb_minimal,
#' annotation = annot_g[[1]])
#' @export


adjusting_HLToInt <-
    function(data,
             Strand = c("+", "-"),
             annotation) {
        # distance HL between two conditions. log2 is applied after dynamic
        # programming to compare HL with intensity
        data[, "distance_HL_log"] <-
            log2(data[, "half_life.cdt1"] / data[, "half_life.cdt2"])
        
        #select unique fragments without outliers
        frag_HL <-
            unique(eliminate_outlier_hl(data)[, "HL_comb_fragment"])
        
        #assign empty vectors
        df <- data.frame()
        Position <- c()
        Region <- c()
        Gene <- c()
        Locus_tag <- c()
        STrand <- c()
        Fragment_HL <- c()
        Fragment_int <- c()
        position_frg_int <- c()
        Mean_HL_fragment <- c()
        Mean_int_fragment <- c()
        Decay_rate <- c()
        Synthesis_rate <- c()
        intensity_FC <- c()
        p_value <- c()
        
        for (j in seq_along(Strand)) {
            #select fragment by strand
            fragments_HL <-
                data[which(data$strand == Strand[j]), "HL_comb_fragment"]
            
            #eliminate outliers
            fragments_HL <-
                unique(fragments_HL[grep(paste0("Dc_\\d+", "$"), fragments_HL)])
            
            #loop into HL fragments
            for (i in seq_along(fragments_HL)) {
                #assign objects
                fg_hl <- NA
                fg_int <- NA
                position <- NA
                region <- NA
                gene <- NA
                locus_tag <- NA
                strand <- NA
                fragment_HL <- NA
                fragment_int <- NA
                mean_HL_fragment <- NA
                mean_int_fragment <- NA
                decay_rate <- NA
                synthesis_rate <- NA
                Description <- NA
                
                #selecting strand
                fg <- data[which(data$strand %in% Strand[j]), ]
                
                #selecting the positions covering HL_comb fragment
                fg_hl_pos <-
                    fg[which(fg$HL_comb_fragment ==
                                 fragments_HL[i]), "position"]
                
                #positions boarders
                pos.1 <- fg_hl_pos[1]
                pos.2 <- last(fg_hl_pos)
                
                #in case of strand is negative, the positions are flipped
                if (j == 2) {
                    pos.1 <- last(fg_hl_pos)
                    pos.2 <- fg_hl_pos[1]
                }
                
                #df with the corresponding positions
                fg_pos <-
                    fg[between(fg$position, pos.1, pos.2),]
                
                # select strand on the annotation data frame
                ann <- as.data.frame(annotation) %>%
                    filter(strand == Strand[j])

                tryCatch({
                    region <-  
                        paste0(annotation_function_df(
                            feature = "region",
                            pos.1 = pos.1, 
                            pos.2 = pos.2,
                            strand = Strand[j],
                            data_annotation = ann
                        ), collapse = "|")
                    
                    gene <- 
                        paste0(annotation_function_df(
                            feature = "gene",
                            pos.1 = pos.1, 
                            pos.2 = pos.2,
                            strand = Strand[j],
                            data_annotation = ann
                        ), collapse = "|")
                    
                    locus_tag <-
                        paste0(annotation_function_df(
                            feature = "locus_tag",
                            pos.1 = pos.1, 
                            pos.2 = pos.2,
                            strand = Strand[j],
                            data_annotation = ann
                        ), collapse = "|")
                    
                    Description <-
                        paste0(annotation_function_df(
                            feature = "Annotation",
                            pos.1 = pos.1, 
                            pos.2 = pos.2,
                            strand = Strand[j],
                            data_annotation = ann
                        ), collapse = "|")
                }, error = function(e) {
                    "No match for the annotation"
                })
                
                #eliminate HL and intensity outliers
                f_withoutOutlier <- eliminate_outlier_hl(fg_pos)
                f_withoutOutlier <-
                    eliminate_outlier_int(f_withoutOutlier)
                f_int <-
                    unique(f_withoutOutlier$intensity_comb_fragment)
                
                #split intensity fragments in case length is more than 1,
                if (length(f_int) > 1) {
                    for (k in seq_along(f_int)) {
                        #delimit the intensity fragment positions
                        hl_int <- f_withoutOutlier %>%
                            filter(intensity_comb_fragment == f_int[k])
                        
                        #extract the mean of distance delimiting int fragment positions
                        fg_hl <- mean_length_hl(hl_int)
                        
                        #delimit the intensity fragment positions
                        fg_int <- mean_length_int(hl_int)
                        
                        # check p_values
                        p_values <- p_value_function(hl_int)
                        
                        Position <-
                            c(Position, f_withoutOutlier$position[1])
                        Region <- c(Region, region)
                        Gene <- c(Gene, gene)
                        Locus_tag <- c(Locus_tag, locus_tag)
                        #Annotation <- c(Annotation, Description)
                        STrand <- c(STrand, Strand[j])
                        
                        #HL fragment
                        Fragment_HL <-
                            c(Fragment_HL, fragments_HL[i])
                        
                        #int fragment
                        Fragment_int <- c(Fragment_int, f_int[k])
                        #border positions of the fragment
                        position_frg_int <-
                            c(position_frg_int,
                              paste(
                                  hl_int$position[1],
                                  last(hl_int$position),
                                  sep = ":"
                              ))
                        
                        if (length(p_values) < 2) {
                            Mean_HL_fragment <- c(Mean_HL_fragment, NA)
                            Mean_int_fragment <-
                                c(Mean_int_fragment, NA)
                            Decay_rate <- c(Decay_rate, NA)
                            Synthesis_rate <- c(Synthesis_rate, NA)
                            intensity_FC <- c(intensity_FC, NA)
                            p_value <- c(p_value, NA)
                        }
                        
                        if (length(p_values) == 2) {
                            #mean of the distance on the fragment hl
                            Mean_HL_fragment <-
                                c(Mean_HL_fragment, fg_hl$mean)
                            #mean of the distance on the fragment int
                            Mean_int_fragment <-
                                c(Mean_int_fragment, fg_int$mean)
                            # decay_rate is calculated dividing log2/mean(HL) of
                            # both condition
                            # please refer to calculating_rate function
                            Decay_rate <- c(Decay_rate,
                                            calculating_rate(hl_int)[[1]])
                            ##ratio between log2FC(decay) and log2FC(intensity)
                            intensity_FC <-
                                c(intensity_FC, log2(
                                    mean(hl_int$intensity.cdt1) /
                                        mean(hl_int$intensity.cdt2)
                                ))
                            # synthesis_rate is calculated from steady-state and
                            # decay rate
                            # multiplication of log2FC(decay) and log2FC(intensity)
                            Synthesis_rate <- c(
                                Synthesis_rate,
                                calculating_rate(hl_int)[[1]] +
                                    calculating_rate(hl_int)[[2]]
                            )
                            if (p_values[1] < .05 |
                                p_values[2] < .05) {
                                p_value <- c(p_value, "*")
                            } else{
                                p_value <- c(p_value, NA)
                            }
                        }
                    }
                }
                if (length(f_int) == 1) {
                    p_values <- p_value_function(fg_pos)
                    fg_hl <- mean_length_hl(f_withoutOutlier)
                    fg_int <- mean_length_int(f_withoutOutlier)
                    Position <-
                        c(Position, f_withoutOutlier$position[1])
                    Region <- c(Region, region)
                    Gene <- c(Gene, gene)
                    Locus_tag <- c(Locus_tag, locus_tag)
                    STrand <- c(STrand, Strand[j])
                    Fragment_HL <- c(Fragment_HL, fragment_HL[i])
                    Fragment_int <- c(Fragment_int, f_int)
                    position_frg_int <- c(
                        position_frg_int,
                        paste(
                            f_withoutOutlier$position[1],
                            last(f_withoutOutlier$position),
                            sep = ":"
                        )
                    )
                    
                    if (length(p_values) < 2) {
                        Mean_HL_fragment <- c(Mean_HL_fragment, NA)
                        Mean_int_fragment <-
                            c(Mean_int_fragment, NA)
                        Decay_rate <- c(Decay_rate, NA)
                        Synthesis_rate <- c(Synthesis_rate, NA)
                        intensity_FC <- c(intensity_FC, NA)
                        p_value <- c(p_value, NA)
                    }
                    if (length(p_values) == 2) {
                        Mean_HL_fragment <- c(Mean_HL_fragment, fg_hl$mean)
                        Mean_int_fragment <-
                            c(Mean_int_fragment, fg_int$mean)
                        Decay_rate <-
                            c(Decay_rate,
                              calculating_rate(f_withoutOutlier)[[1]])
                        Synthesis_rate <- c(
                            Synthesis_rate,
                            calculating_rate(f_withoutOutlier)[[1]] +
                                calculating_rate(f_withoutOutlier)[[2]]
                        )
                        intensity_FC <- c(intensity_FC,
                                          log2(
                                              mean(
                                                  f_withoutOutlier$intensity.cdt1) /
                                              mean(f_withoutOutlier$intensity.cdt2)
                                          ))
                        if (p_values[1] < .05 | p_values[2] < .05) {
                            p_value <- c(p_value, "*")
                        } else{
                            p_value <- c(p_value, NA)
                        }
                    }
                }
            }
            df <-
                cbind.data.frame(
                    Position,
                    Region,
                    Gene,
                    Locus_tag,
                    STrand,
                    Fragment_HL,
                    Fragment_int,
                    position_frg_int,
                    Mean_HL_fragment,
                    Mean_int_fragment,
                    Decay_rate,
                    Synthesis_rate,
                    intensity_FC,
                    p_value
                )
        }
        return(df)
    }
