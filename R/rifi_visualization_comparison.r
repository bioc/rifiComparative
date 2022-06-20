#' rifi_visualization_comparison: plots all the data from two conditions 
#' with fragments and events of both strands.
#' rifi_visualization_comparison: plots the whole genome with genes, 
#' transcription units (TUs), delay, half-life (HL), intensity fragments
#' features, events and annotation from both conditions.
#' rifi_visualization_comparison uses several functions to plot the genes
#' including as-RNA and ncRNA and TUs as segments.
#' The function plots delay, HL and intensity fragments with statistical t-test
#' between the neighboring fragment, significant p-value from t-test and manova 
#' are assigned with '*'.
#' The functions used are:
#' annotation_plot: plots the corresponding annotation
#' strand_comp_function: plots delay, HL, intensity and events as pausing site 
#' and iTSS_I
#' empty_data_positive: plots empty boxes in case no data is available for
#' positive strand
#' empty_data_negative: plots empty boxes in case no data is available for
#' negative strand
#' splitGenome_function: splits the genome into fragments
#' indice_function: assign a new column to the data to distinguish between
#' fragments, outliers from delay or HL or intensity.
#' TU_annotation: designs the segments border for the genes and TUs annotation
#' gene_annot_function: it requires gff3 file, returns a dataframe adjusting
#' each fragment according to its annotation. It allows as well the plot of
#' genes and TUs shared into two pages
#' label_log2_function: used to add log scale to intensity values.
#' label_square_function: used to add square scale to coverage values.
#' coverage_function: this function is used only in case of coverage is
#' available.
#' secondaryAxis: adjusts the half-life or delay to 20 in case of the dataframe
#' row numbers is equal to 1 and
#' the half-life or delay exceed the limit, they are plotted with different
#' shape and color.
#' my_arrow: creates an arrow for the annotation.
#' arrange_byGroup: selects the last row for each segment and add 40 nucleotides
#'  in case of negative strand for a nice plot.


#'
#' @param data dataframe: the probe based dataframe.
#' @param genomeLength integer: genome length output of gff3_preprocess
#' function.
#' @param annot dataframe: the annotation file.
#' @param region dataframe: gff3 features of the genome.
#' @param color_region string vector: vector of colors.
#' @param color_TU string: TU colors
#' @param fontface integer: value assigning labels font
#' @param color_text.1 string: TU color text
#' @param color_text.2 string: genes color text
#' @param size_tu integer: TU size
#' @param size_locusTag integer: locus_tag size
#' @param Limit integer: value for y-axis limit.
#' @param shape integer: value for shape.
#' @param color_TU string. TU color
#' @param limit_intensity integer: intensity limit if applicable.
#' @param face string: label font.
#' @param tick_length integer: value for ticks.
#' @param arrow.color string: arrows color.
#' @param minVelocity integer: threshold to fix the minimum of velocity.
#' @param medianVelocity integer: threshold to fix the maximum of velocity.
#' @param col_above20 string: color for probes/bin above value 20.
#' @param fontface integer: font type
#' @param shape_above20 integer: shape for probes/bins above value 20.
#' @param axis_text_y_size integer: text size for y-axis.
#' @param axis_title_y_size integer: title size for y-axis.
#' @param Alpha integer: color transparency degree.
#' @param size_gene integer: font size for gene annotation.

#'
#' @return The visualization.
#'
#' @examples
#' data(stats_minimal)
#' data(annot_g_minimal)
#' rifi_visualization(data = stats_minimal, genomeLength = annot_g_minimal[[2]],
#' annot = annot_g_minimal[[1]], region = c("CDS","asRNA","5'UTR","ncRNA","3'UTR","tRNA"),
#' color_region = c("grey0", "red", "blue", "orange", "yellow", "green",
#' "white", "darkseagreen1", "grey50", "black"),
#' color_text.1 = "grey0", color_text.2 = "black", color_TU = "blue",
#' size_tu = 1.6, size_locusTag = 1.6, size_gene = 1.6, Limit = 10,
#' shape=22, col_outiler = "grey50", Alpha=0.5, limit_intensity = NA,
#' face="bold", tick_length = .3, arrow.color = "darkseagreen1",
#' minVelocity = 3000, medianVelocity = 6000,
#' col_above20 = "#00FFFF", fontface = "plain", shape_above20 = 14,
#' axis_text_y_size = 3, axis_title_y_size = 6)
#'
#' @export

setwd("/home/loub/rifi_comparison/analysis_window/plot/")
#setwd("/home/lyoussar/rifi_comparison/analysis_window/plot/")
source("annotation_plot_comp.r")
source("functions_visuali_condts.r")
source("strand_comp_Comb_function.r")
load("annot_g_c.rda")
load("annot_g.rda")
setwd("/home/loub/rifi_comparison/outputs")

load("data_combined_se.rdata")
load("df_comb_se.rda")


rifi_visualization_comparison <-
    function(data,
             data_c,
             genomeLength = annot_g[[2]],
             annot = annot_g[[1]],
             condition = c("cdt1", "cdt2"),
             Strand = c("+", "-"),
             region = c("CDS", "asRNA", "5'UTR", "ncRNA", "3'UTR", "tRNA"),
             color_region = c(
                 "grey0",
                 "red",
                 "blue",
                 "orange",
                 "yellow",
                 "green",
                 "white",
                 "darkseagreen1",
                 "grey50",
                 "black"
             ),
             color_TU = c("cyan", "yellow", "orange"),
             scaling_TU = c(0, 3.4, 6.6),
             color_text.1 = "grey0",
             color_text.2 = "black",
             Alpha = 0.5,
             size_tu = 1.6,
             size_locusTag = 1.6,
             size_gene = 1.6,
             Limit = 10,
             shape = 22,
             face = "bold",
             tick_length = .3,
             arrow.color = "darkseagreen1",
             col_above20 = "#00FFFF",
             fontface = "plain",
             shape_above20 = 14,
             axis_text_y_size = 3,
             axis_title_y_size = 6,
             iTSS_threshold = 1.2,
             p_value_manova = 0.05,
             termination_threshold = 0.8
    ) {
        ##########################data preparation#############################
        #II. input for the main features split into 2 data frames according to
        #strand orientation
        tmp.1 <- strand_selection(data, "+")
        tmp.2 <- strand_selection(data, "-")
        tmp.3 <- strand_selection(data_c, "+")
        tmp.4 <- strand_selection(data_c, "-")  
        #replace infinitive in velocity fragment with NA
        tmp.1 <- velo_function(tmp.1)
        tmp.2 <- velo_function(tmp.2)
        #replace slope lower than 0.0009 to 0
        tmp.1 <- slope_function(tmp.1)
        tmp.2 <- slope_function(tmp.2)
        #III. split the genome into fragments
        gLength <- seq_len(genomeLength)
        names(gLength) <- seq_along(gLength)
        fl <- floor(gLength / 10000)
        frag <- splitGenome_function(x = fl, gLength = gLength)
        #################################plot###############################
        #IV. the general plot
        pdf.options(
            onefile = TRUE,
            width = 8,
            height = 5.3
        )
        pdf("genome_fragments_comparison.pdf")
        suppressWarnings(for (i in seq_len(length(frag) - 1)) {#seq_len(length(frag) - 1)
            p <- list()
            p.1 <- list()
            print(i)
            if (i == 1) {
                frag[i] <- 0
            } else if (i == (length(frag) - 1)) {
                #to have homogeneous annotation scaling, 10000 is added to the last frag vector
                frag[i + 1] <- frag[i] + 10000
            }
            #adjust position for genes split on two pages
            pos.1 <- frag[i] - 2000
            pos.2 <- frag[c(i + 1)] + 2000
            ###########################data adjustment###########################
            #define the main dataframe with segments positive strand df1, negative
            #strand df2
            df1 <-
                tmp.1[between(tmp.1$position, frag[i], frag[c(i + 1)]),]
            df2 <-
                tmp.2[between(tmp.2$position, frag[i], frag[c(i + 1)]),]
            df3 <-
                tmp.3[between(tmp.3$position, frag[i], frag[c(i + 1)]),]
            df4 <-
                tmp.4[between(tmp.4$position, frag[i], frag[c(i + 1)]),]
            
            df1_1 <- df1[!is.na(df1$ID),]
            #avoid plot empty pages in case of small data
            if (nrow(df1) == 0 &
                nrow(df2) == 0 & nrow(data) < 10000) {
                next ()
            }
            ##########################annotation section#########################
            #an is the annotation dataframe upon the position on the plot, its used
            # to loop into exactly the number of region contained in the gff3
            an.1 <-
                annot[between(annot$start, frag[i], frag[c(i + 1)]), ]
            an <- annot[between(annot$start, pos.1, pos.2), ]
            an <- an[!duplicated(an), ]
            #in case of no data nor annotation are available
            if (nrow(an.1) == 0 & nrow(df1) == 0 & nrow(df2) == 0) {
                next ()
            }
            p7 <-
                annotation_plot(
                    data_p=df1,
                    data_n=df2,
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
                    pos.2 = pos.2
                )
            
            #########################empty data positive strand##################
            
            if (nrow(df1) == 0) {
                p_positive <- empty_data_positive(
                    data_p = df1,
                    data_n = df2,
                    frag = frag,
                    i = i,
                    axis_title_y_size = axis_title_y_size,
                    axis_text_y_size = axis_text_y_size,
                    Limit = Limit
                )
                p1 <- p_positive[[1]]
                p2 <- p_positive[[2]]
                #p3 <- p_positive[[3]]
            }
            
            ############################positive_strand_plot######################
            if (nrow(df1) != 0) {
                p_positive <- strand_function(
                    data = data,
                    data_p = df1,
                    data_n = df2,
                    data_p_c = df3,
                    data_n_c = df4,
                    Strand = Strand,
                    condition = condition,
                    frag,
                    i,
                    fontface = fontface,
                    axis_title_y_size = axis_title_y_size,
                    axis_text_y_size = axis_text_y_size
                )
                p1 <- p_positive[[1]]
                p2 <- p_positive[[2]]
                #p3 <- p_positive[[3]]
                p4 <- p_positive[[4]]
                p5 <- p_positive[[5]]
                #p6 <- p_positive[[6]]
            }
            #########################empty data reverse strand##################
            if (nrow(df2) == 0) {
                p_negative <- empty_data_negative(
                    data_n = df2,
                    frag = frag,
                    i = i,
                    axis_title_y_size = axis_title_y_size,
                    axis_text_y_size = axis_text_y_size,
                    Limit = Limit
                )
                #p6 <- p_negative[[1]]
                p5 <- p_negative[[2]]
                p4 <- p_negative[[3]]
            }
            ############################Title and plot#############################
            p <- list(p1, p2, p7, p5, p4)
            egg::ggarrange(
                plots = p,
                ncol = 1,
                nrow = 5,
                heights = c(5, 5, 6, 5, 5)
            )
        })
        dev.off()
    }

#rifi_visualization_comparison(data = d, data_c = dc)

rifi_visualization_comparison(
    data = data_combined, #data combined from probe_sta (rifi output)
    data_c = df_comb, #data from DP
    genomeLength = annot_g[[2]],
    annot = annot_g[[1]]
)

