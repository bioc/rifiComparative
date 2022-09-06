#' =============================================================================
#' rifiComparative
#' -----------------------------------------------------------------------------
#' rifiComparative a successor package of rifi. It compares 2 rifi outputs from 
#' 2 different conditions of the same organism.
#'
#' rifiComparative was developed to compare 2 rifi outputs from 2 conditions. 
#' The rifi output may differ significantly from 2 conditions. The complexity of
#' the segments number, position, length and the events make the comparison 
#' between 2 conditions nearly impossible. rifiComparative uses a simple strategy 
#' to generate single segments for both conditions, extract the features and 
#' make them comparable.
#' 
#' Five major steps ate described in rifiComparative:
#' 
#' 1. Joining data
#' 2. penalties
#' 3. fragmentation
#' 4. statistics
#' 5. visualization
#'
#' @import dplyr
#' @import ggplot2
#' @import SummarizedExperiment
#' @import scales
#' @import writexl
#' @import DTA
#' @import ggrepel
#' @import LSD
#' @author Loubna Youssar \email{lyoussar@@gmail.com}
#' @author Jens Georg \email{jens.georg@@biologie.uni-freiburg.de}
#' @name rifiComparative
#' @docType package
#' 
#' 
NULL
