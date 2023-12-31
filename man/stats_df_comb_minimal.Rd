% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.r
\docType{data}
\name{stats_df_comb_minimal}
\alias{stats_df_comb_minimal}
\title{The result of statistics for fragment_int example data
A data frame containing the output of statistics as a data frame}
\format{
A data frame with 500 rows and 26 variables:
\describe{
\item{strand:}{The strand specific}
\item{position:}{The bin/probe specific position}
\item{ID:}{The bin/probe specific ID}
\item{intensity.cdt1:}{The relative intensity at time point 0 for
condition 1}
\item{position_segment:}{The position based segment}
\item{half_life.cdt1:}{The half-life of the bin/probe condition 1}
\item{TI_termination_factor.cdt1:}{The termination factor of the bin/probe
condition 1}
\item{HL_fragment.cdt1:}{The half-life fragment the bin belongs to
condition 1}
\item{intensity_fragment.cdt1:}{The intensity fragment the bin belongs to
condition 1}
\item{TI_termination_fragment.cdt1:}{The TI fragment the bin belongs to
condition 1}
\item{logFC_int:}{The bin/probe log2 fold change of intensity at time 0}
\item{P.Value:}{The bin/probe p_value adjusted }
\item{intensity.cdt2:}{The relative intensity at time point 0
condition 2}
\item{half_life.cdt2:}{The half-life of the bin/probe condition 2}
\item{TI_termination_factor.cdt2:}{The termination factor of the bin/probe
condition 2}
\item{HL_fragment.cdt2:}{The half-life fragment the bin belongs to
condition 2}
\item{intensity_fragment.cdt2:}{The intensity fragment the bin belongs to
condition 2}
\item{TI_termination_fragment.cdt2:}{The TI fragment the bin belongs to
condition 2}
\item{distance_HL:}{The bin/probe difference of half-life from both
conditions}
\item{distance_int:}{The bin/probe log2 fold change of intensity at time 0}
\item{HL_comb_fragment:}{The half-life fragment the bin belongs to both
conditions}
\item{HL_mean_comb_fragment:}{The half-life mean of the fragment the bin
belongs to both conditions}
\item{intensity_comb_fragment:}{The intensity fragment the bin belongs to
both conditions}
\item{intensity_mean_comb_fragment:}{The intensity mean of the fragment the
bin belongs to both conditions}
\item{p_value_distance_HL:}{The p_value adjusted of the half-life fragment
the bin belongs to both conditions}
\item{p_value_distance_intensity:}{The p_value adjusted of the intensity
fragment the bin belongs to both conditions}
}
}
\source{
\url{https://github.com/CyanolabFreiburg/rifiComparative}
}
\usage{
data(stats_df_comb_minimal)
}
\description{
The result of statistics for fragment_int example data
A data frame containing the output of statistics as a data frame
}
\keyword{datasets}
