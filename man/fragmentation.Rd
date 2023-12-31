% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fragmentation.r
\name{fragmentation}
\alias{fragmentation}
\title{fragmentation: Conveniently wraps all fragmentation steps}
\usage{
fragmentation(data, pen_HL, pen_int, cores = 2)
}
\arguments{
\item{data}{data frame: data frame combined data by column}

\item{pen_HL}{list: list of the penalties set optimal for the fragmentation
for half-life}

\item{pen_int}{list: list of the penalties set optimal for the fragmentation
for intensity}

\item{cores}{integer: the number of assigned cores for the task. It needs to
be increased in case of big data.}
}
\value{
Two data frames with half-life and intensity fragments and the mean
of the coefficient fragment based.
}
\description{
fragmentation fragments the half-life and intensity into segments using
the penalties output.
}
\examples{
data(penalties_df)
data(pen_HL)
data(pen_int)
df_comb_minimal <- fragmentation(data = penalties_df, pen_HL, 
pen_int)
}
