\docType{data}
\name{DMR_M3vsWT}
\alias{DMR_M3vsWT}
\title{Transcriptome location and read counts of 1000 candidate DMRs,
and size factors.}
\format{
A list with 5 elements
\describe{
  \item{Regions}{A dataframe of 1000 obs and 5 variables,
  containing the transcriptome location for 1000 candidate
  DMRs.}
  \item{Counts}{A dataframe of 1000 obs and 8 variables,
  caontaining the read counts for candidate DMRs.
  The counts are from both Wild type and METTL3-Knockout
  samples. Each has two paired of IP and input replicates.}
  \item{sf}{A numerical vector,
  containing the size factors of each sample
  estimated from the whole transcriptome}
  \item{Bins}{A dataframe of 737 obs and 5 variables,
  containing the transcriptome location of bins
  overlapping with a subset of 1000 candidate DMRs.}
  \item{BinsCounts}{Read counts (1000-by-8) correspond to ``Bins``.}
  ...
}
}
\usage{
data(DMR_M3vsWT)
}
\description{
A dataset containing the transcriptome location,
read counts of 1000 candidate DMRs from 4 samples.
Each sample contains two paired of IP and input replicates.
It also contains size factors estimated from the whole
transcriptome to normalize sequencing depth of each sample.
In addition, a small proportion of transcriptome bins and
their read counts are also included for the convenience of
visualization.
}
\keyword{datasets}
