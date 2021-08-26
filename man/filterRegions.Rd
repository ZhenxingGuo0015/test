\name{filterRegions}
\alias{filterRegions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter out regions with small marginal variation in methylation.
}
\description{
This function filter out regions who have smaller
(less than quantile of ``quant")
coefficient of variation (CV) in methylation ratio.
}
\usage{
filterRegions(Candidates, quant = 0.25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Candidates}{
A list containing "Counts", ``Regions", "sf" for read counts,
genomic coordinate of each candidate DMR
and size factor of all samples.
}
\item{quant}{A percentage to specify a quantile.
Regions with CV smaller than this quantile would be filtered out.}
}

\value{
The same structure with the input ``candidate" but fewer regions.
}

\examples{
# A toy example
data(DMR_M3vsWT) # data from TRESS
sub.DMR_M3vsWT = filterRegions(DMR_M3vsWT)
}

