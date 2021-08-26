\name{TRESS_DMRCandidates}
\alias{TRESS_DMRCandidates}
\title{
Candidate m6A region identification.
}
\description{
This step performs binomial test and then loosely
identify candidate DMR regions
using an ad hoc algorithm.
}
\usage{
TRESS_DMRCandidates(Counts, bins,
WhichThreshold ="fdr_lfc", pval.cutoff = 1e-5,
fdr.cutoff = 0.05, lfc.cutoff = 0.7,
windlen = 5,lowcount = 30)
}

\arguments{
\item{Counts}{ A data matrix containing
  bin-level (default 50bp) read counts in both
  IP and input samples,
  where the sample order is: input1, ip1, input2, ip2, ...
  }
  \item{bins}{
  A data frame containing the genomic coordinate of
  each bin of fixed length.
  }
  \item{WhichThreshold}{
  A character specifying a threshold for significant bins
  in bump finding using an ad hoc algorithm.
  There are three options: "pval" (only use p-values),
  "fdr" (only use FDR), "lfc" (only use log fold change),
  "pval_lfc" (use both p-values and log fold changes)
  and "fdr_lfc" (use FDR and log fold changes). Default is "fdr_lfc".
  }
  \item{pval.cutoff}{
  A constant indicating the cutoff for p-value. Default is 1e-05.
  }
  \item{fdr.cutoff}{
  A constant indicating the cutoff for FDR. Default is 0.05.
  }
  \item{lfc.cutoff}{
  A constant indicating the cutoff for
  log fold change. Default is 0.7 for fold change of 2.
  }
  \item{windlen}{
  An integer specifying the length of
  consecutive bins used in simple
  moving average smooth of log fold change. Default is 5.
  }
  \item{lowcount}{
  An iteger to filter out candidate regions with lower read counts
  in input. Default is 30.
  }
}

\details{
  The function involves three steps:
  1. Perform binomial test for each bin based bin-level counts.
  2. Use ad hoc algorithm to combine significant bins to form
  candidate regions.
  }

\value{
A list containing
  \item{Regions }{A data frame containng genomic coordinate
  for each candidate region.}
  \item{Counts}{A data matrix containing read counts of
  all samples for each candidate region.}
}



%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
### A toy example, whose results do not have real applications.
data("Basal_binlevel")
Candidates = TRESS_DMRCandidates(Counts = Basal_binlevel$bincounts,
bins = Basal_binlevel$bins)

}
