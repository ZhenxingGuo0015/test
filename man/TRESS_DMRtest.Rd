\name{TRESS_DMRtest}
\alias{TRESS_DMRtest}
\title{
Statistical test for DMR in TRESS.
}
\description{
This function conducts statistical test based on
the user specified contrast of coefficients from the design.
If the contrast is a vector, then TRESS performs Wald test;
if the contrast is a matrix, then TRESS conducts Chi-Squared test.
}
\usage{
TRESS_DMRtest(DMR, contrast, nullModel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DMR}{A list containing IP/input ratio,
  the coefficients estimate,
  variance-covariance estimate.
  This can be obtained from the ouput of ``TRESS_DMRfit".
  }
  \item{contrast}{A contrast for all coefficients in the design.
  It can be either a (p+1) vector or a m-by-(p+1) matrix,
  where p is the number of columns in the design.
  ``m" depends on the number of relationships that users want to test.}

  \item{nullModel}{A character to specify a method to calculate
  p-value based on the statistics. It can be "standN", "2mix",
  "trunN" and "locfdr" for standard normal,
  two-component mixed gaussian, truncated normal, local fdr.
  If it is "locfdr", then there would be only "fdr" and no p-values.
  Defult is "standN".}
}
\value{
This function returned a dataframe containing the testing
results for specified contrast. The columns are
\item{baseMean}{Averaged methylation level cross all samples.}
\item{logOR}{Estimated value of contrast: ``C^T*R".}
\item{lorSE}{Standard error of ``C^T*R".}
\item{stat}{Statistics, basically logOR/lorSE.}
\item{pvalue}{P-values from Wald test.}
\item{padj}{Benjamini-Hochberg adjusted p-values.}
}

\examples{
# A toy example
data(DMR_M3vsWT) # data from TRESS
variable = data.frame(predictor = rep(c("WT", "M3"), c(2, 2)))
model = ~1+predictor
DMR.fit = TRESS_DMRfit(counts = DMR_M3vsWT$Counts,
                       sf = DMR_M3vsWT$sf,
                       variable = variable,
                       model = model)
DMR.test = TRESS_DMRtest(DMR = DMR.fit, contrast = c(0, 1),
                         nullModel = "standN")
head(DMR.test)
head(DMR_M3vsWT$Regions[which(DMR.test$padj < 0.05), ])
}

