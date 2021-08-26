\name{TRESS_DMRfit}
\alias{TRESS_DMRfit}
\title{
Model fitting and parameter estimation by TRESS for each candidate DMR.
}
\description{
TRESS models the read counts in candidate DMR using hierarchical negative
binomial distributions, with methylation level of each
DMR linked to multi-factors in the design by a linear framework.
This function conducts model fitting, parameter estimation,
and the variance-covariance matrix computation.
}
\usage{
TRESS_DMRfit(counts, sf,
             model, variable,
             iniby = "lm",
             shrkPhi = TRUE,
             addsuedo = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{
  A dataframe containing read counts in each candidate DMR
  from all samples.
}
\item{sf}{A vector containing size factor for each sample.}

\item{variable}{A dataframe containing all studied factors.}
\item{model}{ A formula to specify a particular design
based on the variables in ``variable".}

\item{iniby}{A character used for coefficients initialization.
There are three methods.
``lm": linear regression. ``Beta": beta regression.
``Beta-Bin": beta-binomial regression. Defult is ``lm".}

\item{shrkPhi}{A logical value to indicate whether conduct
shringkage estimate for dispersion parameter. Default is ``TRUE".}

\item{addsuedo}{A logical value to indicate whether or not adding a
psuedo count 5 on raw read counts. Default is ``FALSE". }
}

\value{
This function returns a list containing:
\item{Ratio}{A dataframe containing the IP/input ratio from all samples.}
\item{loglik}{A numerical vector containing the log-likelihood of all DMRs.}
\item{Coef}{A matrix containing estimates of coefficients in the design.}
\item{Cov}{A list of variance-covariance matrix estimates for all DMRs.}
\item{shrkPhi}{A vector containing shrinkage estimate for
dispersion of methylation.}
\item{shrkTheta}{A vector containing shrinkage estimate for
scale parameter in gamma model.}
\item{id.rm}{An numberical vector containing the index of
candidate DMRs whoes methylation level are smaller than 0.5
cross all samples.}
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
}
