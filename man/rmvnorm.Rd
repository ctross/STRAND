% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmvnorm.R
\name{rmvnorm}
\alias{rmvnorm}
\title{A function to simulate MV norm data}
\usage{
rmvnorm(
  n,
  mean = rep(0, nrow(sigma)),
  sigma = diag(length(mean)),
  method = c("eigen", "svd", "chol"),
  pre0.9_9994 = FALSE,
  checkSymmetry = TRUE
)
}
\arguments{
\item{n}{number of samples}

\item{mean}{mean vector}

\item{sigma}{vcov matrix}

\item{method}{method to use}

\item{pre0.9_9994}{see mvtnorm}

\item{checkSymmetry}{see mvtnorm}
}
\value{
MV normal samples
}
\description{
See mvtnorm package for details
}
