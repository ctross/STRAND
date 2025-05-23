% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize_bsrm_results.R
\name{summarize_bsrm_results}
\alias{summarize_bsrm_results}
\title{Organize Stan output and provide summaries of model parameters}
\usage{
summarize_bsrm_results(input, include_samples = TRUE, HPDI = 0.9)
}
\arguments{
\item{input}{A STRAND model object, obtained by fitting a combined stochastic block and social relations model.}

\item{include_samples}{An indicator for the user to specify where raw samples, or only the summary statistics should be returned. Samples can take up a lot of space.}

\item{HPDI}{Highest Posterior Density Interval. Ranges in (0,1).}
}
\value{
A STRAND results object including summary table, a summary list, and samples.
}
\description{
This is a function to organize Stan output and provide summaries of key model parameters.
}
\examples{
\dontrun{
res = summarize_bsrm_results(input = fit)
}

}
