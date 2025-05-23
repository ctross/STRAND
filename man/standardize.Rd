% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardize.R
\name{standardize}
\alias{standardize}
\title{A function to standardize data}
\usage{
standardize(input, center = TRUE, type = "sd")
}
\arguments{
\item{input}{A vector of data to standardize.}

\item{center}{Subtract mean?}

\item{type}{Divide by "max" or "sd"?}
}
\value{
A vector of standardize data.
}
\description{
This is a simple helper function to standardize data. We recommmend to standardize covariates as it makes parameter interpretation easier.
}
\examples{
\dontrun{
age_standardize = standardize(input=age)
}

}
