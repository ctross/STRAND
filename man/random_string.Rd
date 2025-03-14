% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/random_string.R
\name{random_string}
\alias{random_string}
\title{A function to make random ID codes}
\usage{
random_string(N = 1, length = 12)
}
\arguments{
\item{N}{Number of ID codes to generate.}

\item{length}{Length of ID codes.}
}
\value{
A vector of unique ID codes.
}
\description{
This is a small helper function to create random ID codes.
}
\examples{
\dontrun{
IDs = random_string(N=100, length=5)
}

}
