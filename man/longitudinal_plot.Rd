% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/longitudinal_plot.R
\name{longitudinal_plot}
\alias{longitudinal_plot}
\title{Plot longitudinal random effects in time-lag form}
\usage{
longitudinal_plot(
  input,
  type = "dyadic",
  mode = "cor",
  parameter_set,
  results,
  normalized = FALSE,
  HPDI = 0.9,
  plot = TRUE,
  export_as_table = FALSE,
  save_plot = NULL,
  height = 6,
  width = 6,
  palette = NULL
)
}
\arguments{
\item{input}{A STRAND model object, obtained by fitting a longitudinal combined stochastic block and social relations model.}

\item{type}{Plot dyadic reciprocity ("dyadic") or generalized reciprocity ("generalized")?. Can use "custom", if the user passes in an appropriate data frame.}

\item{mode}{For dyadic plots, should the dyadic correlation "cor" be plotted, or the dyadic covariance "cov", or the adjusted dyadic+error correlation "adj".}

\item{parameter_set}{A labeled list of paramters to plot. E.g.: parameter_set = list(focal="Age", target="Age", focal="Sex", target="Sex",  dyadic = "Relatedness"). The 
names must be focal, target, or dyadic, and the quoted text must be paramters given in the regression equations. E.g., focal_regression = ~ Age + Sex.}

\item{results}{A results data.frame or matrix with colnames: c("Variable", "Layer", "Target" , "Base" , "Median", "L", "H", "Mean","SD","LayerNumeric")}

\item{normalized}{Should effects be normalized? Only valid for type="coefficient".}

\item{HPDI}{Highest Posterior Density Interval. Ranges in (0,1).}

\item{plot}{Should a plot be displayed?}

\item{export_as_table}{Should the tabular data rather than a ggplot object be returned?}

\item{save_plot}{Should a plot be exported to working directory? If so, set save_plot="desired_filename.pdf".}

\item{height}{Height of exported plot.}

\item{width}{Width of exported plot.}

\item{palette}{Override the default palette with a vector of color codes. "dyadic" need 3 colors, "generalized" needs 5 colors, and "coefficient" needs as many colors as supplied parameters.}
}
\value{
A figure or tabluar data to make a figure.
}
\description{
Plot longitudinal random effects in time-lag form
}
\examples{
\dontrun{
res = longitudinal_plot(input = fit)
}

}
