% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_strand_data.R
\name{make_strand_data}
\alias{make_strand_data}
\title{A function to create a STRAND data object}
\usage{
make_strand_data(
  outcome = NULL,
  self_report = NULL,
  outcome_mode = NULL,
  link_mode = NULL,
  ground_truth = NULL,
  block_covariates = NULL,
  individual_covariates = NULL,
  dyadic_covariates = NULL,
  exposure = NULL,
  m_e_data = NULL,
  mask = NULL,
  diffusion_outcome = NULL,
  diffusion_exposure = NULL,
  diffusion_mask = NULL,
  multiplex = FALSE,
  longitudinal = FALSE,
  directed = TRUE,
  imputation = FALSE,
  check_data_organization = TRUE,
  check_standardization = TRUE,
  standardization_threshold = c(0.5, 2)
)
}
\arguments{
\item{outcome}{A named list of primary network data (e.g., self reports). Each entry in the list must be an adjacency matrix. This will be a list of length 1 for single-sampled networks
and a list of length 2 for double-sampled networks. For douple-sampled networks, data is presumed to be organized such that outcome[[1]][i,j] represents i's reports of transfers from i to j, and outcome[[2]][i,j]
represents i's reports of transfers from j to i. Data should be binary, 0 or 1, unless an alternative outcome_mode is provided. If outcome_mode="poisson", then data can be integer values.
If an exposure variable is provided, outcome can take integer values and outcome_mode="binomial" can be set. For multiplex models, the list can be longer than 2, but all layers must be single sampled, and you must set multiplex = TRUE.}

\item{self_report}{A named list of primary network data (e.g., self reports). This is a deprecated alias for outcome above.}

\item{outcome_mode}{Can be either "bernoulli", "binomial", "poisson", or "gaussian" based on the kind of network data being modeled.}

\item{link_mode}{Can be either "logit", "probit", "log", or "identity"; "log" must be used with "poisson" outcomes; "identity" must be used with "gaussian" outcomes; "logit" is default for Bernoulli and Binomial outcomes; "probit" is basically the same as "logit", but the model fits slower.}

\item{ground_truth}{A list of secondary network data about equivalent latent relationships (i.e., from focal observations). Each entry in the list must be an adjacency matrix. 
Data is presumed to be organized such that ground_truth[[t]][i,j] represents observed transfers from i to j at time-point t.}

\item{block_covariates}{A data.frame of group IDs (e.g., ethnicity, class, religion, etc.) corresponding to the individuals in the 'self_report' network(s). Each variable should be provided as a factor.}

\item{individual_covariates}{An N_id by N_parameters dataframe of all individual-level covariates that are to be included in the model.}

\item{dyadic_covariates}{A named list of N_id by N_id by N_dyadic_parameters matrices.}

\item{exposure}{A named list of matrices matched to the self_report matrices. If self_report is a count data set with binomial outcomes, then this variable holds the sample size information.}

\item{m_e_data}{A list of integer vectors: list(sampled=sampled, sampled_exposure=sampled_exposure, sampled_mask=sampled_mask, detected=detected, detected_exposure=detected_exposure, detected_mask=detected_mask), to be used in measurement error models.}

\item{mask}{A list of matrices matched to the self_report matrices. If mask[i,j,m]==0, then ties between i and j in layer m are detectable. If mask[i,j,m]==1, then i to j ties in layer m are censored (e.g., if i and j were monkeys kept in different enclosures).}

\item{diffusion_outcome}{An N-vector of outcome data for a trait diffusing over a network.}

\item{diffusion_exposure}{An N-vector matched with the diffusion_outcome matrix. If diffusion_outcome is a count data set with binomial outcomes, then this variable holds the sample size information.}

\item{diffusion_mask}{An N-vector of indicators for diffusion outcomes that were masked.}

\item{multiplex}{If TRUE, then all layers in outcome are modeled jointly.}

\item{longitudinal}{If TRUE, then checks for longitudinal data structure are performed.}

\item{directed}{If TRUE, then STRAND will treat the outcomes as directed. If set to FALSE, STRAND will treat the outcomes as undirected; this leads to some addition checks on model definition.}

\item{imputation}{If TRUE, then checks for NAs in data are omitted, and supported models will impute the missings.}

\item{check_data_organization}{If TRUE, then checks that all colnames and rownames match. This will catch missorted data.}

\item{check_standardization}{If TRUE, then checks that all individual and dyadic variables are standardized to have SDs in the range of standardization_threshold. Standardization is important in STRAND so that priors have equal strength across all predictors.}

\item{standardization_threshold}{If check_standardization is TRUE, then individual and dyadic predictors must have SDs in this range.}
}
\value{
A list of data formatted for use by STRAND models.
}
\description{
This function organizes network data and covariates into a form that can be used by STRAND for model fitting. All 
STRAND model fitting functions require their data to be supplied in the format exported here.
}
\examples{
\dontrun{
model_dat = make_strand_data(self_report=LoanData)
}

}
