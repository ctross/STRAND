#' Data on RICH games from Colombia
#'
#' Data on network-structured economic games with additional covariate data. These data have been rescaled, and some missing values were imputed using the mode or median. See the reference to the original study below for full methodological information.
#'
#' @docType data
#'
#' @usage data(RICH_Data)
#'
#' @format An object of class \code{"list"}. N individuals = 93. The list contains three outcome matrices (of dimension 93 x 93) containing data on giving, taking, and reducing behaviours. 
#'
#' @keywords datasets
#'
#' @references Ross and Pisor, 2024.
#' Perceived inequality and variability in the expression of parochial altruism
#'
#' @source \href{http://dx.doi.org/}{Working paper}
#'
#' @examples
#' data(RICH_Data)
"RICH_Data"
