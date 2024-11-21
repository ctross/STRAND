#' Data on Baboon grooming, resting, and playing relationships
#'
#' Data on grooming, resting, and playing paired with additional covariate data. See the reference to the original study below for full methodological information.
#'
#' @docType data
#'
#' @usage data(Baboon_Longitudinal_Data)
#'
#' @format An object of class \code{"list"}. N individuals = 19. The list contains four tensors (of dimension 19x19X14) containing data on grooming, resting, and playing behaviours, hand-coded from focal follows for 28 days.  Each time-step is a two-day period.
#' An exposure tensor gives the count of times in which individual i could have been observed engaging in a given behaviour towards individual j on time-step t. Finally, there are covariate data on sex and age.  
#'
#' @keywords datasets
#'
#' @references Gelardi V, Godard J, Paleressompoulle D, Claidiere N, Barrat A. 2020.
#' Measuring social networks in primates: Wearable sensors versus direct observations.
#' Proc. R. Soc. A 476: 20190737.
#'
#' @source \href{http://dx.doi.org/10.1098/rspa.2019.0737}{RSPA}
#'
#' @examples
#' data(Baboon_Longitudinal_Data)
"Baboon_Longitudinal_Data"
