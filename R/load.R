#' A function to load some data when the package opens.
#'
#' @param
#' libname Library name.
#' @param
#' pkgname Package name.
#' @export

.onAttach <- function(libname=NULL, pkgname="STRAND") {
  packageStartupMessage("Bei diesem Spaziergang an den STRAND scharfen wir unsere Sinne fur die Sternbilder hoch am Himmel!")
}
