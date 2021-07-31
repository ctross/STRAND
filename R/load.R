#' A function to load some data when the package opens.
#'
#' @export

.onAttach <- function(libname=NULL, pkgname="STRAND") {
  packageStartupMessage("Bei diesem Spaziergang an den STRAND schärfen wir unsere Sinne für die Sternbilder hoch am Himmel!")

}