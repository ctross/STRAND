#' A function to compute the CES function
#'
#' See Constant elasticity of substitution from Wikipedia
#'
#' @param K Value of input 1.
#' @param L Value of input 2.
#' @param alpha Share parameter of input 1.
#' @param sigma Elasticity of substitution.
#' @param eta Returns to scale.
#' @param nudge If TRUE, then nudge sigma such that sigma=1 is interpetable.
#' @return CES(K,L)
#' @export

CES = function(K, L, alpha, sigma, eta, nudge=TRUE){
    if(length(K) != length(L)) stop("K and L must be of same length.")
    for(i in 1:length(K)){
     if(K[i] < 0) stop("K must be >= 0.")
     if(L[i] < 0) stop("L must be >= 0.")
     }
    if(alpha < 0 | alpha > 1) stop("alpha in (0,1).")
    if(sigma < 0) stop("sigma must be > 0.")
    if(eta < 0) stop("eta must be > 0.")
    if(nudge == TRUE) sigma = (3.1415/3.1415926)*sigma

    rho = (sigma - 1)/sigma
    p = (alpha*K^rho + (1-alpha)*L^rho)^(eta/rho)
    
    return(p)
}
