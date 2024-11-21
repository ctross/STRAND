#' A function to simulate MV norm data
#'
#' See mvtnorm package for details
#'
#' @param n number of samples
#' @param mean mean vector
#' @param sigma vcov matrix
#' @param method method to use
#' @param pre0.9_9994 see mvtnorm
#' @param checkSymmetry see mvtnorm
#' @return MV normal samples
#' @export

rmvnorm = function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
    method = c("eigen", "svd", "chol"), pre0.9_9994 = FALSE, 
    checkSymmetry = TRUE){
    if (checkSymmetry && !isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
        check.attributes = FALSE)) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != nrow(sigma)) 
        stop("mean and sigma have non-conforming size")
    method = match.arg(method)
    R = if (method == "eigen") {
        ev = eigen(sigma, symmetric = TRUE)
        if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
            warning("sigma is numerically not positive semidefinite")
        }
        t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 
            0))))
    }
    else if (method == "svd") {
        s. = svd(sigma)
        if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {
            warning("sigma is numerically not positive semidefinite")
        }
        t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
    }
    else if (method == "chol") {
        R = chol(sigma, pivot = TRUE)
        R[, order(attr(R, "pivot"))]
    }
    retval = matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*% 
        R
    retval = sweep(retval, 2, mean, "+")
    colnames(retval) = names(mean)
    retval
}
