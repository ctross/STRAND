#' A function to simulate MV Normal data
#'
#' See McElreaths Rethinking package for details
#'
#' @param n Number of samples.
#' @param Mu Mean vector.
#' @param sigma Variance vector.
#' @param Rho Correlation matrix.
#' @param method See details in Rethinking.
#' @return MV normal samples.
#' @export

rmvnorm2 = function (n, Mu = rep(0, length(sigma)), sigma = rep(1, length(Mu)), 
    Rho = diag(length(Mu)), method = "chol"){
    ldim = function(x) {
        z = length(dim(x))
        if (z == 0) 
            z = 1
        return(z)
    }
    if (ldim(Mu) > 1 || ldim(sigma) > 1 || ldim(Rho) > 2) {
        m = 0
        mR = 0
        K_Mu = 1
        if (ldim(Mu) == 2) {
            K_Mu = dim(Mu)[1]
            m = dim(Mu)[2]
        }
        else {
            m = length(Mu)
        }
        K_sigma = 1
        if (ldim(sigma) == 2) {
            K_sigma = dim(sigma)[1]
        }
        K_Rho = 1
        if (ldim(Rho) == 3) {
            K_Rho = dim(Rho)[1]
        }
        mR = dim(Rho)[2]
        K = n
        if (K_Mu < K) {
            Mu_new = matrix(NA, nrow = K, ncol = m)
            row_list = rep(1:K_Mu, length.out = K)
            if (K_Mu == 1) 
                for (i in 1:K) Mu_new[i, ] = Mu
            else for (i in 1:K) Mu_new[i, ] = Mu[row_list[i], 
                ]
            Mu = Mu_new
        }
        if (K_sigma < K) {
            sigma_new = matrix(NA, nrow = K, ncol = m)
            row_list = rep(1:K_sigma, length.out = K)
            for (i in 1:K) sigma_new[i, ] = sigma[row_list[i], 
                ]
            sigma = sigma_new
        }
        if (K_Rho < K) {
            Rho_new = array(NA, dim = c(K, m, m))
            row_list = rep(1:K_Rho, length.out = K)
            for (i in 1:K) Rho_new[i, , ] = Rho[row_list[i], 
                , ]
            Rho = Rho_new
        }
        result = array(NA, dim = c(K, m))
        for (i in 1:n) {
            DS = diag(sigma[i, ])
            SIGMA = DS %*% Rho[i, , ] %*% DS
            result[i, ] = rmvnorm(n = 1, mean = Mu[i, ], sigma = SIGMA, 
                method = method)
        }
        return(result)
    }
    else {
        DS = diag(sigma)
        SIGMA = DS %*% Rho %*% DS
        return(rmvnorm(n = n, mean = Mu, sigma = SIGMA, method = method))
    }
}
