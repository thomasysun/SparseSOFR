#----------------------------------------------------------------------------
# Generate simulated functional data
#----------------------------------------------------------------------------

#' Simulated Functional Data
#'
#' Generate simulated functional data.
#'
#' @param N number of subjects.
#' @param Tn total number of time points.
#' @param K number of basis functions
#' @param sig_noise observation level variance.
#' @param sig_alpha basis coefficients variance.
#' @param seed (optional) sets the seed.
#'
#' @return list with the following elements:
#' * X: simulated functional data with noise
#' * X_true: simulated true functional data
#' * B: true basis used to generate functions
#' @export
#'
#' @importFrom stats poly rnorm
sim_sf_data <- function(N, Tn, K = 4, sig_noise = .01, sig_alpha = 1, seed = NULL){

  tau <- seq(0, 1, by = 1/(Tn-1))

  B <- cbind(1/sqrt(Tn), poly(tau, K-1))

  set.seed(seed)
  alpha_true <- matrix(NA, nrow = K, ncol = N)

  for(i in 1:K){

    a <- rnorm(N, 0, sqrt(sig_alpha)/i)
    alpha_true[i,] <- a

  }
  alphaf_true <- B%*%alpha_true

  noise <- matrix(rnorm(N*Tn, mean = 0, sd = sqrt(sig_noise)), nrow = Tn, ncol = N)

  X_true <- alphaf_true

  X <- X_true + noise


  simdata <- list(X = X,
                  X_true = X_true,
                  B = B)

  return(simdata)

}
