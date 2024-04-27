#----------------------------------------------------------------------------
# Plot output from ssofr object
#----------------------------------------------------------------------------

#' Plot output from ssofr object
#'
#' Plots summary of SOFR model or sparse curve fitting model using the MCMC draws from \code{ssofr()} or \code{sffm()}.
#'
#' @param x \code{ssofr} object
#' @param levels lower and upper quantiles for the credible intervals
#' @param which subset of \code{1:2} specifying which outputs to plot. \code{1} shows the regression coefficient function and \code{2} shows the fitted functional covariate curves.
#' @param ids integer vector specifying index of which fitted functional covariate curves to plot.
#'
#' @return Plots of the posterior mean and pointwise credible interval for the regression coefficient function, as well as the fitted curves if \code{fittedx} is specified.
#' @export plot.ssofr
#' @export
#'
#' @importFrom stats quantile
#' @importFrom graphics matplot polygon abline
#' @importFrom grDevices adjustcolor
#'
plot.ssofr <- function(x, levels = c(0.025, 0.975), which = 1:2, ids = 1:3){

  if('phi_post' %in% names(x)){
  phi_post <- x$phi_post
  }else{
  which = 2
  }
  betaf_post <- x$betaf_post
  Tau <- x$Tau

  if(1 %in% which & 'phi_post' %in% names(x)){
    phi_mean <- apply( phi_post, 1, mean )
    phi_l <- apply( phi_post , 1, function(z) quantile(z, levels[1]) )
    phi_u <- apply( phi_post , 1 , function(z) quantile(z, levels[2]) )

    matplot(x = Tau, cbind(phi_u,  phi_l), type = "l", col = "white", ylab = "Phi")
    polygon(x = c(Tau, rev(Tau)),
            y = c(phi_l,
                  rev(phi_u)),
            col =  adjustcolor(4, alpha.f = .3), border = NA)

    matplot(x = Tau, phi_mean, type = "l", lwd = 4, add=T, col = 4)
    abline(h=0)
  }

  if(2 %in% which){
    beta_mean <- rowMeans(betaf_post, dims = 2)
    beta_l <- apply( betaf_post , 1:2 , function(z) quantile(z, levels[1]) )
    beta_u <- apply( betaf_post , 1:2 , function(z) quantile(z, levels[2]) )

    matplot(x = Tau, beta_mean[,ids], type = "l", lwd = 2, lty = 1, ylab = "x")
    matplot(x = Tau, beta_l[,ids], type = "l", lwd = 1, lty = 2, add = T)
    matplot(x = Tau, beta_u[,ids], type = "l", lwd = 1, lty = 2, add = T)
  }
}
