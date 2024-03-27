
psBasis <- function(x, K = length(x), spline.degree = 3, diff.ord = 2,
                    knots = NULL) {
  if (is.null(knots)) {
    knots.no <- K - spline.degree + 1
    xl <- min(x)
    xr <- max(x)
    xmin <- xl - (xr - xl) / 100
    xmax <- xr + (xr - xl) / 100
    dx <- (xmax - xmin) / (knots.no - 1)
    knots <- seq(xmin - spline.degree * dx, xmax + spline.degree * dx, by = dx)
  }
  X <- splines::spline.des(knots, x, spline.degree + 1, outer.ok = TRUE)$design
  P <- diag(K) # precision
  if (diff.ord > 0) {
    for (d in 1:diff.ord) P <- diff(P)
    #P <- crossprod(P)
  }
  return(list(
    X = X, P = P, knots = knots, K = K, spline.degree = spline.degree,
    diff.ord = diff.ord
  ))
}

fdlm_init = function(Y, tau, K = NULL, use_pace = FALSE){

  # Convert to matrix, if necessary:
  tau = as.matrix(tau)

  # Rescale observation points to [0,1]
  tau01 = apply(tau, 2, function(x) (x - min(x))/(max(x) - min(x)))

  # And the dimensions:
  T = nrow(Y); m = ncol(Y); d = ncol(tau)

  # Compute basic quantities for the FLC splines:
  splineInfo = getSplineInfo_d(tau = tau01,
                               m_eff = floor(median(rowSums(!is.na(Y)))),
                               orthonormalize = TRUE)

  # For initialization: impute
  Y0 = matrix(NA, nrow = T, ncol = m) # Storage for imputed initialization data matrix
  allMissing.t = (rowSums(!is.na(Y))==0)   # Boolean indicator of times at which no points are observed

  # Use PACE, or just impute w/ splines:
  if((d==1) && use_pace && any(is.na(Y)) ){
    # To initialize, use FPCA via PACE:
    fit_fpca = fdapace::FPCA(Ly =  apply(Y, 1, function(y) y[!is.na(y)]),
                    Lt = apply(Y, 1, function(y) tau01[!is.na(y)]),
                    optns = list(dataType = 'Sparse', methodSelectK = K))
    # Fitted FPCA curves:
    Yhat0 = fitted(fit_fpca); t0 = fit_fpca$workGrid

    # For all times at which we observe a curve, impute the full curve (across tau)
    Y0[!allMissing.t,] = t(apply(Yhat0[!allMissing.t,], 1, function(x) splinefun(t0, x, method='natural')(tau01)))

  } else {
    # For all times at which we observe a curve, impute the full curve (across tau)
    Y0[!allMissing.t,] = t(apply(Y[!allMissing.t,], 1, function(x) splinefun(tau01, x, method='natural')(tau01)))
  }

  # Next, impute any times for which no curve is observed (i.e., impute across time)
  Y0 = apply(Y0, 2, function(x){stats::splinefun(1:T, x, method='natural')(1:T)})

  # Compute SVD of the (completed) data matrix:
  # (Complete) data matrix, projected onto basis:
  YB0 = Y0%*%splineInfo$Bmat%*%chol2inv(chol(splineInfo$BtB))
  singVal = svd(YB0)

  # If K is unspecified, select based on cpv
  if(is.null(K)){
    # Cumulative sums of the s^2 proportions (More reasonable when Y has been centered)
    cpv = cumsum(singVal$d^2/sum(singVal$d^2))
    K = max(2, which(cpv >= 0.99)[1])
  }

  # Check to make sure K is less than the number of basis coefficients!
  if(K >= ncol(YB0)){
    warning(paste("K must be less than the number of basis functions used; reducing from K =",K, "to K =", ncol(YB0) - 1))
    K = ncol(YB0) - 1
  }

  # Basis coefficients of FLCs:
  Psi0 = as.matrix(singVal$v[,1:K])

  # Initial FLCs:
  F0 = splineInfo$Bmat%*%Psi0

  # Factors:
  Beta0 = as.matrix((singVal$u%*%diag(singVal$d))[,1:K])

  # Initialize all curves to have positive sums (especially nice for the intercept)
  negSumk = which(colSums(F0) < 0); Psi0[,negSumk] = -Psi0[,negSumk]; Beta0[,negSumk] = -Beta0[,negSumk]

  list(Beta = Beta0, Psi = Psi0, splineInfo = splineInfo, Y0 = Y0)
}



getSplineInfo_d = function(tau, m_eff = NULL, orthonormalize = TRUE){

  # Just in case, reform as matrix
  tau = as.matrix(tau)

  # Number of observation points
  m = nrow(tau)

  # Dimension:
  d = ncol(tau)

  # Order of derivative in penalty:
  m_deriv = 2

  # This is the linear component
  X = cbind(1, tau)

  # Number of effective observation points:
  if(is.null(m_eff)) m_eff = m

  # Number of knots: if more than 25 effective observation points, likely can use fewer knots
  if(m_eff > 25){
    # Guaranteed to be between 20 and 150 knots (but adjust as needed)
    num_knots = max(20, min(ceiling(m_eff/4), 150))
  } else num_knots = max(3, m_eff)

  # Locations of knots:
  if(num_knots < m){
    # Usual case: fewer knots than TOTAL observation points
    if(d == 1){
      # d = 1, just use quantiles of the observed data points:
      knots = as.matrix(quantile(unique(tau), seq(0,1,length=(num_knots+2))[-c(1,(num_knots+2))]))
    } else {
      # d > 1, use space-filling algorithm:
      knots = fields::cover.design(tau, num_knots)$design
    }
  } else knots = tau

  # For the penalty matrix, need to compute distances between obs. points and knots
  dist.mat <- matrix(0, num_knots, num_knots); dist.mat[lower.tri(dist.mat)] <- stats::dist(knots); dist.mat <- dist.mat + t(dist.mat)
  if(d%%2 == 0){
    # Even dim:
    Omega = dist.mat^(2*m_deriv - d)*log(dist.mat)
  } else {
    # Odd dim:
    Omega = dist.mat^(2*m_deriv - d)
  }
  # For numerical stability:
  diag(Omega) = 0

  # Compute the "random effects" matrix
  Zk = matrix(0, nrow=m, ncol=num_knots)
  for (k in 1:num_knots){
    di = sqrt(rowSums((tau - matrix(rep(knots[k,], each = m), nrow=m))^2)) # di = 0; for(j in 1:d) di = di + (tau[,j] - knots[k,j])^2; di = sqrt(di)
    if(d%%2 == 0){# Even dim:
      Zk[,k] = di^(2*m_deriv - d)*log(di)
    } else { # Odd dim:
      Zk[,k] = di^(2*m_deriv - d)
    }
  }
  Zk[is.nan(Zk)] = 0

  # Natural constraints, if necessary:
  if(num_knots > m - 1){Q2 = qr.Q(qr(X), complete=TRUE)[,-(1:2)]; Zk = Zk%*%Q2; Omega = crossprod(Q2, Omega)%*%Q2}

  # SVD of penalty matrix
  # So that the "random effects" have diagonal prior variance
  svd.Omega = svd(Omega)
  sqrt.Omega = t(svd.Omega$v %*%(t(svd.Omega$u)*sqrt(svd.Omega$d)))
  Z = t(solve(sqrt.Omega,t(Zk)))

  # Now combine the linear and nonlinear pieces to obtain the matrix of basis functions evaluated at the obs. points
  Bmat = cbind(X, Z);

  # The penalty matrix:
  Omega = diag(c(rep(0, ncol(X)), rep(1, ncol(Z))))

  if(orthonormalize){
    # QR decomposition:
    qrd = qr(Bmat, complete = TRUE);  R.t = t(qr.R(qrd));
    # Update hte basis and the penalty matrix:
    Bmat = qr.Q(qrd); Omega = forwardsolve(R.t, t(forwardsolve(R.t, Omega, upper.tri = FALSE)), upper.tri = FALSE)

    BtB = diag(1, ncol(Bmat))
  } else BtB = crossprod(Bmat)

  # Return the matrix, the penalty, and the cross product (of the basis)
  list(Bmat = Bmat, Omega = Omega, BtB = BtB)
}


getEffSize = function(postX) {
  if(is.null(dim(postX))) return(coda::effectiveSize(postX))
  summary(coda::effectiveSize(coda::as.mcmc(array(postX, c(dim(postX)[1], prod(dim(postX)[-1]))))))
}

uni.slice <- function (x0, g, w = 1, m = Inf, lower = -Inf, upper = +Inf,
                       gx0 = NULL)
{
  if (!is.numeric(x0) || length(x0) != 1 || !is.function(g) ||
      !is.numeric(w) || length(w) != 1 || w <= 0 || !is.numeric(m) ||
      !is.infinite(m) && (m <= 0 || m > 1e+09 || floor(m) !=
                          m) || !is.numeric(lower) || length(lower) != 1 ||
      x0 < lower || !is.numeric(upper) || length(upper) !=
      1 || x0 > upper || upper <= lower || !is.null(gx0) &&
      (!is.numeric(gx0) || length(gx0) != 1)) {
    stop("Invalid slice sampling argument")
  }
  if (is.null(gx0)) {
    gx0 <- g(x0)
  }
  logy <- gx0 - rexp(1)
  u <- runif(1, 0, w)
  L <- x0 - u
  R <- x0 + (w - u)
  if (is.infinite(m)) {
    repeat {
      if (L <= lower)
        break
      if (g(L) <= logy)
        break
      L <- L - w
    }
    repeat {
      if (R >= upper)
        break
      if (g(R) <= logy)
        break
      R <- R + w
    }
  }
  else if (m > 1) {
    J <- floor(runif(1, 0, m))
    K <- (m - 1) - J
    while (J > 0) {
      if (L <= lower)
        break
      if (g(L) <= logy)
        break
      L <- L - w
      J <- J - 1
    }
    while (K > 0) {
      if (R >= upper)
        break
      if (g(R) <= logy)
        break
      R <- R + w
      K <- K - 1
    }
  }
  if (L < lower) {
    L <- lower
  }
  if (R > upper) {
    R <- upper
  }
  repeat {
    x1 <- runif(1, L, R)
    gx1 <- g(x1)
    if (gx1 >= logy)
      break
    if (x1 > x0) {
      R <- x1
    }
    else {
      L <- x1
    }
  }
  attr(x1, "log.density") <- gx1
  return(x1)
}

fmse <- function(alpha_true, alpha_hat){

  sum((alpha_true - alpha_hat)^2)/length(alpha_true)

}

mciw <- function(u, l){
  mean(u - l)

}

ecp <- function(u, l, alpha){

  mean(alpha < u & alpha > l)


}

rowrep <- function(X, ntimes){
  #as.matrix(as.data.frame(lapply(as.data.frame(X), rep, ntimes)))

  X[rep(seq_along(ntimes), ntimes), ]
}


sampleMGP = function(theta.jh, delta.h, a1 = 2, a2 = 3){

  # Just in case:
  theta.jh = as.matrix(theta.jh)

  # Store the dimensions locally
  p = nrow(theta.jh); K = ncol(theta.jh)

  # Sum over the (squared) replicates:
  sum.theta.l = colSums(theta.jh^2)

  # h = 1 case is separate:
  tau.not.1 = cumprod(delta.h)/delta.h[1]
  delta.h[1] = rgamma(n = 1, shape = a1 + p*K/2,
                      rate = 1 + 1/2*sum(tau.not.1*sum.theta.l))
  # h > 1:
  if(K > 1){for(h in 2:K){
    tau.not.h = cumprod(delta.h)/delta.h[h]
    delta.h[h] = rgamma(n = 1, shape = a2 + p/2*(K - h + 1),
                        rate = 1 + 1/2*sum(tau.not.h[h:K]*sum.theta.l[h:K]))
  }}
  delta.h #list(tau.h = cumprod(delta.h), delta.h = delta.h)
}

sample_lambda = function(lambda, Psi, Omega = NULL, d = 1, uniformPrior = TRUE, orderLambdas = TRUE){
  Lm = nrow(Psi); K = ncol(Psi)

  if(uniformPrior){shape0 = (Lm - d + 1 + 1)/2} else shape0 = (Lm - d - 1)/2 + 0.001; # for Gamma(0.001, 0.001) prior
  #if(uniformPrior){shape0 = (Lm + 1)/2} else shape0 = (Lm - 2)/2 + 0.001; # for Gamma(0.001, 0.001) prior

  for(k in 1:K){
    if(is.null(Omega)){rate0 = crossprod(Psi[-(1:(d+1)),k])/2} else rate0 = crossprod(Psi[,k], Omega)%*%Psi[,k]/2
    #if(is.null(Omega)){rate0 = crossprod(Psi[-(1:2),k])/2} else rate0 = crossprod(Psi[,k], Omega)%*%Psi[,k]/2

    if(!uniformPrior) rate0 = rate0 + 0.001  # for Gamma(0.001, 0.001) prior

    # Lower and upper bounds, w/ ordering constraints (if specified):
    if(orderLambdas){
      lam.l = 10^-8; lam.u = Inf; if(k != 1) lam.u = lambda[k-1];  # if(k != K) lam.l = lambda[k+1];
      lambda[k] = truncdist::rtrunc(1, 'gamma', a=lam.l, b=lam.u, shape=shape0, rate=rate0) # more stable, possibly faster
    } else lambda[k] = rgamma(1, shape = shape0, rate = rate0)
  }
  lambda
}

#' Compute Simultaneous Credible Bands
#'
#' Compute (1-alpha)\% credible BANDS for a function based on MCMC samples using Crainiceanu et al. (2007)
#'
#' @param sampFuns \code{Nsims x m} matrix of \code{Nsims} MCMC samples and \code{m} points along the curve
#' @param alpha confidence level
#'
#' @return \code{m x 2} matrix of credible bands; the first column is the lower band, the second is the upper band
#'
#' @note The input needs not be curves: the simultaneous credible "bands" may be computed
#' for vectors. The resulting credible intervals will provide joint coverage at the (1-alpha)%
#' level across all components of the vector.
#'
#' @export

credBands = function(sampFuns, alpha = .05){

  N = nrow(sampFuns); m = ncol(sampFuns)

  # Compute pointwise mean and SD of f(x):
  Efx = colMeans(sampFuns); SDfx = apply(sampFuns, 2, sd)

  # Compute standardized absolute deviation:
  Standfx = abs(sampFuns - tcrossprod(rep(1, N), Efx))/tcrossprod(rep(1, N), SDfx)

  # And the maximum:
  Maxfx = apply(Standfx, 1, max)

  # Compute the (1-alpha) sample quantile:
  Malpha = quantile(Maxfx, 1-alpha)

  # Finally, store the bands in a (m x 2) matrix of (lower, upper)
  t(cbind(Efx - Malpha*SDfx, Efx + Malpha*SDfx))
}


sampleFLCs_sofr <- function(bzzx, ywa, btil, Beta, Psi, bzzb, Omega, lambda, sig_ex, sig_ey, mu, theta){

  p = nrow(Psi)
  n = nrow(Beta)
  K = ncol(Psi)
  psi_k_0 <- matrix(NA, p, K)
  psi_k <- matrix(NA, p, K)
  nc <- rep(NA, K)

  for(j in 1:K){

    Q_psi <- tcrossprod(btil*t(Beta[,j]^2)[rep_len(1,p),])/(sig_ey) + Reduce('+', Map('*', bzzb, Beta[,j]^2))/sig_ex + lambda[j]*Omega

    psibeta <- matrix(rowSums((t(matrix(rep(t(Psi), n), K, n*p))*rowrep(Beta, rep(p, n)))[,-j]), p, n)

    l_psi0 <- matrix(NA, p, n)

    for(i in 1:n){

      l_psi0[,i] <- Beta[i,j]*btil*(ywa[i] - sum(btil*psibeta[,i]))*(1/(sig_ey)) + Beta[i,j]*(bzzx[[i]] - bzzb[[i]]%*%psibeta[,i])*(1/sig_ex)

    }

    l_psi <- rowSums(l_psi0)

    ch_Q <- chol(matrix(Q_psi, nrow=p, ncol=p))
    psi_k_0[, j] <- backsolve(ch_Q,
                              forwardsolve(t(ch_Q), l_psi) +
                                rnorm(p))


    Lcon <- Psi[,-j]

    Con <- backsolve(ch_Q,
                     forwardsolve(t(ch_Q), Lcon))

    psi_k[, j] <- psi_k_0[, j] - Con%*%solve(crossprod(Lcon, Con))%*%crossprod(Lcon, psi_k_0[, j])

    nc[j] <- sqrt(sum(psi_k[, j]^2))
    psi_k[, j] <- psi_k[, j]/nc[j]
    Psi[, j] <- psi_k[, j]
    Beta[, j] <- Beta[, j]*nc[j]
    mu[j] <- mu[j]*nc[j]
    theta[j,] <- theta[j,]*nc[j]

  }

  return(list(Psi, nc, Beta, mu, theta))

}

sampleFLCs <- function(BtY, Beta, Psi, BtB, Omega, lambda, sigma_e){

  Lm = nrow(Psi)
  n = nrow(Beta)
  K = ncol(Psi)
  psi_k_0 <- matrix(NA, Lm, K)
  psi_k <- matrix(NA, Lm, K)
  nc <- rep(NA, K)

  m <- list()
  m <- lapply(seq_len(n), function(X) Psi)

  for(j in 1:K){

    Q_psi <- Reduce('+', Map('*', BtB, Beta[,j]^2))/sigma_e + lambda[j]*Omega


    psibeta <- matrix(rowSums((t(matrix(rep(t(Psi), n), K, n*Lm))*rowrep(Beta, rep(Lm, n)))[,-j]), Lm, n)

    l_psi0 <- matrix(NA, Lm, n)

    for(i in 1:n){

      l_psi0[,i] <- Beta[i,j]*(BtY[[i]] - BtB[[i]]%*%psibeta[,i])

    }

    l_psi <- rowSums(l_psi0)*(1/sigma_e)

    ch_Q <- chol(matrix(Q_psi, nrow=Lm, ncol=Lm))
    psi_k_0[, j] <- backsolve(ch_Q,
                              forwardsolve(t(ch_Q), l_psi) +
                                rnorm(Lm))


    Lcon <- Psi[,-j]

    Con <- backsolve(ch_Q,
                     forwardsolve(t(ch_Q), Lcon))

    psi_k[, j] <- psi_k_0[, j] - Con%*%solve(crossprod(Lcon, Con))%*%crossprod(Lcon, psi_k_0[, j])

    nc[j] <- sqrt(sum(psi_k[, j]^2))
    psi_k[, j] <- psi_k[, j]/nc[j]
    Psi[, j] <- psi_k[, j]
    Beta[, j] <- Beta[, j]*nc[j]

  }

  return(list(Psi, nc, Beta))

}
