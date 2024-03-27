#----------------------------------------------------------------------------
# Sparse scalar-on-function regression
#----------------------------------------------------------------------------

#' Sparse scalar-on-function regression
#'
#' Performs MCMC estimation for Bayesian scalar-on-function regression with sparse, noisy functional covariate.
#'
#' @param y vector of scalar responses.
#' @param x matrix of functional responses. Each column corresponds to one function and each row is a measurement at a single time point on a common grid. May contain NAs.
#' @param w design matrix of scalar covariates using `model.matrix`. It should include column of 1s for the intercept.
#' @param Tau vector specifying common grid of time points on which function is observed.
#' @param K number of basis functions in factor model.
#' @param Kphi number of basis functions for functional coefficient.
#' @param S number of total MCMC iterations.
#' @param S_burn number of initial MCMC iterations to discard for warm-up.
#' @param showprogress show progress of MCMC algorithm. Defaults to TRUE.
#' @param sparse set to TRUE if functions considered sparsely observed.
#'
#' @return A list with the following elements:
#' * X: the inputted design matrix
#' @export
#'
#' @importFrom stats poly rgamma rnorm dgamma dunif fitted median quantile rexp runif sd splinefun
ssofr <- function(y, x, w = NULL, Tau, K = 10, Kphi = 15, S = 2000, S_burn = 1000, sparse = TRUE, showprogress = TRUE){

  n <- length(y)
  Tn <- nrow(x)
  Tni <- colSums(!is.na(x))

  if(is.null(w) == TRUE){
    w = rep(1, n)
    L = 0
  }else{
    L <- ncol(w) - 1
  }

  Z <- list()
  for(i in 1:n){
    Z[[i]] <- t(
      sapply(which(!is.na(x[,i])), function(z) {
        m <- rep(0, Tn)
        m[z] <- 1
        m
      })
    )
  }

  tau <- seq(0, 1, by = 1/(Tn-1))
  tau <- as.matrix(tau)
  d <- ncol(tau)

  inits = fdlm_init(t(x), tau, K)
  Beta = inits$Beta
  Psi = inits$Psi
  splineInfo = inits$splineInfo
  B <- splineInfo$Bmat
  Omega <- splineInfo$Omega

  bzz <- lapply(Z, function(q) crossprod(B, crossprod(q)))
  bzzb <- lapply(Z, function(q) crossprod(q%*%B))

  x2 <- x
  x2[is.na(x)] <- 0
  xl <- lapply(seq_len(ncol(x)), function(i) x2[,i])
  bzzx <- Map(function(w,v) w%*%v, bzz, xl)
  xlobs <-  lapply(seq_len(ncol(x)), function(i) x[!is.na(x[,i]),i])

  p = ncol(B)

  K = ncol(Beta)

  G = psBasis(1:Tn, K = Kphi, diff.ord=2)$X
  D = psBasis(1:Tn, K = Kphi, diff.ord=2)$P
  P = crossprod(D)

  J <- matrix(NA, nrow = Kphi, ncol = p)
  for(j in 1:Kphi){
    for(k in 1:p){
      J[j,k] <- fdapace::trapzRcpp(Tau, G[, j] * B[, k])
      # J[j,k] <- sum(G[, j] * B[, k])
    }

  }


  Fmat = B %*% Psi
  sig_theta <- 20
  sig_mu <- 20

  lambda = apply(Psi, 2, function(v) (ncol(B) -
                                        (d + 1))/crossprod(v, Omega) %*% v)


  H <- list()
  l_mu0 <- matrix(NA, nrow = K, ncol = n)

  # BtCon = NULL
  # X = NULL
  # X = cbind(rep(1, T), X)
  alpha <- solve(crossprod(w))%*%t(w)%*%(y)
  BetaPsiJ <- J%*%(tcrossprod(Psi, Beta))
  xi <- solve(tcrossprod(BetaPsiJ) + P)%*%(BetaPsiJ)%*%y
  lambda_xi <- 1

  btil <- colSums(c(xi)*(J))
  ywa <- y - w%*%alpha


  mu <- rowMeans(solve(crossprod(Fmat))%*%crossprod(Fmat,t(inits$Y0)))

  # xinit <- t(fpca.face(t(x))$Yhat)
  # Beta <- t(solve(crossprod(Fmat))%*%crossprod(Fmat,xinit))
  # mu <- colMeans(Beta)


  a1_mu = 2
  a2_mu = 3
  delta_mu = sampleMGP(matrix(mu, ncol = K), rep(1, K),
                       a1 = a1_mu, a2 = a2_mu)
  sig_mu = 1/sqrt(cumprod(delta_mu))


  theta = t(Beta - c(mu))
  zeta_theta = t(1/theta^2)
  nu = 3
  a1_eta = 2
  a2_eta = 3
  delta_theta = rep(1, K)
  sig_delta = 1/sqrt(cumprod(delta_theta))
  sig_theta = rep(sig_delta, each = n)/sqrt(zeta_theta)

  sig_ey <- .1
  sig_ex <- .1
  sig_alpha <- .01

  a_alpha <- 1
  b_alpha <- 1
  a_lambda_xi <- 1
  b_lambda_xi <- 1

  nas <- is.na(x)
  nobs <- sum(!(is.na(x)))
  pm <- mean(apply(x, 2, function(zz){
    with(rle(is.na(zz)), max(lengths[values],0))
  }))/Tn

  S <- S
  S_burn <- S_burn
  Sk <- S - S_burn

  theta_post <- list()
  mu_post <- list()
  Y_post <- list()
  Fmat_post <- list()
  Psi_post <- matrix(NA, Sk, K)
  alpha_post <- matrix(NA, Sk, L + 1)
  xi_post <- matrix(NA, Sk, Kphi)
  beta_post <- list()
  lambda_xi_post <- rep(NA, Sk)
  BetaPsiJ_post <- list()
  sig_ey_post <- rep(NA, Sk)

  yt <- rep(NA, S)

  progress <- floor(seq(1, S, length.out= 11))

  for(s in 1:S){
    # print(sum((ywa - t(BetaPsiJ)%*%(xi))^2))
    lambda <- sample_lambda(lambda, Psi, Omega = Omega,
                            d = d, uniformPrior = TRUE, orderLambdas = TRUE)

    sflc <- sampleFLCs_sofr(bzzx, ywa, btil, Beta, Psi, bzzb, Omega, lambda, sig_ex, sig_ey, mu, theta)

    Psi <- sflc[[1]]
    mu <- sflc[[4]]
    theta <- sflc[[5]]

    Fmat = B %*% Psi

    Fmatl <- lapply(1:n, function(zz) Fmat[!nas[,zz], ] )

    btil <- colSums(c(xi)*(J))

    btilpsi <- crossprod(Psi, btil)

    ywa <- y - w%*%alpha

    fma <- sapply(btilpsi^2, function(zz) rep(zz, n))
    if(pm <= .6){
      xft <- list()
      xfm <- list()
      for(i in 1:n){
        xft[[i]] <- xlobs[[i]] -  Fmatl[[i]]*rep(theta[,i], each = Tni[i])
      }

      for(jj in 1:K){

        fmkka <- colSums((theta[-jj,] + mu[-jj])*c(btilpsi[-jj]))

        lmuy <- (1/sig_ey) * btilpsi[jj]*(ywa - theta[jj,]*btilpsi[jj] - fmkka)

        Q_mu <- sum(fma[,jj])*(1/sig_ey)

        l_mu <- rep(NA, n)
        fmb <- rep(NA, n)
        for(i in 1:n){

          fmb[i] <- sum(Fmatl[[i]][,jj]^2)

          fmkkb <- Fmatl[[i]][,-jj]%*%(theta[-jj,i] + mu[-jj])

          l_mu[i] <- lmuy[i] + (1/sig_ex)*crossprod(Fmatl[[i]][,jj],
                                                    xft[[i]][,jj] - fmkkb)

        }
        l_mu <- sum(l_mu)

        Q_mu <- 1/(Q_mu + sum(fmb)*(1/sig_ex) + (1/sig_mu[jj]^2))

        mu[jj] <- rnorm(1, Q_mu*l_mu, sqrt(Q_mu))
      }


      for(i in 1:n){
        xfm[[i]] <- xlobs[[i]] -  Fmatl[[i]]*rep(mu, each = Tni[i])
      }

      for(jj in 1:K){

        fmkka <- colSums((theta[-jj,] + mu[-jj])*c(btilpsi[-jj]))

        lthetay <- (1/sig_ey) * btilpsi[jj]*(ywa - mu[jj]*btilpsi[jj] - fmkka)

        l_theta <- rep(NA, n)

        fmb <- rep(NA, n)
        for(i in 1:n){

          fmb[i] <- sum(Fmatl[[i]][,jj]^2)

          fmkkb <- Fmatl[[i]][,-jj]%*%(theta[-jj,i] + mu[-jj])

          l_theta[i] <- lthetay[i] + (1/sig_ex)*crossprod(Fmatl[[i]][,jj],
                                                          xfm[[i]][,jj] - fmkkb)

        }

        Q_theta <- 1/(fma[,jj]*(1/sig_ey) + fmb*(1/sig_ex) + (1/sig_theta[,jj]^2))

        theta[jj,] <- rnorm(n, Q_theta*l_theta, sqrt(Q_theta))


      }
    }else{
      zf <- lapply(Z, function(q) q%*%Fmat)
      fzzf <- lapply(Z, function(q) crossprod(q%*%Fmat))

      Q <- list()
      l <- list()
      for(i in 1:n){
        Q[[i]] <- solve(diag(1/sig_theta[i,]^2, K) + (1/sig_ex)*fzzf[[i]])

        H[[i]] <- (1/sig_ex)*fzzf[[i]]%*%(diag(K) - (1/sig_ex)*Q[[i]]%*%fzzf[[i]])

        l_mu0[,i] <- (t(zf[[i]]) - (1/sig_ex)*fzzf[[i]]%*%Q[[i]]%*%t(Z[[i]]%*%Fmat))%*%(x[,i])[!nas[,i]]

      }

      Q_mu <- solve(diag(1/sig_mu^2, K) + Reduce("+", H))
      l_mu <- (1/sig_ex)*rowSums(l_mu0)

      mu <- t(Rfast::rmvnorm(1, Q_mu%*%l_mu, Q_mu))
      muf <- Fmat%*%mu
      mu <- c(mu)
      for(i in 1:n){
        l[[i]] <- (1/sig_ex)*crossprod(Fmat[!nas[,i],], (x[,i] - muf)[!nas[,i]])
      }

      theta <- sapply(1:n, function(z)  {
        Rfast::rmvnorm(1, (Q[[z]])%*%(l[[z]]), Q[[z]])
      }
      )

      thetaf <- Fmat%*%theta
    }
    Beta <- t(theta + mu)

    BetaPsiJ <- J%*%(tcrossprod(Psi, Beta))

    Q_xi <- (1/sig_ey)*tcrossprod(BetaPsiJ) + lambda_xi*P
    l_xi <- (1/sig_ey)*crossprod(t(BetaPsiJ), ywa)

    ch_Q <- chol(matrix(Q_xi, nrow=Kphi, ncol=Kphi))
    xi <- backsolve(ch_Q,
                    forwardsolve(t(ch_Q), l_xi) +
                      rnorm((Kphi)))

    Q_alpha <- (1/sig_ey)*crossprod(w) + (1/sig_alpha)*diag(L+1)
    l_alpha <- (1/sig_ey)*crossprod(w, y - crossprod(BetaPsiJ, xi))

    ch_Q <- chol(matrix(Q_alpha, nrow=L+1, ncol=L+1))
    alpha <- backsolve(ch_Q,
                       forwardsolve(t(ch_Q), l_alpha) +
                         rnorm((L+1)))

    # if(s < 500){
    # alpha <- alpha_true
    # }
    lambda_xi <- rgamma(1, a_lambda_xi + Kphi/2, b_lambda_xi + sum((xi^2)/2))

    # crossprod(Fmat[!nas[,i],], x[!nas[,i],i] - outer(Fmat[!nas[,i],1], theta[1,]))
    #
    # aa <- lapply(1:K, function(z) outer(Fmat[,z],theta[z,]) + Fmat[,-z]%*%(theta[-z,] + mu[-z]))


    sig_alpha <- 1/rgamma(1, a_alpha + 1/2, b_alpha + ((alpha^2)/2)  )

    sig_ey <- 1/rgamma(1, .1 +  n/2, .1 +  .5*sum((ywa - crossprod(BetaPsiJ, xi))^2))

    sig_ex <- 1/rgamma(1, .1 + nobs/2, rate = .1 + sum((x - tcrossprod(Fmat, Beta) )^2, na.rm = TRUE)/2)
    #
    # sig_ex <- 1/rtrunc(1, "gamma", a=0.001, b=1000, shape = .1 + nobs/2,
    #                    rate = .1 + sum((x - tcrossprod(Fmat, Beta) )^2, na.rm = TRUE)/2)

    delta_mu = sampleMGP(theta.jh = matrix(mu, ncol = K),
                         delta.h = delta_mu, a1 = a1_mu, a2 = a2_mu)
    sig_mu = 1/sqrt(cumprod(delta_mu))

    a1_mu = uni.slice(a1_mu, g = function(a) {
      dgamma(delta_mu[1], shape = a, rate = 1, log = TRUE) +
        dgamma(a, shape = 2, rate = 1, log = TRUE)
    }, lower = 0, upper = Inf)
    a2_mu = uni.slice(a2_mu, g = function(a) {
      sum(dgamma(delta_mu[-1], shape = a, rate = 1,
                 log = TRUE)) + dgamma(a, shape = 2, rate = 1,
                                       log = TRUE)
    }, lower = 0, upper = Inf)

    delta_theta = sampleMGP(theta.jh = matrix(t(theta) * sqrt(zeta_theta),
                                              ncol = K), delta.h = delta_theta, a1 = a1_eta, a2 = a2_eta)
    sig_delta = 1/sqrt(cumprod(delta_theta))

    a1_eta = uni.slice(a1_eta, g = function(a) {
      dgamma(delta_theta[1], shape = a, rate = 1,
             log = TRUE) + dgamma(a, shape = 2, rate = 1,
                                  log = TRUE)
    }, lower = 0, upper = Inf)
    a2_eta = uni.slice(a2_eta, g = function(a) {
      sum(dgamma(delta_theta[-1], shape = a, rate = 1,
                 log = TRUE)) + dgamma(a, shape = 2, rate = 1,
                                       log = TRUE)
    }, lower = 0, upper = Inf)

    zeta_theta = matrix(rgamma(n = n * K, shape = nu/2 +
                                 1/2, rate = nu/2 + (t(theta)/rep(sig_delta, each = n))^2/2), nrow = n)
    sig_theta = rep(sig_delta, each = n)/sqrt(zeta_theta)

    nu = uni.slice(nu, g = function(nu) {
      sum(dgamma(zeta_theta, shape = nu/2, rate = nu/2,
                 log = TRUE)) + dunif(nu, min = 2, max = 128,
                                      log = TRUE)
    }, lower = 2, upper = 128)
    # matplot(Fmat%*%t(Beta), type = "l")

    if(s > S_burn){
      sk <- s - S_burn

      alpha_post[sk,] <- alpha
      xi_post[sk,] <- xi
      mu_post[[sk]] <- mu
      theta_post[[sk]] <- theta
      Fmat_post[[sk]] <- Fmat
      lambda_xi_post[sk] <- lambda_xi
      beta_post[[sk]] <- theta + c(mu)
      BetaPsiJ_post[[sk]] <- BetaPsiJ
      sig_ey_post[sk] <- sig_ey

    }

    if(s %in% progress & showprogress == TRUE){
      cat(paste0("MCMC draws: [", s, "/", S, "]\n"))
    }

  }

  betaf_post <- array( NA , c(Tn,n,Sk) )
  muf_post <- matrix(NA, Tn, Sk)
  phi_post <- matrix(NA, Tn, Sk)

  for(i in (S_burn+1):S){
    betaf_post[,,i-S_burn] <- (Fmat_post[[i-S_burn]])%*%(beta_post[[i-S_burn]])

    muf_post[,i-S_burn] <- (Fmat_post[[i-S_burn]])%*%(mu_post[[i-S_burn]])

    phi_post[,i-S_burn] <- G%*%(xi_post[i-S_burn,])
  }

  output <- list(betaf_post = betaf_post,
                 muf_post = muf_post,
                 Fmat_post = Fmat_post,
                 phi_post = phi_post,
                 alpha_post = alpha_post,
                 BetaPsiJ_post = BetaPsiJ_post,
                 sig_ey_post = sig_ey_post
  )
  return(output)

}
