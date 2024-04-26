#----------------------------------------------------------------------------
# Function factor model for sparsely observed curves
#----------------------------------------------------------------------------

#' Sparse functional factor model
#'
#' Performs MCMC estimation for Bayesian functional factor model.
#'
#' @param Ym matrix of functional responses. Each column corresponds to one function and each row is a measurement at a single time point on a common grid. May contain NAs.
#' @param K number of basis functions.
#' @param S number of total MCMC iterations
#' @param S_burn number of initial MCMC iterations to discard for warm-up
#' @param sparse set to TRUE if functions considered sparsely observed.
#'
#' @return A list with the following elements:
#' * X: the inputted design matrix
#' @export
#'
#' @importFrom stats poly rgamma rnorm dgamma dunif fitted median quantile rexp runif sd splinefun
sffm <- function(Ym, K = 10, S = 2000, S_burn = S/2, sparse = TRUE){

n <- ncol(Ym)
Tn <- nrow(Ym)
K <- K

Z <- list()
for(i in 1:n){
  Z[[i]] <- t(
    sapply(which(!is.na(Ym[,i])), function(z) {
      m <- rep(0, Tn)
      m[z] <- 1
      m
    })
  )
}

tau <- seq(0, 1, by = 1/(Tn-1))
tau <- as.matrix(tau)
d <- ncol(tau)

inits = fdlm_init(t(Ym), tau, K)
Beta = inits$Beta
Psi = inits$Psi
splineInfo = inits$splineInfo
B <- splineInfo$Bmat

Lm <- ncol(B)

K = ncol(Beta)

Y = inits$Y0

# Initialize MCMC
Fmat = B %*% Psi
BetaPsit = tcrossprod(Beta, Psi)
Btheta = tcrossprod(BetaPsit, B)
sig_nu = 0
sig_e = sd(Y - Btheta, na.rm = TRUE)
sig_theta <- 20
sig_mu <- 20

tau_f_k = apply(Psi, 2, function(x) (ncol(B) -
                                       (d + 1))/crossprod(x, splineInfo$Omega) %*% x)

H <- list()
l_mu0 <- matrix(NA, nrow = K, ncol = n)


mu <- rowMeans(solve(crossprod(Fmat))%*%crossprod(Fmat,t(Y)))

a1_mu = 2
a2_mu = 3
delta_mu = sampleMGP(matrix(mu, ncol = K), rep(1, K),
                     a1 = a1_mu, a2 = a2_mu)
sig_mu = 1/sqrt(cumprod(delta_mu))


theta = t(Beta - c(mu))
xi_theta = t(1/theta^2)
nu = 3
a1_eta = 2
a2_eta = 3
delta_theta = rep(1, K)
sig_delta = 1/sqrt(cumprod(delta_theta))
sig_theta = rep(sig_delta, each = n)/sqrt(xi_theta)

nas <- is.na(Ym)
nobs <- sum(!(is.na(Ym)))

bzzb <- lapply(Z, function(q) crossprod(q%*%B))

bzzy <- lapply(1:n, function(i) crossprod(B[!nas[,i],], Ym[!nas[,i],i]))

S <- S
S_burn <- S_burn
Sk <- S - S_burn

beta_post <- list()
mu_post <- list()
Y_post <- list()
Fmat_post <- list()
Psi_post <- matrix(NA, Sk, Lm)

yt <- rep(NA, S)

progress <- floor(seq(1, S, length.out= 11))

s1t1 <- Sys.time()
for(s in 1:S){

  # Sample loading curve parameters
  tau_f_k <- sample_lambda(tau_f_k, Psi, Omega = splineInfo$Omega,
                           d = d, uniformPrior = FALSE, orderLambdas = FALSE)

  sflc <- sampleFLCs(bzzy, Beta = Beta, Psi, bzzb, splineInfo$Omega, tau_f_k, sigma_e = sig_e)

  Psi <- sflc[[1]]

  mu <- mu*sflc[[2]]
  theta <- theta*sflc[[2]]

  Fmat = B %*% Psi

  fzzf <- lapply(1:n, function(i) crossprod(Fmat[!nas[,i],]))


  # Sample factor model coefficients
  Q <- list()
  l <- list()
  for(i in 1:n){
    Q[[i]] <- solve(diag(1/sig_theta[i,]^2, K) + (1/sig_e)*fzzf[[i]])

    H[[i]] <- (1/sig_e)*fzzf[[i]]%*%(diag(K) - (1/sig_e)*Q[[i]]%*%fzzf[[i]])

    l_mu0[,i] <- (t(Fmat[!nas[,i],]) - (1/sig_e)*fzzf[[i]]%*%tcrossprod(Q[[i]], Fmat[!nas[,i],]))%*%(Ym[,i][!nas[,i]])

  }

  Q_mu <- solve(diag(1/sig_mu^2, K) + Reduce("+", H))
  l_mu <- (1/sig_e)*rowSums(l_mu0)

  mu <- t(Rfast::rmvnorm(1, Q_mu%*%l_mu, Q_mu))
  muf <- Fmat%*%mu

  for(i in 1:n){
    l[[i]] <- (1/sig_e)*crossprod(Fmat[!nas[,i],], (Ym[,i] - muf)[!nas[,i]])
  }

  theta <- sapply(1:n, function(zz)  {
    Rfast::rmvnorm(1, (Q[[zz]])%*%(l[[zz]]), Q[[zz]])
  }
  )


  thetaf <- Fmat%*%theta

  Beta <- t(theta + c(mu))

  Y <- matrix(rnorm(length(Ym), (thetaf + c(muf)), sqrt(sig_e)), nrow = Tn, ncol = n)


  # Sample variance parameters and multiplicative gamma process
  sig_e <- 1/rgamma(1, nobs/2, rate = sum((Ym - thetaf - c(muf) )^2, na.rm = TRUE)/2)

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

  delta_theta = sampleMGP(theta.jh = matrix(t(theta) * sqrt(xi_theta),
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

  xi_theta = matrix(rgamma(n = n * K, shape = nu/2 +
                             1/2, rate = nu/2 + (t(theta)/rep(sig_delta, each = n))^2/2), nrow = n)
  sig_theta = rep(sig_delta, each = n)/sqrt(xi_theta)

  nu = uni.slice(nu, g = function(nu) {
    sum(dgamma(xi_theta, shape = nu/2, rate = nu/2,
               log = TRUE)) + dunif(nu, min = 2, max = 128,
                                    log = TRUE)
  }, lower = 2, upper = 128)

  if(s > S_burn){
    sk <- s - S_burn

    beta_post[[sk]] <- theta + c(mu)

    mu_post[[sk]] <- mu

    Y_post[[sk]] <- Y

    Fmat_post[[sk]] <- Fmat

  }

  if(s %in% progress){
    cat(paste0("MCMC draws: [", s, "/", S, "]\n"))
  }

}
s1t2 <- difftime(Sys.time(), s1t1, units = "secs")

betaf_post <- array( NA , c(Tn,n,Sk) )
for(i in (S_burn+1):S){
  betaf_post[,,i-S_burn] <- (Fmat_post[[i-S_burn]])%*%(beta_post[[i-S_burn]])
}


muf_post <- matrix(NA, Tn, Sk)

for(i in (S_burn+1):S){
  muf_post[,i-S_burn] <- (Fmat_post[[i-S_burn]])%*%(mu_post[[i-S_burn]])
}

output <- list(betaf_post = betaf_post,
               muf_post = muf_post,
               # Y_post = Y_post,
               # Fmat_post = Fmat_post
               )
class(output) <- c('list', 'sffm')
return(output)

}


