################################################################################
##  Heckman selection model with nonparametric Gaussian scale-mixture errors
##  (HeckmanEM  x  Seo, Noh, Lee & Yoon (2017, JKSS) 결합)
##
##  Model
##    Y*_i = x_i' beta  + e1_i        (outcome, observed only if cc_i = 1)
##    S*_i = w_i' gamma + e2_i        (selection, cc_i = 1{S*_i > 0})
##    (e1_i, e2_i) | tau_i  ~  N2(0, tau_i * Sigma),   tau_i ~ G  (G unspecified)
##    Sigma = [ sigma^2   rho*sigma ]
##            [ rho*sigma     1     ]
##
##    Normal : G = delta_1
##    t_nu   : G = law of 1/U, U ~ Gamma(nu/2, nu/2)       (HeckmanEM "T")
##    CN     : G = (1-nu1) delta_1 + nu1 delta_{1/nu2}       (HeckmanEM "CN")
##    여기서는 G를 NPMLE로 추정 (유한 개의 support point를 갖는 이산분포)
##
##  Identification : sigma 와 gamma 의 scale 이 G 의 scale 과 교락되므로
##                   E_G(tau) = 1 로 정규화 (=> Var(e) = Sigma, Var(e2) = 1)
##
##  Algorithm (Seo et al. 2017, Sec 2.2 와 동일한 교대 최적화)
##    (1) theta 고정 -> G 의 NPMLE : CNM (Wang, 2007) = gradient 로 support 추가 +
##                                   NNLS 로 weight 갱신 + line search
##    (2) G 고정     -> theta      : HeckmanEM 의 ECM 을 K-성분 scale mixture 로 확장
##                                   (절단정규 적률이 닫힌형이라 MomTrunc/mvtnorm 불필요)
##  두 단계 모두 관측 로그우도를 감소시키지 않으므로 전체가 단조 증가.
##
##  Dependencies: base R (stats) only
################################################################################

## ---------------------------------------------------------------- utilities
.lse_rows <- function(M) {                     # row-wise log-sum-exp
  m <- apply(M, 1, max)
  m + log(rowSums(exp(M - m)))
}

.mills <- function(a) {                        # phi(a) / Phi(a), stable
  exp(dnorm(a, log = TRUE) - pnorm(a, log.p = TRUE))
}

## Lawson-Hanson non-negative least squares (base R)
.nnls <- function(A, b, tol = 1e-10, maxit = 3 * ncol(A) + 50) {
  p <- ncol(A); x <- numeric(p); P <- rep(FALSE, p)
  w <- drop(crossprod(A, b - A %*% x))
  it <- 0
  ls_P <- function(P) {
    z <- numeric(p)
    co <- qr.coef(qr(A[, P, drop = FALSE]), b)
    co[is.na(co)] <- 0
    z[P] <- co
    z
  }
  while (any(!P & w > tol) && it < maxit) {
    it <- it + 1
    cand <- which(!P)
    P[cand[which.max(w[cand])]] <- TRUE
    z <- ls_P(P)
    inner <- 0
    while (any(z[P] <= 0) && inner < maxit) {
      inner <- inner + 1
      idx <- P & z <= 0
      den <- x[idx] - z[idx]
      ratio <- ifelse(den > 0, x[idx] / den, 0)
      alpha <- min(ratio)
      x <- x + alpha * (z - x)
      P <- P & x > tol
      x[!P] <- 0
      z <- ls_P(P)
    }
    x <- z
    w <- drop(crossprod(A, b - A %*% x))
  }
  pmax(x, 0)
}

## ------------------------------------------------- component log-densities
## n x K matrix: log f_i(tau_k ; theta)
.logf_mat <- function(tau, th, y, x, w, cc) {
  sigma <- th$sigma; rho <- th$rho
  mu1 <- drop(x %*% th$beta); mu2 <- drop(w %*% th$gamma)
  n <- length(cc); K <- length(tau)
  L <- matrix(0, n, K)
  o <- cc == 1; c0 <- !o
  e <- (y - mu1)[o]
  for (k in seq_len(K)) {
    st <- sqrt(tau[k])
    if (any(o))
      L[o, k] <- dnorm(e, 0, sigma * st, log = TRUE) +
        pnorm((mu2[o] + rho * e / sigma) / (st * sqrt(1 - rho^2)), log.p = TRUE)
    if (any(c0))
      L[c0, k] <- pnorm(-mu2[c0] / st, log.p = TRUE)
  }
  L
}

loglik.npsm <- function(th, G, y, x, w, cc) {
  L <- .logf_mat(G$tau, th, y, x, w, cc)
  sum(.lse_rows(sweep(L, 2, log(G$prob), "+")))
}

## directional derivative D(t) = sum_i f_i(t)/f_i(G) - n
.grad_G <- function(t, lfG, th, y, x, w, cc) {
  Lt <- .logf_mat(t, th, y, x, w, cc)
  colSums(exp(Lt - lfG)) - length(lfG)
}

## ------------------------------------------------------- G-step : CNM
cnm_G <- function(th, G, y, x, w, cc, tau.min = 0.01, maxit = 100,
                  tol = 1e-6, ngrid = 100, verbose = FALSE) {
  n <- length(cc)
  o <- cc == 1
  r2 <- c(((y[o] - drop(x[o, , drop = FALSE] %*% th$beta)) / th$sigma)^2,
          drop(w[!o, , drop = FALSE] %*% th$gamma)^2)
  tau.max <- max(50, 3 * max(r2))
  grid <- exp(seq(log(tau.min), log(tau.max), length.out = ngrid))

  ll <- loglik.npsm(th, G, y, x, w, cc)
  for (it in seq_len(maxit)) {
    L   <- .logf_mat(G$tau, th, y, x, w, cc)
    lfG <- .lse_rows(sweep(L, 2, log(G$prob), "+"))

    ## (a) new support points: local maxima of D(t) with D > 0
    dg <- .grad_G(grid, lfG, th, y, x, w, cc)
    K  <- length(dg)
    loc <- which(dg > tol &
                 c(TRUE, dg[-1] >= dg[-K]) & c(dg[-K] >= dg[-1], TRUE))
    newt <- numeric(0); newd <- numeric(0)
    for (j in loc) {
      lo <- log(grid[max(1, j - 1)]); hi <- log(grid[min(K, j + 1)])
      if (hi > lo) {
        op <- optimize(function(s) .grad_G(exp(s), lfG, th, y, x, w, cc),
                       c(lo, hi), maximum = TRUE)
        newt <- c(newt, exp(op$maximum)); newd <- c(newd, op$objective)
      } else { newt <- c(newt, grid[j]); newd <- c(newd, dg[j]) }
    }
    maxD <- max(c(newd, .grad_G(G$tau, lfG, th, y, x, w, cc)))
    if (maxD < tol) break

    tau  <- c(G$tau, newt)
    prob <- c(G$prob, rep(0, length(newt)))

    ## (b) weights: min ||S a - 2||, a >= 0, sum(a)=1  (Wang 2007), then line search
    Lall <- cbind(L, if (length(newt)) .logf_mat(newt, th, y, x, w, cc))
    S    <- exp(Lall - lfG)
    big  <- 1e3 * sqrt(n)
    a    <- .nnls(rbind(S, rep(big, ncol(S))), c(rep(2, n), big))
    if (sum(a) <= 0) break
    a    <- a / sum(a)

    alpha <- 1; accepted <- FALSE
    for (h in 1:40) {
      pn <- prob + alpha * (a - prob)
      keep <- pn > 1e-12
      Gn <- list(tau = tau[keep], prob = pn[keep] / sum(pn[keep]))
      lln <- loglik.npsm(th, Gn, y, x, w, cc)
      if (lln >= ll - 1e-12) { accepted <- TRUE; break }
      alpha <- alpha / 2
    }
    if (!accepted) break
    ## merge numerically identical support points
    ord <- order(Gn$tau); Gn$tau <- Gn$tau[ord]; Gn$prob <- Gn$prob[ord]
    dup <- c(FALSE, diff(log(Gn$tau)) < 1e-8)
    if (any(dup)) {
      grp <- cumsum(!dup)
      Gn <- list(tau = tapply(Gn$tau, grp, `[`, 1), prob = tapply(Gn$prob, grp, sum))
      Gn <- lapply(Gn, as.numeric)
    }
    dll <- lln - ll
    G <- Gn; ll <- lln
    if (verbose) cat(sprintf("  CNM %3d  K=%d  maxD=%.3e  ll=%.6f\n", it, length(G$tau), maxD, ll))
    if (dll < 1e-10 * abs(ll)) break
  }
  list(G = G, ll = ll)
}

## ----------------------------------------------------- theta-step : ECM
## HeckmanEM 의 ECM 을 G = sum_k pi_k delta_{tau_k} 로 일반화
ecm_theta <- function(th, G, y, x, w, cc, rho.max = 0.999) {
  n <- length(cc); p <- ncol(x); q <- ncol(w)
  beta <- th$beta; gama <- th$gamma; sigma <- th$sigma; rho <- th$rho
  rhoa <- rho * sigma
  Sigma <- matrix(c(sigma^2, rhoa, rhoa, 1), 2, 2)
  mu1 <- drop(x %*% beta); mu2 <- drop(w %*% gama)

  L   <- sweep(.logf_mat(G$tau, th, y, x, w, cc), 2, log(G$prob), "+")
  Pik <- exp(L - .lse_rows(L))                     # posterior of component k
  K <- length(G$tau)

  u0 <- numeric(n); u1 <- matrix(0, n, 2); u2 <- array(0, c(n, 2, 2))
  o <- cc == 1; c0 <- !o
  yo <- y; yo[c0] <- 0
  for (k in seq_len(K)) {
    tk <- G$tau[k]; kap <- Pik[, k] / tk           # E[1/tau * 1{k}]
    ES <- ES2 <- EY <- EY2 <- EYS <- numeric(n)
    if (any(o)) {                                  # S* | y  ~ N(m, v), S* > 0
      m <- mu2[o] + rho * (y[o] - mu1[o]) / sigma
      s <- sqrt(tk * (1 - rho^2)); a <- m / s; lam <- .mills(a)
      ES[o]  <- m + s * lam
      ES2[o] <- s^2 * (1 - lam * (lam + a)) + ES[o]^2
      EY[o]  <- y[o]; EY2[o] <- y[o]^2; EYS[o] <- y[o] * ES[o]
    }
    if (any(c0)) {                                 # (Y*,S*) ~ N2(mu, tk Sigma), S* < 0
      s <- sqrt(tk); b <- -mu2[c0] / s; lam <- .mills(b)
      ES[c0]  <- mu2[c0] - s * lam
      VS      <- tk * (1 - b * lam - lam^2)
      ES2[c0] <- VS + ES[c0]^2
      D  <- ES[c0] - mu2[c0]; D2 <- ES2[c0] - 2 * mu2[c0] * ES[c0] + mu2[c0]^2
      EY[c0]  <- mu1[c0] + rhoa * D
      EYS[c0] <- mu1[c0] * ES[c0] + rhoa * (ES2[c0] - mu2[c0] * ES[c0])
      EY2[c0] <- tk * sigma^2 * (1 - rho^2) + mu1[c0]^2 + 2 * mu1[c0] * rhoa * D + rhoa^2 * D2
    }
    u0 <- u0 + kap
    u1[, 1] <- u1[, 1] + kap * EY;  u1[, 2] <- u1[, 2] + kap * ES
    u2[, 1, 1] <- u2[, 1, 1] + kap * EY2; u2[, 2, 2] <- u2[, 2, 2] + kap * ES2
    u2[, 1, 2] <- u2[, 1, 2] + kap * EYS
  }
  u2[, 2, 1] <- u2[, 1, 2]

  ## CM-step 1 : (beta, gamma) | Sigma   (GLS with weights u0)
  SigI <- solve(Sigma)
  A11 <- crossprod(x, u0 * x) * SigI[1, 1]; A22 <- crossprod(w, u0 * w) * SigI[2, 2]
  A12 <- crossprod(x, u0 * w) * SigI[1, 2]
  Amat <- rbind(cbind(A11, A12), cbind(t(A12), A22))
  bvec <- c(crossprod(x, SigI[1, 1] * u1[, 1] + SigI[1, 2] * u1[, 2]),
            crossprod(w, SigI[2, 1] * u1[, 1] + SigI[2, 2] * u1[, 2]))
  bC <- solve((Amat + t(Amat)) / 2, bvec)
  beta <- bC[1:p]; gama <- bC[(p + 1):(p + q)]

  ## CM-step 2 : Sigma | (beta, gamma), with Sigma22 = 1
  a1 <- drop(x %*% beta); a2 <- drop(w %*% gama)
  O11 <- u2[, 1, 1] - 2 * a1 * u1[, 1] + u0 * a1^2
  O22 <- u2[, 2, 2] - 2 * a2 * u1[, 2] + u0 * a2^2
  O12 <- u2[, 1, 2] - a1 * u1[, 2] - a2 * u1[, 1] + u0 * a1 * a2
  rhoa <- sum(O12) / sum(O22)
  psi  <- mean(O11 - 2 * rhoa * O12 + rhoa^2 * O22)
  sigma <- sqrt(psi + rhoa^2)
  rho <- max(min(rhoa / sigma, rho.max), -rho.max)

  list(beta = beta, gamma = gama, sigma = sigma, rho = rho)
}

## ---------------------------------------------------- initial values
heckit2step <- function(y, x, w, cc) {
  pr <- suppressWarnings(glm.fit(w, cc, family = binomial(link = "probit")))
  gama <- pr$coefficients
  o <- cc == 1
  lp <- drop(w %*% gama)[o]; imr <- .mills(lp)
  fit <- lm.fit(cbind(x[o, , drop = FALSE], imr), y[o])
  p <- ncol(x); beta <- fit$coefficients[1:p]; bl <- fit$coefficients[p + 1]
  delta <- imr * (imr + lp)
  sigma <- sqrt(mean(fit$residuals^2) + bl^2 * mean(delta))
  rho <- max(min(bl / sigma, 0.9), -0.9)
  list(beta = unname(beta), gamma = unname(gama), sigma = unname(sigma), rho = unname(rho))
}

## ========================================================= main function
##  y  : outcome (NA 또는 아무 값이어도 cc==0 이면 무시)
##  x  : outcome design matrix (intercept 포함),  w : selection design matrix
##  cc : 1 = observed, 0 = censored/not selected
##  fixG = TRUE 이면 G = delta_1 로 고정 → HeckmanEM(family="Normal") 과 같은 모형
HeckmanNPSM <- function(y, x, w, cc, tau.min = 0.01, error = 1e-8,
                        iter.max = 500, fixG = FALSE, se = TRUE,
                        init = NULL, verbose = FALSE) {
  x <- as.matrix(x); w <- as.matrix(w); cc <- as.numeric(cc)
  y <- as.numeric(y); y[cc == 0] <- 0
  n <- length(cc)

  th <- if (is.null(init)) heckit2step(y, x, w, cc) else init
  G  <- list(tau = 1, prob = 1)
  ll <- loglik.npsm(th, G, y, x, w, cc)
  lls <- ll

  for (it in seq_len(iter.max)) {
    ## (1) G-step
    if (!fixG) {
      G <- cnm_G(th, G, y, x, w, cc, tau.min = tau.min)$G
      ## normalisation E_G(tau) = 1  (likelihood invariant)
      m <- sum(G$tau * G$prob)
      G$tau <- G$tau / m
      th$gamma <- th$gamma / sqrt(m)
      th$sigma <- th$sigma * sqrt(m)
    }
    ## (2) theta-step
    th <- ecm_theta(th, G, y, x, w, cc)
    ll.new <- loglik.npsm(th, G, y, x, w, cc)
    lls <- c(lls, ll.new)
    if (verbose) cat(sprintf("iter %3d  logL = %.8f  K = %d\n", it, ll.new, length(G$tau)))
    if (ll.new < ll - 1e-6 * abs(ll)) warning("log-likelihood decreased at iteration ", it)
    conv <- abs(ll.new - ll) < error * abs(ll)
    ll <- ll.new
    if (conv) break
  }

  ## robustness weights  E[1/tau | data]  (paper's w_i)
  L   <- sweep(.logf_mat(G$tau, th, y, x, w, cc), 2, log(G$prob), "+")
  Pik <- exp(L - .lse_rows(L))
  wt  <- drop(Pik %*% (1 / G$tau))

  out <- list(beta = th$beta, gamma = th$gamma, sigma = th$sigma, rho = th$rho,
              G = data.frame(tau = G$tau, prob = G$prob),
              logL = ll, loglik.path = lls, iter = it, converged = it < iter.max,
              weights = wt, n = n, p = ncol(x), q = ncol(w),
              y = y, x = x, w = w, cc = cc, fixG = fixG)
  ## AIC/BIC: G 의 자유모수는 2K-2 (정규화 1개 제외) 로 계산 (참고용)
  k.par <- ncol(x) + ncol(w) + 2 + if (fixG) 0 else 2 * nrow(out$G) - 2
  out$npar <- k.par
  out$AIC <- -2 * ll + 2 * k.par
  out$BIC <- -2 * ll + log(n) * k.par
  if (se) out$se <- HeckmanNPSM.se(out)
  class(out) <- "HeckmanNPSM"
  out
}

## ----------------------------------------------------- standard errors
## G 를 추정값 Ĝ 에 고정한 관측정보행렬 (Seo et al. 의 Î_n 에 해당; G 추정의
## 불확실성은 반영하지 않으므로 보수적으로는 bootstrap 권장)
HeckmanNPSM.se <- function(fit) {
  p <- fit$p; q <- fit$q
  G <- list(tau = fit$G$tau, prob = fit$G$prob)
  unpack <- function(v) list(beta = v[1:p], gamma = v[(p + 1):(p + q)],
                             sigma = exp(v[p + q + 1]), rho = tanh(v[p + q + 2]))
  f <- function(v) loglik.npsm(unpack(v), G, fit$y, fit$x, fit$w, fit$cc)
  v0 <- c(fit$beta, fit$gamma, log(fit$sigma), atanh(fit$rho))
  H <- optimHess(v0, f)
  V <- tryCatch(solve(-H), error = function(e) MASS_ginv(-H))
  J <- diag(c(rep(1, p + q), fit$sigma, 1 - fit$rho^2))    # delta method
  Vt <- J %*% V %*% J
  s <- sqrt(pmax(diag(Vt), 0))
  names(s) <- c(paste0("beta", seq_len(p) - 1), paste0("gamma", seq_len(q) - 1), "sigma", "rho")
  s
}
MASS_ginv <- function(A) { s <- svd(A); d <- ifelse(s$d > 1e-10 * max(s$d), 1 / s$d, 0)
  s$v %*% (d * t(s$u)) }

print.HeckmanNPSM <- function(x, digits = 4, ...) {
  est <- c(x$beta, x$gamma, x$sigma, x$rho)
  nm  <- c(paste0("beta", seq_len(x$p) - 1), paste0("gamma", seq_len(x$q) - 1), "sigma", "rho")
  tab <- cbind(Estimate = est)
  if (!is.null(x$se)) tab <- cbind(tab, SE = x$se, z = est / x$se)
  rownames(tab) <- nm
  cat("Heckman selection model with", if (x$fixG) "normal errors (G = delta_1)\n"
      else "nonparametric Gaussian scale-mixture errors (NPMLE of G)\n")
  print(round(tab, digits))
  cat("\nMixing distribution G (normalised E[tau] = 1):\n")
  print(round(x$G, digits), row.names = FALSE)
  cat(sprintf("\nlogL = %.4f   AIC = %.2f   BIC = %.2f   iter = %d   converged = %s\n",
              x$logL, x$AIC, x$BIC, x$iter, x$converged))
  invisible(x)
}

## ------------------------------------------------ fitted error density
## outcome error e1 의 추정 밀도 (논문 Fig.3 과 같은 용도)
dens.e1 <- function(fit, e) {
  sapply(e, function(z) sum(fit$G$prob * dnorm(z, 0, fit$sigma * sqrt(fit$G$tau))))
}

## ------------------------------------------------ data generator
##  rtau : function(n) -> tau  (예: t_nu : function(n) 1/rgamma(n, nu/2, nu/2))
rHeckmanSMN <- function(n, beta, gamma, sigma, rho, x, w, rtau = function(n) rep(1, n)) {
  tau <- rtau(n)
  z1 <- rnorm(n); z2 <- rnorm(n)
  e2 <- sqrt(tau) * z2
  e1 <- sqrt(tau) * sigma * (rho * z2 + sqrt(1 - rho^2) * z1)
  ys <- drop(x %*% beta) + e1
  cc <- as.numeric(drop(w %*% gamma) + e2 > 0)
  y <- ifelse(cc == 1, ys, NA)
  list(y = y, cc = cc, tau = tau)
}
