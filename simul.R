source("C:\\Users\\heeju\\Downloads\\HeckmanNPSM.R")
set.seed(2026)
n <- 500; R <- 25
beta <- c(1, 1, -1); gamma <- c(0.5, 1, -0.5, 1); rho <- 0.6
res <- NULL
for (sc in c("normal", "t3", "CN")) {
  rtau <- switch(sc, normal = function(n) rep(1, n),
                 t3 = function(n) 1 / rgamma(n, 1.5, 1.5),
                 CN = function(n) ifelse(runif(n) < 0.1, 10, 1))
  for (r in 1:R) {
    x1 <- rnorm(n); x2 <- rnorm(n); z <- rnorm(n)
    x <- cbind(1, x1, x2); w <- cbind(1, x1, x2, z)
    d <- rHeckmanSMN(n, beta, gamma, 1, rho, x, w, rtau)
    fn <- HeckmanNPSM(d$y, x, w, d$cc, fixG = TRUE, se = FALSE)
    fm <- HeckmanNPSM(d$y, x, w, d$cc, se = FALSE)
    res <- rbind(res,
                 data.frame(sc = sc, method = "Normal", b1 = fn$beta[2], b2 = fn$beta[3], rho = fn$rho, dl = 0, K = 1),
                 data.frame(sc = sc, method = "NPSM",   b1 = fm$beta[2], b2 = fm$beta[3], rho = fm$rho,
                            dl = fm$logL - fn$logL, K = nrow(fm$G)))
  }
}
agg <- function(v, t) sprintf("%.4f (%+.4f)", mean((v - t)^2), mean(v - t))
out <- do.call(rbind, lapply(split(res, list(res$method, res$sc)), function(d)
  data.frame(scenario = d$sc[1], method = d$method[1],
             `beta1 MSE(bias)` = agg(d$b1, 1), `beta2 MSE(bias)` = agg(d$b2, -1),
             `rho MSE(bias)` = agg(d$rho, rho), `mean dlogL` = round(mean(d$dl), 2),
             `mean K` = mean(d$K), check.names = FALSE)))
rownames(out) <- NULL
print(out[order(out$scenario), ])

fm$G
