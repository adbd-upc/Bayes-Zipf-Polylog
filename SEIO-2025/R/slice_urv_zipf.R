
if (!requireNamespace("gsl", quietly=TRUE))
  install.packages("gsl", repos="https://cloud.r-project.org")
library(gsl)

slice1_zipf = function(x0, logf, w=1, m=50, lower=-Inf, upper=Inf) {
  f0 = logf(x0)
  if (!is.finite(f0)) return(x0)
  logy = f0 - rexp(1)
  # stepping out
  u = runif(1, 0, w)
  L = x0 - u; R = L + w
  J = floor(runif(1, 0, m)); K = m - 1 - J
  while (J > 0 && L > lower && logf(L) > logy) { L = L - w; J = J - 1 }
  while (K > 0 && R < upper && logf(R) > logy) { R = R + w; K = K - 1 }
  # shrinkage
  repeat {
    x1 = runif(1, L, R)
    f1 = logf(x1)
    if (f1 >= logy) return(x1)
    if (x1 < x0) L = x1 else R = x1
  }
}


logf_alpha_zipf = function(alpha, T2, n) {
  if (alpha <= 1) return(-Inf)
  z = zeta(alpha)    # gsl::zeta
  if (!is.finite(z) || z <= 0) return(-Inf)
  -alpha * T2 - n * log(z)
}

zipf_slice = function(y,
                       n_iter  = 5000,
                       burn_in = 1000,
                       w       = 1.0,
                       m       = 50,
                       alpha_0 = 2.0,
                       progress = TRUE) {
  T2 = sum(log(y))
  n  = length(y)
  alpha = numeric(n_iter)
  alpha[1] = alpha_0
  for (t in 2:n_iter) {
    alpha[t] = slice1_zipf(alpha[t-1],
                       function(a) logf_alpha_zipf(a, T2, n),
                       w = w, m = m,
                       lower = 1 + 1e-8, upper = Inf)
    if (progress && t %% 1000 == 0)
      cat(sprintf("iter %4d: α = %.4f\n", t, alpha[t]))
  }
  alpha[-seq_len(burn_in)]
}
