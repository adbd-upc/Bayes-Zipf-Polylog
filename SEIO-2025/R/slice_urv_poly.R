Li = function(alpha, beta, tol=1e-12, max_terms=1e7) {
  beta = as.numeric(beta)
  k = length(beta)
  s = numeric(k)
  below  = !is.na(beta) & beta <= 0
  above  = !is.na(beta) & beta >= 1
  inside = !is.na(beta) & beta > 0 & beta < 1
  s[below] = 0
  s[above] = Inf
  if (any(inside)) {
    b_vec   = beta[inside]
    m       = length(b_vec)
    partial = numeric(m)
    term    = b_vec
    active  = rep(TRUE, m)
    y       = 1L
    while (any(active)) {
      t = term[active] / (y^alpha)
      partial[active] = partial[active] + t
      active[active]  = abs(t) >= tol
      y = y + 1L
      if (y > max_terms) stop("Li(): series failed to converge")
      term[active] = term[active] * b_vec[active]
    }
    s[inside] = partial
  }
  s
}

safe_Li = function(alpha, beta) {
  out = try(Li(alpha, beta), silent=TRUE)
  if (inherits(out, "try-error")) return(rep(NA_real_, length(beta)))
  out
}

stats_T = function(y) list(T1 = sum(y), T2 = sum(log(y)), n = length(y))

slice1 = function(x0, logf, w=1, m=100, lower=-Inf, upper=Inf) {
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

logf_alpha = function(a, beta, T2, n) {
  li = safe_Li(a, beta)
  if (!is.finite(li) || li <= 0) return(-Inf)
  -(a * T2) - n * log(li)
}

logf_u = function(u, alpha, T1, n) {
  # transform back
  beta = 1 / (1 + exp(-u))
  # jacobian: 
  J = log(beta) + log1p(-beta)
  li = safe_Li(alpha, beta)
  if (!is.finite(li) || li <= 0) return(-Inf)
  (T1 - 1) * log(beta) - n * log(li) + J
}

polylog_slice = function(y,
                          n_iter=12000,
                          burn_in=2000,
                          w_alpha=1,
                          w_beta=1,
                          m=50,
                          alpha_0=0,
                          beta_0=0.5,
                          progress=TRUE) {
  S  = stats_T(y); T1 = S$T1; T2 = S$T2; n = S$n
  alpha = numeric(n_iter); beta = numeric(n_iter)
  alpha[1] = alpha_0
  beta[1]  = beta_0
  u        = log(beta_0 / (1 - beta_0))
  for (t in 2:n_iter) {
    # slice‐sample alpha ∈ ℝ
    alpha[t] = slice1(alpha[t-1],
                       function(a) logf_alpha(a, beta[t-1], T2, n),
                       w=w_alpha, m=m, lower=-Inf, upper=Inf)
    # slice‐sample u = logit(β)
    u        = slice1(u,
                       function(u0) logf_u(u0, alpha[t], T1, n),
                       w=w_beta, m=m, lower=-Inf, upper=Inf)
    beta[t]  = 1 / (1 + exp(-u))
    if (progress && t %% 1000 == 0)
      cat(sprintf("iter %4d: α=%.3f  β=%.3f\n", t, alpha[t], beta[t]))
  }
  draws = cbind(alpha, beta)
  draws[-seq_len(burn_in), ]
}

