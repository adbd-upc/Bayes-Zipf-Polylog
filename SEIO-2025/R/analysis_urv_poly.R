# source sampler
source(url("https://raw.githubusercontent.com/adbd-upc/Bayes-Zipf-Polylog/main/SEIO-2025/R/slice_urv_poly.R"))

# load data
load(url("https://raw.githubusercontent.com/adbd-upc/Bayes-Zipf-Polylog/main/SEIO-2025/R/URV_data.RData"))
urv = with(degreeSeq, rep(Var1, Freq))




urv = with(degreeSeq, rep(Var1, Freq))


set.seed(2025)
y = urv


draws = polylog_slice(y,
                      n_iter=18000,
                      burn_in=3000,
                      w_alpha=1,
                      w_beta=1,
                      m=100,
                      alpha_0=0.1774,
                      beta_0=0.9108,
                      progress=TRUE)


###############################
# posterior summaries Polylog #
###############################

print(apply(draws,2,quantile,probs=c(0.025,0.5,0.975)))


library(ggplot2)
library(viridis)

df_draws = as.data.frame(draws)
colnames(df_draws) = c("alpha","beta")

ggplot(df_draws, aes(x = alpha, y = beta)) +
  stat_density2d(
    aes(fill = ..density..),
    geom = "raster",
    contour = FALSE,
    n = 200
  ) +
  scale_fill_viridis_c(option = "plasma", name = "Density") +
  geom_point(alpha = 0.2, size = 0.5, color = "black") +
  geom_point(aes(x = 0.1774, y = 0.9108),
             color = "red", size = 3, shape = 4, stroke = 1.5) +
  labs(
    x = expression(alpha),
    y = expression(beta),
    title = "Posterior samples of (α,β) with density heatmap",
    subtitle = "Red cross is the MLE"
  ) +
  theme_minimal(base_size = 14)


sim_polylog_fast = function(N, alpha, beta, ymax) {
  ks  = seq_len(ymax)
  pmf = beta^ks / (ks^alpha)
  pmf = pmf / sum(pmf)
  sample(ks, size = N, replace = TRUE, prob = pmf)
}

###############################
# posterior predictive checks #
###############################

set.seed(123)
n_sims = 15000
idx     = sample(nrow(draws), n_sims)
ymax    = max(y)
n_obs   = length(y)

pp_mat = matrix(0, nrow = n_sims, ncol = ymax)
for (i in seq_len(n_sims)) {
  a     = draws[idx[i], "alpha"]
  b     = draws[idx[i], "beta"]
  y_rep = sim_polylog_fast(n_obs, a, b, ymax)
  pp_mat[i, ] = tabulate(y_rep, nbins = ymax) / n_obs
}

mean_rel = colMeans(pp_mat)
lwr      = apply(pp_mat, 2, quantile, probs = 0.025)
upr      = apply(pp_mat, 2, quantile, probs = 0.975)

obs_tab = tabulate(y, nbins = ymax)
obs_rel = obs_tab / sum(obs_tab)

df = data.frame(
  y        = 1:ymax,
  mean_rel = mean_rel,
  lwr      = lwr,
  upr      = upr,
  obs_rel  = obs_rel
)
coverage = mean(df$obs_rel >= df$lwr & df$obs_rel <= df$upr)

ggplot(df, aes(x = y)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = "95% predictive interval"),
              alpha = 0.2) +
  geom_line(aes(y = mean_rel, color = "Mean predicted rel. freq"), size = 1) +
  geom_point(aes(y = obs_rel, color = "Observed rel. freq"), size = 2) +
  scale_fill_manual(name = NULL,
                    values = c("95% predictive interval" = "steelblue")) +
  scale_color_manual(name = NULL,
                     values = c("Mean predicted rel. freq" = "steelblue",
                                "Observed rel. freq"       = "red")) +
  labs(
    title   = "Posterior–Predictive Check: Relative Frequencies",
    subtitle= paste0("Coverage: ",
                     scales::percent(coverage, accuracy = 0.1),
                     " of observed values within 95% bands"),
    x       = "Value of y",
    y       = "Relative frequency"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 12)
  )
