# source sampler
source(url("https://raw.githubusercontent.com/adbd-upc/Bayes-Zipf-Polylog/main/SEIO-2025/R/slice_urv_zipf.R"))

# load data
load(url("https://raw.githubusercontent.com/adbd-upc/Bayes-Zipf-Polylog/main/SEIO-2025/R/URV_data.RData"))
urv = with(degreeSeq, rep(Var1, Freq))


set.seed(2025)
y = urv

alpha_post = zipf_slice(y,
                         n_iter  = 18000,
                         burn_in = 3000,
                         w       = 1.0,
                         m       = 50,
                         alpha_0 = 2.0,
                         progress = FALSE)

############################
# posterior summaries Zipf #
############################


library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)

print(quantile(alpha_post, probs = c(0.05, 0.5, 0.95)))

df_alpha = data.frame(alpha = alpha_post)
ggplot(df_alpha, aes(x = alpha)) +
  geom_histogram(aes(y = ..density..), bins = 40, fill = "grey80", color = "black") +
  geom_density(color = "steelblue", size = 1) +
  labs(
    x     = expression(alpha),
    y     = "Density",
    title = "Posterior Density of α (Zipf Model)"
  ) +
  theme_minimal(base_size = 14)


sim_zipf_fast = function(N, alpha, ymax) {
  ks  = seq_len(ymax)
  pmf = ks^(-alpha)
  pmf = pmf / sum(pmf)
  sample(ks, size = N, replace = TRUE, prob = pmf)
}

set.seed(123)
n_sims = 15000
idx     = sample(length(alpha_post), n_sims)
ymax    = max(y)
n_obs   = length(y)

pp_mat = matrix(0, nrow = n_sims, ncol = ymax)
for (i in seq_len(n_sims)) {
  a     = alpha_post[idx[i]]
  y_rep = sim_zipf_fast(n_obs, a, ymax)
  pp_mat[i, ] = tabulate(y_rep, nbins = ymax) / n_obs
}

mean_rel = colMeans(pp_mat)
lwr = apply(pp_mat, 2, quantile, probs = 0.025)
upr = apply(pp_mat, 2, quantile, probs = 0.975)

obs_tab = tabulate(y, nbins = ymax)
obs_rel = obs_tab / sum(obs_tab)

df_ppc = data.frame(
  y = 1:ymax,
  mean_rel = mean_rel,
  lwr = lwr,
  upr = upr,
  obs_rel = obs_rel
)
coverage = mean(df_ppc$obs_rel >= df_ppc$lwr & df_ppc$obs_rel <= df_ppc$upr)

ggplot(df_ppc, aes(x = y)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = "95% predictive interval"),
              alpha = 0.2) +
  geom_line(aes(y = mean_rel, color = "Mean predicted"), size = 1) +
  geom_point(aes(y = obs_rel, color = "Observed"), size = 2) +
  scale_fill_manual(name = NULL,
                    values = c("95% predictive interval" = "steelblue")) +
  scale_color_manual(name = NULL,
                     values = c("Mean predicted" = "steelblue",
                                "Observed"       = "red")) +
  labs(
    title    = "Posterior–Predictive Check: Zipf Relative Frequencies",
    subtitle = paste0("Coverage: ", scales::percent(coverage, accuracy=0.1),
                      " of observed values within 95% bands"),
    x        = "Value of y",
    y        = "Relative frequency"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

