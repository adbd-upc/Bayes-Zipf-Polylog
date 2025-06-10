# load posterior draws and data
load(url("https://raw.githubusercontent.com/adbd-upc/Bayes-Zipf-Polylog/main/SEIO-2025/R/draws_poly.RData"))
load(url("https://raw.githubusercontent.com/adbd-upc/Bayes-Zipf-Polylog/main/SEIO-2025/R/URV_data.RData"))
urv = with(degreeSeq, rep(Var1, Freq))
y = urv

library(mvtnorm)
library(ggplot2)

draws = as.data.frame(draws)
n = length(y)
# MLE as in Valero et al
alpha_hat = 0.1774
beta_hat  = 0.9108

K   = max(2000, ceiling(max(y)*5))
ks  = seq_len(K)
logks = log(ks)

# approx
w = beta_hat^ks / (ks^alpha_hat)
Li = sum(w)              
pks = w / Li              

E_Y = sum( ks      * pks )
E_Y2 = sum((ks^2)   * pks )
E_logY = sum( logks    * pks )
E_logY2 = sum((logks^2)* pks )
E_Y_logY = sum( ks*logks * pks )

Var_Y = E_Y2 - E_Y^2
Var_logY = E_logY2 - E_logY^2
Cov_Y_logY = E_Y_logY - E_Y * E_logY


I11 = Var_logY
I22 = Var_Y / beta_hat^2
I12 = - Cov_Y_logY / beta_hat

I_mat = matrix(c(I11, I12,
                  I12, I22),
                nrow = 2, byrow = TRUE)
Sigma_bvm = solve(n * I_mat)

alpha_seq = seq(min(draws$alpha), max(draws$alpha), length = 100)
beta_seq  = seq(min(draws$beta),  max(draws$beta),  length = 100)
grid = expand.grid(alpha = alpha_seq, beta = beta_seq)

grid_mat = as.matrix(grid[, c("alpha", "beta")])
grid$dens = dmvnorm(grid_mat,
                     mean  = c(alpha_hat, beta_hat),
                     sigma = Sigma_bvm)

p_joint = ggplot() +
  geom_point(aes(alpha, beta), data = draws, alpha = 0.3) +
  geom_contour(aes(alpha, beta, z = dens),
               data = grid, color = "red") +
  geom_point(aes(x = 0.1774, y = 0.9108),
             color = "red", size = 3, shape = 4, stroke = 1.5) +
  labs(
    title = "Slice-Sampler Draws vs BvM Gaussian Contours",
    subtitle = "Red cross is the MLE",
    x     = expression(alpha),
    y     = expression(beta)
  ) +
  theme_minimal()

# marignals
p_alpha = ggplot(draws, aes(x = alpha)) +
  geom_histogram(aes(y = ..density..),
                 bins = 30, fill = "grey80", color = "white") +
  stat_function(fun = dnorm,
                args = list(mean = alpha_hat,
                            sd   = sqrt(Sigma_bvm[1,1])),
                color = "blue", size = 1) +
  labs(title = "Marginal α: sample vs. BvM normal") +
  theme_minimal()

p_beta = ggplot(draws, aes(x = beta)) +
  geom_histogram(aes(y = ..density..),
                 bins = 30, fill = "grey80", color = "white") +
  stat_function(fun = dnorm,
                args = list(mean = beta_hat,
                            sd   = sqrt(Sigma_bvm[2,2])),
                color = "blue", size = 1) +
  labs(title = "Marginal β: sample vs. BvM normal") +
  theme_minimal()

print(p_joint)
print(p_alpha)
print(p_beta)
