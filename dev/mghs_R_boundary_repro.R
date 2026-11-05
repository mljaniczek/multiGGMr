# Minimal reproduction: mGHS chains with few observations per group run the
# correlation of edge scales R to +-1 under the paper's uniform prior, and stop
# at the global-scale guard. See NEWS for multiGGM 0.8.0 and ?multiggm_mcmc.
library(multiGGM)

sim <- simulate_multiggm(K = 2, p = 28, n = 100, seed = 7)

# 1. Uniform prior on R (default): stops, and the error reports 1 - R12^2.
try(multiggm_mcmc(data_list = sim$data_list, method = "mghs", burnin = 2000,
                  nsave = 2000, seed = 123, max_attempts = 1, save_draws = FALSE))

# 2. With n = 1000 per group the same model keeps R inside.
big <- simulate_multiggm(K = 2, p = 28, n = 1000, seed = 7)
f_big <- multiggm_mcmc(data_list = big$data_list, method = "mghs", burnin = 2000,
                       nsave = 2000, seed = 123, max_attempts = 1, save_draws = FALSE)
summary(f_big$R_save[1, 2, ])

# 3. An LKJ(10) prior on R keeps the small-n chain away from the boundary.
f_lkj <- multiggm_mcmc(data_list = sim$data_list, method = "mghs", burnin = 2000,
                       nsave = 2000, seed = 123, max_attempts = 1, save_draws = FALSE,
                       hyper = list(R_eta = 10))
summary(f_lkj$R_save[1, 2, ])
