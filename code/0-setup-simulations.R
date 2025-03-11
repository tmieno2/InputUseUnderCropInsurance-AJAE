## -------------------------------------------------
library("dplyr")
library("parallel")
library("data.table")
library("mgcv")
library("mc2d")
library("Rcpp")
library("quadprog")


## -------------------------------------------------
#--- crop insurance parameters and functions ---#
source("code/R/utility.R")

#--- sourcing cpp simulation function ---#
Rcpp::sourceCpp(here::here("code/R/calc_premium.cpp"))


## -------------------------------------------------
#--- county specific crop insurance parameters ---#
source("code/R/Sioux_IA_parameters.R")


## -------------------------------------------------
N_min <- 0
N_max <- 240
N_seq <- seq(N_min, N_max, by = 2)

search_increment <- 0.1
N_seq_for_search <- seq(N_min, N_max, by = search_increment)

N_data <-
  list(
    N_min = N_min,
    N_max = N_max,
    N_seq = N_seq,
    N_seq_for_search = N_seq_for_search
  )


## -------------------------------------------------
APH_min <- 100
APH_max <- 160
APH_list <- seq(APH_min, APH_max, by = 1)

APH_data <-
  list(
    APH_min = APH_min,
    APH_max = APH_max,
    APH_list = APH_list
  )


## -------------------------------------------------
prod_data <-
  list(
    #--- production acres (does not affect optimal N) ---#
    acres = 100,
    #--- crop function parameters ---#
    min_y = 48,
    max_y = 202,
    p_pars = c(3.14, -0.0921, 0.00603),
    q_pars = c(12.30, -1.353, 0.0456)
  )


## -------------------------------------------------
#--- total number of periods ---#
T <- 50

#--- Number of knots for Bernstein smoothing ---#
Nk <- 5

#--- Number of years to be used to calculate APH ---#
n_window <- 10

#--- discount rate ---#
disc <- 0.04


## -------------------------------------------------
#--- mean and variance of the corn price ---#
# do not use ln_var and ln_mean as they are used in premium calculation
# parameters defined here are used for simulating the actual realization
# of harvest price
p_price <- 2.20 # projected price
pvf <- 0.17535 # price volatility factor
ln_var_p <- log(pvf^2 + 1)
ln_mean_p <- log(p_price) - ln_var_p / 2
hp_prop <- c(ln_mean_p, sqrt(ln_var_p))

price_data <-
  list(
    p_price = p_price, # projected price
    pvf = pvf, # price volatility factor
    ln_var_p = ln_var_p,
    ln_mean_p = ln_mean_p,
    hp_prop = hp_prop,
    p_N = 0.15, # N price ($/lb)
    other_cost = 0 # ($/acre)
  )


## -------------------------------------------------
rho <- -0.3 # correlation coefficient
mu <- c(0, 0) # mean
sigma <- matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2) # covariance mat


## -------------------------------------------------
ci_data <-
  list(
    c_rate_now =
      county_rate[year == 2016, list(reference_amount, exponent_value, reference_rate, fixed_rate)] %>%
        unlist(),
    c_rate_prev =
      county_rate[year == 2015, list(reference_amount, exponent_value, reference_rate, fixed_rate)] %>%
        unlist(),
    c_pars = c(
      rate_differential_factor, py_rate_differential_factor,
      unit_residual_factor, py_unit_residual_factor
    ),
    beta_RP = c(
      beta_RP_0, beta_RP_1, beta_RP_2, beta_RP_3, beta_RP_4,
      beta_RP_5, beta_RP_6, beta_RP_7, beta_RP_8, beta_RP_9,
      beta_RP_10, beta_RP_11, beta_RP_12, beta_RP_13, beta_RP_14
    ),
    beta_RPHPE = c(
      beta_RPHPE_0, beta_RPHPE_1, beta_RPHPE_2, beta_RPHPE_3, beta_RPHPE_4,
      beta_RPHPE_5, beta_RPHPE_6, beta_RPHPE_7, beta_RPHPE_8, beta_RPHPE_9,
      beta_RPHPE_10, beta_RPHPE_11, beta_RPHPE_12, beta_RPHPE_13, beta_RPHPE_14
    ),
    sec5 = c(mean_quantity, sd_quantity),
    sec6 = c(
      capping_reference_yield, py_capping_reference_yield,
      capping_exponent_value, py_capping_exponent_value,
      capping_reference_rate, py_capping_reference_rate,
      capping_fixed_rate, py_capping_fixed_rate,
      capping_year
    ),
    county_TY = 170,
    #--- current year ---#
    c_year = 2016
  )

