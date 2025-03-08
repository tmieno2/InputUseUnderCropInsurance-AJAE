# ===================================
# Set parameters that are likely to stay across simulations
# ===================================
#--------------------------
# N
#--------------------------
N_min <- 0
N_max <- 240
N_seq <- seq(N_min, N_max, by = 2)
N_len <- length(N_seq)

search_increment <- 0.1
search_N <- seq(min(N_seq), max(N_seq), by = search_increment)

#--------------------------
# APH
#--------------------------
APH_min <- 100
APH_max <- 160
APH_list <- seq(APH_min, APH_max, by = 1)
APH_len <- length(APH_list)

#--------------------------
# Total number of periods
#--------------------------
T <- 50

#--------------------------
# Number of knots for Bernstein
#--------------------------
Nk <- 5

#--------------------------
# Production parameters
#--------------------------
acres <- 100

#--- crop function parameters ---#
min_y <- 48
max_y <- 202
p_pars <- c(3.14, -0.0921, 0.00603)
q_pars <- c(12.30, -1.353, 0.0456)

#--------------------------
# Current year
#--------------------------
c_year <- 2016

#--------------------------
# Price
#--------------------------
#--- mean and variance of the corn price ---#
# do not use ln_var and ln_mean as they are used in premium calculation
# parameters defined here are used for simulating the actual realization
# of harvest price
p_price <- 2.20 # projected price
pvf <- 0.17535 # price volatility factor
ln_var_p <- log(pvf^2 + 1)
ln_mean_p <- log(p_price) - ln_var_p / 2
hp_prop <- c(ln_mean_p, sqrt(ln_var_p))

#--- N price ---#
p_N <- 0.15 # ($/lb)

#--- non-N costs (other costs) ---#
other_cost <- 0 # ($/acre)

#--------------------------
# parameters of joint price-yield distribution
#--------------------------
rho <- -0.3 # correlation coefficient
mu <- c(0, 0) # mean
sigma <- matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2) # covariance mat

#--------------------------
# Insurance parameters
#--------------------------
c_rate_now <- county_rate[year == c_year, list(reference_amount, exponent_value, reference_rate, fixed_rate)] %>%
  unlist()
c_rate_prev <- county_rate[year == c_year - 1, list(reference_amount, exponent_value, reference_rate, fixed_rate)] %>%
  unlist()
c_pars <- c(
  rate_differential_factor, py_rate_differential_factor,
  unit_residual_factor, py_unit_residual_factor
)
beta_RP <- c(
  beta_RP_0, beta_RP_1, beta_RP_2, beta_RP_3, beta_RP_4,
  beta_RP_5, beta_RP_6, beta_RP_7, beta_RP_8, beta_RP_9,
  beta_RP_10, beta_RP_11, beta_RP_12, beta_RP_13, beta_RP_14
)
beta_RPHPE <- c(
  beta_RPHPE_0, beta_RPHPE_1, beta_RPHPE_2, beta_RPHPE_3, beta_RPHPE_4,
  beta_RPHPE_5, beta_RPHPE_6, beta_RPHPE_7, beta_RPHPE_8, beta_RPHPE_9,
  beta_RPHPE_10, beta_RPHPE_11, beta_RPHPE_12, beta_RPHPE_13, beta_RPHPE_14
)
sec5 <- c(mean_quantity, sd_quantity)
sec6 <- c(
  capping_reference_yield, py_capping_reference_yield,
  capping_exponent_value, py_capping_exponent_value,
  capping_reference_rate, py_capping_reference_rate,
  capping_fixed_rate, py_capping_fixed_rate,
  capping_year
)
county_TY <- 170

#--- list of insurance types ---#
insurance_type <- c("YP", "RPHPE", "RP", "non")
ins_len <- length(insurance_type)

#--------------------------
# Utility types
#--------------------------
utility_type <- c("RN", "CA1", "CA2", "CR")
uti_len <- length(utility_type)

#--------------------------
# Subsidy Percent Lookup
#--------------------------
subsidy_percent <- data.table(
  coverage_level = c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.90),
  subsidy_percent = c(0.67, 0.64, 0.64, 0.59, 0.59, 0.55, 0.48, 0.38, 0.28)
)

#--------------------------
# Number of years to be used to calculate APH
#--------------------------
n_window <- 10
