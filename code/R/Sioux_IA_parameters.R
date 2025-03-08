# ===================================
# Exogenous Parameters
# ===================================
#--------------------------
# section 1
#--------------------------
#--- Projected Price ---#
# vary by year, set nationally
# proj_price <- 4.15

# This is farmers' choice. But, it seems okay to
# assume farmers insure the whole acres
insured_share_percent <- 1

#--------------------------
# section 2
#--------------------------
# ask Cory
optional_unit_discount_factor <- 1
unit_structure_discount_factor <- 1

#--------------------------
# section 3
#--------------------------
# This information is available at ADM "A01040"
# It should be okay to assume these are fixed over time
rate_differential_factor <- 2.057
py_rate_differential_factor <- 1.545
unit_residual_factor <- 1.00
py_unit_residual_factor <- 1.00

#--------------------------
# section 4
#--------------------------
# Where can we obtain there information?
# It should be okay to assume these are fixed over time
additive_optional_rate_adj_factor <- 0
multiplicative_optional_rate_adj_factor <- 1

#--------------------------
# section 5
#--------------------------
revenue_lookup_adj_factor <- 1

#--- mean quantity and sd quantity ---#
# information available at ADM
mean_quantity <- 100
sd_quantity <- 25.275654190

#--- price volatility factor ---#
# vary by year, set nationally
# pvf <- 0.21

#--------------------------
# section 6
#--------------------------
capping_reference_yield <- 157
py_capping_reference_yield <- 157
capping_exponent_value <- -2.051
py_capping_exponent_value <- -2.050
capping_reference_rate <- 0.0120
py_capping_reference_rate <- 0.0120
capping_fixed_rate <- 0.008
py_capping_fixed_rate <- 0.008
capping_year <- 2011

#--- beta RPHPE ---#
# Information available at ADM
# Assume fixed?
beta_RPHPE_0 <- -0.079815390
beta_RPHPE_1 <- 1.219518
beta_RPHPE_2 <- -0.4361823
beta_RPHPE_3 <- -0.062436830
beta_RPHPE_4 <- 0.2539943
beta_RPHPE_5 <- 0.07154643
beta_RPHPE_6 <- -0.00367864
beta_RPHPE_7 <- -0.2409182
beta_RPHPE_8 <- 0.384326
beta_RPHPE_9 <- -0.01062107
beta_RPHPE_10 <- -0.04992516
beta_RPHPE_11 <- -0.1323247
beta_RPHPE_12 <- -0.09314001
beta_RPHPE_13 <- 0.2016753
beta_RPHPE_14 <- 0.00120584

#--- beta RP ---#
# Information available at ADM
# Assume fixed?
beta_RP_0 <- -0.096524510
beta_RP_1 <- 1.393955000
beta_RP_2 <- -0.653385200
beta_RP_3 <- -0.052425310
beta_RP_4 <- 0.273246100
beta_RP_5 <- 0.074885210
beta_RP_6 <- 0.001166760
beta_RP_7 <- -0.312272900
beta_RP_8 <- 0.269245800
beta_RP_9 <- -0.226560800
beta_RP_10 <- 0.043531740
beta_RP_11 <- 0.503836700
beta_RP_12 <- -0.110972300
beta_RP_13 <- 0.515274900
beta_RP_14 <- -0.032281610

#--------------------------
# section 9
#--------------------------
# assume fixed
experience_factor <- 1
premium_surcharge_percent <- 1
total_premium_multiplicative_optional_rate_adj_factor <- 1
multiple_commodity_adj_factor <- 1
beg_fr_subsidy_percent <- 0
native_sod_subsidy_amount <- 0

# ===================================
# County-specific parameters
# ===================================
county_rate <- data.table(
  year = c(2015, 2016),
  reference_amount = c(79, 80),
  exponent_value = c(-1.104, -0.983),
  reference_rate = c(0.022, 0.018),
  fixed_rate = c(0.006, 0.006)
)
