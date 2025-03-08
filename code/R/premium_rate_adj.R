######################################
# Premium Rate Calculation
######################################
# written by Taro Mieno on 1/07/2016

# ===================================
# Premium calculation
# ===================================
premium_calc_ya <- function(APH_true, APH_adj, cov_level, proj_price, pvf, acres, county_rate, c_year) {
  #--------------------------
  # Section 1: Liability Calculation
  #--------------------------
  # maximum potential cost to insurer
  premium_liability_amount <- APH_adj * cov_level * proj_price * acres * insured_share_percent

  #--------------------------
  # Section 2: Unit Discount Calculation
  #--------------------------
  # see county specific rate structure.
  # ex. Chase_rate.xls for Chase county
  # optional_unit_discount_factor <- 1
  # unit_structure_discount_factor <- 1

  #--------------------------
  # Section 3: Bare Premium Rate Calculation
  #--------------------------
  # rate yield: basically the same as APH (see the handbook for the definition)
  # reference amount: county specific (look up the rate information)
  # prior year reference amount (py_reference_amount): county specific (look up the rate information)

  rate_yield <- APH_true
  cy_yield_ratio <- round(rate_yield / county_rate[year == c_year, reference_amount], digits = 2)
  py_yield_ratio <- round(rate_yield / county_rate[year == c_year - 1, reference_amount], digits = 2)

  # exponent value: county specific (look up the rate information)
  # prior year exponent value: county specific (look up the rate information)
  cy_rate_multiplier <- cy_yield_ratio^county_rate[year == c_year, exponent_value]
  py_rate_multiplier <- py_yield_ratio^county_rate[year == c_year - 1, exponent_value]

  # reference rate: county specific (look up the rate information)
  # prior year reference rate: county specific (look up the rate information)
  # fixed rate: county specific (look up the rate information)
  # prior year fixed rate: county specific (look up the rate information)
  cy_base_rate <- cy_rate_multiplier * county_rate[year == c_year, reference_rate] + county_rate[year == c_year, fixed_rate]
  py_base_rate <- py_rate_multiplier * county_rate[year == c_year - 1, reference_rate] + county_rate[year == c_year - 1, fixed_rate]

  # rate differential factor: what the hell is this?
  # unit residual factor: 1 for operational unit?
  cy_base_premium_rate <- cy_base_rate * rate_differential_factor * unit_residual_factor
  py_base_premium_rate <- py_base_rate * py_rate_differential_factor * py_unit_residual_factor

  base_premium_rate <- min(cy_base_premium_rate, py_base_premium_rate * 1.2)
  revenue_lookup_rate <- min(cy_base_rate, py_base_rate * 1.2)

  #--------------------------
  # Section 4: Optional Coverage Calculation
  #--------------------------
  # additive_optional_rate_adj_factor
  # multiplicative_optional_rate_adj_factor

  #--------------------------
  # Section 5: Revenue Coverage Add On Rates
  #--------------------------
  # mean quantity: ?
  # standard deviation quantity: ?
  # revenue_lookup_adj_factor: ?
  # price volatility factor (pvf): set nationally

  lookup_rate <- round(revenue_lookup_rate * revenue_lookup_adj_factor, digits = 4)
  adj_mean_quantity <- APH_adj * mean_quantity / 100
  adj_sd_quantity <- APH_adj * sd_quantity / 100
  ln_var <- log(pvf^2 + 1)
  ln_mean <- log(proj_price) - ln_var / 2

  #--- calculate simulated yield and revenue loss ---#
  loss <- loss_simu(100000, proj_price, APH_adj, cov_level, adj_mean_quantity, adj_sd_quantity, ln_mean, ln_var)

  simulated_yield_protection_base_premium_rate <- loss$sypl / (APH_adj * cov_level)
  simulated_revenue_protection_base_premium_rate <- loss$srpl / (APH_adj * cov_level * proj_price)

  preliminary_revenue_protection_add_on_rate <-
    max(
      simulated_revenue_protection_base_premium_rate - simulated_yield_protection_base_premium_rate,
      0.01 * base_premium_rate
    )

  #--------------------------
  # Section 6: Historical Revenue Capping
  #--------------------------
  # Note: capping by historical rate does not happen for
  # coverage level under 0.65.
  # capping reference yield?
  # prior capping reference yield?

  if (cov_level >= 0.65) {
    #--- 6.01 ---#
    capping_yield_ratio <- rate_yield / capping_reference_yield
    py_capping_yield_ratio <- rate_yield / py_capping_reference_yield

    #--- 6.02 ---#
    capping_rate_multiplier <- capping_yield_ratio^capping_exponent_value
    py_capping_rate_multiplier <- py_capping_yield_ratio^py_capping_exponent_value

    #--- 6.03 ---#
    historical_capping_base_rate <- capping_rate_multiplier * capping_reference_rate + capping_fixed_rate
    historical_py_capping_base_rate <- py_capping_rate_multiplier * py_capping_reference_rate + py_capping_fixed_rate

    #--- 6.04 ---#
    hisotircal_basic_unit_base_rate <- 0.9 * min(0.999, historical_py_capping_base_rate * 1.2, historical_capping_base_rate)

    #--- 6.05 ---#
    historical_revenue_protection_base_premium_rate <-
      (beta_0
      + round(beta_1 * hisotircal_basic_unit_base_rate, 8)
        + round(beta_2 * hisotircal_basic_unit_base_rate^2, 8)
        + round(beta_3 * cov_level, 8)
        + round(beta_4 * cov_level^2, 8)
        + round(beta_5 * (APH_adj / capping_reference_yield), 8)
        + round(beta_6 * (APH_adj / capping_reference_yield)^2, 8)
        + round(beta_7 * pvf, 8)
        + round(beta_8 * pvf^2, 8)
        + round(beta_9 * hisotircal_basic_unit_base_rate * cov_level, 8)
        + round(beta_10 * hisotircal_basic_unit_base_rate * (APH_adj / capping_reference_yield), 8)
        + round(beta_11 * hisotircal_basic_unit_base_rate * pvf, 8)
        + round(beta_12 * cov_level * (APH_adj / capping_reference_yield), 8)
        + round(beta_13 * cov_level * pvf, 8)
        + round(beta_14 * (APH_adj / capping_reference_yield) * pvf, 8)
      ) * unit_residual_factor * 1.1

    #--- 6.06 ---#
    a1 <- base_premium_rate + preliminary_revenue_protection_add_on_rate
    a2 <- historical_revenue_protection_base_premium_rate * 1.2^(c_year - capping_year)
    capped_revenue_protection_add_on_rate <- min(a1, a2) - base_premium_rate
  } else {
    capped_revenue_protection_add_on_rate <- preliminary_revenue_protection_add_on_rate
  }

  #--------------------------
  # Section 7: Whole Farm Unit Discount
  #--------------------------

  #--------------------------
  # Section 8: Premium Rate
  #--------------------------
  #--- 8.01 ---#
  capped_revenue_protection_add_on_factor <- capped_revenue_protection_add_on_rate

  #--- 8.02 ---#
  premium_rate <- min(base_premium_rate * unit_structure_discount_factor * multiplicative_optional_rate_adj_factor + additive_optional_rate_adj_factor + capped_revenue_protection_add_on_factor, 0.999)

  #--------------------------
  # Section 9: Total premium, subsidy, and producer premium calculation
  #--------------------------
  #--- 9.01 ---#
  preliminary_total_premium_amount <- premium_liability_amount * premium_rate * experience_factor * premium_surcharge_percent * total_premium_multiplicative_optional_rate_adj_factor

  #--- 9.02 ---#
  total_premium_amount <- preliminary_total_premium_amount * multiple_commodity_adj_factor

  #--- 9.03 ---#
  sub_per <- subsidy_percent[coverage_level == cov_level, subsidy_percent]
  base_subsidy_amount <- total_premium_amount * sub_per

  #--- 9.04 ---#
  beg_farmer_rancher_subsidy_amount <- total_premium_amount * beg_fr_subsidy_percent

  #--- 9.05 (Native Sod Subsidy)---#
  native_sod_subsidy_amount <- 0

  #--- 9.06 (Subsidy amount)---#
  subsidy_amount <- base_subsidy_amount + beg_farmer_rancher_subsidy_amount + native_sod_subsidy_amount

  #--- 9.07 (Producer premium amount) ---#
  premium <- total_premium_amount - subsidy_amount

  #--------------------------
  # Section 10: Miscellaneous Calculations
  #--------------------------
  # This is only for insurance agents.
  # AO_expense_subsidy_amount <- total_premium_amount*AO_expense_subsidy_percent

  #--------------------------
  # return value
  #--------------------------
  return(data.table(
    APH = APH_true,
    rate_yield = APH_adj,
    cov_level = cov_level,
    premium = premium,
    premium_rate = premium_rate,
    sub_rate = premium_rate * (1 - sub_per),
    subsidy = subsidy_amount,
    total_premium = total_premium_amount
  ))
}

# ===================================
# Sub functions
# ===================================
loss_simu <- function(N, proj_price, APH_adj, cov_level, adj_mean, adj_sd, ln_mean, ln_var) {
  sim_yield <- round(pmax(0, rnorm(N, mean = adj_mean, sd = adj_sd)), digits = 12)

  #--- simulated yield protection loss ---#
  sypl <- pmax(0, APH_adj * cov_level - sim_yield) %>% mean()

  #--- simulated harvest price ---#
  sim_price <- exp(rnorm(N) * sqrt(ln_var) + ln_mean)

  #--- simulated revenue protection loss ---#
  min_price <- pmin(2 * proj_price, sim_price)
  max_price <- pmax(proj_price, min_price)
  srpl <- pmax(0, APH_adj * cov_level * max_price - sim_yield * min_price) %>% mean()

  return(list(sypl = sypl, srpl = srpl))
}

# loss_simu(500,APH,cov_level,adj_mean,adj_sd,ln_mean,ln_var)
