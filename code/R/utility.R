#++++++++++++++++++++++++++++++++++++
#+ Generate utility function
#++++++++++++++++++++++++++++++++++++
utility_gen <- function(utility_type = "RN") {
  if (utility_type == "CA1") { # CAA 1
    U <- function(x) {
      u <- 1 - exp(-0.0046 * x / 100)
      return(u)
    }
  } else if (utility_type == "CA2") { # CAA2
    U <- function(x) {
      u <- 1 - exp(-0.0100 * x / 100)
      return(u)
    }
  } else if (utility_type == "CR") { # CRA
    U <- function(x) {
      u <- (x^(1 - 0.6) - 1) / (1 - 0.6)
      return(u)
    }
  } else { # risk neutral
    U <- function(x) {
      u <- x
      return(u)
    }
  }
  return(U)
}

#++++++++++++++++++++++++++++++++++++
#+ Find optimal N by GAM smoothing
#++++++++++++++++++++++++++++++++++++
find_opt_N_by_gam <- function(data, dep_var = "value", N_seq_for_search_ls) {
  w_data <- copy(data)
  setnames(w_data, dep_var, "value")

  gam_e_u <- gam(value ~ s(N, k = 5), data = w_data)

  #--- optimal N ---#
  pi_data <-
    data.table(N = N_seq_for_search) %>%
    .[, e_u := predict.gam(gam_e_u, newdata = .)] %>%
    .[e_u == max(e_u), ] %>%
    setnames(c("N", "e_u"), c("N_star", "max_eu"))

  return(pi_data)
}

#++++++++++++++++++++++++++++++++++++
#+ Smooth premium payment and premium rate
#++++++++++++++++++++++++++++++++++++
smooth_premium <- function(data) {
  #---------------------
  #- Premium
  #---------------------
  premium_smooth_gam <- gam(premium ~ s(APH, k = 3), data = data)

  premium_smoothed <- predict(premium_smooth_gam, data.frame(APH = data$APH))

  data$premium <- premium_smoothed

  #---------------------
  #- Premium rate
  #---------------------
  premium_rate_smooth_gam <- gam(premium_rate ~ s(APH, k = 3), data = data)

  premium_rate_smoothed <- predict(premium_rate_smooth_gam, data.frame(APH = data$APH))

  data$premium_rate <- premium_rate_smoothed

  return(data)
}

# !===========================================================
# ! Generate yield-hprice data
# !===========================================================
generate_yield_hprice_data <- function(B, sigma, price_data) {
  #--- generate price and yield in U() ---#
  half_U <- rmvnorm(mean = c(0, 0), sig = sigma, n = B / 2) %>% pnorm()

  U <- rbind(half_U, 1 - half_U)

  h_price <- qlnorm(U[, 1], meanlog = price_data$ln_mean_p, sdlog = sqrt(price_data$ln_var_p))

  yield_unif <- U[, 2]

  #--- generate yield and harvest price ---#
  yield_hprice_data <-
    data.table(
      yield_u = yield_unif,
      h_price = h_price
    )

  return(yield_hprice_data)
}

expand_grid_df <- function(data_1, data_2) {
  data_1_ex <-
    data_1[rep(1:nrow(data_1), each = nrow(data_2)), ] %>%
    data.table() %>%
    .[, rowid := 1:nrow(.)]

  data_2_ex <-
    data_2[rep(1:nrow(data_2), nrow(data_1)), ] %>%
    data.table() %>%
    .[, rowid := 1:nrow(.)]

  expanded_data <-
    data_1_ex[data_2_ex, on = "rowid"] %>%
    .[, rowid := NULL]

  if ("tbl" %in% class(data_1)) {
    expanded_data <- as_tibble(expanded_data)
  }

  if ("rowwise_df" %in% class(data_1)) {
    expanded_data <- rowwise(expanded_data)
  }

  return(expanded_data)
}

#++++++++++++++++++++++++++++++++++++
#+ Calculate premium
#++++++++++++++++++++++++++++++++++++
get_premium <- function(insurance_type, APH, cov_level, subsidy_percent, price_data, ci_data, acres) {
  rate_yield <- APH

  county_TY <- ci_data$county_TY
  c_year <- ci_data$c_year
  c_rate_now <- ci_data$c_rate_now
  c_rate_prev <- ci_data$c_rate_prev
  c_pars <- ci_data$c_pars
  beta_RPHPE <- ci_data$beta_RPHPE
  beta_RP <- ci_data$beta_RP
  sec5 <- ci_data$sec5
  sec6 <- ci_data$sec6
  pp <- price_data$p_price
  pvf <- price_data$pvf

  if (insurance_type == "YP") {
    premium <- premium_calc_YP(
      rate_yield, APH, cov_level, subsidy_percent, pp,
      pvf, acres, county_TY, c_year, c_rate_now, c_rate_prev,
      c_pars
    )
  } else if (insurance_type == "RPHPE") {
    premium <- premium_calc_RPHPE(
      rate_yield, APH, cov_level, subsidy_percent, pp,
      pvf, acres, county_TY, c_year, c_rate_now, c_rate_prev,
      c_pars, beta_RPHPE, sec5, sec6
    )
  } else if (insurance_type == "RP") {
    premium <- premium_calc_RP(
      rate_yield, APH, cov_level, subsidy_percent, pp,
      pvf, acres, county_TY, c_year, c_rate_now, c_rate_prev,
      c_pars, beta_RP, sec5, sec6
    )
  }

  return(premium)
}

#++++++++++++++++++++++++++++++++++++
#+ Generate function to calculate revenue by insurance type
#++++++++++++++++++++++++++++++++++++
gen_revenue_function <- function(insurance_type) {
  if (insurance_type == "YP") {
    get_revenue <- function(raw_revenue, price_data, prod_data, APH, cov_level, yield, h_price) {
      raw_revenue + pmax(0, price_data$p_price * (APH * cov_level - yield) * prod_data$acres)
    }
  } else if (insurance_type == "RPHPE") {
    get_revenue <- function(raw_revenue, price_data, prod_data, APH, cov_level, yield, h_price) {
      revenue <- pmax(raw_revenue, APH * cov_level * prod_data$acres * price_data$p_price)
    }
  } else if (insurance_type == "RP") {
    get_revenue <- function(raw_revenue, price_data, prod_data, APH, cov_level, yield, h_price) {
      revenue <- pmax(raw_revenue, pmax(h_price, price_data$p_price) * APH * cov_level * prod_data$acres)
    }
  }
  return(get_revenue)
}


# !===========================================================
# ! Bernstenin bases creation and regression
# !===========================================================
#++++++++++++++++++++++++++++++++++++
#+ Bernstein Basis
#++++++++++++++++++++++++++++++++++++
BernBasis <- function(x, N) {
  # ptm <- proc.time()
  B <- do.call(cbind, lapply(0:N, function(t) factorial(N) / factorial(N - t) / factorial(t) * x^t * (1 - x)^(N - t)))
  # proc.time() - ptm

  return(B)
}

#++++++++++++++++++++++++++++++++++++
#+ Bernstein Basis
#++++++++++++++++++++++++++++++++++++
BernBasis_fast <- function(x, N, core = 5) {
  num_each <- length(x) / core

  BernBasis_by_part <- function(i) {
    start <- 1 + (i - 1) * num_each
    end <- i * num_each
    y <- x[start:end]
    B_temp <- do.call(cbind, lapply(0:N, function(t) factorial(N) / factorial(N - t) / factorial(t) * y^t * (1 - y)^(N - t)))
    return(B_temp)
  }

  B <- do.call(rbind, mclapply(1:core, BernBasis_by_part, mc.cores = core))
  return(B)
}

#++++++++++++++++++++++++++++++++++++
#+ Semi-parametric regression using Bernstein Basis
#++++++++++++++++++++++++++++++++++++

semi_B <- function(y, x, Nk, decrease = TRUE) {
  X <- BernBasis(x, Nk)

  XX <- t(X) %*% X
  L <- chol(XX)
  z <- solve(t(L)) %*% t(X) %*% y
  D <- diag(dim(XX)[1]) * 2
  dvec <- 2 * z

  #--- constraint matrix ---#
  if (decrease) {
    S <- cbind(diag(Nk), rep(0, Nk)) + cbind(rep(0, Nk), -diag(Nk))
  } else {
    S <- cbind(-diag(Nk), rep(0, Nk)) + cbind(rep(0, Nk), diag(Nk))
  }
  A <- S %*% solve(L)

  sol <- solve.QP(Dmat = D, dvec, Amat = t(A))
  beta <- solve(L) %*% sol$solution
  return(beta)
}

# !===========================================================
# ! Backward induction
# !===========================================================
BI <- function(VF_T, production_data, APH_data) {
  #++++++++++++++++++++++++++++++++++++
  #+ Setup (Create place holders)
  #++++++++++++++++++++++++++++++++++++
  #--- coefficients of semi-parametric regressions with Bernstein bases that allow for mapping APH to V()  ---#
  VF <- vector(mode = "list", length = T)

  #--- GAM regression results that map APH to optimal N at each t (optimal decision) ---#
  N_star_gam <- vector(mode = "list", length = T)

  #--- highest APH level to consider at each t ---#
  APH_threshold <- vector(mode = "list", length = T)

  APH_mat <- get_APH_Bernstein_mat(production_data, APH_data)

  #++++++++++++++++++++++++++++++++++++
  #+ t = T (final period)
  #++++++++++++++++++++++++++++++++++++
  # Run semi-parametric regression using Bernstein polynomials. The resulting coefficient estimates on the bases are stored.
  #---------------------
  #- Find the mapping of APH to profit/utility at t = T
  #---------------------
  VF[[1]] <- semi_B(VF_T[, max_eu], VF_T[, (APH - APH_data$APH_min) / (APH_data$APH_max - APH_data$APH_min)], Nk, decrease = FALSE)

  #--------------------------
  #- Find pptimal N conditional on APH at t = T
  #--------------------------
  #--- identify the highest value of APH at which optimal N is NOT 0 ---#
  if (any(VF_T[, N_star == 0])) {
    APH_threshold[[1]] <- max(VF_T[N_star != 0, APH])
  } else {
    APH_threshold[[1]] <- 999
  }

  #--- semi-parameteric regresion ---#
  # the resulting gam object can map (predict) APH level to optimal N
  N_star_gam[[1]] <- gam(N_star ~ s(APH, k = 5), data = VF_T[APH < APH_threshold[[1]], ])

  #++++++++++++++++++++++++++++++++++++
  #+ t < T
  #++++++++++++++++++++++++++++++++++++
  for (s in 1:(T - 1)) { # s=1 is T-1

    print(paste0(s, "/", T - 1))

    #---------------------
    #- Find optimal N conditional on APH at t = T - s
    #---------------------
    opt_N_on_APH <-
      production_data %>%
      .[, sV := utility + APH_mat %*% VF[[s]] / (1 + disc)] %>%
      #--- find mean VF by N and APH ---#
      .[, .(sV = mean(sV)), by = .(APH, N)] %>%
      #--- find optimal N by APH ---#
      dplyr::nest_by(APH) %>%
      dplyr::mutate(opt_N_data = list(
        find_opt_N_by_gam(data = data, dep_var = "sV", N_seq_for_search_ls = N_data$N_seq_for_search)
      )) %>%
      dplyr::select(APH, opt_N_data) %>%
      tidyr::unnest(opt_N_data) %>%
      data.table()

    #--- update value function (continuous) mapping ---#
    VF[[s + 1]] <- semi_B(opt_N_on_APH[, max_eu], opt_N_on_APH[, (APH - APH_data$APH_min) / (APH_data$APH_max - APH_data$APH_min)], Nk, decrease = FALSE)

    #--- APH threshold ---#
    if (any(opt_N_on_APH[, N_star == 0])) {
      APH_threshold[[s + 1]] <- max(opt_N_on_APH[N_star != 0, APH])
    } else {
      APH_threshold[[s + 1]] <- 999
    }

    #--- find optimal N function ---#
    N_star_gam[[s + 1]] <- gam(N_star ~ s(APH, k = 5), data = opt_N_on_APH[APH < APH_threshold[[s + 1]], ])
  }

  gc(APH_mat, production_data)

  #--- return the results in the storage ---#
  return(list(VF, N_star_gam, APH_threshold))
}


get_APH_Bernstein_mat <- function(data, APH_data) {
  data_for_mat <-
    data %>%
    #--- bound APH_next between APH_min and APH_max ---#
    .[APH_next <= APH_data$APH_min, APH_next := APH_data$APH_min] %>%
    .[APH_next >= APH_data$APH_max, APH_next := APH_data$APH_max] %>%
    #--- normalize for Bernstein basis evaluation ---#
    .[, APH_norm := (APH_next - APH_data$APH_min) / (APH_data$APH_max - APH_data$APH_min)]

  #--- evaluate APH_norm using Bernstein bases ---#
  APH_mat <- BernBasis_fast(data_for_mat[, APH_norm], Nk)

  return(APH_mat)
}


# !===========================================================
# ! Simulate N-APH path
# !===========================================================
simulate_path <- function(BI_results, prod_data) {
  #--- pre-allocate values ---#
  sim_data_N <- matrix(0, B_sim, T)
  sim_data_APH <- matrix(0, B_sim, T)
  sim_data_APH[, 1] <- APH_0

  APH_cutoff <- BI_results[[3]]
  opt_N_path <- BI_results[[2]]

  for (t in 1:T) {
    print(paste(t, "/", T, sep = ""))

    cutoff_index <- APH_cutoff[[T + 1 - t]] < sim_data_APH[, t]
    sim_data_N[cutoff_index, t] <- 0

    #--- find optimal N ---#
    sim_data_N[!cutoff_index, t] <-
      predict.gam(
        opt_N_path[[T + 1 - t]],
        newdata = data.table(APH = sim_data_APH[!cutoff_index, t])
      )

    #--- simulate yield and APH---#
    yield <-
      lapply(
        sim_data_N[, t],
        \(x) {
          rbetagen(
            1, prod_data$p_pars[1] + prod_data$p_pars[2] * sqrt(x) + prod_data$p_pars[3] * x,
            prod_data$q_pars[1] + prod_data$q_pars[2] * sqrt(x) + prod_data$q_pars[3] * x, prod_data$min_y, prod_data$max_y
          )
        }
      ) %>%
      unlist()

    if (t < T) {
      sim_data_APH[, t + 1] <- sim_data_APH[, t] * (n_window - 1) / n_window + yield / n_window
    }
  }

  #--------------------------
  # summarize N path
  #--------------------------
  colnames(sim_data_N) <- paste("N_", seq(1, T), sep = "")

  N_sim <-
    data.table(sim_data_N) %>%
    melt() %>%
    mutate(variable = as.numeric(gsub("N_", "", variable))) %>%
    data.table() %>%
    setnames(c("variable", "value"), c("t", "N"))

  EN <- N_sim[, list(N = mean(N)), by = t]
  EN_gam <- gam(N ~ s(t, k = 5), data = EN)
  EN[, EN_sm := predict.gam(EN_gam, newdata = data.table(t = seq(1, T)))]

  #--------------------------
  # APH summary
  #--------------------------
  colnames(sim_data_APH) <- paste("APH_", seq(1, T), sep = "")

  APH_sim <-
    data.table(sim_data_APH) %>%
    melt() %>%
    mutate(variable = as.numeric(gsub("APH_", "", variable))) %>%
    data.table() %>%
    setnames(c("variable", "value"), c("t", "APH"))

  E_APH <- APH_sim[, list(APH = mean(APH)), by = t]
  E_APH_gam <- gam(APH ~ s(t, k = 5), data = E_APH)
  E_APH[, APH_sm := predict.gam(E_APH_gam, newdata = data.table(t = seq(1, T)))]

  #--------------------------
  # combine N and APH path
  #--------------------------
  path <- cbind(EN, E_APH[, .(APH, APH_sm)])

  #--------------------------
  # return the results
  #--------------------------
  return(path)
}
