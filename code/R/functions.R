######################################
# Collection of functions
######################################
# written by Taro Mieno on 01/17/2016

# ===================================
# Define function that generate simulation function
# ===================================
sim_gen <- function(min_y, max_y, p_price, p_N, pvf, APH, rate_yield) {
  opt_N_aph <- function(i) {
    #--- N ---#
    N <- N_seq[i]
    pN <- p_pars[1] + p_pars[2] * sqrt(N) + p_pars[3] * N
    qN <- q_pars[1] + q_pars[2] * sqrt(N) + q_pars[3] * N
    y_prop <- c(pN, qN, min_y, max_y)
    N_cost <- p_N * N * acres
    cost <- N_cost + other_cost * acres

    #--- coverage level ---#
    data_temp <- pi_insured(
      APH, rate_yield, B, mu, sigma, hp_prop, y_prop, acres, cost, cov_level, p_price,
      pvf, sub_per, c_year, county_TY, c_rate_now, c_rate_prev,
      c_pars, beta_RPHPE, beta_RP, sec5, sec6
    )

    data <- data.table(
      pi_YP = data_temp$pi_YP,
      pi_RPHPE = data_temp$pi_RPHPE,
      pi_RP = data_temp$pi_RP,
      pi_non = data_temp$pi_non,
      yield = data_temp$yield,
      # indemnity_YP=data_temp$indemnity_YP,
      # indemnity_RPHPE=data_temp$indemnity_RPHPE,
      # indemnity_RP=data_temp$indemnity_RP,
      # revenue_YP=data_temp$revenue_YP,
      # revenue_RPHPE=data_temp$revenue_RPHPE,
      # revenue_RP=data_temp$revenue_RP,
      # raw_revenue=data_temp$raw_revenue,
      h_price = data_temp$h_price,
      premium_YP = data_temp$premium_YP,
      premium_RPHPE = data_temp$premium_RPHPE,
      premium_RP = data_temp$premium_RP,
      APH = data_temp$APH,
      rate_yield = data_temp$rate_yield,
      cost = cost
    )

    data[, `:=`(N = N, cov_level = cov_level)]

    return(data)
  }
}

# ===================================
# Generate Yield
# ===================================
APH_update <- function(APH, N, B) {
  #--- generate shape parameters ---#
  pN <- 3.14 - 0.0921 * sqrt(N) + 0.00603 * N
  qN <- 12.30 - 1.353 * sqrt(N) + 0.0456 * N

  #--- yield ---#
  half_U <- runif(B / 2)
  yield <- qbetagen(c(half_U, -half_U), pN, qN, min_y, max_y)
  # yield <- rbetagen(B,p_N,q_N,min_y,max_y)

  #--- update APH ---#
  APH_next <- APH * 9 / 10 + yield / 10

  return(data.table(APH = APH_next, N = N))
}


#--------------------------
# Bernstein Basis
#--------------------------
BernBasis <- function(x, N) {
  # ptm <- proc.time()
  B <- do.call(cbind, lapply(0:N, function(t) factorial(N) / factorial(N - t) / factorial(t) * x^t * (1 - x)^(N - t)))
  # proc.time() - ptm

  return(B)
}

#--------------------------
# Bernstein Basis
#--------------------------
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

#--------------------------
# semi-parametric regression using Bernstein Basis
#--------------------------

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

# ===================================
# Utility functions
# ===================================
#--------------------------
# CRRA utility function
#--------------------------
crra_u <- function(rho, x) {
  u <- (x^(1 - rho) - 1) / (1 - rho)
  return(u)
}

#--------------------------
# CARA utility function
#--------------------------
cara_u <- function(alpha, x) {
  u <- 1 - exp(-alpha * x / 100)
  return(u)
}

utility_gen <- function(u_type) {
  if (u_type == "CA1") { # CAA 1
    U <- function(x) {
      u <- 1 - exp(-0.0046 * x / 100)
      return(u)
    }
  } else if (u_type == "CA2") { # CAA2
    U <- function(x) {
      u <- 1 - exp(-0.0100 * x / 100)
      return(u)
    }
  } else if (u_type == "CR") { # CRA
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

# ===================================
# Generate value function
# ===================================
# sim_for_BI is defined in BI.cpp
gen_yield_profit_APHnext <- function(APH) {
  rate_yield <- APH

  sim_BI <- sim_for_BI(
    yield_U, h_price, APH, N_seq, N_len, p_N, rate_yield, B_bi, max_y,
    min_y, p_pars, q_pars, acres, cov_level, p_price, pvf, sub_per, c_year, county_TY, c_rate_now, c_rate_prev, c_pars, beta_RPHPE, beta_RP, sec5, sec6, n_window
  )

  return(data.table(
    pi_YP = sim_BI$pi_YP,
    pi_RPHPE = sim_BI$pi_RPHPE,
    pi_RP = sim_BI$pi_RP,
    APH_next = sim_BI$APH_next,
    APH = APH,
    N = rep(N_seq, each = B_bi)
  ))
}


find_opt_N_t <- function(APH, VF_c) {
  APH_temp <- APH

  #--- find optimal N by APH ---#
  VF_gam <- gam(sV ~ s(N, k = 10), data = VF_c[APH == APH_temp, ])

  #--- find the profit by N and APH ---#
  pi_series <- predict.gam(VF_gam, newdata = data.table(N = search_N))

  #--- maximum profit ---#
  max_pi <- max(pi_series)

  #--- optimal N ---#
  opt_N <- search_N[which.max(pi_series)]

  return(
    data.table(
      APH = APH_temp,
      max_pi = max_pi,
      opt_N = opt_N
    )
  )
}

# ===================================
# Add breaks to figures
# ===================================
x_axis_break_add <- function(gplot, breaks) {
  x.axis.break <- breaks
  x.axis.label <- vector()
  for (i in 1:length(x.axis.break)) {
    gplot <- gplot + geom_vline(xintercept = x.axis.break[i], size = 0.2, linetype = 3)
    x.axis.label <- c(x.axis.label, as.character(x.axis.break[i]))
  }
  return(gplot)
}

y_axis_break_add <- function(gplot, breaks) {
  y.axis.break <- breaks
  y.axis.label <- vector()
  for (i in 1:length(y.axis.break)) {
    gplot <- gplot + geom_hline(yintercept = y.axis.break[i], size = 0.2, linetype = 3)
    y.axis.label <- c(y.axis.label, as.character(y.axis.break[i]))
  }
  return(gplot)
}
