######################################
# Find the optimal N rate for the static problems
######################################

find_opt_N_static <- function(APH) {
  print(paste("working on APH = ", APH, sep = ""))

  #--------------------------
  # set parameters
  #--------------------------
  rate_yield <- APH

  #--------------------------
  # simulation (steady-state with APH taken into account)
  #--------------------------
  sim <- sim_gen(min_y, max_y, p_price, p_N, pvf, APH, rate_yield)
  # sim_results <- mclapply(1:N_len,sim,mc.cores=2) %>% rbindlist()
  sim_results <- lapply(1:N_len, function(x) sim(x)) %>% rbindlist()
  # saveRDS(sim_results,'./Data/sim_with_aph.rds')

  #--------------------------
  # Expected Profit
  #--------------------------
  pi_sim <- sim_results[, list(
    pi_YP_RN = mean(pi_YP),
    pi_RPHPE_RN = mean(pi_RPHPE),
    pi_RP_RN = mean(pi_RP),
    pi_non_RN = mean(pi_non),
    pi_YP_CA1 = mean(cara_u(0.0046, pi_YP)),
    pi_RPHPE_CA1 = mean(cara_u(0.0046, pi_RPHPE)),
    pi_RP_CA1 = mean(cara_u(0.0046, pi_RP)),
    pi_non_CA1 = mean(cara_u(0.0046, pi_non)),
    pi_YP_CA2 = mean(cara_u(0.0100, pi_YP)),
    pi_RPHPE_CA2 = mean(cara_u(0.0100, pi_RPHPE)),
    pi_RP_CA2 = mean(cara_u(0.0100, pi_RP)),
    pi_non_CA2 = mean(cara_u(0.0100, pi_non)),
    pi_YP_CR = mean(crra_u(0.6, pi_YP)),
    pi_RPHPE_CR = mean(crra_u(0.6, pi_RPHPE)),
    pi_RP_CR = mean(crra_u(0.6, pi_RP)),
    pi_non_CR = mean(crra_u(0.6, pi_non))
  ), by = list(N, cov_level)] %>%
    melt(id.vars = c("N", "cov_level")) %>%
    data.table() %>%
    setnames("variable", "case")

  pi_sim[, N_norm := (N - N_min) / (N_max - N_min)]
  pi_sim[, APH := APH]

  smoothing <- function(k) {
    temp_case <- paste0("pi_", full_ui_pairs_ls[k, 1], "_", full_ui_pairs_ls[k, 2])
    temp_data <- pi_sim[case == temp_case, ]

    e_pi <- gam(value ~ s(N, k = 5), data = temp_data)

    #--- optimal N ---#
    fitted_non <- predict.gam(e_pi, newdata = data.table(N = search_N))
    opt_N_non <- min(search_N) + search_increment * (which.max(fitted_non) - 1)
    temp_data <- data.table(
      APH = APH,
      uti_type = full_ui_pairs_ls[k, utility_type],
      ins_type = full_ui_pairs_ls[k, ins_type],
      N_star = opt_N_non,
      max_pi = max(fitted_non)
    )
    return(temp_data)
  }

  VF <- lapply(1:(ins_len * uti_len), smoothing) %>% rbindlist()
  VF[, coverage := cov_level]

  return(VF)
}
