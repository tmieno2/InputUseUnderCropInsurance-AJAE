######################################
# Backward Induction
######################################

BI <- function(k) { # loop over insurance-utility cases

  #++++++++++++++++++++++++++++++++++++
  #+ Setup
  #++++++++++++++++++++++++++++++++++++
  #---------------------
  #- Define utility and insurance types to work on
  #---------------------
  U_type <- ui_pairs_ls[k, utility_type]
  I_type <- ui_pairs_ls[k, ins_type]

  #---------------------
  #- Create place holders
  #---------------------
  #--- coefficients of semi-parametric regressions that allows for mapping APH to V()  ---#
  VF <- vector(mode = "list", length = T)
  #--- gam regression results that map APH to optimal N at each t (optimal decision) ---#
  N_star_gam <- vector(mode = "list", length = T)
  #--- highest APH level to consider at each t ---#
  APH_threshold <- vector(mode = "list", length = T)

  #++++++++++++++++++++++++++++++++++++
  #+ t = T (final period)
  #++++++++++++++++++++++++++++++++++++
  #--- Extract the relevant data from the value function at t = T ---#
  VT <- VF_T[uti_type == U_type & ins_type == I_type, ]

  #--- Define utility function ---#
  Util <- utility_gen(U_type)

  #--- Calculate utility level ---#
  c_U <- Util(V_c[, paste("pi_", I_type, sep = ""), with = FALSE])

  #--------------------------
  #- VF_t as a function of APH
  #--------------------------
  # Run semi-parametric regression using Bernstein polynomials. The resulting coefficient estimates on the bases are stored.
  VF[[1]] <- semi_B(VT[, max_pi], VT[, (APH - APH_min) / (APH_max - APH_min)], Nk, decrease = FALSE)

  #--------------------------
  #- Optimal N conditional on APH
  #--------------------------
  #--- identify the highest value of APH at which optimal N is NOT 0 ---#
  if (any(VT[, N_star == 0])) {
    APH_threshold[[1]] <- max(VT[N_star != 0, APH])
  } else {
    APH_threshold[[1]] <- 999
  }
  #--- semi-parameteric regresion ---#
  # the resulting gam object can map (predict) APH level to optimal N
  N_star_gam[[1]] <- gam(N_star ~ s(APH, k = 5), data = VT[APH < APH_threshold[[1]], ])

  #++++++++++++++++++++++++++++++++++++
  #+ t < T
  #++++++++++++++++++++++++++++++++++++
  for (s in 1:(T - 1)) { # k=1 is T-1

    print(paste0(s, "/", T - 1, " of ", k, "/", nrow(ui_pairs_ls), " cases"))

    #--- find value function ---#
    summed_V <- c_U + APH_mat %*% VF[[s]] / (1 + disc)
    # summed_V <- Util(V_c[,paste('pi_',I_type,sep=''),with=FALSE]) + predict.gam(VF[[s]],newdata=V_c[,.(APH=APH_next)])/(1+disc)
    V_c[, sV := summed_V]

    #--- find mean VF by N and APH ---#
    VF_c <- V_c[, list(sV = mean(sV)), by = .(APH, N)]
    VF_data <-
      lapply(
        1:APH_len,
        \(x) find_opt_N_t(APH_list[x], VF_c)
      ) %>%
      rbindlist()
    # VF_data <- lapply(1:APH_len, function(x) VF_gen_par(x,VF_c)) %>% rbindlist()

    #--- update value function (continuous) ---#
    VF[[s + 1]] <- semi_B(VF_data[, max_pi], VF_data[, (APH - APH_min) / (APH_max - APH_min)], Nk, decrease = FALSE)

    #--- APH threshold ---#
    if (any(VF_data[, opt_N == 0])) {
      APH_threshold[[s + 1]] <- max(VF_data[opt_N != 0, APH])
    } else {
      APH_threshold[[s + 1]] <- 999
    }

    #--- find optimal N function ---#
    N_star_gam[[s + 1]] <- gam(opt_N ~ s(APH, k = 5), data = VF_data[APH < APH_threshold[[s + 1]], ])
  }

  #--- return the results in the storage ---#
  return(list(VF, N_star_gam, APH_threshold))
}
