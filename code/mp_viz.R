######################################
# Simulation with APH
######################################
# written by Taro Mieno on 02/24/2016

# ===================================
# 0. preliminary operations
# ===================================
setwd("~/Dropbox/MyResearch/CropInsuranceProgram/MoralHazard/Simulation")
source("~/Dropbox/R_libraries_load/library.R")
source("~/Dropbox/Teaching/R_functions/functions.R")
library("mc2d")
library("quadprog")

#--- crop insurance parameters and functions ---#
source("./Codes/Base/premium_rate_adj.R")
source("./Codes/Base/functions.R")


#--- county specific crop insurance parameters ---#
source("./Codes/Base/Sioux_IA_parameters.R")

#--- other parameters ---#
source("./Codes/Simulation2/parameters2.R")

#--- sourcing cpp simulation function ---#
sourceCpp("./Codes/Simulation2/BI.cpp")

#--------------------------
# simulation-related preparation
#--------------------------
#--- seed ---#
set.seed(48723)

#--- create lists of insurance and utility types ---#
insurance_type <- c("YP", "RPHPE", "RP", "non")
ins_len <- length(insurance_type)
utility_type <- c("RN", "CA1", "CA2", "CR")
uti_len <- length(utility_type)
full_ui_pairs_ls <- expand.grid(utility_type, insurance_type)
ui_pairs_list <- expand.grid(utility_type, insurance_type[-4])
ui_pairs_len <- dim(ui_pairs_list)[1]

#--- parameters for iterations ---#
coverage_list <- c(0.70, 0.80, 0.85, 0.90)
subisty_list <- c(1, 0.59, 0.48, 0.38, 0)
paired_list <- expand.grid(coverage_list, subisty_list)
pairs_len <- nrow(paired_list)


#--------------------------
# marginal value breakdown
#--------------------------
B_bi <- 50000
disc <- 0.04 # discount rate

half_U <- rmvnorm(mean = c(0, 0), sig = sigma, n = B_bi) %>% pnorm()
U <- rbind(half_U, 1 - half_U)
h_price <- qlnorm(U[, 1], meanlog = ln_mean_p, sdlog = sqrt(ln_var_p))
yield_U <- U[, 2]

#--- generate APH_next mat ---#
cov_level <- 0.8 # temporary
sub_per <- 0.38 # temporary
V_c <- VF_sim(37) # temporary
V_c[APH_next <= APH_min, APH_next := APH_min]
V_c[APH_next >= APH_max, APH_next := APH_max]
V_c[, APH_norm := (APH_next - APH_min) / (APH_max - APH_min)]
APH_mat <- BernBasis_fast(V_c[, APH_norm], Nk)

MP_gen <- function(i) {
  print(i)

  #--- coverage level ---#
  cov_level <- paired_list[i, 1]

  #--- subsidy level ---#
  sub_per <- paired_list[i, 2]

  #--- VF_path ---#
  VF_path <- readRDS(paste("./Results2/VF_path_", min_y, "_", max_y, "_cov_", cov_level * 100, "_rho_", abs(rho) * 10,
    "_sub_", sub_per * 100, "_pCorn_", p_price * 100, "_pN_", p_N * 100, ".rds",
    sep = ""
  ))

  V_c <- VF_sim(37)

  temp_MP_store <- list()
  for (k in 1:ui_pairs_len) {
    U_type <- ui_pairs_list[k, 1]
    I_type <- ui_pairs_list[k, 2]

    #--- current profit ---#
    temp_pi_c <- eval(parse(text = paste("temp_data <- V_c[,pi_", I_type, "]", sep = "")))

    #--- VF ---#
    V_next <- APH_mat %*% VF_path[[k]][[40]]

    #--- merge ---#
    temp_data <- V_c[, .(APH_next, N)]
    temp_data[, pi_c := temp_pi_c]
    temp_data[, VF := V_next]

    #--- find mean by N ---#
    temp_MP <- temp_data[, .(pi_c = mean(pi_c), VF = mean(VF)), by = N]
    pi_c_gam <- gam(pi_c ~ s(N, k = 5), data = temp_MP)
    VF_gam <- gam(VF ~ s(N, k = 5), data = temp_MP)
    pi_c <- predict.gam(pi_c_gam, newdata = data.table(N = N_seq_for_search))
    VF <- predict.gam(VF_gam, newdata = data.table(N = N_seq_for_search))

    MP_data <- data.table(
      pi_c = pi_c,
      VF = VF,
      N = N_seq_for_search
    )

    MP_data[, rev_c := pi_c + N * p_N * acres] # revenue
    MP_data[, rev_next := c(MP_data[, rev_c][-1], NA)]
    MP_data[, rev_dif := rev_next - rev_c]
    MP_data[, VF_next := c(MP_data[, VF][-1], NA)]
    MP_data[, VF_dif := VF_next - VF]
    MP_data[, tot_dif := rev_dif + VF_dif / (1 + disc)]
    MP_data[, tot_prof := pi_c + VF / (1 + disc)]
    MP_data[, insurance := I_type]
    MP_data[, utility := U_type]

    # N_seq_for_search[which.max(MP_data[,tot_prof])]

    # ggplot(data=MP_data) +
    # 	geom_line(aes(x=N,y=tot_dif/acres/search_increment),size=0.3) +
    # 	h_line(xmin=150,xmax=220,y=0.15,size=0.2) +
    # 	xlim(150,220) +
    # 	theme(
    # 		legend.position='bottom'
    # 	)

    #--- save ---#
    temp_MP_store[[k]] <- MP_data
  }
  MP_store <- rbindlist(temp_MP_store)
  MP_store[, coverage := cov_level]
  MP_store[, subsidy := sub_per]
  return(MP_store)
}

MP_data <- mclapply(1:pairs_len, MP_gen, mc.cores = 5) %>% rbindlist()


#--------------------------
# find marginal revenue and marginal VF
#--------------------------
N_seq_for_search[which.max(MP_data[coverage == 0.85 & insurance == "RP" & utility == "RN" & subsidy == 0.38, tot_prof])]
N_seq_for_search[which.max(MP_data[coverage == 0.70 & insurance == "RP" & utility == "RN" & subsidy == 0.59, tot_prof])]


#--------------------------
# How subsidy affects MB
#--------------------------
g_sub <- ggplot(data = MP_data[coverage == 0.85 & insurance == "RP" & utility == "RN" & subsidy %in% c(0, 1), ]) +
  geom_line(aes(x = N, y = tot_dif / acres / search_increment, color = factor(subsidy)), size = 0.3) +
  h_line(xmin = 150, xmax = 220, y = 0.15, size = 0.2) +
  scale_y_continuous(breaks = seq(0.12, 0.18, by = 0.01), limits = c(0.12, 0.18)) +
  scale_x_continuous(breaks = seq(190, 200, by = 1), limits = c(190, 200)) +
  scale_color_discrete(name = "Subsidy Rate") +
  ylab("$/acre") +
  xlab("N (lb/acre)") +
  theme(
    legend.position = "bottom"
  )
y_axis_break_add(g_sub, breaks = seq(0.12, 0.18, by = 0.01))
ggsave("./Graphs/MR_MC_RN_85_subsidy.pdf", width = 6, height = 3.5)

#--------------------------
# Static vs Dynamic
#--------------------------
# col <- c('MR_t + VF_t'='red','MR_t'='blue','N price'='black')
lt <- c("MR_t + VF_t" = "dashed", "MR_t" = "dotted", "N price" = "solid")
g_mr <- ggplot(data = MP_data[coverage == 0.85 & insurance == "RP" & utility == "RN" & subsidy == 0.38, ]) +
  geom_line(aes(x = N, y = tot_dif / acres / search_increment, linetype = "MR_t + VF_t"), size = 0.4) +
  geom_line(aes(x = N, y = rev_dif / acres / search_increment, linetype = "MR_t"), size = 0.4) +
  geom_line(data = data.table(x = N_seq_for_search, y = p_N), aes(x = x, y = y, linetype = "N price"), size = 0.4) +
  # scale_color_manual(values=col,name='',labels=expression(MR[t],MR[t]+MVF[t],np)) +
  scale_linetype_manual(values = lt, name = "", labels = expression(MR[t], MR[t] + MVF[t], np)) +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1), limits = c(0, 0.5)) +
  scale_x_continuous(breaks = seq(0, 240, by = 20)) +
  ylab("$/acre") +
  xlab("N (lb/acre)") +
  theme(
    legend.position = "bottom"
  )
y_axis_break_add(g_mr, breaks = seq(0, 0.5, by = 0.1))
ggsave("./Graphs/MR_MC_RN_85_38.pdf", width = 6, height = 3.5)
ggsave("~/Dropbox/MyResearch/CropInsuranceProgram/MoralHazard/APH_on_MH_dynamic/Revision_1/RevisedVersion/MR_MC_RN_85_38.pdf", width = 6, height = 3.5)
