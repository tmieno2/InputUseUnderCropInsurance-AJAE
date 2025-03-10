######################################
# Simulation with APH
######################################
# written by Taro Mieno on 02/24/2016

# ===================================
# 0. preliminary operations
# ===================================
setwd("~/Dropbox/MyResearch/CropInsuranceProgram/MoralHazard/Simulation")
source("~/Dropbox/R_libraries_load/library.R")
library("mc2d")
library("quadprog")

#--- crop insurance parameters and functions ---#
source("./Codes/Base/premium_rate_adj.R")
source("./Codes/Base/functions.R")

#--- county specific crop insurance parameters ---#
source("./Codes/Base/Sioux_IA_parameters.R")

#--- other parameters ---#
source("./Codes/Simulation/parameters.R")

#--------------------------
# source functions
#--------------------------
#--- sourcing cpp simulation function ---#
sourceCpp("./Codes/Simulation/profit_insured_2.cpp")
sourceCpp("./Codes/Simulation/BI.cpp")

#--- 1st stage ---#
source("./Codes/Simulation/static_simulate_par.R")

#--- 2nd stage ---#
source("./Codes/Simulation/backward_induction_par.R")

#--- 3rd stage ---#
source("./Codes/Simulation/simulate_path_par.R")

#--------------------------
# simulation-related preparation
#--------------------------
#--- seed ---#
set.seed(783473)

#--- coverage level ---#
cov_level <- 0.85

#--- subsidy percent ---#
sub_per <- 0.38

#--- read BI results ---#
opt_N_path <- readRDS(paste("./Results/opt_N_path_", min_y, "_", max_y, "_cov_", cov_level * 100, "_rho_", abs(rho) * 10, "_sub_", sub_per * 100, "_pCorn_", p_price * 100, "_pN_", p_N * 100, ".rds", sep = ""))

#--- insurance and utility type ---#
U_type <- "CR"
I_type <- "RP"

#--- specify the starting point ---#
APH_0 <- 136.5 # starting APH

#--- number of iterations for simulation ---#
B_sim <- 30000

# ===================================
# Simulate
# ===================================
#--- pre-allocate values ---#
sim_data_N <- matrix(0, B_sim, T)
sim_data_APH <- matrix(0, B_sim, T)
sim_data_APH[, 1] <- APH_0

for (t in 1:T) {
  print(paste(t, "/", T - 1, sep = ""))
  #--- find optimal N ---#
  # 12 is for CR and RP
  sim_data_N[, t] <- predict.gam(opt_N_path[[12]][[T + 1 - t]], newdata = data.table(APH = sim_data_APH[, t]))

  #--- simulate yield and APH---#
  yield <- mclapply(
    sim_data_N[, t],
    function(x) {
      rbetagen(
        1, p_pars[1] + p_pars[2] * sqrt(x) + p_pars[3] * x,
        q_pars[1] + q_pars[2] * sqrt(x) + q_pars[3] * x, min_y, max_y
      )
    },
    mc.cores = 2
  ) %>%
    unlist()

  if (t < T) {
    sim_data_APH[, t + 1] <- sim_data_APH[, t] * (n_window - 1) / n_window + yield / n_window
  }
}

colnames(sim_data_N) <- paste("N_", seq(1, T), sep = "")
N_sim <- data.table(sim_data_N) %>%
  mutate(id = 1:B_sim) %>%
  melt(id.var = "id") %>%
  mutate(variable = as.numeric(gsub("N_", "", variable))) %>%
  data.table() %>%
  setnames(c("variable", "value"), c("t", "N"))


colnames(sim_data_APH) <- paste("APH_", seq(1, T), sep = "")
APH_sim <- data.table(sim_data_APH) %>%
  mutate(id = 1:B_sim) %>%
  melt(id.var = "id") %>%
  mutate(variable = as.numeric(gsub("APH_", "", variable))) %>%
  data.table() %>%
  setnames(c("variable", "value"), c("t", "APH"))


# ===================================
# Example N and APH path
# ===================================
# identify ids that show very high and low N rates and mean N rates
mean_N_by_id <- N_sim[, mean(N), by = id]

high_id <- which.max(mean_N_by_id[, V1])
mean_id <- 1
low_id <- which.min(mean_N_by_id[, V1])

#--------------------------
# N path
#--------------------------
plot_N <- N_sim[id %in% c(high_id, mean_id, low_id), ]
plot_N[id == 1, type := "Average N"]
plot_N[id == high_id, type := "High N"]
plot_N[id == low_id, type := "Low N"]
plot_N[, value := N]
plot_N[, N := NULL]
plot_N[, value_type := "N (lb/acre)"]

# ggplot(data=plot_N) +
# 	geom_line(aes(x=t,y=N,color=factor(type))) +
# 	scale_color_discrete(name='') +
# 	theme(
# 		legend.position='bottom'
# 		)
# ggsave('./Graphs/N_path_example')


#--------------------------
# APH path
#--------------------------
plot_APH <- APH_sim[id %in% c(high_id, mean_id, low_id), ]
plot_APH[id == 1, type := "Average N"]
plot_APH[id == high_id, type := "High N"]
plot_APH[id == low_id, type := "Low N"]
plot_APH[, value := APH]
plot_APH[, APH := NULL]
plot_APH[, value_type := "APH (bu/acre)"]

# ggplot(data=plot_APH) +
# 	geom_line(aes(x=t,y=APH,color=factor(type))) +
# 	scale_color_discrete(name='') +
# 	theme(
# 		legend.position='bottom'
# 		)
# ggsave('./Graphs/N_path_example')

#--------------------------
# combine
#--------------------------
plot_data <- rbind(plot_N, plot_APH)

ggplot(data = plot_data) +
  geom_line(aes(x = t, y = value, color = factor(type))) +
  facet_grid(value_type ~ ., scales = "free") +
  scale_color_discrete(name = "") +
  ylab("") +
  theme(
    legend.position = "bottom"
  )
ggsave("./Graphs/N_path_example.pdf", width = 6, height = 3.5)
ggsave("./Graphs/N_path_example.png", width = 6, height = 3.5)
