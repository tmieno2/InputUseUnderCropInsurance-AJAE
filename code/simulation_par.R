######################################
# Simulation with APH 
######################################
# written by Taro Mieno on 02/24/2016

#===================================
# 0. Preliminary Operations
#===================================
setwd('~/Dropbox/MyResearch/CropInsuranceProgram/MoralHazard/Simulation')
source('~/Dropbox/R_libraries_load/library.R')
library('mc2d')
library('quadprog')

#--- crop insurance parameters and functions ---#
source('./Codes/Base/premium_rate_adj.R')
source('./Codes/Base/functions.R')

#--- county specific crop insurance parameters ---#
source('./Codes/Base/Sioux_IA_parameters.R')

#--- other parameters ---#
source('./Codes/BackwardInduction/parameters.R')

#--- sourcing cpp simulation function ---#
sourceCpp('./Codes/BackwardInduction/profit_insured.cpp')
sourceCpp('./Codes/BackwardInduction/BI.cpp')

#===================================
# Optimal N rate simulation for the static problems 
#===================================
#--- seed ---#
set.seed(48723)

#--------------------------
# specify coverage and subsidy level
#--------------------------
#--- coverage level ---#
cov_level <- 0.85

#--- subsidy level ---#
# sub_per <- 0.99
sub_per <- subsidy_percent[coverage_level==cov_level,subsidy_percent]

#--------------------------
# run simulations
#--------------------------
#--- # of iterations ---#
B <- 50000
# B <- 50

#--- run simulations ---#
source('./Codes/BackwardInduction/static_simulate_par.R')
VF_T_all <- mclapply(1:APH_len,static_sim,mc.cores=26) %>% rbindlist() %>% setkey(uti_type,ins_type,APH)

#--- save optimal N (conditional on APH) ---#
saveRDS(VF_T_all,paste('./Results/VF_T_all_',min_y,'_',max_y,'_cov_',cov_level*100,'_rho_',abs(rho)*10,'_sub_',sub_per*100,'.rds',sep=''))

#===================================
# Step 2: Backward Induction
#===================================
#--- set parameters for simulations ---#
B_bi <- 5000 # number of iterations  
# B_bi <- 50 # number of iterations  
disc <- 0.04 # discount rate

#--- run BI ---#
source('./Codes/BackwardInduction/backward_induction_par.R')
BI_results <- mclapply(1:pairs_len,BI,mc.cores=12)

#--- arrange the results ---#
VF_path <- sapply(BI_results,'[',1) 
opt_N_path <- sapply(BI_results,'[',2) 

#--- save results ---#
saveRDS(VF_path,paste('./Results/VF_path_',min_y,'_',max_y,'_cov_',cov_level*100,'_rho_',abs(rho)*10,'_sub_',sub_per*100,'.rds',sep=''))
saveRDS(opt_N_path,paste('./Results/opt_N_path_',min_y,'_',max_y,'_cov_',cov_level*100,'_rho_',abs(rho)*10,'_sub_',sub_per*100,'.rds',sep=''))

#===================================
# Simulate Expected Path of N and APH
#===================================
#--- specify the starting point ---#
APH_0 <- 136.5 # starting APH
B_sim <- 50000 # number of iterations for simulation
# B_sim <- 50 # number of iterations for simulation

#--- run simulation ---#
source('./Codes/BackwardInduction/simulate_path_par.R')
path <- mclapply(1:pairs_len,sim_path,mc.cores=12) %>% rbindlist()

#--- save the results ---#
saveRDS(path,paste('./Results/path_',min_y,'_',max_y,'_cov_',cov_level*100,'_rho_',abs(rho)*10,'_sub_',sub_per*100,'.rds',sep=''))


