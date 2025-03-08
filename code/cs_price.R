######################################
# Simulation with APH 
######################################
# written by Taro Mieno on 02/24/2016

#===================================
# 0. preliminary operations
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
source('./Codes/Simulation2/parameters2.R')

#--------------------------
# source functions 
#--------------------------
#--- sourcing cpp simulation function ---#
sourceCpp('./Codes/Simulation2/profit_insured_2.cpp')
sourceCpp('./Codes/Simulation2/BI.cpp')

#--- 1st stage ---#
source('./Codes/Simulation2/static_simulate_par.R')

#--- 2nd stage ---#
source('./Codes/Simulation2/backward_induction_par.R')

#--- 3rd stage ---#
source('./Codes/Simulation2/simulate_path_par.R')

#--------------------------
# simulation-related preparation
#--------------------------
#--- seed ---#
set.seed(48723)

#--- create lists of insurance and utility types ---#
insurance_type <- c('YP','RPHPE','RP','non')
ins_len <- length(insurance_type)
utility_type <- c('RN','CA1','CA2','CR')
uti_len <- length(utility_type)
full_ui_pairs_ls <- expand.grid(utility_type,insurance_type)
ui_pairs_list <- expand.grid(utility_type,insurance_type[-4]) 
ui_pairs_len <- dim(ui_pairs_list)[1]

#--- parameters for iterations ---#
coverage_ls  <- c(0.70,0.80,0.85,0.90)
cov_len <- length(coverage_ls)
subsidy_ls <- c(0.59,0.48,0.38,0.28)

#--- price lists ---#
# p_corn_ls <- c(3.90,6.01,2.32,6.01,2.32,5.65)
# p_N_ls <- c(0.33,0.51,0.15,0.45,0.25,0.51)
p_corn_ls <- c(2.32,2.59)
p_N_ls <- c(0.15,0.31)
p_len <- length(p_corn_ls)

#===================================
# 1st stage
#===================================
#--- number of simulation ---#
B <- 50000 # (1st stage)
# B <- 50 # (1st stage)

for (i in 1:p_len){
	#--- corn price ---#
	p_price <- p_corn_ls[i]
	ln_mean_p <- log(p_price)-ln_var_p/2
	hp_prop <- c(ln_mean_p,sqrt(ln_var_p))

	#--- N price ---#
	p_N <- p_N_ls[i]

	for (j in 1:cov_len){
		#--- coverage level ---#
		cov_level <- coverage_ls[j]

		#--- subsidy percent ---#
		sub_per <- subsidy_ls[j] 

		#--- simulate ---#
		VF_T_all <- mclapply(1:APH_len,static_sim,mc.cores=21) %>% 
			rbindlist() %>% 
			setkey(uti_type,ins_type,APH)  

		#--- save optimal N (conditional on APH) ---#
		saveRDS(VF_T_all,paste('./Results2/VF_T_all_',min_y,'_',max_y,'_cov_',cov_level*100,'_rho_',abs(rho)*10,
			'_sub_',sub_per*100,'_pCorn_',p_price*100,'_pN_',p_N*100,'.rds',sep=''))
	}
	
}

#===================================
# 2nd stage: Backward Induction
#===================================
#--- number of simulation ---#
B_bi <- 5000 # (2nd stage)

#--- discount rate (2nd stage) ---#
disc <- 0.04 # discount rate

#--------------------------
# Generate concurrent profit and APH in the next period
#--------------------------
# This part finds
# 1. the profits in the current year for each combination of N and current APH
# 2. APH in the next period for each combination of N and current APH
# The results are used for all the periods during the backward induction procedure

#--- generate price and yield in U() ---#
half_U <- rmvnorm(mean=c(0,0),sig=sigma,n=B_bi) %>% pnorm() 
U <- rbind(half_U,1-half_U)  
yield_U <- U[,2] 

#--------------------------
# Simulate
#--------------------------

for (i in 1:p_len){
	#--- corn price ---#
	p_price <- p_corn_ls[i]
	ln_mean_p <- log(p_price)-ln_var_p/2
	hp_prop <- c(ln_mean_p,sqrt(ln_var_p))

	#--- simulated harvest price ---#
	h_price  <- qlnorm(U[,1],meanlog=ln_mean_p,sdlog=sqrt(ln_var_p))

	#--- N price ---#
	p_N <- p_N_ls[i]

	for (j in 1:cov_len){
		#--- coverage level ---#
		cov_level <- coverage_ls[j]

		#--- subsidy percent ---#
		sub_per <- subsidy_ls[j]

		#--- import the 1st stage results ---#
		VF_T <- readRDS(paste('./Results2/VF_T_all_',min_y,'_',max_y,'_cov_',cov_level*100,'_rho_',abs(rho)*10,
			'_sub_',sub_per*100,'_pCorn_',p_price*100,'_pN_',p_N*100,'.rds',sep=''))
		VF_T[,ins_type:=factor(ins_type,levels=c('YP','RPHPE','RP'))]

		#--- generate yield, APH, and profit ---#
		V_c <- mclapply(1:APH_len,VF_sim,mc.cores=31) %>% 
			rbindlist() 

		V_c[APH_next<=APH_min,APH_next:=APH_min]
		V_c[APH_next>=APH_max,APH_next:=APH_max]
		V_c[,APH_norm:=(APH_next-APH_min)/(APH_max-APH_min)]
		
		#--- APH evaluated at Bernstein bases ---#
		APH_mat <- BernBasis_fast(V_c[,APH_norm],Nk)
	
		#--- simulate ---#
		BI_results <- mclapply(1:ui_pairs_len, function(x) BI(x), mc.cores=6)

		#--- arrange the results ---#
		VF_path <- sapply(BI_results,'[',1) 
		opt_N_path <- sapply(BI_results,'[',2) 
		APH_cutoff <- sapply(BI_results,'[',3) 

		#--- save results ---#
		saveRDS(VF_path,paste('./Results2/VF_path_',min_y,'_',max_y,'_cov_',cov_level*100,'_rho_',abs(rho)*10,
			'_sub_',sub_per*100,'_pCorn_',p_price*100,'_pN_',p_N*100,'.rds',sep=''))
		saveRDS(opt_N_path,paste('./Results2/opt_N_path_',min_y,'_',max_y,'_cov_',cov_level*100,'_rho_',abs(rho)*10,
			'_sub_',sub_per*100,'_pCorn_',p_price*100,'_pN_',p_N*100,'.rds',sep=''))
		saveRDS(APH_cutoff,paste('./Results2/APH_cutoff_',min_y,'_',max_y,'_cov_',cov_level*100,'_rho_',abs(rho)*10,
			'_sub_',sub_per*100,'_pCorn_',p_price*100,'_pN_',p_N*100,'.rds',sep=''))
	}
}

#===================================
# 3rd stage: Simulate Expected Path of N and APH
#===================================
#--- specify the starting point ---#
APH_0 <- 136.5 # starting APH
	
#--- number of iterations for simulation ---#
B_sim <- 30000

for (i in 1:p_len){
	#--- corn price ---#
	# no need to redefine hp_prop as prices are not
	# simulated
	p_price <- p_corn_ls[i]

	#--- N price ---#
	p_N <- p_N_ls[i]

	for (j in 1:cov_len){
		#--- coverage level ---#
		cov_level <- coverage_ls[j]

		#--- subsidy percent ---#
		sub_per <- subsidy_ls[j]

		#--- read BI results ---#
		opt_N_path <- readRDS(paste('./Results2/opt_N_path_',min_y,'_',max_y,'_cov_',cov_level*100,'_rho_',abs(rho)*10,
			'_sub_',sub_per*100,'_pCorn_',p_price*100,'_pN_',p_N*100,'.rds',sep=''))
		APH_cutoff <- readRDS(paste('./Results2/APH_cutoff_',min_y,'_',max_y,'_cov_',cov_level*100,'_rho_',abs(rho)*10,
			'_sub_',sub_per*100,'_pCorn_',p_price*100,'_pN_',p_N*100,'.rds',sep=''))		#--- run simulation ---#

		# path <- mclapply(1:ui_pairs_len,sim_path,mc.cores=12) %>% rbindlist()
		path <- mclapply(1:ui_pairs_len,sim_path,mc.cores=12) %>% rbindlist()

		#--- save the results ---#
		saveRDS(path,paste('./Results2/path_',min_y,'_',max_y,'_cov_',cov_level*100,'_rho_',abs(rho)*10,
			'_sub_',sub_per*100,'_pCorn_',p_price*100,'_pN_',p_N*100,'.rds',sep=''))
	}
}


