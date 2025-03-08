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
actual_pairs <-	paste(c(0.70,0.80,0.85,0.90),'_',c(0.59,0.48,0.38,0.28),sep='')
coverage_list  <- c(0.70,0.80,0.85,0.90)
subisty_list <- c(1,0.59,0.48,0.38,0)
paired_list <- expand.grid(coverage_list,subisty_list)
index <- paste(paired_list[,1],'_',paired_list[,2],sep='') %in% actual_pairs
paired_list <- paired_list[!index,]
pairs_len <- nrow(paired_list)

#===================================
# 1st stage
#===================================
# 1 st stage is not necessary as the solutions to the static problems
# depends on coverage level, but not on subsidy level. Therefore,
# you can use VF_T_all that were obtained by running main.R

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
h_price  <- qlnorm(U[,1],meanlog=ln_mean_p,sdlog=sqrt(ln_var_p))
yield_U <- U[,2] 

#--------------------------
# Simulate
#--------------------------
for (i in 1:pairs_len){
	#--- coverage level ---#
	cov_level <- paired_list[i,1]

	#--- subsidy percent to import the appropriate VF_T_all ---#
	sub_file <- subsidy_percent[coverage_level==cov_level,subsidy_percent]

	#--- subsidy level used to calculate premium ---#
	sub_per <- paired_list[i,2]

	#--- import the 1st stage results ---#
	VF_T <- readRDS(paste('./Results2/VF_T_all_',min_y,'_',max_y,'_cov_',cov_level*100,'_rho_',abs(rho)*10,
		'_sub_',sub_file*100,'_pCorn_',p_price*100,'_pN_',p_N*100,'.rds',sep=''))[ins_type!='non',]
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
	BI_results <- mclapply(1:ui_pairs_len,BI,mc.cores=6)

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

#===================================
# 3rd stage: Simulate Expected Path of N and APH
#===================================
#--- specify the starting point ---#
APH_0 <- 136.5 # starting APH
	
#--- number of iterations for simulation ---#
B_sim <- 50000

for (i in 1:pairs_len){
	#--- coverage level ---#
	cov_level <- paired_list[i,1]
	
	#--- subsidy level ---#
	sub_per <- paired_list[i,2]

		#--- read BI results ---#
	opt_N_path <- readRDS(paste('./Results2/opt_N_path_',min_y,'_',max_y,'_cov_',cov_level*100,'_rho_',abs(rho)*10,
		'_sub_',sub_per*100,'_pCorn_',p_price*100,'_pN_',p_N*100,'.rds',sep=''))
	APH_cutoff <- readRDS(paste('./Results2/APH_cutoff_',min_y,'_',max_y,'_cov_',cov_level*100,'_rho_',abs(rho)*10,
		'_sub_',sub_per*100,'_pCorn_',p_price*100,'_pN_',p_N*100,'.rds',sep=''))

	#--- run simulation ---#
	# path <- mclapply(1:ui_pairs_len,sim_path,mc.cores=12) %>% rbindlist()
	path <- mclapply(1:ui_pairs_len,sim_path,mc.cores=12) %>% rbindlist()

	#--- save the results ---#
	saveRDS(path,paste('./Results2/path_',min_y,'_',max_y,'_cov_',cov_level*100,'_rho_',abs(rho)*10,
		'_sub_',sub_per*100,'_pCorn_',p_price*100,'_pN_',p_N*100,'.rds',sep=''))
}


