######################################
# The impact of APH adjustment on N use 
######################################
# written by Taro Mieno on 07/20/2016

#===================================
# 0. preliminary operations
#===================================
setwd('~/Dropbox/MyResearch/CropInsuranceProgram/MoralHazard/Simulation')
source('~/Dropbox/R_libraries_load/library.R')
library('mc2d')

#--- crop insurance parameters and functions ---#
source('./Codes/Base/functions.R')

#--- county specific crop insurance parameters ---#
source('./Codes/Base/Sioux_IA_parameters.R')

#--- other parameters ---#
source('./Codes/Simulation2/parameters.R')

#--------------------------
# source functions 
#--------------------------
#--- import gam that calculates premium ---#
premium_YP_gam <- readRDS('./Results/premium_YP_gam.rds') 
premium_RPHPE_gam <- readRDS('./Results/premium_RPHPE_gam.rds') 
premium_RP_gam <- readRDS('./Results/premium_RP_gam.rds') 

#--- source the simulation function ---#
source('./Codes/Simulation2/simulate_path_APH_adj_par.R')

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
ui_pairs_list <- expand.grid(utility_type,insurance_type[-4]) 
ui_pairs_len <- dim(ui_pairs_list)[1]

#--- parameters for iterations ---#
coverage_ls  <- c(0.70,0.80,0.85)
cov_len <- length(coverage_ls)
subsidy_ls <- c(0.59,0.48,0.38)

#--- APH list ---#
APH_0_list <- c(110,120,122.8,130,132.8,140,142.8,150,152.8)
APH_len <- length(APH_0_list)

#===================================
# Simulate 
#===================================
#--- number of iterations for simulation ---#
B_sim <- 30000
# B_sim <- 10

for (i in 1:APH_len){

	#--- starting APH ---#
	APH_0 <- APH_0_list[i]
	
	for (j in 1:cov_len){
		#--- coverage level ---#
		cov_level <- coverage_ls[j]

		#--- subsidy percent ---#
		sub_per <- subsidy_ls[j]

		#--- read BI results ---#
		opt_N_path <- readRDS(paste('./Results2/opt_N_path_',min_y,'_',max_y,'_cov_',cov_level*100,'_rho_',abs(rho)*10,
			'_sub_',sub_per*100,'_pCorn_',p_price*100,'_pN_',p_N*100,'.rds',sep=''))

		APH_cutoff <- readRDS(paste('./Results2/APH_cutoff_',min_y,'_',max_y,'_cov_',cov_level*100,'_rho_',abs(rho)*10,
			'_sub_',sub_per*100,'_pCorn_',p_price*100,'_pN_',p_N*100,'.rds',sep=''))
		
		#--- run simulation ---#
		# path <- mclapply(1:ui_pairs_len,sim_path,mc.cores=12) %>% rbindlist()
		path <- mclapply(1:ui_pairs_len,sim_path,mc.cores=12) %>% rbindlist()
		path[,p_corn:=p_price]
		path[,p_N:=p_N]
		path[,APH_0:=APH_0]

		#--- save the results ---#
		saveRDS(path,paste('./Results2/path_APH_adj_',min_y,'_',max_y,'_cov_',cov_level*100,'_rho_',abs(rho)*10,
			'_sub_',sub_per*100,'_pCorn_',p_price*100,'_pN_',p_N*100,'_APH0_',APH_0,'.rds',sep=''))
	}
}



