#===================================
# Define functions that simulate N and APH paths
#===================================
#--------------------------
# preparation
#--------------------------
sim_path <- function(i){

	U_type <- ui_pairs_list[i,1]	
	I_type <- ui_pairs_list[i,2]
	
	#--- pre-allocate values ---#
	sim_data_N <- matrix(0,B_sim,T)
	sim_data_APH <- matrix(0,B_sim,T)
	sim_data_APH[,1] <- APH_0

	for (t in 1:T){
	
		print(paste(t,'/',T,sep=''))

		cutoff_index <- APH_cutoff[[i]][[T+1-t]] < sim_data_APH[,t]
		sim_data_N[cutoff_index,t] <- 0

		#--- find optimal N ---#
		sim_data_N[!cutoff_index,t] <- predict.gam(
			opt_N_path[[i]][[T+1-t]],
			newdata=data.table(APH=sim_data_APH[!cutoff_index,t])
			) 
	
		#--- simulate yield and APH---#	
		yield <- mclapply(
			sim_data_N[,t], 
			function(x) rbetagen(
				1,p_pars[1]+p_pars[2]*sqrt(x)+p_pars[3]*x,
				q_pars[1]+q_pars[2]*sqrt(x)+q_pars[3]*x,min_y,max_y)
			,mc.cores=2) %>% 
			unlist()
	
		if(t<T){
			sim_data_APH[,t+1] <- sim_data_APH[,t]*(n_window-1)/n_window + yield/n_window
		}
	
	}

	#--------------------------
	# summarize N path
	#--------------------------
	colnames(sim_data_N) <- paste('N_',seq(1,T),sep='')
	N_sim <- data.table(sim_data_N) %>% 
		melt() %>% 
		mutate(variable=as.numeric(gsub('N_','',variable))) %>%  
		data.table() %>% 
		setnames(c('variable','value'),c('t','N')) 
	
	EN <- N_sim[,list(N=mean(N)),by=t]
	EN_gam <- gam(N~s(t,k=5),data=EN)
	EN[,EN_sm:=predict.gam(EN_gam,newdata=data.table(t=seq(1,T)))] 

	#--------------------------
	# APH summary
	#--------------------------
	colnames(sim_data_APH) <- paste('APH_',seq(1,T),sep='')
	APH_sim <- data.table(sim_data_APH) %>% 
		melt() %>% 
		mutate(variable=as.numeric(gsub('APH_','',variable))) %>% 
		data.table() %>%
		setnames(c('variable','value'),c('t','APH')) 
	
	E_APH <- APH_sim[,list(APH=mean(APH)),by=t]
	E_APH_gam <- gam(APH~s(t,k=5),data=E_APH)
	E_APH[,APH_sm:=predict.gam(E_APH_gam,newdata=data.table(t=seq(1,T)))]

	#--------------------------
	# combine N and APH path
	#--------------------------
	path <- cbind(EN,E_APH[,.(APH,APH_sm)]) 

	#--------------------------
	# add other variables
	#--------------------------
	path[,coverage:=cov_level]
	path[,subsidy:=sub_per]
	path[,rho:=rho]
	path[,min_y:=min_y]
	path[,max_y:=max_y]
	path[,utility:=U_type]
	path[,insurance:=I_type]

	#--------------------------
	# return the results
	#--------------------------
	return(path)
}


