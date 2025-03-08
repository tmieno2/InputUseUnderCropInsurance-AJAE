#===================================
# Define functions that simulate N and APH paths
#===================================
#--------------------------
# profit calculation
#--------------------------
pi_calc <- function(hp,yield,N,APH,premium,I_type,U_type){
 	raw_revenue <- hp*yield*acres;
 	N_cost <- p_N*N*acres

	if (I_type=='YP'){
	 	revenue_YP <- raw_revenue + pmax(0,p_price*(APH*cov_level-yield)*acres);
		profit <- revenue_YP - N_cost - premium

	} else if (I_type=='RPHPE'){
		revenue_RPHPE <- pmax(raw_revenue,APH*cov_level*acres*p_price);
		profit <- revenue_RPHPE - N_cost - premium

	} else if (I_type=='RP'){
		revenue_RP <- pmax(raw_revenue,APH*cov_level*acres*pmax(hp,p_price));
		profit <- revenue_RP - N_cost - premium

	} else{
		print('Insurance type is not correct')
	}

	Util <- utility_gen(U_type)
	return(Util(profit))
}

#--------------------------
# preparation
#--------------------------
sim_path_APH_adj <- function(k){

	U_type <- ui_pairs_list[k,1]	
	I_type <- ui_pairs_list[k,2]
	
	#--- pre-allocate values ---#
	sim_data_N <- matrix(0,B_sim,T)
	sim_data_APH <- matrix(0,B_sim,T)
	sim_data_APH[,1] <- APH_0
	sim_data_yield <- matrix(0,B_sim,T)
	sim_data_hp <- matrix(0,B_sim,T)

	for (t in 1:T){
	
		print(paste(t,'/',T,sep=''))

		#--- find optimal N ---#
		sim_data_N[,t] <- predict.gam(opt_N_path[[k]][[T+1-t]],newdata=data.table(APH=sim_data_APH[,t])) 

		#--- simulate yield and harvest price ---#	
		half_U <- rmvnorm(mean=c(0,0),sig=sigma,n=B_sim/2) %>% pnorm() 
		U <- rbind(half_U,1-half_U)  
		price_U <- U[,1]
		yield_U <- U[,2] 

		#--- harvest price ---#
		sim_data_hp[,t] <- qlnorm(price_U,meanlog=ln_mean_p,sdlog=sqrt(ln_var_p))

		#--- simulate yield and APH---#	
		sim_data_yield[,t] <- mclapply(
			1:B_sim, 
			function(x) qbetagen(
				yield_U[x],p_pars[1]+p_pars[2]*sqrt(sim_data_N[,t][x])+p_pars[3]*sim_data_N[,t][x],
				q_pars[1]+q_pars[2]*sqrt(sim_data_N[,t][x])+q_pars[3]*sim_data_N[,t][x],min_y,max_y)
			,mc.cores=5) %>% 
			unlist()

		if(t<T){
			sim_data_APH[,t+1] <- sim_data_APH[,t]*(n_window-1)/n_window + sim_data_yield[,t]/n_window
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
	
	#--------------------------
	# APH summary
	#--------------------------
	colnames(sim_data_APH) <- paste('APH_',seq(1,T),sep='')
	APH_sim <- data.table(sim_data_APH) %>% 
		melt() %>% 
		mutate(variable=as.numeric(gsub('APH_','',variable))) %>% 
		data.table() %>%
		setnames(c('variable','value'),c('t','APH')) 
	
	#--------------------------
	# yield summary
	#--------------------------
	colnames(sim_data_yield) <- paste('yield_',seq(1,T),sep='')
	yield_sim <- data.table(sim_data_yield) %>% 
		melt() %>% 
		mutate(variable=as.numeric(gsub('yield_','',variable))) %>% 
		data.table() %>%
		setnames(c('variable','value'),c('t','yield')) 
	
	#--------------------------
	# pi summary
	#--------------------------
	colnames(sim_data_hp) <- paste('hp_',seq(1,T),sep='')
	hp_sim <- data.table(sim_data_hp) %>% 
		melt() %>% 
		mutate(variable=as.numeric(gsub('hp_','',variable))) %>% 
		data.table() %>%
		setnames(c('variable','value'),c('t','hp')) 

	#--------------------------
	# merge all the data
	#--------------------------	
	temp_path <- cbind(N_sim,APH_sim[,.(APH)],yield_sim[,.(yield)],hp_sim[,.(hp)]) 
	temp_path[,id:=1:B_sim]

	#--------------------------
	# premium calculation
	#--------------------------
	cov_text <- paste('cov_',cov_level,sep='')
	if (I_type=='YP'){

		temp_path[,premium:=predict.gam(premium_YP_gam[[cov_text]],newdata=data.table(APH))]

	} else if (I_type=='RPHPE'){

		temp_path[,premium:=predict.gam(premium_RPHPE_gam[[cov_text]],newdata=data.table(APH))]

	} else if (I_type=='RP'){

		temp_path[,premium:=predict.gam(premium_RP_gam[[cov_text]],newdata=data.table(APH))]
		
	}
	
	#--------------------------
	# Utility calculation
	#--------------------------	
	temp_path[,u:=pi_calc(hp,yield,N,APH,premium,I_type,U_type)]

	path <- temp_path[,.(N=mean(N),yield=mean(yield),APH=mean(APH),u=mean(u)),by=t]
	path[,coverage:=cov_level]
	path[,subsidy:=sub_per]
	path[,rho:=rho]
	path[,min_y:=min_y]
	path[,max_y:=max_y]
	path[,utility:=U_type]
	path[,insurance:=I_type]

	return(path)

}


