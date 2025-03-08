######################################
# Compare how different actual APH is from the approximated APH
######################################
# written by Taro Mieno on 07/18/2016

######################################
# Simulation with APH 
######################################
# written by Taro Mieno on 02/24/2016

#===================================
# 0. preliminary operations
#===================================
setwd('~/Dropbox/MyResearch/CropInsuranceProgram/MoralHazard/Simulation')
source('~/Dropbox/R_libraries_load/library.R')
library('zoo')
library('mc2d')

#--- crop insurance parameters and functions ---#
source('./Codes/Base/premium_rate_adj.R')
source('./Codes/Base/functions.R')

#--- county specific crop insurance parameters ---#
source('./Codes/Base/Sioux_IA_parameters.R')

#--- other parameters ---#
source('./Codes/Simulation/parameters.R')

#===================================
# 1. 
#===================================
set.seed(398574895)
N <- 180
p_N <- 3.14-0.0921*sqrt(N)+0.00603*N
q_N <- 12.30-1.353*sqrt(N)+0.0456*N
Y_mean <- qbetagen(runif(10000),p_N,q_N,48,202) %>% mean()

Y_0 <- rep(Y_mean,9)
U <- runif(50)
Y_gen <- qbetagen(U,p_N,q_N,48,202)
Y_true <- c(Y_0,Y_gen) 

#--- true APH ---#
APH_true <- data.table(
	APH = rollmean(Y_true,10),
	t=1:50,
	type='True'
	)

#--- approximate APH ---#
APH_app <- vector()
APH_app[1] <- Y_mean
for (i in 1:length(Y_gen)){
	APH_app[i+1] <- 0.9*APH_app[i]+0.1*Y_gen[i]
}

APH_app <- data.table(
	APH = APH_app[-1],
	t=1:50,
	type='Approximated'
	)
	
data <- rbind(APH_true,APH_app)
ggplot(data=data) +
	geom_line(aes(x=t,y=APH,color=type)) +
	ylim(100,150) +
	scale_color_discrete(name='') +
	theme(
		legend.position='bottom'	
		)
ggsave('./Graphs/true_vs_approximate_APH.pdf',width=6,height=3.5)
ggsave('./Graphs/true_vs_approximate_APH.png',width=6,height=3.5)


