#include <RcppArmadillo.h> 
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//------------------
// multivariate normal 
//------------------
// mu: mean
// sigma: covariance matrix
arma::mat mvrnormArma(int B, arma::vec mu, arma::mat sigma){
	int ncols = sigma.n_cols;
	arma::mat Y = arma::randn(B, ncols); 
	return arma::repmat(mu, 1, B).t() + Y * arma::chol(sigma);
}

arma::vec loss_simu(
	int B, double proj_price, double APH, double cov_level, 
	double adj_mean, double adj_sd, double ln_mean, double ln_var)
{
	//--- declare the storage of simulation results ---//
	arma::vec sim_results(2);

	//--- simulated yield ---//
	NumericVector sim_yield = pmax(0,rnorm(B,adj_mean,adj_sd));
	// arma::vec sim_yield = round(,digits=12);

	//--- simulated yield protection loss ---//
	sim_results(0) = mean(pmax(0,APH*cov_level-sim_yield));

	//--- simulated harvest price ---//	
	NumericVector sim_price = exp(rnorm(B)*sqrt(ln_var)+ln_mean);

	//--- simulated revenue protection loss ---//
	NumericVector min_price = pmin(2*proj_price,sim_price);
	NumericVector max_price = pmax(proj_price,min_price);
	sim_results(1) = mean(pmax(0,APH*cov_level*max_price-sim_yield*min_price));

	return sim_results;
}

double premium_calc(
	double rate_yield, double APH, double cov_level, double subsidy_percent, 
	double proj_price, double pvf, double acres, double county_rate, int year, 
	arma::vec c_rate_now, arma::vec c_rate_prev, arma::vec c_pars,
	arma::vec beta, arma::vec sec5, arma::vec sec6){

	//===================================
	// Glossary
	//===================================
	// c_rate_now: county-specific rate in the current year
	// c_rate_now(0): reference amount 
	// c_rate_now(1): exponent value 
	// c_rate_now(2): reference rate
	// c_rate_now(3): fixed rate 

	// c_rate_prev: county-specific rate in the previous year
	// c_rate_prev(0): reference amount
	// c_rate_prev(1): exponent value 
	// c_rate_prev(2): reference rate
	// c_rate_prev(3): fixed rate

	//--------------------------
	// County specific parameters
	//--------------------------
	// c_pars
	// c_pars(0): rate_differential_factor 
	// c_pars(1): py_rate_differential_factor 
	// c_pars(2): unit_residual_factor	
	// c_pars(3): py_unit_residual_factor
	
	//--------------------------
	// section 4
	//--------------------------	
	// sec4(0): additive_optional_rate_adj_factor 
	// sec4(1): multiplicative_optional_rate_adj_factor

	//--------------------------
	// Section 5
	//--------------------------
	// sec5(0): revenue_lookup_adj_factor
	// sec5(1): mean_quantity
	// sec5(2): sd_quantity

	//--------------------------
	// Section 6 (capping )
	//--------------------------	
	// sec6(0): capping_reference_yield
	// sec6(1): py_capping_reference_yield
	// sec6(2): capping_exponent_value
	// sec6(3): py_capping_exponent_value
	// sec6(4): capping_reference_rate
	// sec6(5): py_capping_reference_rate
	// sec6(6): capping_fixed_rate
	// sec6(7): py_capping_fixed_rate
	// sec6(8): capping_year
	
	//===================================
	// Main
	//===================================
	//--------------------------
	// Section 1: Liability Calculation
	//--------------------------
	// maximum potential cost to insurer
	double insured_share_percent = 1;
	double premium_liability_amount = APH*cov_level*proj_price*acres*insured_share_percent;

	//--------------------------
	// Section 2: Unit Discount Calculation
	//--------------------------
	// see county specific rate structure. 
	// ex. Chase_rate.xls for Chase county
	double optional_unit_discount_factor = 1;
	double unit_structure_discount_factor = optional_unit_discount_factor;
	
	//--------------------------
	// Section 3: Bare Premium Rate Calculation
	//--------------------------
	// rate yield: basically the same as APH (see the handbook for the definition)
	// reference amount: county specific (look up the rate information)
	// prior year reference amount (py_reference_amount): county specific (look up the rate information)

	double cy_reference_amount = c_rate_now(0);
	double py_reference_amount = c_rate_prev(0);
	double cy_yield_ratio = round(rate_yield/cy_reference_amount*1e2)/1e2; 
	double py_yield_ratio = round(rate_yield/py_reference_amount*1e2)/1e2;
	
	// exponent value: county specific (look up the rate information) 
	// prior year exponent value: county specific (look up the rate information) 

	double cy_exponent_value = c_rate_now(1);
	double py_exponent_value = c_rate_prev(1);
	double cy_rate_multiplier = pow(cy_yield_ratio,cy_exponent_value);
	double py_rate_multiplier = pow(py_yield_ratio,py_exponent_value);	

	// reference rate: county specific (look up the rate information) 
	// prior year reference rate: county specific (look up the rate information)
	// fixed rate: county specific (look up the rate information) 
	// prior year fixed rate: county specific (look up the rate information)

	double cy_ref_rate =	c_rate_now(2); 
	double cy_fixed_rate = c_rate_now(3); 
	double py_ref_rate = c_rate_prev(2); 
	double py_fixed_rate = c_rate_prev(3);
	double cy_base_rate = cy_rate_multiplier*cy_ref_rate+cy_fixed_rate;
	double py_base_rate = py_rate_multiplier*py_ref_rate+py_fixed_rate;
	
	// rate differential factor: what the hell is this?
	// unit residual factor: 1 for operational unit?
	double rate_differential_factor = c_pars(0);
	double py_rate_differential_factor = c_pars(1);
	double unit_residual_factor = c_pars(2);
	double py_unit_residual_factor = c_pars(3);
	double cy_base_premium_rate = cy_base_rate*rate_differential_factor*unit_residual_factor;
	double py_base_premium_rate = py_base_rate*py_rate_differential_factor*py_unit_residual_factor;
	
	double base_premium_rate = std::min(cy_base_premium_rate,py_base_premium_rate*1.2);
	// double revenue_lookup_rate = std::min(cy_base_rate,py_base_rate*1.2);

	//--------------------------
	// Section 4: Optional Coverage Calculation
	//--------------------------
	// additive_optional_rate_adj_factor
	// multiplicative_optional_rate_adj_factor
	double additive_optional_rate_adj_factor = 0;
	double multiplicative_optional_rate_adj_factor = 1;

	//--------------------------
	// Section 5: Revenue Coverage Add On Rates
	//--------------------------
	// mean quantity: ?
	// standard deviation quantity: ? 
	// revenue_lookup_adj_factor: ?
	// price volatility factor (pvf): set nationally
	double mean_quantity = sec5(0);
	double sd_quantity = sec5(1);
	// double revenue_lookup_adj_factor = 1;
	// double lookup_rate = round(revenue_lookup_rate*revenue_lookup_adj_factor*1e4)/1e4 ;
	double adj_mean_quantity = APH*mean_quantity/100;
	double adj_sd_quantity = APH*sd_quantity/100;
	double ln_var = log(pow(pvf,2)+1);
	double ln_mean = log(proj_price)-ln_var/2;
	
	//--- calculate simulated yield and revenue loss --//#
	arma::vec loss = loss_simu(10000,proj_price,APH,cov_level,adj_mean_quantity,adj_sd_quantity,ln_mean,ln_var);
	
	double simulated_yield_protection_base_premium_rate = loss(0)/(APH*cov_level);
	double simulated_revenue_protection_base_premium_rate = loss(1)/(APH*cov_level*proj_price);
	
	double preliminary_revenue_protection_add_on_rate = 
		std::max(simulated_revenue_protection_base_premium_rate-simulated_yield_protection_base_premium_rate,
			0.01*base_premium_rate);

	//--------------------------
	// Section 6: Historical Revenue Capping
	//--------------------------
	// Note: capping by historical rate does not happen for 
	// coverage level under 0.65.
	// capping reference yield?
	// prior capping reference yield?

	double capped_revenue_protection_add_on_rate;

	if (cov_level>=0.65){
			
		double capping_reference_yield = sec6(0);
		double py_capping_reference_yield = sec6(1);
		double capping_exponent_value = sec6(2);	
		double py_capping_exponent_value = sec6(3);
		double capping_reference_rate = sec6(4);
		double py_capping_reference_rate = sec6(5);
		double capping_fixed_rate = sec6(6);
		double py_capping_fixed_rate = sec6(7);
		double capping_year = sec6(8);

		//--- 6.01 --//#
		double capping_yield_ratio = rate_yield/capping_reference_yield;
		double py_capping_yield_ratio = rate_yield/py_capping_reference_yield;
	
		//--- 6.02 --//#
		double capping_rate_multiplier = pow(capping_yield_ratio,capping_exponent_value);
		double py_capping_rate_multiplier = pow(py_capping_yield_ratio,py_capping_exponent_value);
	
		//--- 6.03 --//#
		double historical_capping_base_rate = capping_rate_multiplier*capping_reference_rate+capping_fixed_rate;
		double historical_py_capping_base_rate = py_capping_rate_multiplier*py_capping_reference_rate+py_capping_fixed_rate;
	
		//--- 6.04 --//#
		arma::vec hist_rate_vec(3);
		hist_rate_vec(0) = 0.999;
		hist_rate_vec(1) = historical_py_capping_base_rate*1.2;
		hist_rate_vec(2) = historical_capping_base_rate;
		double hisotircal_basic_unit_base_rate = 0.9*arma::min(hist_rate_vec);
	
		//--- 6.05 --//#
		double historical_revenue_protection_base_premium_rate = 
		(beta(0) 
			+ round(beta(1)*hisotircal_basic_unit_base_rate*1e8)/1e8 
			+ round(beta(2)*pow(hisotircal_basic_unit_base_rate,2)*1e8)/1e8
			+ round(beta(3)*cov_level*1e8)/1e8
			+ round(beta(4)*pow(cov_level,2)*1e8)/1e8
			+ round(beta(5)*(APH/capping_reference_yield)*1e8)/1e8
			+ round(beta(6)*pow((APH/capping_reference_yield),2)*1e8)/1e8
			+ round(beta(7)*pvf*1e8)/1e8
			+ round(beta(8)*pow(pvf,2)*1e8)/1e8
			+ round(beta(9)*hisotircal_basic_unit_base_rate*cov_level*1e8)/1e8
			+ round(beta(10)*hisotircal_basic_unit_base_rate*(APH/capping_reference_yield)*1e8)/1e8
			+ round(beta(11)*hisotircal_basic_unit_base_rate*pvf*1e8)/1e8
			+ round(beta(12)*cov_level*(APH/capping_reference_yield)*1e8)/1e8
			+ round(beta(13)*cov_level*pvf*1e8)/1e8
			+ round(beta(14)*(APH/capping_reference_yield)*pvf*1e8)/1e8
		)*c_pars(2)*1.1;
	
		//--- 6.06 --//#
		double a1 = base_premium_rate + preliminary_revenue_protection_add_on_rate;
		double a2 = historical_revenue_protection_base_premium_rate*pow(1.2,(year-capping_year));
		capped_revenue_protection_add_on_rate = std::min(a1,a2)-base_premium_rate;
	} 
	else{
		capped_revenue_protection_add_on_rate = preliminary_revenue_protection_add_on_rate;
	}
	
	//--------------------------
	// Section 7: Whole Farm Unit Discount 
	//--------------------------

	//--------------------------
	// Section 8: Premium Rate
	//--------------------------
	//--- 8.01 --//#
	double capped_revenue_protection_add_on_factor = capped_revenue_protection_add_on_rate;
	
	//--- 8.02 --//#
	double premium_rate = std::min(base_premium_rate*unit_structure_discount_factor*multiplicative_optional_rate_adj_factor
		+additive_optional_rate_adj_factor+capped_revenue_protection_add_on_factor,0.999);
	
	//--------------------------
	// Section 9: Total premium, subsidy, and producer premium calculation
	//--------------------------
	double experience_factor =  1;
	double premium_surcharge_percent = 1;
	double total_premium_multiplicative_optional_rate_adj_factor = 1;
	double multiple_commodity_adj_factor = 1;
	double beg_fr_subsidy_percent = 0;
	
	//--- 9.01 --//#
	double preliminary_total_premium_amount = premium_liability_amount*premium_rate*experience_factor*premium_surcharge_percent*total_premium_multiplicative_optional_rate_adj_factor;
	
	//--- 9.02 --//#
	double total_premium_amount = preliminary_total_premium_amount*multiple_commodity_adj_factor;
	
	//--- 9.03 --//#
	double base_subsidy_amount = total_premium_amount*subsidy_percent;
	
	//--- 9.04 ---// 
	double beg_farmer_rancher_subsidy_amount = total_premium_amount*beg_fr_subsidy_percent;
	
	//--- 9.05 (Native Sod Subsidy)--//#
	double native_sod_subsidy_amount = 0;
	
	//--- 9.06 (Subsidy amount)--//#
	double subsidy_amount = base_subsidy_amount + beg_farmer_rancher_subsidy_amount + native_sod_subsidy_amount;
	
	//--- 9.07 (Producer premium amount) --//#
	double premium = total_premium_amount-subsidy_amount;

	return premium;
}

//===================================
// Constant Relative Risk Aversion 
//===================================
arma::vec crra(double rho, arma::vec x){
	arma::vec y = (pow(x,(1-rho))-1)/(1-rho);
	return y;
}

//===================================
// Simulation 
//===================================
// [[Rcpp::export]]
List sim_ci(
	int B_aph, int B, arma::vec mu, arma::mat sigma, 
	arma::vec hp_prop, arma::vec y_prop, double acres, double cost, 
	double cov_level, double pp, double pvf, double subsidy_percent, 
	int year, double county_TY, arma::vec c_rate_now, arma::vec c_rate_prev, 
	arma::vec c_pars, arma::vec beta, arma::vec sec5, arma::vec sec6
	)
{
	//------------------
	// Glossary
	//------------------
	// hp_prop: statistical property of harvest price
	// hp_prop(0) = log mean
	// hp_prop(1) = sqrt(log var)
	// y_prop: statistical property of harvest price
	// y_prop(0) = 1st shaper parameter
	// y_prop(1) = 2nd shaper parameter
	// y_prop(2) = min
	// y_prop(3) = max
	// pp: projected price
	// acres: acres
	// cov_level: coverage level

	//===================================
	// Calculate APH
	//===================================
	arma::mat u_aph_1 = runif(B_aph/2); 
	arma::mat u_aph_2 = 1-u_aph_1; 
	arma::mat u_aph = join_cols(u_aph_1,u_aph_2);
	arma::vec raw_yield_temp = qbeta(NumericVector(u_aph.begin(),u_aph.end()),y_prop(0),y_prop(1));
	arma::vec raw_yield = raw_yield_temp*(y_prop(3)-y_prop(2)) + y_prop(2);
	double rate_yield = mean(raw_yield);

	//------------------
	// YE, TA, YS
	//------------------
	//--- to be edited ---//
	double APH = rate_yield;

	//===================================
	// Premium
	//===================================
	double premium = premium_calc(
		rate_yield, APH, cov_level, subsidy_percent, pp, 
		pvf, acres, county_TY, year, c_rate_now, c_rate_prev, 
		c_pars, beta, sec5, sec6
		);

	//===================================
	// Generate harvest price and yield
	//===================================
	//--- generate multivariate normal ---//
	arma::mat mn = mvrnormArma(B/2,mu,sigma); 
	arma::vec v1 = mn.col(0);
	arma::vec v2 = mn.col(1);

	//--- transform two variables into uniform ---//
	arma::vec vec1_1 = pnorm(NumericVector(v1.begin(),v1.end()));
	arma::vec vec1_2 = 1-vec1_1;
	arma::vec vec2_1 = pnorm(NumericVector(v2.begin(),v2.end()));
	arma::vec vec2_2 = 1-vec2_1;
	arma::vec u1 = join_cols(vec1_1,vec1_2);
	arma::vec u2 = join_cols(vec2_1,vec2_2);

	//------------------
	// harvest price generation 
	//------------------
	NumericVector hp = qlnorm(NumericVector(u1.begin(),u1.end()),hp_prop(0),hp_prop(1));

	//------------------
	// yield generation
	//------------------
	NumericVector yield_temp = qbeta(NumericVector(u2.begin(),u2.end()),y_prop(0),y_prop(1));
	NumericVector yield = yield_temp*(y_prop(3)-y_prop(2)) + y_prop(2);

	//===================================
	// Revenue and Profit (with or without insurance)
	//===================================
	//--- raw revenue ---//
	NumericVector raw_revenue = hp*yield*acres;

	//--- revenue guarantee ---//
	NumericVector revenue_guarantee = APH*cov_level*acres*pmax(hp,pp);

	//--- revenue with indemnity payment ---//
	NumericVector revenue = pmax(raw_revenue,revenue_guarantee);

	//--- indemnity ---//
	NumericVector indemnity = pmax(revenue_guarantee-raw_revenue,0);

	//--- profit ---//
	NumericVector profit_insured = revenue - cost - premium;
	NumericVector profit_non = raw_revenue - cost;

	// //--- utility ---//
	// arma::vec u_insured = crra(rho,profit_insured-min_profit);
	// arma::vec u_non = crra(rho,profit_non-min_profit);

	return Rcpp::List::create(
		Rcpp::Named("pi_ins")=profit_insured,
		Rcpp::Named("pi_non")=profit_non,
		Rcpp::Named("yield")=yield,
		Rcpp::Named("indemnity")=indemnity,
		Rcpp::Named("revenue")=revenue,
		Rcpp::Named("revenue_guarantee")=revenue_guarantee,
		Rcpp::Named("raw_revenue")=raw_revenue,
		Rcpp::Named("h_price")=hp,
		Rcpp::Named("premium")=premium,
		Rcpp::Named("rate_yield")=rate_yield,
		Rcpp::Named("APH")=APH
		);
}


//--- 1. at T,   ---//
// find the optimal N  
// plug in the optimal N and simulate 
// value function for a series of values of 
// APH





