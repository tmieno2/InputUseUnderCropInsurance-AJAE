#include <RcppArmadillo.h> 
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
List premium_calc_YP(
	double rate_yield, double APH, double cov_level, double subsidy_percent, 
	double proj_price, double pvf, double acres, double county_rate, int year, 
	arma::vec c_rate_now, arma::vec c_rate_prev, arma::vec c_pars){

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
	// Section 8: Premium Rate
	//--------------------------
	//--- 8.01 --//#
	double capped_revenue_protection_add_on_factor = 0;
	
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

	return Rcpp::List::create(
		Rcpp::Named("premium")=premium,
		Rcpp::Named("premium_rate")=premium_rate
		);
}
