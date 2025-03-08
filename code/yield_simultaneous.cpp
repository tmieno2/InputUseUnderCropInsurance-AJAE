#include <RcppArmadillo.h> 
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//===================================
// Simulation 
//===================================
// [[Rcpp::export]]
List yield(
	arma::vec y_prop, 
		)
{
	//------------------
	// Glossary
	//------------------
	// y_prop(0) = 1st shaper parameter
	// y_prop(1) = 2nd shaper parameter
	// y_prop(2) = min
	// y_prop(3) = max
	
	//===================================
	// Calculate APH
	//===================================
	arma::mat u_aph = runif(B_aph); 
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





