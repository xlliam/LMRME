//******************************************************************************************************//
//                                                                                                      //
//                                                                                                      //
//  This document includes two functions:                                                               //
//  1) double func(double x, double beta1, double Lambda) is used to  calcualte                         //
//  \phi_{K} / \phi_{U} in corrected kernel K^{*}.                                                      //
//  2）double FFT_AP(NumericVector& parms, NumericVector& x_cont, NumericVector& y_o,                   //
//                  const NumericVector& input, double Lambda,double beta, const NumericVector&         //
//                  mconst)                                                                             //
//  is used to compute the objective funciton in corrected kernel method by using fast                  //
//  fourier transformation.                                                                             //
//                                                                                                      //
//******************************************************************************************************//

//******************************************************************************************************//
//required library
//******************************************************************************************************//
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

//******************************************************************************************************//
// func function is used to compute \phi_{K} / \phi_{U} in corrected kenrel K^{*}
// Arguments:
// beta1  - slope parameter in mode regressioni
// Lambda - bandwidth in corrected kernel
//******************************************************************************************************//
// [[Rcpp::export]]
double func(double x, double beta1, double Lambda)
{
    if (  x <=1 && x >=-1 ) {
        return pow(1-pow(x,2),3) / (2*pow(Lambda,2) / ( 2*pow(Lambda,2) + pow(0.3390,2)*pow(beta1,2)*pow(x,2) ));
    }else{
        return 0;
    }
    
}

//******************************************************************************************************//
// FFT_AP function is is used to compute the objective funciton in corrected kernel method by
// using fast fourier transformation(FFT).
// Arguments:
// parms  - slope parameter in mode regressioni.
// x_cont - contaminated covariate W.
// y_o    - response variable Y.
// Lambda - bandwidth in corrected kernel.
// input  - m input values of K^{*}(t), m = 2^16 in our simulation.
// mconst - positive and negative sign corresponding to m input values in FFT.
// beta   - the interval in t for the m input values of K^{*}(t).
//******************************************************************************************************//
// [[Rcpp::export]]
double FFT_AP(NumericVector& parms, NumericVector& x_cont, NumericVector& y_o, const NumericVector& input, double Lambda,double beta, const NumericVector& mconst)
{
    int m = input.size();
    NumericVector FKoutput(m);
    NumericVector re = (y_o-parms[0]-parms[1]*x_cont)/Lambda;

    int m_mid = m/2;
    FKoutput[m_mid] = func(input[m_mid],parms[1],Lambda);
    int i_d = 1;
    double indicator1 = 1;
    double indicator2 = 1;
    while ( ( abs(indicator1) > 1e-30 || abs(indicator2) > 1e-30 ) &&  (i_d < m/2) ) {
        indicator1 = func(input[m_mid-i_d],parms[1],Lambda);
        indicator2 = func(input[m_mid+i_d],parms[1],Lambda);
        FKoutput[m_mid-i_d] = indicator1;
        FKoutput[m_mid+i_d] = indicator2;
        i_d=i_d+1;
    }

    arma::vec FK_op(FKoutput.begin(), m, false);
    arma::vec FK_i_f(m,fill::zeros);
    arma::cx_vec FK_f(FK_op,FK_i_f);
    arma::vec mcon(const_cast<NumericVector&>(mconst).begin(), m, false);
   
    arma::cx_vec fXF = mcon*beta%arma::ifft(mcon%FK_f)/2/PI*m;
    
    arma::vec fhat = arma::real(fXF);

    int n_re = re.size();
    NumericVector f_estimate(n_re);
    
    for (int i=0; i<n_re; i++) {
        int ind = (int)(round(re[i]/beta+m_mid));
        f_estimate[i] = fhat[ind];
    }
    
    double f_e = Rcpp::sum(f_estimate)/n_re/-Lambda ;
    
    return f_e;
    
}


