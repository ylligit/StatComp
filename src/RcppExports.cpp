// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rwMetropolisC
NumericVector rwMetropolisC(int N, double sigma, double x0);
RcppExport SEXP _StatComp20006_rwMetropolisC(SEXP NSEXP, SEXP sigmaSEXP, SEXP x0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    rcpp_result_gen = Rcpp::wrap(rwMetropolisC(N, sigma, x0));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_StatComp20006_rwMetropolisC", (DL_FUNC) &_StatComp20006_rwMetropolisC, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_StatComp20006(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
