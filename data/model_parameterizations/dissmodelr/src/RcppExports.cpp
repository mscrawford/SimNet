// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// simulate
List simulate(std::string expParFilename);
RcppExport SEXP _dissmodelr_simulate(SEXP expParFilenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type expParFilename(expParFilenameSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate(expParFilename));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dissmodelr_simulate", (DL_FUNC) &_dissmodelr_simulate, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_dissmodelr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
