// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Delta_fast
Rcpp::List Delta_fast(const arma::vec& Avec, const arma::mat& A_test, const arma::mat& Z_matrix, int iStar, int g, int h, int k, int l, const arma::uvec& ind_iStar, int n_kl);
RcppExport SEXP _noisysbmGGM_Delta_fast(SEXP AvecSEXP, SEXP A_testSEXP, SEXP Z_matrixSEXP, SEXP iStarSEXP, SEXP gSEXP, SEXP hSEXP, SEXP kSEXP, SEXP lSEXP, SEXP ind_iStarSEXP, SEXP n_klSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type Avec(AvecSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A_test(A_testSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z_matrix(Z_matrixSEXP);
    Rcpp::traits::input_parameter< int >::type iStar(iStarSEXP);
    Rcpp::traits::input_parameter< int >::type g(gSEXP);
    Rcpp::traits::input_parameter< int >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type ind_iStar(ind_iStarSEXP);
    Rcpp::traits::input_parameter< int >::type n_kl(n_klSEXP);
    rcpp_result_gen = Rcpp::wrap(Delta_fast(Avec, A_test, Z_matrix, iStar, g, h, k, l, ind_iStar, n_kl));
    return rcpp_result_gen;
END_RCPP
}
// Delta_orange_5_fast
arma::vec Delta_orange_5_fast(const arma::vec& dataVec, const arma::mat& Z_matrix, int iStar, int g, int h, int k, int l, const arma::vec& I_kl, int n_kl, int n_kl_test, const arma::vec& ind_iStar, const arma::vec& k_add, const arma::vec& l_remove);
RcppExport SEXP _noisysbmGGM_Delta_orange_5_fast(SEXP dataVecSEXP, SEXP Z_matrixSEXP, SEXP iStarSEXP, SEXP gSEXP, SEXP hSEXP, SEXP kSEXP, SEXP lSEXP, SEXP I_klSEXP, SEXP n_klSEXP, SEXP n_kl_testSEXP, SEXP ind_iStarSEXP, SEXP k_addSEXP, SEXP l_removeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type dataVec(dataVecSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z_matrix(Z_matrixSEXP);
    Rcpp::traits::input_parameter< int >::type iStar(iStarSEXP);
    Rcpp::traits::input_parameter< int >::type g(gSEXP);
    Rcpp::traits::input_parameter< int >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type I_kl(I_klSEXP);
    Rcpp::traits::input_parameter< int >::type n_kl(n_klSEXP);
    Rcpp::traits::input_parameter< int >::type n_kl_test(n_kl_testSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ind_iStar(ind_iStarSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type k_add(k_addSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type l_remove(l_removeSEXP);
    rcpp_result_gen = Rcpp::wrap(Delta_orange_5_fast(dataVec, Z_matrix, iStar, g, h, k, l, I_kl, n_kl, n_kl_test, ind_iStar, k_add, l_remove));
    return rcpp_result_gen;
END_RCPP
}
// Delta_SBM_2_fast
arma::vec Delta_SBM_2_fast(const arma::vec& Z, const arma::mat& Z_matrix, int iStar, int g, int h, int k, int l, int n0, double eta0, double zeta0, int n_kl, int n_kl_test);
RcppExport SEXP _noisysbmGGM_Delta_SBM_2_fast(SEXP ZSEXP, SEXP Z_matrixSEXP, SEXP iStarSEXP, SEXP gSEXP, SEXP hSEXP, SEXP kSEXP, SEXP lSEXP, SEXP n0SEXP, SEXP eta0SEXP, SEXP zeta0SEXP, SEXP n_klSEXP, SEXP n_kl_testSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z_matrix(Z_matrixSEXP);
    Rcpp::traits::input_parameter< int >::type iStar(iStarSEXP);
    Rcpp::traits::input_parameter< int >::type g(gSEXP);
    Rcpp::traits::input_parameter< int >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type n0(n0SEXP);
    Rcpp::traits::input_parameter< double >::type eta0(eta0SEXP);
    Rcpp::traits::input_parameter< double >::type zeta0(zeta0SEXP);
    Rcpp::traits::input_parameter< int >::type n_kl(n_klSEXP);
    Rcpp::traits::input_parameter< int >::type n_kl_test(n_kl_testSEXP);
    rcpp_result_gen = Rcpp::wrap(Delta_SBM_2_fast(Z, Z_matrix, iStar, g, h, k, l, n0, eta0, zeta0, n_kl, n_kl_test));
    return rcpp_result_gen;
END_RCPP
}
// I_fast
Rcpp::List I_fast(const arma::umat& ind_all, const arma::vec& Avec, const arma::mat& Z_matrix, int k, int l);
RcppExport SEXP _noisysbmGGM_I_fast(SEXP ind_allSEXP, SEXP AvecSEXP, SEXP Z_matrixSEXP, SEXP kSEXP, SEXP lSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::umat& >::type ind_all(ind_allSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Avec(AvecSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z_matrix(Z_matrixSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    rcpp_result_gen = Rcpp::wrap(I_fast(ind_all, Avec, Z_matrix, k, l));
    return rcpp_result_gen;
END_RCPP
}
// convertNodePair_fast
arma::vec convertNodePair_fast(int i, const arma::vec& j, int p, bool directed);
RcppExport SEXP _noisysbmGGM_convertNodePair_fast(SEXP iSEXP, SEXP jSEXP, SEXP pSEXP, SEXP directedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type j(jSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< bool >::type directed(directedSEXP);
    rcpp_result_gen = Rcpp::wrap(convertNodePair_fast(i, j, p, directed));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_noisysbmGGM_Delta_fast", (DL_FUNC) &_noisysbmGGM_Delta_fast, 10},
    {"_noisysbmGGM_Delta_orange_5_fast", (DL_FUNC) &_noisysbmGGM_Delta_orange_5_fast, 13},
    {"_noisysbmGGM_Delta_SBM_2_fast", (DL_FUNC) &_noisysbmGGM_Delta_SBM_2_fast, 12},
    {"_noisysbmGGM_I_fast", (DL_FUNC) &_noisysbmGGM_I_fast, 5},
    {"_noisysbmGGM_convertNodePair_fast", (DL_FUNC) &_noisysbmGGM_convertNodePair_fast, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_noisysbmGGM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
