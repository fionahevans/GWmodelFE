// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// coordinate_rotate
arma::mat coordinate_rotate(arma::mat coords, double theta);
RcppExport SEXP _GWmodelFE_coordinate_rotate(SEXP coordsSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type coords(coordsSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(coordinate_rotate(coords, theta));
    return rcpp_result_gen;
END_RCPP
}
// eu_dist_mat
arma::mat eu_dist_mat(arma::mat in_locs, arma::mat out_locs);
RcppExport SEXP _GWmodelFE_eu_dist_mat(SEXP in_locsSEXP, SEXP out_locsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type in_locs(in_locsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type out_locs(out_locsSEXP);
    rcpp_result_gen = Rcpp::wrap(eu_dist_mat(in_locs, out_locs));
    return rcpp_result_gen;
END_RCPP
}
// eu_dist_smat
arma::mat eu_dist_smat(arma::mat in_locs);
RcppExport SEXP _GWmodelFE_eu_dist_smat(SEXP in_locsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type in_locs(in_locsSEXP);
    rcpp_result_gen = Rcpp::wrap(eu_dist_smat(in_locs));
    return rcpp_result_gen;
END_RCPP
}
// eu_dist_vec
arma::vec eu_dist_vec(arma::mat in_locs, arma::vec out_loc);
RcppExport SEXP _GWmodelFE_eu_dist_vec(SEXP in_locsSEXP, SEXP out_locSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type in_locs(in_locsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type out_loc(out_locSEXP);
    rcpp_result_gen = Rcpp::wrap(eu_dist_vec(in_locs, out_loc));
    return rcpp_result_gen;
END_RCPP
}
// md_dist_mat
arma::mat md_dist_mat(arma::mat in_locs, arma::mat out_locs);
RcppExport SEXP _GWmodelFE_md_dist_mat(SEXP in_locsSEXP, SEXP out_locsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type in_locs(in_locsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type out_locs(out_locsSEXP);
    rcpp_result_gen = Rcpp::wrap(md_dist_mat(in_locs, out_locs));
    return rcpp_result_gen;
END_RCPP
}
// md_dist_smat
arma::mat md_dist_smat(arma::mat in_locs);
RcppExport SEXP _GWmodelFE_md_dist_smat(SEXP in_locsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type in_locs(in_locsSEXP);
    rcpp_result_gen = Rcpp::wrap(md_dist_smat(in_locs));
    return rcpp_result_gen;
END_RCPP
}
// md_dist_vec
arma::vec md_dist_vec(arma::mat in_locs, arma::vec out_loc);
RcppExport SEXP _GWmodelFE_md_dist_vec(SEXP in_locsSEXP, SEXP out_locSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type in_locs(in_locsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type out_loc(out_locSEXP);
    rcpp_result_gen = Rcpp::wrap(md_dist_vec(in_locs, out_loc));
    return rcpp_result_gen;
END_RCPP
}
// cd_dist_mat
arma::mat cd_dist_mat(arma::mat in_locs, arma::mat out_locs);
RcppExport SEXP _GWmodelFE_cd_dist_mat(SEXP in_locsSEXP, SEXP out_locsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type in_locs(in_locsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type out_locs(out_locsSEXP);
    rcpp_result_gen = Rcpp::wrap(cd_dist_mat(in_locs, out_locs));
    return rcpp_result_gen;
END_RCPP
}
// cd_dist_smat
arma::mat cd_dist_smat(arma::mat in_locs);
RcppExport SEXP _GWmodelFE_cd_dist_smat(SEXP in_locsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type in_locs(in_locsSEXP);
    rcpp_result_gen = Rcpp::wrap(cd_dist_smat(in_locs));
    return rcpp_result_gen;
END_RCPP
}
// cd_dist_vec
arma::vec cd_dist_vec(arma::mat in_locs, arma::vec out_loc);
RcppExport SEXP _GWmodelFE_cd_dist_vec(SEXP in_locsSEXP, SEXP out_locSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type in_locs(in_locsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type out_loc(out_locSEXP);
    rcpp_result_gen = Rcpp::wrap(cd_dist_vec(in_locs, out_loc));
    return rcpp_result_gen;
END_RCPP
}
// mk_dist_mat
arma::mat mk_dist_mat(arma::mat in_locs, arma::mat out_locs, double p);
RcppExport SEXP _GWmodelFE_mk_dist_mat(SEXP in_locsSEXP, SEXP out_locsSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type in_locs(in_locsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type out_locs(out_locsSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(mk_dist_mat(in_locs, out_locs, p));
    return rcpp_result_gen;
END_RCPP
}
// mk_dist_smat
arma::mat mk_dist_smat(arma::mat in_locs, double p);
RcppExport SEXP _GWmodelFE_mk_dist_smat(SEXP in_locsSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type in_locs(in_locsSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(mk_dist_smat(in_locs, p));
    return rcpp_result_gen;
END_RCPP
}
// mk_dist_vec
arma::vec mk_dist_vec(arma::mat in_locs, arma::vec out_loc, double p);
RcppExport SEXP _GWmodelFE_mk_dist_vec(SEXP in_locsSEXP, SEXP out_locSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type in_locs(in_locsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type out_loc(out_locSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(mk_dist_vec(in_locs, out_loc, p));
    return rcpp_result_gen;
END_RCPP
}
// bisq_wt_vec
arma::vec bisq_wt_vec(arma::vec distv, double bw);
RcppExport SEXP _GWmodelFE_bisq_wt_vec(SEXP distvSEXP, SEXP bwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type distv(distvSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    rcpp_result_gen = Rcpp::wrap(bisq_wt_vec(distv, bw));
    return rcpp_result_gen;
END_RCPP
}
// bisq_wt_mat
arma::mat bisq_wt_mat(arma::mat distm, arma::vec bw);
RcppExport SEXP _GWmodelFE_bisq_wt_mat(SEXP distmSEXP, SEXP bwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type distm(distmSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bw(bwSEXP);
    rcpp_result_gen = Rcpp::wrap(bisq_wt_mat(distm, bw));
    return rcpp_result_gen;
END_RCPP
}
// gauss_wt_vec
arma::vec gauss_wt_vec(arma::vec distv, double bw);
RcppExport SEXP _GWmodelFE_gauss_wt_vec(SEXP distvSEXP, SEXP bwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type distv(distvSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    rcpp_result_gen = Rcpp::wrap(gauss_wt_vec(distv, bw));
    return rcpp_result_gen;
END_RCPP
}
// gauss_wt_mat
arma::mat gauss_wt_mat(arma::mat distm, arma::vec bw);
RcppExport SEXP _GWmodelFE_gauss_wt_mat(SEXP distmSEXP, SEXP bwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type distm(distmSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bw(bwSEXP);
    rcpp_result_gen = Rcpp::wrap(gauss_wt_mat(distm, bw));
    return rcpp_result_gen;
END_RCPP
}
// tri_wt_vec
arma::vec tri_wt_vec(arma::vec distv, double bw);
RcppExport SEXP _GWmodelFE_tri_wt_vec(SEXP distvSEXP, SEXP bwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type distv(distvSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    rcpp_result_gen = Rcpp::wrap(tri_wt_vec(distv, bw));
    return rcpp_result_gen;
END_RCPP
}
// tri_wt_mat
arma::mat tri_wt_mat(arma::mat distm, arma::vec bw);
RcppExport SEXP _GWmodelFE_tri_wt_mat(SEXP distmSEXP, SEXP bwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type distm(distmSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bw(bwSEXP);
    rcpp_result_gen = Rcpp::wrap(tri_wt_mat(distm, bw));
    return rcpp_result_gen;
END_RCPP
}
// exp_wt_vec
arma::vec exp_wt_vec(arma::vec distv, double bw);
RcppExport SEXP _GWmodelFE_exp_wt_vec(SEXP distvSEXP, SEXP bwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type distv(distvSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    rcpp_result_gen = Rcpp::wrap(exp_wt_vec(distv, bw));
    return rcpp_result_gen;
END_RCPP
}
// exp_wt_mat
arma::mat exp_wt_mat(arma::mat distm, arma::vec bw);
RcppExport SEXP _GWmodelFE_exp_wt_mat(SEXP distmSEXP, SEXP bwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type distm(distmSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bw(bwSEXP);
    rcpp_result_gen = Rcpp::wrap(exp_wt_mat(distm, bw));
    return rcpp_result_gen;
END_RCPP
}
// gw_reg
List gw_reg(arma::mat x, arma::vec y, arma::vec w, bool hatmatrix, int focus);
RcppExport SEXP _GWmodelFE_gw_reg(SEXP xSEXP, SEXP ySEXP, SEXP wSEXP, SEXP hatmatrixSEXP, SEXP focusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< bool >::type hatmatrix(hatmatrixSEXP);
    Rcpp::traits::input_parameter< int >::type focus(focusSEXP);
    rcpp_result_gen = Rcpp::wrap(gw_reg(x, y, w, hatmatrix, focus));
    return rcpp_result_gen;
END_RCPP
}
// trhat2
arma::vec trhat2(arma::mat S);
RcppExport SEXP _GWmodelFE_trhat2(SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(trhat2(S));
    return rcpp_result_gen;
END_RCPP
}
// fitted
arma::vec fitted(arma::mat X, arma::mat beta);
RcppExport SEXP _GWmodelFE_fitted(SEXP XSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(fitted(X, beta));
    return rcpp_result_gen;
END_RCPP
}
// ehat
arma::vec ehat(arma::vec y, arma::mat X, arma::mat beta);
RcppExport SEXP _GWmodelFE_ehat(SEXP ySEXP, SEXP XSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(ehat(y, X, beta));
    return rcpp_result_gen;
END_RCPP
}
// rss
double rss(arma::vec y, arma::mat X, arma::mat beta);
RcppExport SEXP _GWmodelFE_rss(SEXP ySEXP, SEXP XSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(rss(y, X, beta));
    return rcpp_result_gen;
END_RCPP
}
// gwr_diag
arma::vec gwr_diag(arma::vec y, arma::mat x, arma::mat beta, arma::mat S);
RcppExport SEXP _GWmodelFE_gwr_diag(SEXP ySEXP, SEXP xSEXP, SEXP betaSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(gwr_diag(y, x, beta, S));
    return rcpp_result_gen;
END_RCPP
}
// gwr_diag1
arma::vec gwr_diag1(arma::vec y, arma::mat x, arma::mat beta, arma::vec s_hat);
RcppExport SEXP _GWmodelFE_gwr_diag1(SEXP ySEXP, SEXP xSEXP, SEXP betaSEXP, SEXP s_hatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type s_hat(s_hatSEXP);
    rcpp_result_gen = Rcpp::wrap(gwr_diag1(y, x, beta, s_hat));
    return rcpp_result_gen;
END_RCPP
}
// AICc
double AICc(arma::vec y, arma::mat x, arma::mat beta, arma::mat S);
RcppExport SEXP _GWmodelFE_AICc(SEXP ySEXP, SEXP xSEXP, SEXP betaSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(AICc(y, x, beta, S));
    return rcpp_result_gen;
END_RCPP
}
// AICc1
double AICc1(arma::vec y, arma::mat x, arma::mat beta, arma::vec s_hat);
RcppExport SEXP _GWmodelFE_AICc1(SEXP ySEXP, SEXP xSEXP, SEXP betaSEXP, SEXP s_hatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type s_hat(s_hatSEXP);
    rcpp_result_gen = Rcpp::wrap(AICc1(y, x, beta, s_hat));
    return rcpp_result_gen;
END_RCPP
}
// AICc_rss
arma::vec AICc_rss(arma::vec y, arma::mat x, arma::mat beta, arma::mat S);
RcppExport SEXP _GWmodelFE_AICc_rss(SEXP ySEXP, SEXP xSEXP, SEXP betaSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(AICc_rss(y, x, beta, S));
    return rcpp_result_gen;
END_RCPP
}
// Ci_mat
arma::mat Ci_mat(arma::mat x, arma::vec w);
RcppExport SEXP _GWmodelFE_Ci_mat(SEXP xSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(Ci_mat(x, w));
    return rcpp_result_gen;
END_RCPP
}
// scgwr_pre
List scgwr_pre(arma::mat x, arma::vec y, int bw, int poly, double b0, arma::mat g0, arma::mat neighbour);
RcppExport SEXP _GWmodelFE_scgwr_pre(SEXP xSEXP, SEXP ySEXP, SEXP bwSEXP, SEXP polySEXP, SEXP b0SEXP, SEXP g0SEXP, SEXP neighbourSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< int >::type poly(polySEXP);
    Rcpp::traits::input_parameter< double >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type g0(g0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type neighbour(neighbourSEXP);
    rcpp_result_gen = Rcpp::wrap(scgwr_pre(x, y, bw, poly, b0, g0, neighbour));
    return rcpp_result_gen;
END_RCPP
}
// scgwr_loocv
double scgwr_loocv(arma::vec target, arma::mat x, arma::vec y, int bw, int poly, arma::mat Mx0, arma::mat My0, arma::mat XtX, arma::mat XtY);
RcppExport SEXP _GWmodelFE_scgwr_loocv(SEXP targetSEXP, SEXP xSEXP, SEXP ySEXP, SEXP bwSEXP, SEXP polySEXP, SEXP Mx0SEXP, SEXP My0SEXP, SEXP XtXSEXP, SEXP XtYSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type target(targetSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< int >::type poly(polySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Mx0(Mx0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type My0(My0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XtX(XtXSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XtY(XtYSEXP);
    rcpp_result_gen = Rcpp::wrap(scgwr_loocv(target, x, y, bw, poly, Mx0, My0, XtX, XtY));
    return rcpp_result_gen;
END_RCPP
}
// scgwr_reg
List scgwr_reg(arma::mat x, arma::vec y, int bw, int poly, arma::mat G0, arma::mat Mx0, arma::mat My0, arma::mat XtX, arma::mat XtY, arma::mat neighbour, arma::vec parameters);
RcppExport SEXP _GWmodelFE_scgwr_reg(SEXP xSEXP, SEXP ySEXP, SEXP bwSEXP, SEXP polySEXP, SEXP G0SEXP, SEXP Mx0SEXP, SEXP My0SEXP, SEXP XtXSEXP, SEXP XtYSEXP, SEXP neighbourSEXP, SEXP parametersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< int >::type poly(polySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type G0(G0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Mx0(Mx0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type My0(My0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XtX(XtXSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XtY(XtYSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type neighbour(neighbourSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type parameters(parametersSEXP);
    rcpp_result_gen = Rcpp::wrap(scgwr_reg(x, y, bw, poly, G0, Mx0, My0, XtX, XtY, neighbour, parameters));
    return rcpp_result_gen;
END_RCPP
}
// box_wt_vec
arma::vec box_wt_vec(arma::vec distv, double bw);
RcppExport SEXP _GWmodelFE_box_wt_vec(SEXP distvSEXP, SEXP bwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type distv(distvSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    rcpp_result_gen = Rcpp::wrap(box_wt_vec(distv, bw));
    return rcpp_result_gen;
END_RCPP
}
// box_wt_vec_ad
arma::vec box_wt_vec_ad(arma::vec distv, double bw);
RcppExport SEXP _GWmodelFE_box_wt_vec_ad(SEXP distvSEXP, SEXP bwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type distv(distvSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    rcpp_result_gen = Rcpp::wrap(box_wt_vec_ad(distv, bw));
    return rcpp_result_gen;
END_RCPP
}
// gau_wt_vec_ad
arma::vec gau_wt_vec_ad(arma::vec distv, double bw);
RcppExport SEXP _GWmodelFE_gau_wt_vec_ad(SEXP distvSEXP, SEXP bwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type distv(distvSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    rcpp_result_gen = Rcpp::wrap(gau_wt_vec_ad(distv, bw));
    return rcpp_result_gen;
END_RCPP
}
// bis_wt_vec_ad
arma::vec bis_wt_vec_ad(arma::vec distv, double bw);
RcppExport SEXP _GWmodelFE_bis_wt_vec_ad(SEXP distvSEXP, SEXP bwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type distv(distvSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    rcpp_result_gen = Rcpp::wrap(bis_wt_vec_ad(distv, bw));
    return rcpp_result_gen;
END_RCPP
}
// tri_wt_vec_ad
arma::vec tri_wt_vec_ad(arma::vec distv, double bw);
RcppExport SEXP _GWmodelFE_tri_wt_vec_ad(SEXP distvSEXP, SEXP bwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type distv(distvSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    rcpp_result_gen = Rcpp::wrap(tri_wt_vec_ad(distv, bw));
    return rcpp_result_gen;
END_RCPP
}
// exp_wt_vec_ad
arma::vec exp_wt_vec_ad(arma::vec distv, double bw);
RcppExport SEXP _GWmodelFE_exp_wt_vec_ad(SEXP distvSEXP, SEXP bwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type distv(distvSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    rcpp_result_gen = Rcpp::wrap(exp_wt_vec_ad(distv, bw));
    return rcpp_result_gen;
END_RCPP
}
// gw_weight
arma::vec gw_weight(arma::vec vdist, double bw, std::string kernel, bool adaptive);
RcppExport SEXP _GWmodelFE_gw_weight(SEXP vdistSEXP, SEXP bwSEXP, SEXP kernelSEXP, SEXP adaptiveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type vdist(vdistSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< bool >::type adaptive(adaptiveSEXP);
    rcpp_result_gen = Rcpp::wrap(gw_weight(vdist, bw, kernel, adaptive));
    return rcpp_result_gen;
END_RCPP
}
// gw_reg_2
arma::vec gw_reg_2(arma::mat x, arma::vec y, arma::vec w);
RcppExport SEXP _GWmodelFE_gw_reg_2(SEXP xSEXP, SEXP ySEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(gw_reg_2(x, y, w));
    return rcpp_result_gen;
END_RCPP
}
// gwr_q
arma::mat gwr_q(arma::mat x, arma::vec y, arma::mat dMat, double bw, std::string kernel, bool adaptive);
RcppExport SEXP _GWmodelFE_gwr_q(SEXP xSEXP, SEXP ySEXP, SEXP dMatSEXP, SEXP bwSEXP, SEXP kernelSEXP, SEXP adaptiveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dMat(dMatSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< bool >::type adaptive(adaptiveSEXP);
    rcpp_result_gen = Rcpp::wrap(gwr_q(x, y, dMat, bw, kernel, adaptive));
    return rcpp_result_gen;
END_RCPP
}
// e_vec
arma::vec e_vec(int m, int n);
RcppExport SEXP _GWmodelFE_e_vec(SEXP mSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(e_vec(m, n));
    return rcpp_result_gen;
END_RCPP
}
// gwr_mixed_trace
double gwr_mixed_trace(arma::mat x1, arma::mat x2, arma::vec y, arma::mat dMat, double bw, std::string kernel, bool adaptive);
RcppExport SEXP _GWmodelFE_gwr_mixed_trace(SEXP x1SEXP, SEXP x2SEXP, SEXP ySEXP, SEXP dMatSEXP, SEXP bwSEXP, SEXP kernelSEXP, SEXP adaptiveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dMat(dMatSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< bool >::type adaptive(adaptiveSEXP);
    rcpp_result_gen = Rcpp::wrap(gwr_mixed_trace(x1, x2, y, dMat, bw, kernel, adaptive));
    return rcpp_result_gen;
END_RCPP
}
// gwr_mixed_2
List gwr_mixed_2(arma::mat x1, arma::mat x2, arma::vec y, arma::mat dMat, double bw, std::string kernel, bool adaptive);
RcppExport SEXP _GWmodelFE_gwr_mixed_2(SEXP x1SEXP, SEXP x2SEXP, SEXP ySEXP, SEXP dMatSEXP, SEXP bwSEXP, SEXP kernelSEXP, SEXP adaptiveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dMat(dMatSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< bool >::type adaptive(adaptiveSEXP);
    rcpp_result_gen = Rcpp::wrap(gwr_mixed_2(x1, x2, y, dMat, bw, kernel, adaptive));
    return rcpp_result_gen;
END_RCPP
}
