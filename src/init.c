#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


extern SEXP _GWmodelFE_AICc(SEXP, SEXP, SEXP, SEXP);
extern SEXP _GWmodelFE_AICc_rss(SEXP, SEXP, SEXP, SEXP);
extern SEXP _GWmodelFE_bisq_wt_mat(SEXP, SEXP);
extern SEXP _GWmodelFE_bisq_wt_vec(SEXP, SEXP);
extern SEXP _GWmodelFE_cd_dist_mat(SEXP, SEXP);
extern SEXP _GWmodelFE_cd_dist_smat(SEXP);
extern SEXP _GWmodelFE_cd_dist_vec(SEXP, SEXP);
extern SEXP _GWmodelFE_Ci_mat(SEXP, SEXP);
extern SEXP _GWmodelFE_coordinate_rotate(SEXP, SEXP);
extern SEXP _GWmodelFE_ehat(SEXP, SEXP, SEXP);
extern SEXP _GWmodelFE_eu_dist_mat(SEXP, SEXP);
extern SEXP _GWmodelFE_eu_dist_smat(SEXP);
extern SEXP _GWmodelFE_eu_dist_vec(SEXP, SEXP);
extern SEXP _GWmodelFE_exp_wt_mat(SEXP, SEXP);
extern SEXP _GWmodelFE_exp_wt_vec(SEXP, SEXP);
extern SEXP _GWmodelFE_fitted(SEXP, SEXP);
extern SEXP _GWmodelFE_gauss_wt_mat(SEXP, SEXP);
extern SEXP _GWmodelFE_gauss_wt_vec(SEXP, SEXP);
extern SEXP _GWmodelFE_gw_reg(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GWmodelFE_gwr_diag(SEXP, SEXP, SEXP, SEXP);
extern SEXP _GWmodelFE_md_dist_mat(SEXP, SEXP);
extern SEXP _GWmodelFE_md_dist_smat(SEXP);
extern SEXP _GWmodelFE_md_dist_vec(SEXP, SEXP);
extern SEXP _GWmodelFE_mk_dist_mat(SEXP, SEXP, SEXP);
extern SEXP _GWmodelFE_mk_dist_smat(SEXP, SEXP);
extern SEXP _GWmodelFE_mk_dist_vec(SEXP, SEXP, SEXP);
extern SEXP _GWmodelFE_rss(SEXP, SEXP, SEXP);
extern SEXP _GWmodelFE_tri_wt_mat(SEXP, SEXP);
extern SEXP _GWmodelFE_tri_wt_vec(SEXP, SEXP);
extern SEXP _GWmodelFE_scgwr_pre(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GWmodelFE_scgwr_loocv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GWmodelFE_scgwr_reg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GWmodelFE_AICc1(SEXP, SEXP, SEXP, SEXP);
extern SEXP _GWmodelFE_gwr_diag1(SEXP, SEXP, SEXP, SEXP);
extern SEXP _GWmodelFE_BIC(SEXP, SEXP, SEXP, SEXP);
extern SEXP _GWmodelFE_box_wt_vec(SEXP, SEXP);
extern SEXP _GWmodelFE_box_wt_vec_ad(SEXP, SEXP);
extern SEXP _GWmodelFE_gau_wt_vec_ad(SEXP, SEXP);
extern SEXP _GWmodelFE_bis_wt_vec_ad(SEXP, SEXP);
extern SEXP _GWmodelFE_tri_wt_vec_ad(SEXP, SEXP);
extern SEXP _GWmodelFE_exp_wt_vec_ad(SEXP, SEXP);
extern SEXP _GWmodelFE_gw_reg_2(SEXP, SEXP, SEXP);
extern SEXP _GWmodelFE_gwr_q(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GWmodelFE_e_vec(SEXP, SEXP);
extern SEXP _GWmodelFE_gwr_mixed_trace(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _GWmodelFE_gwr_mixed_2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_GWmodelFE_AICc",              (DL_FUNC) &_GWmodelFE_AICc,              4},
    {"_GWmodelFE_AICc_rss",          (DL_FUNC) &_GWmodelFE_AICc_rss,          4},
    {"_GWmodelFE_bisq_wt_mat",       (DL_FUNC) &_GWmodelFE_bisq_wt_mat,       2},
    {"_GWmodelFE_bisq_wt_vec",       (DL_FUNC) &_GWmodelFE_bisq_wt_vec,       2},
    {"_GWmodelFE_cd_dist_mat",       (DL_FUNC) &_GWmodelFE_cd_dist_mat,       2},
    {"_GWmodelFE_cd_dist_smat",      (DL_FUNC) &_GWmodelFE_cd_dist_smat,      1},
    {"_GWmodelFE_cd_dist_vec",       (DL_FUNC) &_GWmodelFE_cd_dist_vec,       2},
    {"_GWmodelFE_Ci_mat",            (DL_FUNC) &_GWmodelFE_Ci_mat,            2},
    {"_GWmodelFE_coordinate_rotate", (DL_FUNC) &_GWmodelFE_coordinate_rotate, 2},
    {"_GWmodelFE_ehat",              (DL_FUNC) &_GWmodelFE_ehat,              3},
    {"_GWmodelFE_eu_dist_mat",       (DL_FUNC) &_GWmodelFE_eu_dist_mat,       2},
    {"_GWmodelFE_eu_dist_smat",      (DL_FUNC) &_GWmodelFE_eu_dist_smat,      1},
    {"_GWmodelFE_eu_dist_vec",       (DL_FUNC) &_GWmodelFE_eu_dist_vec,       2},
    {"_GWmodelFE_exp_wt_mat",        (DL_FUNC) &_GWmodelFE_exp_wt_mat,        2},
    {"_GWmodelFE_exp_wt_vec",        (DL_FUNC) &_GWmodelFE_exp_wt_vec,        2},
    {"_GWmodelFE_fitted",            (DL_FUNC) &_GWmodelFE_fitted,            2},
    {"_GWmodelFE_gauss_wt_mat",      (DL_FUNC) &_GWmodelFE_gauss_wt_mat,      2},
    {"_GWmodelFE_gauss_wt_vec",      (DL_FUNC) &_GWmodelFE_gauss_wt_vec,      2},
    {"_GWmodelFE_gw_reg",            (DL_FUNC) &_GWmodelFE_gw_reg,            5},
    {"_GWmodelFE_gwr_diag",          (DL_FUNC) &_GWmodelFE_gwr_diag,          4},
    {"_GWmodelFE_md_dist_mat",       (DL_FUNC) &_GWmodelFE_md_dist_mat,       2},
    {"_GWmodelFE_md_dist_smat",      (DL_FUNC) &_GWmodelFE_md_dist_smat,      1},
    {"_GWmodelFE_md_dist_vec",       (DL_FUNC) &_GWmodelFE_md_dist_vec,       2},
    {"_GWmodelFE_mk_dist_mat",       (DL_FUNC) &_GWmodelFE_mk_dist_mat,       3},
	{"_GWmodelFE_mk_dist_smat",      (DL_FUNC) &_GWmodelFE_mk_dist_smat,      2},
    {"_GWmodelFE_mk_dist_vec",       (DL_FUNC) &_GWmodelFE_mk_dist_vec,       3},
    {"_GWmodelFE_rss",               (DL_FUNC) &_GWmodelFE_rss,               3},
    {"_GWmodelFE_tri_wt_mat",        (DL_FUNC) &_GWmodelFE_tri_wt_mat,        2},
    {"_GWmodelFE_tri_wt_vec",        (DL_FUNC) &_GWmodelFE_tri_wt_vec,        2},
    {"_GWmodelFE_scgwr_pre",         (DL_FUNC) &_GWmodelFE_scgwr_pre,         7},
    {"_GWmodelFE_scgwr_loocv",       (DL_FUNC) &_GWmodelFE_scgwr_loocv,       9},
    {"_GWmodelFE_scgwr_reg",         (DL_FUNC) &_GWmodelFE_scgwr_reg,         11},
	{"_GWmodelFE_AICc1",              (DL_FUNC) &_GWmodelFE_AICc1,             4},
	{"_GWmodelFE_gwr_diag1",          (DL_FUNC) &_GWmodelFE_gwr_diag1,         4},
	{"_GWmodelFE_BIC",               (DL_FUNC) &_GWmodelFE_BIC,                4},
	{"_GWmodelFE_box_wt_vec",        (DL_FUNC) &_GWmodelFE_box_wt_vec,         2},
	{"_GWmodelFE_box_wt_vec_ad",     (DL_FUNC) &_GWmodelFE_box_wt_vec_ad,      2},
	{"_GWmodelFE_gau_wt_vec_ad",     (DL_FUNC) &_GWmodelFE_gau_wt_vec_ad,      2},
	{"_GWmodelFE_bis_wt_vec_ad",     (DL_FUNC) &_GWmodelFE_bis_wt_vec_ad,      2},
	{"_GWmodelFE_tri_wt_vec_ad",     (DL_FUNC) &_GWmodelFE_tri_wt_vec_ad,      2},
	{"_GWmodelFE_exp_wt_vec_ad",     (DL_FUNC) &_GWmodelFE_exp_wt_vec_ad,      2},
	{"_GWmodelFE_gw_reg_2",          (DL_FUNC) &_GWmodelFE_gw_reg_2,           3},
	{"_GWmodelFE_gwr_q",             (DL_FUNC) &_GWmodelFE_gwr_q,              6},
	{"_GWmodelFE_e_vec",             (DL_FUNC) &_GWmodelFE_e_vec,              2},
	{"_GWmodelFE_gwr_mixed_trace",   (DL_FUNC) &_GWmodelFE_gwr_mixed_trace,    7},
	{"_GWmodelFE_gwr_mixed_2",       (DL_FUNC) &_GWmodelFE_gwr_mixed_2,        8},
    {NULL, NULL, 0}
};

void R_init_GWmodelFE(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
