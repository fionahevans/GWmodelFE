// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
#include <iostream>
using namespace Rcpp;
using namespace arma;

//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------


const arma::mat mSum(2, 1, fill::ones);

//distance matrix calculation
//coords must be a matrix with 2 columns
// [[Rcpp::export]]
arma::mat coordinate_rotate(arma::mat coords, double theta)
{
  int n = coords.n_rows;
  arma::mat rotated_coords(n, 2);
  rotated_coords.col(0) = coords.col(0) * cos(theta) - coords.col(1) * sin(theta);
  rotated_coords.col(1) = coords.col(0) * sin(theta) + coords.col(1) * cos(theta);
  return rotated_coords;
}

//Eudclidean distance matrix
// [[Rcpp::export]]
arma::mat eu_dist_mat(arma::mat in_locs, arma::mat out_locs)
{
  int n_in = in_locs.n_rows;
  int n_out = out_locs.n_rows;
  arma::mat eu_dist(n_in, n_out);
  int i = 0, j = 0;
  for (i = 0; i < n_in; i++)
  {
    for (j = 0; j < n_out; j++)
    {
      eu_dist(i,j) = sum(pow(in_locs.row(i) - out_locs.row(j),2));
    }
  }
  return sqrt(eu_dist);
}
//symmetrical distance matrix
// [[Rcpp::export]]
arma::mat eu_dist_smat(arma::mat in_locs)
{
  int n = in_locs.n_rows;
  arma::mat eu_dist(n, n);
  for (int k = 0; k < n * n; k++)
  {
    int i = k / n, j = k % n;
    eu_dist(i, j) = sum(pow(in_locs.row(i) - in_locs.row(j), 2));
    eu_dist(j, i) = eu_dist(i, j);
  }
  return sqrt(eu_dist);
}

// [[Rcpp::export]]
arma::vec eu_dist_vec(arma::mat in_locs, arma::vec out_loc)
{
  int n_in = in_locs.n_rows;
  arma::vec eu_dist(n_in);
  for (int i = 0; i < n_in; i++)
  {
    eu_dist(i) = sum(pow(in_locs.row(i) - trans(out_loc), 2));
  }
  return sqrt(eu_dist);
  // arma::mat v_span(n_in, 1, fill::ones);
  // arma::mat m_diff = in_locs - v_span * trans(out_loc);
  // return sqrt(m_diff % m_diff * mSum);
}

//Manhattan distance matrix
// [[Rcpp::export]]
arma::mat md_dist_mat(arma::mat in_locs, arma::mat out_locs)
{
  int n_in = in_locs.n_rows;
  int n_out = out_locs.n_rows;
  arma::mat md_dist(n_in, n_out);
  for (int i = 0; i < n_in; i++)
  {
    for (int j = 0; j < n_out; j++)
    {
      md_dist(i, j) = sum(abs(in_locs.row(i) - out_locs.row(j)));
    }
  }
  return md_dist;
}

//symmetrical distance matrix
// [[Rcpp::export]]
arma::mat md_dist_smat(arma::mat in_locs)
{
  int n = in_locs.n_rows;
  arma::mat md_dist(n, n);
  for (int i = 0; i < n; i++)
  {
    for (int j = i; j < n; j++)
    {
      md_dist(i, j) = sum(abs(in_locs.row(i) - in_locs.row(j)));
      md_dist(j, i) = md_dist(i, j);
    }
  }
  return md_dist;
}
// [[Rcpp::export]]
arma::vec md_dist_vec(arma::mat in_locs, arma::vec out_loc)
{
  int n_in = in_locs.n_rows;
  arma::vec md_dist(n_in);
  for (int i = 0; i < n_in; i++)
  {
    md_dist(i) = sum(abs(in_locs.row(i) - trans(out_loc)));
  }
  return md_dist;
}

//Chebyshev distance matrix
// [[Rcpp::export]]
arma::mat cd_dist_mat(arma::mat in_locs, arma::mat out_locs)
{
  int n_in = in_locs.n_rows;
  int n_out = out_locs.n_rows;
  arma::mat cd_dist(n_in, n_out);
  for (int i = 0; i < n_in; i++)
  {
    for (int j = i; j < n_out; j++)
    {
      cd_dist(i, j) = max(abs(in_locs.row(i) - out_locs.row(j)));
      cd_dist(j, i) = cd_dist(i, j);
    }
  }
  return cd_dist;
}

//symmetrical distance matrix
// [[Rcpp::export]]
arma::mat cd_dist_smat(arma::mat in_locs)
{
  int n = in_locs.n_rows;
  arma::mat cd_dist(n, n);
  for (int i = 0; i < n; i++)
  {
    for (int j = i; j < n; j++)
    {
      cd_dist(i, j) = max(abs(in_locs.row(i) - in_locs.row(j)));
      cd_dist(j, i) = cd_dist(i, j);
    }
  }
  return cd_dist;
}
// [[Rcpp::export]]
arma::vec cd_dist_vec(arma::mat in_locs, arma::vec out_loc)
{
  int n_in = in_locs.n_rows;
  arma::vec cd_dist(n_in);
  for (int i = 0; i < n_in; i++)
  {
    cd_dist(i) = max(abs(in_locs.row(i) - trans(out_loc)));
  }
  return cd_dist;
}

//Minkowski distance matrix
// [[Rcpp::export]]
arma::mat mk_dist_mat(arma::mat in_locs, arma::mat out_locs, double p)
{
  int n_in = in_locs.n_rows;
  int n_out = out_locs.n_rows;
  arma::mat mk_dist(n_in, n_out);
  for (int i = 0; i < n_in; i++)
  {
    for (int j = 0; j < n_out; j++)
    {
      mk_dist(i, j) = pow(sum(pow(abs(in_locs.row(i) - out_locs.row(j)), p)), 1.0 / p);
    }
  }
  return mk_dist;
}
//sqrt(sum(pow(in_locs.row(i) - trans(out_loc),2)))
//symmetrical distance matrix
// [[Rcpp::export]]
arma::mat mk_dist_smat(arma::mat in_locs, double p)
{
  int n = in_locs.n_rows;
  arma::mat mk_dist(n, n);
  for (int i = 0; i < n; i++)
  {
    for (int j = i; j < n; j++)
    {
      mk_dist(i, j) = pow(sum(pow(abs(in_locs.row(i) - in_locs.row(j)), p)), 1.0 / p);
      mk_dist(j, i) = mk_dist(i, j);
    }
  }
  return mk_dist;
}

// [[Rcpp::export]]
arma::vec mk_dist_vec(arma::mat in_locs, arma::vec out_loc, double p)
{
  int n_in = in_locs.n_rows;
  arma::vec mk_dist(n_in);
  for (int i = 0; i < n_in; i++)
  {
    mk_dist(i) = pow(sum(pow(abs(in_locs.row(i) - trans(out_loc)), p)), 1.0 / p);
  }
  return mk_dist;
}
//Weight matrix
//Bisuqare weight
// [[Rcpp::export]]
arma::vec bisq_wt_vec(arma::vec distv, double bw)
{
  int n = distv.n_elem;
  arma::vec wtv;
  wtv.zeros(n);
  // #pragma omp parallel for num_threads(omp_get_num_procs() - 1)
  for (int i = 0; i < n; i++)
  {
    if (distv(i) <= bw)
      wtv(i) = pow(1 - pow(distv(i) / bw, 2), 2);
  }
  return wtv;
}
//Calculated by column, the length of bw must be the same the number of columns of distm
// [[Rcpp::export]]
arma::mat bisq_wt_mat(arma::mat distm, arma::vec bw)
{
  int m = distm.n_cols;
  int n = distm.n_rows;
  arma::mat wtm;
  wtm.zeros(n, m);
  // #pragma omp parallel for num_threads(omp_get_num_procs() - 1)
  for (int k = 0; k < m * n; k++)
  {
    int i = k / n, j = k % n;
    if (distm(j, i) <= bw(i))
      wtm(j, i) = (1 - distm(j, i) / bw(i) * distm(j, i) / bw(i)) * (1 - distm(j, i) / bw(i) * distm(j, i) / bw(i));
  }
  return wtm;
}

//Gaussian weight
// [[Rcpp::export]]
arma::vec gauss_wt_vec(arma::vec distv, double bw)
{
  int n = distv.n_elem;
  arma::vec wtv;
  wtv.zeros(n);
  // #pragma omp parallel for num_threads(omp_get_num_procs() - 1)
  for (int i = 0; i < n; i++)
  {
    wtv(i) = exp(pow(distv(i), 2) / ((-2) * pow(bw, 2)));
  }
  return wtv;
}
// [[Rcpp::export]]
arma::mat gauss_wt_mat(arma::mat distm, arma::vec bw)
{
  int m = distm.n_cols;
  int n = distm.n_rows;
  arma::mat wtm;
  wtm.zeros(n, m);
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      wtm(j, i) = exp(pow(distm(j, i), 2) / ((-2) * pow(bw(i), 2)));
    }
  }
  return wtm;
}


//Tricube weight
// [[Rcpp::export]]
arma::vec tri_wt_vec(arma::vec distv, double bw)
{
  int n = distv.n_elem;
  arma::vec wtv;
  wtv.zeros(n);
  for (int i = 0; i < n; i++)
  {
    if (distv(i) <= bw)
      wtv(i) = pow(1 - pow(distv(i), 3) / pow(bw, 3), 3);
  }
  return wtv;
}
// [[Rcpp::export]]
arma::mat tri_wt_mat(arma::mat distm, arma::vec bw)
{
  int m = distm.n_cols;
  int n = distm.n_rows;
  arma::mat wtm;
  wtm.zeros(n, m);
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      if (distm(j, i) <= bw(i))
        wtm(j, i) = pow(1 - pow(distm(j, i), 3) / pow(bw(i), 3), 3);
    }
  }
  return wtm;
}
//exponential kernel weight adaptive=F
// [[Rcpp::export]]
arma::vec exp_wt_vec(arma::vec distv, double bw)
{
  int n = distv.n_elem;
  arma::vec wtv;
  wtv.zeros(n);
  for (int i = 0; i < n; i++)
  {
    wtv(i) = exp(-distv(i) / bw);
  }
  return wtv;
}
// [[Rcpp::export]]
arma::mat exp_wt_mat(arma::mat distm, arma::vec bw)
{
  int m = distm.n_cols;
  int n = distm.n_rows;
  arma::mat wtm;
  wtm.zeros(n, m);
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      wtm(j, i) = exp(-distm(j, i) / bw(i));
    }
  }
  return wtm;
}
//GWR calibration
// [[Rcpp::export]]
List gw_reg(arma::mat x, arma::vec y, arma::vec w, bool hatmatrix, int focus)
{
  arma::mat wspan(1, x.n_cols, fill::ones);
  arma::mat xtw = trans(x % (w * wspan));
  arma::mat xtwx = xtw * x;
  arma::mat xtwy = trans(x) * (w % y);
  arma::mat xtwx_inv = inv(xtwx);
  arma::vec beta = xtwx_inv * xtwy;
  if (hatmatrix)
  {
    arma::mat ci = xtwx_inv * xtw;
    arma::mat s_ri = x.row(focus - 1) * ci;
    return List::create(
      Named("beta") = beta,
      Named("S_ri") = s_ri,
      Named("Ci") = ci);
  }
  else
  {
    return List::create(
      Named("beta") = beta);
  }
}



//Trace of hatmatrix + trace of HH' in one function. Used in beta_se

// [[Rcpp::export]]
arma::vec trhat2(arma::mat S)
{
  int n_obs = S.n_rows;
  double htr = 0.0;
  double htr2 = 0.0;
  arma::vec result(2);
  for (int i = 0; i < n_obs; i++)
  {
    htr += S(i, i);
    htr2 += sum(S.row(i) % S.row(i));
  }
  result(0) = htr;
  result(1) = htr2;
  return result;
}
// Fited values

// [[Rcpp::export]]
arma::vec fitted(arma::mat X, arma::mat beta)
{
  arma::vec fitted = sum(beta % X, 1);
  return fitted;
}

// Residuals

// [[Rcpp::export]]
arma::vec ehat(arma::vec y, arma::mat X, arma::mat beta)
{
  arma::vec fitted = sum(beta % X, 1);
  return y - fitted;
}

// Residual sum-of squares
// [[Rcpp::export]]
double rss(arma::vec y, arma::mat X, arma::mat beta)
{
  arma::vec r = ehat(y, X, beta);
  return sum(r % r);
}

// [[Rcpp::export]]
arma::vec gwr_diag(arma::vec y,arma::mat x, arma::mat beta, arma::mat S) {
  double ss = rss(y,x,beta);
  arma::vec s_hat = trhat2(S);
  int n = S.n_rows;
  arma::vec result(9);
  double AIC = n*log(ss/n)+n*log(2*datum::pi)+n+s_hat(0); //AIC
  double AICc = n*log(ss/n)+n*log(2*datum::pi)+n*((n+s_hat(0))/(n-2-s_hat(0))); //AICc
  double edf = n-2*s_hat(0) + s_hat(1); //edf
  double enp = 2*s_hat(0) - s_hat(1); // enp
  double yss = sum(pow(y-mean(y),2)); //yss.g
  double r2 = 1 - ss/yss;
  double r2_adj = 1-(1-r2)*(n-1)/(edf-1);
  result(0) = AIC;
  result(1) = AICc;
  result(2) = edf;
  result(3) = enp;
  result(4) = ss;
  result(5) = r2;
  result(6) = r2_adj;
  result(7) = s_hat(0);
  result(8) = s_hat(1);
  return result;
  //return 2*enp + 2*n*log(ss/n) + 2*enp*(enp+1)/(n - enp - 1);
}

// [[Rcpp::export]]
arma::vec gwr_diag1(arma::vec y, arma::mat x, arma::mat beta, arma::vec s_hat)
{
  double ss = rss(y, x, beta);
  // arma::vec s_hat = trhat2(S);
  int n = x.n_rows;
  // arma::vec result(9);
  arma::vec result(7);
  double AIC = n * log(ss / n) + n * log(2 * datum::pi) + n + s_hat(0);																//AIC
  double AICc = n * log(ss / n) + n * log(2 * datum::pi) + n * ((n + s_hat(0)) / (n - 2 - s_hat(0))); //AICc
  double edf = n - 2 * s_hat(0) + s_hat(1);																														//edf
  double enp = 2 * s_hat(0) - s_hat(1);																																// enp
  double yss = sum(pow(y - mean(y), 2));																															//yss.g
  double r2 = 1 - ss / yss;
  double r2_adj = 1 - (1 - r2) * (n - 1) / (edf - 1);
  result(0) = AIC;
  result(1) = AICc;
  result(2) = edf;
  result(3) = enp;
  result(4) = ss;
  result(5) = r2;
  result(6) = r2_adj;
  // result(7) = s_hat(0);
  // result(8) = s_hat(1);
  return result;
  //return 2*enp + 2*n*log(ss/n) + 2*enp*(enp+1)/(n - enp - 1);
}

// return the AICc - nice tool to give to 'optimise' to find 'best' bandwidth
// [[Rcpp::export]]
double AICc(arma::vec y, arma::mat x, arma::mat beta, arma::mat S)
{
  double ss = rss(y, x, beta);
  arma::vec s_hat = trhat2(S);
  int n = S.n_rows;
  double AICc = n * log(ss / n) + n * log(2 * datum::pi) + n * ((n + s_hat(0)) / (n - 2 - s_hat(0))); //AICc
  return AICc;
  //return 2*enp + 2*n*log(ss/n) + 2*enp*(enp+1)/(n - enp - 1);
}

// [[Rcpp::export]]
double AICc1(arma::vec y, arma::mat x, arma::mat beta, arma::vec s_hat)
{
  double ss = rss(y, x, beta);
  // arma::vec s_hat = trhat2(S);
  int n = x.n_rows;
  double AICc = n * log(ss / n) + n * log(2 * datum::pi) + n * ((n + s_hat(0)) / (n - 2 - s_hat(0))); //AICc
  return AICc;
  //return 2*enp + 2*n*log(ss/n) + 2*enp*(enp+1)/(n - enp - 1);
}

// return the AICc and RSS , used for function model.selection
// [[Rcpp::export]]
arma::vec AICc_rss(arma::vec y, arma::mat x, arma::mat beta, arma::mat S)
{
  arma::vec result(3);
  double ss = rss(y, x, beta);
  result[0] = ss;
  arma::vec s_hat = trhat2(S);
  int n = S.n_rows;
  double AIC = n * log(ss / n) + n * log(2 * datum::pi) + n + s_hat(0);
  double AICc = n * log(ss / n) + n * log(2 * datum::pi) + n * ((n + s_hat(0)) / (n - 2 - s_hat(0))); //AICc
  result[1] = AIC;
  result[2] = AICc;
  return result;
  //return 2*enp + 2*n*log(ss/n) + 2*enp*(enp+1)/(n - enp - 1);
}

//Caculate the i row of
// [[Rcpp::export]]
arma::mat Ci_mat(arma::mat x, arma::vec w)
{
  return inv(trans(x) * diagmat(w) * x) * trans(x) * diagmat(w);
}
//Scalable GWR C++ functions
// [[Rcpp::export]]
List scgwr_pre(arma::mat x, arma::vec y, int bw, int poly, double b0, arma::mat g0, arma::mat neighbour) {
  int n = x.n_rows;
  int k = x.n_cols;
  arma::mat g0s(g0.n_cols + 1, g0.n_rows, fill::ones);
  arma::mat g0t = trans(g0);
  for (int i = 0; i < bw; i++) {
    g0s.row(i + 1) = g0t.row(i);
  }
  g0s = trans(g0s);
  arma::mat Mx0((poly + 1)*k*k, n, fill::zeros);
  arma::mat My0((poly + 1)*k, n, fill::zeros);
  arma::mat spanXnei(1, poly + 1, fill::ones);
  arma::mat spanXtG(1, k, fill::ones);
  for (int i = 0; i < n; i++) {
    arma::mat g(poly + 1, bw + 1, fill::ones);
    for (int p = 0; p < poly; p++) {
      g.row(p + 1) = pow(g0s.row(i), pow(2.0, poly/2.0)/pow(2.0, p + 1));
    }
    g = trans(g);
    g = g.rows(1, bw);
    arma::mat xnei(bw, k, fill::zeros);
    arma::vec ynei(bw, fill::zeros);
    for (int j = 0; j < bw; j++) {
      int inei = int(neighbour(i, j) - 1);
      xnei.row(j) = x.row(inei);
      ynei.row(j) = y(inei);
    }
    for (int k1 = 0; k1 < k; k1++) {
      arma::mat XtG = xnei.col(k1) * spanXnei % g;
      for (int p = 0; p < (poly + 1); p++) {
        arma::mat XtGX = XtG.col(p) * spanXtG % xnei;
        for (int k2 = 0; k2 < k; k2++) {
          int xindex = (k1 * (poly + 1) + p) * k + k2;
          Mx0(xindex, i) = sum(XtGX.col(k2));
        }
        int yindex = p * k + k1;
        arma::vec XtGY = XtG.col(p) % ynei;
        My0(yindex, i) = sum(XtGY);
      }
    }
  }
  return List::create(
    Named("Mx0") = Mx0,
    Named("My0") = My0
  );
}


// [[Rcpp::export]]
double scgwr_loocv(arma::vec target, arma::mat x, arma::vec y, int bw, int poly, arma::mat Mx0, arma::mat My0, arma::mat XtX, arma::mat XtY) {
  int n = x.n_rows, k = x.n_cols, poly1 = poly + 1;
  double b = target(0) * target(0), a = target(1) * target(1);
  arma::vec R0 = arma::vec(poly1) * b;
  for (int p = 1; p < poly1; p++) {
    R0(p) = pow(b, p + 1);
  }
  arma::vec Rx(k*k*poly1, fill::zeros), Ry(k*poly1, fill::zeros);
  for (int p = 0; p < poly1; p++) {
    for (int k2 = 0; k2 < k; k2++) {
      for (int k1 = 0; k1 < k; k1++) {
        int xindex = k1*poly1*k + p*k + k2;
        Rx(xindex) = R0(p);
      }
      int yindex = p*k + k2;
      Ry(yindex) = R0(p);
    }
  }
  arma::mat Mx = Rx * arma::mat(1, n, fill::ones) % Mx0, My = Ry * arma::mat(1, n, fill::ones) % My0;
  arma::vec yhat(n, 1, fill::zeros);
  for (int i = 0; i < n; i++) {
    arma::mat sumMx(k, k, fill::zeros);
    arma::vec sumMy(k, fill::zeros);
    for (int k2 = 0; k2 < k; k2++) {
      for (int p = 0; p < poly1; p++) {
        for (int k1 = 0; k1 < k; k1++) {
          int xindex = k1*poly1*k + p*k + k2;
          sumMx(k1, k2) += Mx(xindex, i);
        }
        int yindex = p*k + k2;
        sumMy(k2) += My(yindex, i);
      }
    }
    sumMx = sumMx + a * XtX;
    sumMy = sumMy + a * XtY;
    if (det(sumMx) < 1e-10) {
      return 1e6;
    } else {
      arma::mat beta = solve(sumMx, sumMy);
      yhat.row(i) = x.row(i) * beta;
    }
  }
  return sum((y - yhat) % (y - yhat));
}


// [[Rcpp::export]]
List scgwr_reg(arma::mat x, arma::vec y, int bw, int poly, arma::mat G0, arma::mat Mx0, arma::mat My0, arma::mat XtX, arma::mat XtY, arma::mat neighbour, arma::vec parameters) {
  int n = x.n_rows, k = x.n_cols, poly1 = poly + 1;
  double b = parameters(0), a = parameters(1);
  /**
   * Calculate Rx, Ry, and R0.
   */
  // printf("Calculate Rx, Ry, and R0 ");
  arma::vec R0 = arma::vec(poly1, fill::ones) * b;
  for (int p = 1; p < poly1; p++) {
    R0(p) = pow(b, p + 1);
  }
  arma::vec Rx(k*k*poly1, fill::zeros), Ry(k*poly1, fill::zeros);
  for (int p = 0; p < poly1; p++) {
    for (int k2 = 0; k2 < k; k2++) {
      for (int k1 = 0; k1 < k; k1++) {
        Rx(k1*poly1*k + p*k + k2) = R0(p);
      }
      Ry(p*k + k2) = R0(p);
    }
  }
  /**
   * Create G0.
   */
  // printf("Create G0 ");
  arma::mat G0s(G0.n_cols + 1, G0.n_rows, fill::ones);
  arma::mat G0t = trans(G0);
  G0s.rows(1, bw) = G0t.rows(0, bw - 1);
  G0s = trans(G0s);
  arma::mat G2(n, bw + 1, fill::zeros);
  for (int p = 0; p < poly; p++) {
    G2 += pow(G0s, pow(2.0, poly/2.0)/pow(2.0, p + 1)) * R0(poly - 1);
  }
  /**
   * Calculate Mx, My.
   */
  // printf("Calculate Mx, My ");
  // arma::mat Mx00(Mx0), My00(My0);
  for (int i = 0; i < n; i++) {
    for (int k1 = 0; k1 < k; k1++) {
      for (int p = 0; p < poly1; p++) {
        for (int k2 = 0; k2 < k; k2++) {
          Mx0((k1 * (poly + 1) + p) * k + k2, i) += x(i, k1) * x(i, k2);
        }
        My0(p * k + k1, i) += x(i, k1) * y(i);
      }
    }
  }
  arma::mat Mx = (Rx * arma::mat(1, n, fill::ones)) % Mx0, My = (Ry * arma::mat(1, n, fill::ones)) % My0;
  /**
   * Regression.
   */
  // printf("Regression ");
  arma::mat Xp(bw + 1, k * poly1, fill::zeros);
  arma::mat rowsumG(poly1, 1, fill::ones);
  arma::mat colsumXp(1, bw + 1, fill::zeros);
  arma::mat spanG(1, k, fill::ones);
  arma::mat spanX(1, poly1, fill::ones);
  arma::mat betas(n, k, fill::zeros);
  arma::mat betasSE(n, k, fill::zeros);
  double trS = 0.0, trStS = 0.0;
  for (int i = 0; i < n; i++) {
    /**
     * Calculate G.
     */
    arma::mat G = arma::mat(poly1, bw + 1, fill::ones) * R0(0);
    for (int p = 0; p < poly; p++) {
      G.row(p + 1) = pow(G0s.row(i), pow(2.0, poly/2.0)/pow(2.0, p + 1)) * R0(p);
    }
    G = trans(G);
    arma::mat g = G * rowsumG;
    /**
     * Calculate Xp.
     */
    arma::mat xnei(bw + 1, k, fill::zeros);
    arma::vec ynei(bw + 1, fill::zeros);
    xnei.row(0) = x.row(i);
    ynei.row(0) = y.row(i);
    for (int j = 0; j < bw; j++) {
      int inei = int(neighbour(i, j) - 1);
      xnei.row(j+1) = x.row(inei);
      ynei.row(j+1) = y(inei);
    }
    /**
     * Calculate sumMx, sumMy.
     */
    arma::mat sumMx(k, k, fill::zeros);
    arma::vec sumMy(k, fill::zeros);
    for (int k2 = 0; k2 < k; k2++) {
      for (int p = 0; p < poly1; p++) {
        for (int k1 = 0; k1 < k; k1++) {
          int xindex = k1*poly1*k + p*k + k2;
          sumMx(k1, k2) += Mx(xindex, i);
        }
        int yindex = p*k + k2;
        sumMy(k2) += My(yindex, i);
      }
    }
    sumMx += a * XtX;
    sumMy += a * XtY;
    arma::mat invMx = inv(sumMx);
    betas.row(i) = trans(invMx * sumMy);
    /**
     * Calculate Diagnostic statistics, trS and trStS.
     */
    arma::mat StS = invMx * trans(x.row(i));
    trS += det((x.row(i) * g(0, 0)) * StS);
    arma::mat XG2X(k, k, fill::zeros);
    for (int k1 = 0; k1 < k; k1++) {
      for (int k2 = 0; k2 < k; k2++) {
        arma::mat Gi = G2.row(i);
        XG2X(k1, k2) = sum(xnei.col(k1) % trans(Gi % Gi + 2 * a * Gi) % xnei.col(k2)) + a * a * XtX(k1, k2);
      }
    }
    arma::mat XX = invMx * XG2X * invMx;
    arma::mat xi = x.row(i);
    trStS += det(sum(xi * XX * trans(xi)));
    betasSE.row(i) = trans(sqrt(XX.diag()));
  }
  return List::create(
    Named("betas") = betas,
    Named("tr.S") = trS,
    Named("tr.StS") = trStS,
    Named("betas.SE") = betasSE
  );
}


// FE EDITS BELOW BY FIONA.H.EVANS@GMAIL.COM


// For debugging
void printVec(arma::vec v) {
  int n = v.size();
  n = 10;
  
  for (int i=0; i<n; i++) {
    Rprintf("%f ", v(i));
  }
  Rprintf("\n");
}
void printMat(arma::mat m) {
  int n = m.n_rows;
  int p = m.n_cols;
  
  n = 10;
  
  for (int i=0; i<n; i++) {
    for (int j=0; j<p; j++) {
      Rprintf("%f ", m(i, j));
    }
    Rprintf("\n");
  }
  Rprintf("\n");
}

// Boxcar weights 
// [[Rcpp::export]]
arma::vec box_wt_vec(arma::vec distv, double bw)
{
  int n = distv.n_elem;
  arma::vec wtv(n, fill::zeros);
  
  arma::uvec u = arma::find(distv <= bw);
  wtv.elem(u).fill(1);
  
  return wtv;
}

// Boxcar adaptive weights
// [[Rcpp::export]]
arma::vec box_wt_vec_ad(arma::vec distv, double bw)
{
  int n = distv.n_elem;
  arma::vec wtv;
  wtv.zeros(n);
  double bwd = bw;
  if (bw >= n) bwd = n;
  
  // equivalent to R function rank(distv, ties.method='first')
  arma::uvec rnk1 = arma::sort_index(distv) + 1;
  arma::uvec rnk = arma::sort_index(rnk1) + 1;
  
  arma::uvec u = arma::find(rnk <= bwd);
  wtv.elem(u).fill(1);
  
  return wtv;
}

// Gaussian adaptive weights
// [[Rcpp::export]]
arma::vec gau_wt_vec_ad(arma::vec distv, double bw)
{
  int n = distv.n_elem;
  double bwd = bw/n * distv.max();
  
  if (bw <= n){
    // equivalent to R function rank(distv, ties.method='first')
    arma::uvec rnk1 = arma::sort_index(distv) + 1;
    arma::uvec rnk = arma::sort_index(rnk1) + 1;
    
    arma::uvec u = arma::find(rnk == bw);
    bwd = distv(u(0));
  }
  
  arma::vec wtv = exp(pow(distv, 2) / ((-2) * pow(bwd, 2)));
  
  return wtv;
}

// Bisquare adaptive weights
// [[Rcpp::export]]
arma::vec bis_wt_vec_ad(arma::vec distv, double bw)
{
  int n = distv.n_elem;
  double bwd = bw/n * distv.max();
  
  if (bw <= n){
    // equivalent to R function rank(distv, ties.method='first')
    arma::uvec rnk1 = arma::sort_index(distv) + 1;
    arma::uvec rnk = arma::sort_index(rnk1) + 1;
    
    arma::uvec u = arma::find(rnk == bw);
    bwd = distv(u(0));
  }
  
  arma::vec wtv = bisq_wt_vec(distv, bwd);
  
  return wtv;
}

// Tricube adaptive weight
// [[Rcpp::export]]
arma::vec tri_wt_vec_ad(arma::vec distv, double bw)
{
  int n = distv.n_elem;
  double bwd = bw/n * distv.max();
  
  if (bw <= n){
    // equivalent to R function rank(distv, ties.method='first')
    arma::uvec rnk1 = arma::sort_index(distv) + 1;
    arma::uvec rnk = arma::sort_index(rnk1) + 1;
    
    arma::uvec u = arma::find(rnk == bw);
    bwd = distv(u(0));
  }
  
  arma::vec wtv = tri_wt_vec(distv, bwd);
  
  return wtv;
}

// Exponential adaptive weights
// [[Rcpp::export]]
arma::vec exp_wt_vec_ad(arma::vec distv, double bw)
{
  int n = distv.n_elem;
  double bwd = bw/n * distv.max();
  
  if (bw <= n){
    // equivalent to R function rank(distv, ties.method='first')
    arma::uvec rnk1 = arma::sort_index(distv) + 1;
    arma::uvec rnk = arma::sort_index(rnk1) + 1;
    
    arma::uvec u = arma::find(rnk == bw);
    bwd = distv(u(0));
  }
  
  arma::vec wtv = exp_wt_vec(distv, bwd);
  
  return wtv;
}

// For use in gw_weight
enum string_code{ga, bi, tr, bo, ex};
string_code hashit (std::string const& inString) {
  if (inString == "gaussian") return ga;
  if (inString == "bisquare") return bi;
  if (inString == "tricube") return tr;
  if (inString == "boxcar") return bo;
  if (inString == "exponential") return ex;
}

// Geographic weights (gw.weight equivalent)
// [[Rcpp::export]]
arma::vec gw_weight(arma::vec vdist, double bw, std::string kernel, bool adaptive)
{
  if (adaptive) switch(hashit(kernel)){
  case ga: return gau_wt_vec_ad (vdist, bw);
  case bi: return bis_wt_vec_ad (vdist, bw);
  case tr: return tri_wt_vec_ad (vdist, bw);
  case bo: return box_wt_vec_ad (vdist, bw); 
  case ex: return exp_wt_vec_ad (vdist, bw);
  }
  else switch(hashit(kernel)){
  case ga: return gauss_wt_vec (vdist, bw);
  case bi: return bisq_wt_vec (vdist, bw);
  case tr: return tri_wt_vec (vdist, bw);
  case bo: return box_wt_vec (vdist, bw);
  case ex: return exp_wt_vec (vdist, bw);
  }
}

// GWR calibration, returns betas only
// [[Rcpp::export]]
arma::vec gw_reg_2(arma::mat x, arma::vec y, arma::vec w)
{
  arma::mat wspan(1, x.n_cols, fill::ones);
  arma::mat xtw = trans(x % (w * wspan));
  arma::mat xtwx = xtw * x;
  arma::mat xtwy = trans(x) * (w % y);
  arma::mat xtwx_inv = inv(xtwx);
  arma::vec beta = xtwx_inv * xtwy;
  
  return beta;
}

// C++ version of gwr.q, used in gwr.mixed
// [[Rcpp::export]]
arma::mat gwr_q(arma::mat x, arma::vec y, 
                arma::mat dMat, double bw, std::string kernel, bool adaptive)
{
  // loc and out.loc only used to create distances
  int n =  dMat.n_cols;  // loc.n_rows() 
  int m =  x.n_cols;
  arma::mat beta(n, m);
  arma::vec distv;
  arma::vec w;
  
  for (int i = 0; i < n; i++) {
    distv = dMat.col(i);
    w = gw_weight(distv, bw, kernel, adaptive);
    beta.row(i) = gw_reg_2(x, y, w);
  }
  
  return beta;
}

// [[Rcpp::export]]
arma::vec e_vec(int m, int n){
  arma::vec e = linspace(0, n-1, n);
  arma::vec ret(n, fill::zeros);
  arma::uvec u = arma::find(e == m);
  ret.elem(u).fill(1);

  return ret;
}

// [[Rcpp::export]]
double gwr_mixed_trace(arma::mat x1, arma::mat x2, arma::vec y, 
                       arma::mat dMat, double bw, std::string kernel, bool adaptive){
  int i;
  int n = x1.n_rows;
  int nc2 = x2.n_cols;
  arma::mat mtemp, model1, model2;
  arma::mat x3(n, nc2);
  arma::vec y2(n);
  arma::vec hii(n, fill::zeros);
  arma::mat m(1,n);
  double s1, s2;
  
  for (i = 0; i < nc2; i++) {
    mtemp = gwr_q(x1, x2.col(i), dMat, bw, kernel, adaptive);
    x3.col(i) = x2.col(i) - fitted(x1, mtemp);
  }
  
  for (i = 0; i < n; i++) {
    mtemp = gwr_q(x1, e_vec(i, n), dMat, bw, kernel, adaptive);
    y2 = e_vec(i, n) - fitted(x1, mtemp);
    model2 = gwr_q(x3, y2, dMat, 100000, "boxcar", true);
    m.row(0) = dMat.row(i);
    model1 = gwr_q(x1, e_vec(i, n) - fitted(x2, model2), dMat.col(i), bw, kernel, adaptive);
    model2 = gwr_q(x3, y2, dMat.col(i), 100000, "boxcar", true);
    s1 = fitted(x1.row(i), model1)(0);
    s2 = fitted(x2.row(i), model2)(0);
    hii(i) = s1 + s2;
  }
  
  return sum(hii);
}

// [[Rcpp::export]]
List gwr_mixed_2(arma::mat x1, arma::mat x2, arma::vec y, 
                       arma::mat dMat, arma::mat dMat_rp,
                       double bw, std::string kernel, bool adaptive){
  int i;
  int n = x1.n_rows;
  int nc2 = x2.n_cols;
  arma::mat mtemp, model1, model2;
  arma::mat x3(n, nc2);
  arma::vec y2(n);
  arma::vec hii(n, fill::zeros);
  arma::mat m(1,n);
  double s1, s2;
  
  for (i = 0; i < nc2; i++) {
    mtemp = gwr_q(x1, x2.col(i), dMat, bw, kernel, adaptive);
    x3.col(i) = x2.col(i) - fitted(x1, mtemp);
  }
  
  mtemp = gwr_q(x1, y, dMat, bw, kernel, adaptive);
  y2 = y - fitted(x1, mtemp);
  model2 = gwr_q(x3, y2, dMat, 100000, "boxcar", true);
  
  model1 = gwr_q(x1, y-fitted(x2, model2), dMat_rp, bw, kernel, adaptive);
  model2 = gwr_q(x3, y2, dMat_rp, 100000, "boxcar", true);
  
  return List::create(
    Named("local") = model1,
    Named("global") = model2
  );
}