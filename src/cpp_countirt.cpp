#include <Rcpp.h>
#include <math.h>
#include <Rmath.h>
#include <RcppGSL.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
using namespace Rcpp;
using namespace R;

// [[Rcpp::depends(RcppGSL)]]

// [[Rcpp::export]]
double logFactorial( int n )
{
  if( n <= 0 )  return 0 ;
  double res = log(n) ;
  while( --n>1 ) res += log(n) ;
  return res ;
}
// computes log(n!) insted of n!

// global iter_max
int iter_max = 1e4;

// [[Rcpp::export]]
double computeA(double lambda, double mu, double nu, double log_Z, double min_iter) {

  // double Q = 0.0;

  // Series summation
  int index;              // Current Index
  double dA;           // summation increment of Q;
  double reltol = 1e-12;  // Break if Term_current / Sum_Current < reltol
  // Initialize largest term and sum
  int index_mode = floor(mu);
  double A =  (double) logFactorial(index_mode) * (index_mode - mu) *
    exp(index_mode*log(lambda) - log_Z - nu*lgamma(index_mode+1));
  // berechne nur die Teile mit Exponenten auf log Ebene und exponentiere dann fuer mehr
  // numerische Stabilitaet; der Rest kann nicht erst auf log-Ebene berechnet werden

  // Left tail
  for (int j = 1; j < iter_max; j++) {
    if (j == iter_max - 1) {
      warning("max iter. reached");
      break;
    }
    index = index_mode - j;
    if (index < 0) break;
    dA = (double) logFactorial(index) * (index - mu) *
      exp(index*log(lambda) - log_Z - nu*lgamma(index+1));
    A += dA;
    if ((abs(dA) < reltol) & (j > min_iter)) break;
  }

  // Right tail
  for (int j = 1; j < iter_max; j++) {
    if (j == iter_max - 1) {
      warning("max iter. reached");
      break;
    }
    index = index_mode + j;
    dA = (double) logFactorial(index) * (index - mu) *
      exp(index*log(lambda) - log_Z - nu*lgamma(index+1));
    A += dA;
    if ((abs(dA) < reltol) & (j > min_iter)) break;
  }

  return A;
}

// [[Rcpp::export]]
double computeB(double lambda, double mu, double nu, double log_Z, double min_iter) {

  // double Q = 0.0;

  // Series summation
  int index;              // Current Index
  double dB;           // summation increment of Q;
  double reltol = 1e-12;  // Break if Term_current / Sum_Current < reltol
  // Initialize largest term and sum
  int index_mode = floor(mu);
  double B =  (double) logFactorial(index_mode) *
    exp(index_mode*log(lambda) - log_Z - nu*lgamma(index_mode+1));
  // berechne nur die Teile mit Exponenten auf log Ebene und exponentiere dann fuer mehr
  // numerische Stabilitaet; der Rest kann nicht erst auf log-Ebene berechnet werden

  // Left tail
  for (int j = 1; j < iter_max; j++) {
    if (j == iter_max - 1) {
      warning("max iter. reached");
      break;
    }
    index = index_mode - j;
    if (index < 0) break;
    dB = (double) logFactorial(index) *
      exp(index*log(lambda) - log_Z - nu*lgamma(index+1));
    B += dB;
    if ((abs(dB) < reltol) & (j > min_iter)) break;
  }

  // Right tail
  for (int j = 1; j < iter_max; j++) {
    if (j == iter_max - 1) {
      warning("max iter. reached");
      break;
    }
    index = index_mode + j;
    dB = (double) logFactorial(index) *
      exp(index*log(lambda) - log_Z - nu*lgamma(index+1));
    B += dB;
    if ((abs(dB) < reltol) & (j > min_iter)) break;
  }

  return B;
}

// [[Rcpp::export]]
double computeQ(double lambda, double mu, double nu, double min_iter) {

  // double Q = 0.0;
  double loglambda = log(lambda);

  // Series summation
  int index;              // Current Index
  double dQ;           // summation increment of Q;
  double reltol = 1e-12;  // Break if Term_current / Sum_Current < reltol
  // Initialize largest term and sum
  int index_mode = floor(mu);
  double Q = (index_mode - mu) * logFactorial(index_mode) *
    exp(loglambda * (double) index_mode - (double) logFactorial(index_mode) * nu);
  // berechne nur die Teile mit Exponenten auf log Ebene und exponentiere dann fuer mehr
  // numerische Stabilitaet; der Rest kann nicht erst auf log-Ebene berechnet werden

  // Left tail
  for (int j = 1; j < iter_max; j++) {
    if (j == iter_max - 1) {
      warning("max iter. reached");
      break;
    }
    index = index_mode - j;
    if (index < 0) break;
    dQ = (index - mu) * logFactorial(index) *
      exp(loglambda * index - logFactorial(index) * nu);
    Q += dQ;
    if ((abs(dQ) < reltol) & (j > min_iter)) break;
  }

  // Right tail
  for (int j = 1; j < iter_max; j++) {
    if (j == iter_max - 1) {
      warning("max iter. reached");
      break;
    }
    index = index_mode + j;
    dQ = (index - mu) * logFactorial(index) *
      exp(loglambda * index - logFactorial(index) * nu);
    Q += dQ;
    if ((abs(dQ) < reltol) & (j > min_iter)) break;
  }

  return Q;
}

// [[Rcpp::export]]
double computeW(double lambda, double mu, double nu, double min_iter) {

  // double W = 0.0;
  double loglambda = log(lambda);

  // Series summation
  int index;              // Current Index
  double dW;           // summation increment of W;
  double reltol = 1e-12;  // Break if Term_current / Sum_Current < reltol
  // Initialize largest term and sum
  int index_mode = floor(mu);
  double Q = computeQ(lambda, mu, nu, min_iter);
  double log_Q = log(Q);
  double W = pow(index_mode - mu, 2) *
    exp(loglambda * index_mode - logFactorial(index_mode) * nu - log_Q);
  // berechne nur die Teile mit Exponenten auf log Ebene und exponentiere dann fuer mehr
  // numerische Stabilitaet; der Rest kann nicht erst auf log-Ebene berechnet werden

  // Left tail
  for (int j = 1; j < iter_max; j++) {
    if (j == iter_max - 1) {
      warning("max iter. reached");
      break;
    }
    index = index_mode - j;
    if (index < 0) break;
    dW = pow(index - mu, 2) * exp(loglambda * index - logFactorial(index) * nu - log_Q);
    W += dW;
    if ((abs(dW) < reltol) & (j > min_iter)) break;
  }

  // Right tail
  for (int j = 1; j < iter_max; j++) {
    if (j == iter_max - 1) {
      warning("max iter. reached");
      break;
    }
    index = index_mode + j;
    dW = pow(index - mu, 2) * exp(loglambda * index - logFactorial(index) * nu - log_Q);
    W += dW;
    if ((abs(dW) < reltol) & (j > min_iter)) break;
  }

  return W;
  // Q als Konstante kann man auch vor die Reihe ziehen und dann erst die Reihe berechnen
  // und dann erst die Reihe berechnen und die fertige Summe einmal durch Q teilen
}

// [[Rcpp::export]]
double computeR(double lambda, double mu, double nu,
                double log_Z, double min_iter) {

  // double W = 0.0;
  double loglambda = log(lambda);

  // Series summation
  int index;              // Current Index
  double dR;           // summation increment of W;
  double reltol = 1e-12;  // Break if Term_current / Sum_Current < reltol
  // Initialize largest term and sum
  int index_mode = floor(mu);
  double W = computeW(lambda, mu, nu, min_iter);
  double R = exp(loglambda * index_mode - logFactorial(index_mode) * nu - log_Z) *
    (index_mode / W - logFactorial(index_mode));
  // berechne nur die Teile mit Exponenten auf log Ebene und exponentiere dann fuer mehr
  // numerische Stabilitaet; der Rest kann nicht erst auf log-Ebene berechnet werden

  // Left tail
  for (int j = 1; j < iter_max; j++) {
    if (j == iter_max - 1) {
      warning("max iter. reached");
      break;
    }
    index = index_mode - j;
    if (index < 0) break;
    dR = exp(loglambda * index - logFactorial(index) * nu - log_Z) *
      (index / W - logFactorial(index));
    R += dR;
    if ((abs(dR) < reltol) & (j > min_iter)) break;
  }

  // Right tail
  for (int j = 1; j < iter_max; j++) {
    if (j == iter_max - 1) {
      warning("max iter. reached");
      break;
    }
    index = index_mode + j;
    dR = exp(loglambda * index - logFactorial(index) * nu - log_Z) *
      (index / W - logFactorial(index));
    R += dR;
    if ((abs(dR) < reltol) & (j > min_iter)) break;
  }

  return R;
}


// [[Rcpp::export]]
NumericVector computepp(NumericVector dens, double quad_weight, NumericVector marg_prob){
  // computes posterior probability at one node and for one item for each persons, i.e., returns
  // a vector of length n_persons with posterior probabilities for one node
  // dens = densities for response vector to one item j (of all N participants)
  // quad_weight = Quadrature weight for node k
  // marg_prob = marginal probability for data for each participant, so is a vector of length N
  // (is marginalized over nodes, so not node-specific)

  int N = dens.size();
  NumericVector pp(N);

  for(int i=0;i<N;i++){
    // in C++ indexing starts at 0, no of elements is length of array minus 1
    pp[i] = (dens[i] * quad_weight) / marg_prob[i];
  }

  return pp;
}

// [[Rcpp::export]]
NumericMatrix computepp_allnodes(NumericMatrix dens, NumericVector quad_weight, NumericVector marg_prob){
  // in this all nodes version of computepp, we output a matrix of posterior probabilities
  // where we do one column per node and each column per node contains the output from computepp
  // it thus also need dens and marg_prob in matrix format where each node gets one column
  // that contains entries for all N persons; in general, the matrices need to be of the shape
  // persons (rows) + nodes (columns)
  // dens = densities for response vector to one item j (of all N participants); matrix: persons (rows) + nodes (columns)
  // quad_weight = all quadrature weights (in a vector); this vector should have the same length as ncol(dens)
  // marg_prob = marginal probability for data for each participant, so is a vector of length N
  // (is marginalized over nodes, so not node-specific)

  int N = dens.nrow();
  int K = dens.ncol();
  NumericMatrix pp_an(N, K);

  for(int j=0;j<K;j++){
    // loop through all columns and use computepp in each column
    pp_an( _ , j) = computepp(dens( _ , j), quad_weight[j], marg_prob);
    // we don't need an index over marg_prob because it is not node-specific
  }

  return pp_an;
}

// [[Rcpp::export]]
double computer(NumericVector resp, NumericVector pp) {
  // computes sufficient statistic r at one node and for one item
  // resp = a vector of length n_persons with responses of all participants to one item
  // pp = a vector of length n_persons with post. probs of all participants to one item at one node

  double r = 0;
  int N = resp.size();

  for(int i=0;i<N;i++){
    // in C++ indexing starts at 0, no of elements is length of array minus 1
    r = r + resp[i] * pp[i];
  }

  return r;
}

// [[Rcpp::export]]
NumericVector computer_allnodes(NumericVector resp, NumericMatrix pp){
  // computes sufficient statistic r for all nodes and for one item
  // outputs a vector with r values (for one item) for all nodes, has length no. of nodes
  // resp = a vector of length n_persons with responses of all participants to one item
  // pp = matrix as returned by computepp_allnodes of dimensions persons (rows) + nodes (columns)

  int K = pp.ncol();
  NumericVector r(K);

  for(int j=0;j<K;j++){
    // we loop over the nodes and produce a vector as the output which is as long as no. of nodes
    r[j] = computer(resp, pp( _ , j));
  }

  return r;
}

// [[Rcpp::export]]
double computef(NumericVector pp){
  // computes sufficient statistic f at one node and for one item
  // pp = a vector of length n_persons with post. probs of all participants to one item at one node

  double f = 0;
  int N = pp.size();

  for(int i=0;i<N;i++){
    // in C++ indexing starts at 0, no of elements is length of array minus 1
    f = f + pp[i];
  }

  return f;
}

// [[Rcpp::export]]
NumericVector computef_allnodes(NumericMatrix pp){
  // computes sufficient statistic f for all nodes and for one item
  // pp = matrix as returned by computepp_allnodes of dimensions persons (rows) + nodes (columns)

  int K = pp.ncol();
  NumericVector f(K);

  for(int j=0;j<K;j++){
    // we loop over the nodes and produce a vector as the output which is as long as no. of nodes
    f[j] = computef(pp( _ , j));
  }

  return f;
}


// [[Rcpp::export]]
double computeh(NumericVector resp, NumericVector pp){
  // computes sufficient statistic h at one node and for one item
  // resp = a vector of length n_persons with responses of all participants to one item
  // pp = a vector of length n_persons with post. probs of all participants to one item at one node

  double h = 0;
  int N = pp.size();

  for(int i=0;i<N;i++){
    // in C++ indexing starts at 0, no of elements is length of array minus 1
    h = h + logFactorial(resp[i]) * pp[i];
  }

  return h;
}

// [[Rcpp::export]]
NumericVector computeh_allnodes(NumericVector resp, NumericMatrix pp){
  // computes sufficient statistic h for all nodes and for one item
  // resp = a vector of length n_persons with responses of all participants to one item
  // pp = matrix as returned by computepp_allnodes of dimensions persons (rows) + nodes (columns)

  int K = pp.ncol();
  NumericVector h(K);

  for(int j=0;j<K;j++){
    // we loop over the nodes and produce a vector as the output which is as long as no. of nodes
    h[j] = computeh(resp, pp( _ , j));
  }

  return h;
}


// [[Rcpp::export]]
double interp_from_grid(NumericVector mu, NumericVector nu,
                   NumericVector grid_long, const double mu0, const double nu0){

  // get input into right format for interpolation function
  int K = mu.size();
  double mu_d[K];
  for(int i=0;i<K;i++){
    mu_d[i] = mu[i];
  }

  int L = nu.size();
  double nu_d[L];
  for(int i=0;i<L;i++){
    nu_d[i] = nu[i];
  }

  int M = grid_long.size();
  double grid_long_d[M];
  for(int i=0;i<M;i++){
    grid_long_d[i] = grid_long[i];
  }

  const gsl_interp2d_type *T = gsl_interp2d_bicubic;
  const size_t nx = sizeof(mu_d) / sizeof(double);
  const size_t ny = sizeof(nu_d) / sizeof(double);
  //gsl_spline2d *spline = gsl_spline2d_alloc(T, nx, ny);
  gsl_interp2d *interp = gsl_interp2d_alloc(T, nx, ny);
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();

  gsl_interp2d_init(interp, mu_d, nu_d, grid_long_d, nx, ny);
  //gsl_spline2d_init(spline, mu_d, nu_d, grid_long_d, nx, ny);

  double out = gsl_interp2d_eval_extrap(interp, mu_d, nu_d, grid_long_d,
                                        mu0, nu0, xacc, yacc);
  //double out = gsl_spline2d_eval_extrap(spline, mu0, nu0, xacc, yacc);

  gsl_interp2d_free (interp);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);

  return out;
}

// [[Rcpp::export]]
NumericVector interp_from_grid_v(NumericVector mu, NumericVector nu,
                                  NumericVector grid_long,
                                  NumericVector mu0_v,
                                  NumericVector nu0_v){
  // mu0_v and nu0_v must have the same length

  int m = mu0_v.size();
  NumericVector out(m);

  // get input into right format for interpolation function
  int K = mu.size();
  double mu_d[K];
  for(int i=0;i<K;i++){
    mu_d[i] = mu[i];
  }

  int L = nu.size();
  double nu_d[L];
  for(int i=0;i<L;i++){
    nu_d[i] = nu[i];
  }

  int M = grid_long.size();
  double grid_long_d[M];
  for(int i=0;i<M;i++){
    grid_long_d[i] = grid_long[i];
  }

  const gsl_interp2d_type *T = gsl_interp2d_bicubic;
  const size_t nx = sizeof(mu_d) / sizeof(double);
  const size_t ny = sizeof(nu_d) / sizeof(double);
  //gsl_spline2d *spline = gsl_spline2d_alloc(T, nx, ny);
  gsl_interp2d *interp = gsl_interp2d_alloc(T, nx, ny);
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();

  gsl_interp2d_init(interp, mu_d, nu_d, grid_long_d, nx, ny);
  //gsl_spline2d_init(spline, mu_d, nu_d, grid_long_d, nx, ny);

  for(int i=0;i<m;i++){
    double mu0 = mu0_v[i];
    double nu0 = nu0_v[i];
    out[i] = gsl_interp2d_eval_extrap(interp, mu_d, nu_d, grid_long_d,
                                      mu0, nu0, xacc, yacc);
    //out[i] = gsl_spline2d_eval(spline, mu0, nu0, xacc, yacc);
  }

  gsl_interp2d_free (interp);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);

  return(out);
}

// [[Rcpp::export]]
NumericMatrix interp_from_grid_m(NumericVector mu, NumericVector nu,
                                 NumericVector grid_long,
                                 NumericMatrix mu0_m,
                                 NumericMatrix nu0_m){
  // mu0_m and nu0_m must have the same dimensions
  // they contain the values we are trying to interpolate

  int n = mu0_m.nrow(); // no. of rows
  int m = mu0_m.ncol(); // no. of columns
  NumericMatrix out(n, m);

  // get input into right format for interpolation function
  int K = mu.size();
  double mu_d[K];
  for(int i=0;i<K;i++){
    mu_d[i] = mu[i];
  }

  int L = nu.size();
  double nu_d[L];
  for(int i=0;i<L;i++){
    nu_d[i] = nu[i];
  }

  int M = grid_long.size();
  double grid_long_d[M];
  for(int i=0;i<M;i++){
    grid_long_d[i] = grid_long[i];
  }

  const gsl_interp2d_type *T = gsl_interp2d_bicubic;
  const size_t nx = sizeof(mu_d) / sizeof(double);
  const size_t ny = sizeof(nu_d) / sizeof(double);
  //gsl_spline2d *spline = gsl_spline2d_alloc(T, nx, ny);
  gsl_interp2d *interp = gsl_interp2d_alloc(T, nx, ny);
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();

  gsl_interp2d_init(interp, mu_d, nu_d, grid_long_d, nx, ny);
  //gsl_spline2d_init(spline, mu_d, nu_d, grid_long_d, nx, ny);

  for(int i=0;i<n;i++){
    // over rows
    for(int j=0;j<m;j++){
      // over columns
    double mu0 = mu0_m(i,j);
    double nu0 = nu0_m(i,j);
    out(i,j) = gsl_interp2d_eval_extrap(interp, mu_d, nu_d, grid_long_d,
        mu0, nu0, xacc, yacc);
    //out[i] = gsl_spline2d_eval(spline, mu0, nu0, xacc, yacc);
    }
  }

  gsl_interp2d_free (interp);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);

  return(out);
}

// [[Rcpp::export]]
NumericMatrix interp_from_grid_lin_m(NumericVector mu, NumericVector nu,
                                     NumericVector grid_long,
                                     NumericMatrix mu0_m,
                                     NumericMatrix nu0_m){
  // mu0_m and nu0_m must have the same dimensions
  // they contain the values we are trying to interpolate

  int n = mu0_m.nrow(); // no. of rows
  int m = mu0_m.ncol(); // no. of columns
  NumericMatrix out(n, m);

  // get input into right format for interpolation function
  int K = mu.size();
  double mu_d[K];
  for(int i=0;i<K;i++){
    mu_d[i] = mu[i];
  }

  int L = nu.size();
  double nu_d[L];
  for(int i=0;i<L;i++){
    nu_d[i] = nu[i];
  }

  int M = grid_long.size();
  double grid_long_d[M];
  for(int i=0;i<M;i++){
    grid_long_d[i] = grid_long[i];
  }

  const gsl_interp2d_type *T = gsl_interp2d_bilinear;
  const size_t nx = sizeof(mu_d) / sizeof(double);
  const size_t ny = sizeof(nu_d) / sizeof(double);
  //gsl_spline2d *spline = gsl_spline2d_alloc(T, nx, ny);
  gsl_interp2d *interp = gsl_interp2d_alloc(T, nx, ny);
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();

  gsl_interp2d_init(interp, mu_d, nu_d, grid_long_d, nx, ny);
  //gsl_spline2d_init(spline, mu_d, nu_d, grid_long_d, nx, ny);

  for(int i=0;i<n;i++){
    // over rows
    for(int j=0;j<m;j++){
      // over columns
      double mu0 = mu0_m(i,j);
      double nu0 = nu0_m(i,j);
      out(i,j) = gsl_interp2d_eval_extrap(interp, mu_d, nu_d, grid_long_d,
          mu0, nu0, xacc, yacc);
      //out[i] = gsl_spline2d_eval(spline, mu0, nu0, xacc, yacc);
    }
  }

  gsl_interp2d_free (interp);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);

  return(out);
}

// [[Rcpp::export]]
NumericVector dcmp_cpp(NumericVector data, NumericVector mu, NumericVector nu,
                       bool logprob,
                       NumericVector all_mus_lambda,
                       NumericVector all_nus_lambda,
                       NumericVector all_mus_logZ,
                       NumericVector all_nus_logZ,
                       NumericVector grid_long_log_lambda,
                       NumericVector grid_long_logZ){
  // length data = length mu = length nu (all checks in R)
  // handle cases (e.g., same nu for all mus, same data for all nus) in R
  // for numerical stability we use log lambda table

  int m = mu.size();
  NumericVector out(m);
  NumericVector log_lambdas = interp_from_grid_v(all_mus_lambda, all_nus_lambda,
                                                 grid_long_log_lambda,
                                                 mu, nu);
  NumericVector logZ = interp_from_grid_v(all_mus_logZ, all_nus_logZ,
                                          grid_long_logZ,
                                          mu, nu);

  for(int i=0;i<m;i++){
    out[i] = data[i]*log_lambdas[i] - logZ[i] - nu[i]*lgamma(data[i]+1);
  }

  if (logprob) {
    return(out);
  }
  else {
    return(exp(out));
  }

}


// [[Rcpp::export]]
NumericVector grad_cmp_cpp(NumericVector alphas,
                           NumericVector deltas,
                           NumericVector disps,
                           NumericVector nodes,
                           NumericVector weights,
                           NumericMatrix r,
                           NumericMatrix f,
                           NumericMatrix h,
                           NumericVector grid_mus,
                           NumericVector grid_nus,
                           NumericVector grid_cmp_var_long,
                           NumericVector grid_log_lambda_long,
                           NumericVector grid_logZ_long,
                           NumericVector grid_W_long,
                           NumericVector grid_R_long,
                           double max_mu,
                           double min_mu) {

  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h

  int m = alphas.size();
  int n_nodes = nodes.size();
  NumericVector grad_alphas(m);
  NumericVector grad_deltas(m);
  NumericVector grad_disps(m);
  NumericVector out(3*m);

  // set up mu's and nu's for interpolation function to be computed all in one

  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int i=0;i<m;i++){
    // loop over items (columns)
    for(int j=0;j<n_nodes;j++) {
      // loop over nodes (rows)
      mu(j,i) = exp(alphas[i] * nodes[j] + deltas[i]);
      mu_interp(j,i) = mu(j,i);
      if (mu(j,i) > max_mu) { mu_interp(j,i) = max_mu; }
      if (mu(j,i) < min_mu) { mu_interp(j,i) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(j,i) = disps[i];
    }
  }

  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  // NEW
  NumericMatrix W(n_nodes, m);
  NumericMatrix R(n_nodes, m);
  // NEW
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // NEW
  W = interp_from_grid_m(grid_mus, grid_nus,
                         grid_W_long,
                         mu_interp, disp_interp);
  R = interp_from_grid_m(grid_mus, grid_nus,
                         grid_R_long,
                         mu_interp, disp_interp);
  // NEW
  // V and log_lambda are matrices with as many columns as we have nodes and
  // as many columns as we have items


  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    grad_alphas[i] = 0;
    grad_deltas[i] = 0;
    grad_disps[i] = 0;

    for(int j=0;j<n_nodes;j++) {
      // over nodes (rows in my matrices)
      //double term = 0;
      //if ((r(j,i) > 1e-8) | (f(j,i) > 1e-8) | (h(j,i) > 1e-8)) {
      //if (weights[j] > 1e-4) {
        // NEW
        //double lambda = exp(log_lambda(j,i));
        //double W = computeW(lambda, mu_interp(j,i), disps[i], 10);
        //double R = computeR(lambda, mu_interp(j,i), disps[i], log_Z(j,i), 10);
        // NEW
        double frac_r_W = r(j,i) / W(j,i);
        double f_R = f(j,i)*R(j,i);
        double term = frac_r_W - h(j,i) - f_R;
      //}
      //double frac_r_W = 0;
      //if (r(j,i) > 1e-8) {
      //  double frac_r_W = r(j,i) / W;
      //}
      //double frac_muf_W = 1e-8;
      //if (mu(j,i)*f(j,i) > 1e-8) {
      //  frac_muf_W = mu(j,i)*f(j,i) / W;
      //}
      //double e_logfac = 1e-8;
      //if (f(j,i) > 1e-8) {
      //  double lambda = exp(log_lambda(j,i));
      //  e_logfac = tmbElogFactorial2(lambda, disps[i]);
      //}
      //double f_R = 0;
      //if (f(j,i) > 1e-8) {
      //  double f_R = f(j,i)*R;
      //}

      grad_alphas[i] = grad_alphas[i] +
        (nodes[j]*mu_interp(j,i) / V(j,i))*(r(j,i) - mu_interp(j,i)*f(j,i));
      grad_deltas[i] = grad_deltas[i] +
        (mu_interp(j,i) / V(j,i))*(r(j,i) - mu_interp(j,i)*f(j,i));
      //grad_disps[i] = grad_disps[i] +
      //  disps[i]*(frac_r_W - frac_muf_W - h(j,i) + e_logfac*f(j,i));
      grad_disps[i] = grad_disps[i] + disps[i]*term;
    }
  }

  // fill up output vector
  for(int i=0;i<m;i++){
    out[i] = grad_alphas[i];
    out[i + m] = grad_deltas[i];
    out[i + 2*m] = grad_disps[i];
  }

  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_fixdisps_cpp(NumericVector alphas,
                                    NumericVector deltas,
                                    NumericVector disps,
                                    NumericVector nodes,
                                    NumericMatrix r,
                                    NumericMatrix f,
                                    NumericVector grid_mus,
                                    NumericVector grid_nus,
                                    NumericVector grid_cmp_var_long,
                                    NumericVector grid_log_lambda_long,
                                    double max_mu) {

  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h

  int m = alphas.size();
  int n_nodes = nodes.size();
  NumericVector grad_alphas(m);
  NumericVector grad_deltas(m);
  NumericVector out(2*m);

  // set up mu's and nu's for interpolation function to be computed all in one

  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int i=0;i<m;i++){
    // loop over items (columns)
    for(int j=0;j<n_nodes;j++) {
      // loop over nodes (rows)
      mu(j,i) = exp(alphas[i] * nodes[j] + deltas[i]);
      mu_interp(j,i) = mu(j,i);
      if (mu(j,i) > max_mu) { mu_interp(j,i) = max_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(j,i) = disps[i];
    }
  }

  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many columns as we have nodes and
  // as many columns as we have items


  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    grad_alphas[i] = 0;
    grad_deltas[i] = 0;

    for(int j=0;j<n_nodes;j++) {
      // over nodes (rows in my matrices)
      grad_alphas[i] = grad_alphas[i] +
        (nodes[j]*mu_interp(j,i) / V(j,i))*(r(j,i) - mu_interp(j,i)*f(j,i));
      grad_deltas[i] = grad_deltas[i] +
        (mu_interp(j,i) / V(j,i))*(r(j,i) - mu_interp(j,i)*f(j,i));
    }
  }

  // fill up output vector
  for(int i=0;i<m;i++){
    out[i] = grad_alphas[i];
    out[i + m] = grad_deltas[i];
  }

  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_fixalphas_cpp(NumericVector alphas,
                                    NumericVector deltas,
                                    NumericVector disps,
                                    NumericVector nodes,
                                    NumericMatrix r,
                                    NumericMatrix f,
                                    NumericMatrix h,
                                    NumericVector grid_mus,
                                    NumericVector grid_nus,
                                    NumericVector grid_cmp_var_long,
                                    NumericVector grid_log_lambda_long,
                                    NumericVector grid_logZ_long,
                                    double max_mu) {

  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h

  int m = alphas.size();
  int n_nodes = nodes.size();
  NumericVector grad_disps(m);
  NumericVector grad_deltas(m);
  NumericVector out(2*m);

  // set up mu's and nu's for interpolation function to be computed all in one

  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int i=0;i<m;i++){
    // loop over items (columns)
    for(int j=0;j<n_nodes;j++) {
      // loop over nodes (rows)
      mu(j,i) = exp(alphas[i] * nodes[j] + deltas[i]);
      mu_interp(j,i) = mu(j,i);
      if (mu(j,i) > max_mu) { mu_interp(j,i) = max_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(j,i) = disps[i];
    }
  }

  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many columns as we have nodes and
  // as many columns as we have items


  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    grad_deltas[i] = 0;
    grad_disps[i] = 0;

    for(int j=0;j<n_nodes;j++) {
      // over nodes (rows in my matrices)
      double lambda = exp(log_lambda(j,i));
      double W = computeW(lambda, mu_interp(j,i), disps[i], 10);
      double R = computeR(lambda, mu_interp(j,i), disps[i], log_Z(j,i), 10);
      double frac_r_W = r(j,i) / W;
      double f_R = f(j,i)*R;
      double term = frac_r_W - h(j,i) - f_R;

      grad_deltas[i] = grad_deltas[i] +
        (mu_interp(j,i) / V(j,i))*(r(j,i) - mu_interp(j,i)*f(j,i));
      grad_disps[i] = grad_disps[i] + disps[i]*term;
    }
  }

  // fill up output vector
  for(int i=0;i<m;i++){
    out[i] = grad_deltas[i];
    out[i + m] = grad_disps[i];
  }

  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_logdisps_cpp(NumericVector alphas,
                                    NumericVector deltas,
                                    NumericVector disps,
                                    NumericVector nodes,
                                    NumericVector weights,
                                    NumericMatrix r,
                                    NumericMatrix f,
                                    NumericMatrix h,
                                    NumericVector grid_mus,
                                    NumericVector grid_nus,
                                    NumericVector grid_cmp_var_long,
                                    NumericVector grid_log_lambda_long,
                                    NumericVector grid_logZ_long,
                                    double max_mu) {

  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h

  int m = alphas.size();
  int n_nodes = nodes.size();
  NumericVector grad_disps(m);


  // set up mu's and nu's for interpolation function to be computed all in one

  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int i=0;i<m;i++){
    // loop over items (columns)
    for(int j=0;j<n_nodes;j++) {
      // loop over nodes (rows)
      mu(j,i) = exp(alphas[i] * nodes[j] + deltas[i]);
      mu_interp(j,i) = mu(j,i);
      if (mu(j,i) > max_mu) { mu_interp(j,i) = max_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(j,i) = disps[i];
    }
  }

  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many columns as we have nodes and
  // as many columns as we have items


  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    grad_disps[i] = 0;

    for(int j=0;j<n_nodes;j++) {
      // over nodes (rows in my matrices)
      double term = 0;
      //if ((r(j,i) > 1e-8) | (f(j,i) > 1e-8) | (h(j,i) > 1e-8)) {
      if (weights[j] > 1e-4) {
        double lambda = exp(log_lambda(j,i));
        double W = computeW(lambda, mu_interp(j,i), disps[i], 10);
        double R = computeR(lambda, mu_interp(j,i), disps[i], log_Z(j,i), 10);
        double frac_r_W = r(j,i) / W;
        double f_R = f(j,i)*R;
        term = frac_r_W - h(j,i) - f_R;
      }

      grad_disps[i] = grad_disps[i] + disps[i]*term;
    }
  }

  return(grad_disps);
}

// [[Rcpp::export]]
NumericVector grad_cmp_samedisp_cpp(NumericVector alphas,
                                    NumericVector deltas,
                                    double disp,
                                    NumericVector nodes,
                                    NumericVector weights,
                                    NumericMatrix r,
                                    NumericMatrix f,
                                    NumericMatrix h,
                                    NumericVector grid_mus,
                                    NumericVector grid_nus,
                                    NumericVector grid_cmp_var_long,
                                    NumericVector grid_log_lambda_long,
                                    NumericVector grid_logZ_long,
                                    NumericVector grid_W_long,
                                    NumericVector grid_R_long,
                                    double max_mu,
                                    double min_mu) {

  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h

  int m = alphas.size();
  int n_nodes = nodes.size();
  NumericVector grad_alphas(m);
  NumericVector grad_deltas(m);
  NumericVector grad_disps(m);
  NumericVector out(2*m+1);

  // set up mu's and nu's for interpolation function to be computed all in one

  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int i=0;i<m;i++){
    // loop over items (columns)
    for(int j=0;j<n_nodes;j++) {
      // loop over nodes (rows)
      mu(j,i) = exp(alphas[i] * nodes[j] + deltas[i]);
      mu_interp(j,i) = mu(j,i);
      if (mu(j,i) > max_mu) { mu_interp(j,i) = max_mu; }
      if (mu(j,i) < min_mu) { mu_interp(j,i) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(j,i) = disp;
    }
  }

  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  // NEW
  NumericMatrix W(n_nodes, m);
  NumericMatrix R(n_nodes, m);
  // NEW
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // NEW
  W = interp_from_grid_m(grid_mus, grid_nus,
                             grid_W_long,
                             mu_interp, disp_interp);
  R = interp_from_grid_m(grid_mus, grid_nus,
                         grid_R_long,
                         mu_interp, disp_interp);
  // NEW
  // V and log_lambda are matrices with as many columns as we have nodes and
  // as many columns as we have items

  double grad_disp_sum_across_items = 0;
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    grad_alphas[i] = 0;
    grad_deltas[i] = 0;
    grad_disps[i] = 0;

    for(int j=0;j<n_nodes;j++) {
      // over nodes (rows in my matrices)
      //double term = (-1)*h(j,i);
      //if ((r(j,i) > 1e-8) | (f(j,i) > 1e-8) | (h(j,i) > 1e-8)) {
      //if (weights[j] > 1e-32) {
        double lambda = exp(log_lambda(j,i));
        // NEW
        //double W = computeW(lambda, mu_interp(j,i), disp, 10);
        //double R = computeR(lambda, mu_interp(j,i), disp, log_Z(j,i), 10);
        double frac_r_W = r(j,i) / W(j,i);
        double f_R = f(j,i)*R(j,i);
        // NEW
        double term = frac_r_W - h(j,i) - f_R;
      //}

      grad_alphas[i] = grad_alphas[i] +
        (nodes[j]*mu_interp(j,i) / V(j,i))*(r(j,i) - mu_interp(j,i)*f(j,i));
      grad_deltas[i] = grad_deltas[i] +
        (mu_interp(j,i) / V(j,i))*(r(j,i) - mu_interp(j,i)*f(j,i));
      //grad_disps[i] = grad_disps[i] +
      //  disps[i]*(frac_r_W - frac_muf_W - h(j,i) + e_logfac*f(j,i));
      grad_disps[i] = grad_disps[i] + disp*term;
    }
    grad_disp_sum_across_items = grad_disp_sum_across_items + grad_disps[i];
  }

  // fill up output vector
  for(int i=0;i<m;i++){
    out[i] = grad_alphas[i];
    out[i + m] = grad_deltas[i];
  }
  out[2*m] = grad_disp_sum_across_items;

  return(out);
}


// [[Rcpp::export]]
NumericVector grad_cmp_samedisp_cpp_lininterp(NumericVector alphas,
                                      NumericVector deltas,
                                      double disp,
                                      NumericVector nodes,
                                      NumericVector weights,
                                      NumericMatrix r,
                                      NumericMatrix f,
                                      NumericMatrix h,
                                      NumericVector grid_mus,
                                      NumericVector grid_nus,
                                      NumericVector grid_cmp_var_long,
                                      NumericVector grid_log_lambda_long,
                                      NumericVector grid_W_long,
                                      NumericVector grid_R_long,
                                      double max_mu) {

  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h

  int m = alphas.size();
  int n_nodes = nodes.size();
  NumericVector grad_alphas(m);
  NumericVector grad_deltas(m);
  NumericVector grad_disps(m);
  NumericVector out(2*m+1);

  // set up mu's and nu's for interpolation function to be computed all in one

  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int i=0;i<m;i++){
    // loop over items (columns)
    for(int j=0;j<n_nodes;j++) {
      // loop over nodes (rows)
      mu(j,i) = exp(alphas[i] * nodes[j] + deltas[i]);
      mu_interp(j,i) = mu(j,i);
      if (mu(j,i) > max_mu) { mu_interp(j,i) = max_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(j,i) = disp;
    }
  }

  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  NumericMatrix W(n_nodes, m);
  NumericMatrix R(n_nodes, m);
  V = interp_from_grid_lin_m(grid_mus, grid_nus,
                           grid_cmp_var_long,
                           mu_interp, disp_interp);
  log_lambda = interp_from_grid_lin_m(grid_mus, grid_nus,
                                    grid_log_lambda_long,
                                    mu_interp, disp_interp);
  W = interp_from_grid_lin_m(grid_mus, grid_nus,
                             grid_W_long,
                             mu_interp, disp_interp);
  R = interp_from_grid_lin_m(grid_mus, grid_nus,
                             grid_R_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many columns as we have nodes and
  // as many columns as we have items

  double grad_disp_sum_across_items = 0;
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    grad_alphas[i] = 0;
    grad_deltas[i] = 0;
    grad_disps[i] = 0;

    for(int j=0;j<n_nodes;j++) {
      // over nodes (rows in my matrices)
      //double term = (-1)*h(j,i);
      //if ((r(j,i) > 1e-8) | (f(j,i) > 1e-8) | (h(j,i) > 1e-8)) {
      //if (weights[j] > 1e-32) {
      //double lambda = exp(log_lambda(j,i));
      double frac_r_W = r(j,i) / W(j,i);
      double f_R = f(j,i)*R(j,i);
      double term = frac_r_W - h(j,i) - f_R;
      //}

      grad_alphas[i] = grad_alphas[i] +
        (nodes[j]*mu_interp(j,i) / V(j,i))*(r(j,i) - mu_interp(j,i)*f(j,i));
      grad_deltas[i] = grad_deltas[i] +
        (mu_interp(j,i) / V(j,i))*(r(j,i) - mu_interp(j,i)*f(j,i));
        //grad_disps[i] = grad_disps[i] +
        //  disps[i]*(frac_r_W - frac_muf_W - h(j,i) + e_logfac*f(j,i));
      grad_disps[i] = grad_disps[i] + disp*term;
    }
    grad_disp_sum_across_items = grad_disp_sum_across_items + grad_disps[i];
  }

  // fill up output vector
  for(int i=0;i<m;i++){
    out[i] = grad_alphas[i];
    out[i + m] = grad_deltas[i];
  }
  out[2*m] = grad_disp_sum_across_items;

  return(out);
}



// [[Rcpp::export]]
NumericVector grad_cmp_samealpha_cpp(double alpha,
                                    NumericVector deltas,
                                    NumericVector disps,
                                    NumericVector nodes,
                                    NumericVector weights,
                                    NumericMatrix r,
                                    NumericMatrix f,
                                    NumericMatrix h,
                                    NumericVector grid_mus,
                                    NumericVector grid_nus,
                                    NumericVector grid_cmp_var_long,
                                    NumericVector grid_log_lambda_long,
                                    NumericVector grid_logZ_long,
                                    double max_mu) {

  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h

  int m = deltas.size();
  int n_nodes = nodes.size();
  NumericVector grad_alphas(m);
  NumericVector grad_deltas(m);
  NumericVector grad_disps(m);
  NumericVector out(2*m+1);

  // set up mu's and nu's for interpolation function to be computed all in one

  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int i=0;i<m;i++){
    // loop over items (columns)
    for(int j=0;j<n_nodes;j++) {
      // loop over nodes (rows)
      mu(j,i) = exp(alpha * nodes[j] + deltas[i]);
      mu_interp(j,i) = mu(j,i);
      if (mu(j,i) > max_mu) { mu_interp(j,i) = max_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(j,i) = disps[i];
    }
  }

  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many columns as we have nodes and
  // as many columns as we have items

  double grad_alpha_sum_across_items = 0;
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    grad_alphas[i] = 0;
    grad_deltas[i] = 0;
    grad_disps[i] = 0;

    for(int j=0;j<n_nodes;j++) {
      // over nodes (rows in my matrices)
      //double term = (-1)*h(j,i);
      //if ((r(j,i) > 1e-8) | (f(j,i) > 1e-8) | (h(j,i) > 1e-8)) {
      //if (weights[j] > 1e-32) {
      double lambda = exp(log_lambda(j,i));
      double W = computeW(lambda, mu_interp(j,i), disps[i], 10);
      double R = computeR(lambda, mu_interp(j,i), disps[i], log_Z(j,i), 10);
      double frac_r_W = r(j,i) / W;
      double f_R = f(j,i)*R;
      double term = frac_r_W - h(j,i) - f_R;
      //}

      grad_alphas[i] = grad_alphas[i] +
        (nodes[j]*mu_interp(j,i) / V(j,i))*(r(j,i) - mu_interp(j,i)*f(j,i));
      grad_deltas[i] = grad_deltas[i] +
        (mu_interp(j,i) / V(j,i))*(r(j,i) - mu_interp(j,i)*f(j,i));
      //grad_disps[i] = grad_disps[i] +
      //  disps[i]*(frac_r_W - frac_muf_W - h(j,i) + e_logfac*f(j,i));
      grad_disps[i] = grad_disps[i] + disps[i]*term;
    }
    grad_alpha_sum_across_items = grad_alpha_sum_across_items + grad_alphas[i];
  }

  // fill up output vector
  out[0] = grad_alpha_sum_across_items;
  for(int i=1;i<m+1;i++){
    out[i] = grad_deltas[i-1];
    out[i + m] = grad_disps[i-1];
  }

  return(out);
}


// [[Rcpp::export]]
double ell_cmp_cpp (NumericVector alphas,
                    NumericVector deltas,
                    NumericVector disps,
                    NumericVector nodes,
                    NumericMatrix r,
                    NumericMatrix f,
                    NumericMatrix h,
                    NumericVector grid_mus,
                    NumericVector grid_nus,
                    NumericVector grid_logZ_long,
                    NumericVector grid_log_lambda_long,
                    double max_mu) {

  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h

  int m = alphas.size();
  int n_nodes = nodes.size();
  double out = 0;

  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int i=0;i<m;i++){
    // loop over items (columns)
    for(int j=0;j<n_nodes;j++) {
      // loop over nodes (rows)
      mu(j,i) = exp(alphas[i] * nodes[j] + deltas[i]);
      mu_interp(j,i) = mu(j,i);
      if (mu(j,i) > max_mu) { mu_interp(j,i) = max_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(j,i) = disps[i];
    }
  }

  NumericMatrix log_Z(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many columns as we have nodes and
  // as many columns as we have items

  for(int i=0;i<m;i++){
    double sum_for_item = 0;
    for(int j=0;j<n_nodes;j++) {
      sum_for_item = sum_for_item + (r(j,i)*log_lambda(j,i) - disps[i]*h(j,i) - f(j,i)*log_Z(j,i));
    }
    out = out + sum_for_item;
  }

  return(out);
}


// [[Rcpp::export]]
NumericMatrix e_values_cpp (NumericMatrix data,
                            NumericVector alphas,
                            NumericVector deltas,
                            NumericVector disps,
                            NumericVector nodes,
                            NumericVector weights,
                            NumericVector grid_mus,
                            NumericVector grid_nus,
                            NumericVector grid_logZ_long,
                            NumericVector grid_log_lambda_long,
                            double max_mu,
                            double min_mu) {

  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h

  int m = alphas.size();
  int n_nodes = nodes.size();
  int N = data.nrow();

  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int i=0;i<m;i++){
    // loop over items (columns)
    for(int j=0;j<n_nodes;j++) {
      // loop over nodes (rows)
      mu(j,i) = exp(alphas[i] * nodes[j] + deltas[i]);
      mu_interp(j,i) = mu(j,i);
      if (mu(j,i) > max_mu) { mu_interp(j,i) = max_mu; }
      if (mu(j,i) < min_mu) { mu_interp(j,i) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(j,i) = disps[i];
    }
  }

  NumericMatrix log_Z(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                               grid_logZ_long,
                               mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                    grid_log_lambda_long,
                                    mu_interp, disp_interp);
  // V and log_lambda are matrices with as many columns as we have nodes and
  // as many columns as we have items

  NumericMatrix r(n_nodes, m);
  NumericMatrix f(n_nodes, m);
  NumericMatrix h(n_nodes, m);
  for (int k=0;k<n_nodes;k++) {
    for (int j=0;j<m;j++) {
      r(k, j) = 0;
      f(k, j) = 0;
      h(k, j) = 0;
    }
  }

  // cmp log density
  for (int j=0;j<m;j++){
    NumericMatrix item_pp(N, n_nodes);
    NumericMatrix item_log_dens(N, n_nodes);
    NumericMatrix item_pp_num(N, n_nodes);
    NumericVector item_marg_prob(N);
    // compute densities for all persons and nodes (for the item we're currently on)
    // and compute (log) marginal probability for all persons (for the item we're currently on)
    for(int i=0;i<N;i++){
      item_marg_prob[i] = 0;
      for (int k=0;k<n_nodes;k++){
        item_log_dens(i,k) = data(i,j)*log_lambda(k,j) - log_Z(k,j) - disps[j]*lgamma(data(i,j)+1);
        item_pp_num(i,k) = exp(item_log_dens(i,k) + log(weights[k]));
        item_marg_prob[i] = item_marg_prob[i] + item_pp_num(i,k);
      }
      // then compute pp for all persons and nodes (for the item we're currently on)
      // because pp are person and node specific, we need to again loop over nodes, but
      // we needed to complete first loop over nodes first because we need the marg probs
      for (int k=0;k<n_nodes;k++){
        item_pp(i,k) = item_pp_num(i,k) / item_marg_prob[i];
        r(k, j) = r(k, j) + data(i,j) * item_pp(i,k);
        f(k, j) = f(k, j) + item_pp(i,k);
        h(k, j) = h(k, j) + logFactorial(data(i,j)) * item_pp(i,k);
      }
    }
  }
  // pp and density are node, item and person specific
  // marginal prob is only item and person specific (summed over nodes)

  NumericMatrix out(3*n_nodes, m);
  // fill up output vector
  for(int k=0;k<n_nodes;k++){
    for(int j=0;j<m;j++) {
      out(k, j ) = r(k, j);
      out(k + n_nodes, j) = f(k, j);
      out(k + 2*n_nodes, j) = h(k, j);
    }
  }

  return(out);
}

// [[Rcpp::export]]
NumericMatrix e_values_cpp_lininterp (NumericMatrix data,
                            NumericVector alphas,
                            NumericVector deltas,
                            NumericVector disps,
                            NumericVector nodes,
                            NumericVector weights,
                            NumericVector grid_mus,
                            NumericVector grid_nus,
                            NumericVector grid_logZ_long,
                            NumericVector grid_log_lambda_long,
                            double max_mu) {

  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h

  int m = alphas.size();
  int n_nodes = nodes.size();
  int N = data.nrow();

  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int i=0;i<m;i++){
    // loop over items (columns)
    for(int j=0;j<n_nodes;j++) {
      // loop over nodes (rows)
      mu(j,i) = exp(alphas[i] * nodes[j] + deltas[i]);
      mu_interp(j,i) = mu(j,i);
      if (mu(j,i) > max_mu) { mu_interp(j,i) = max_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(j,i) = disps[i];
    }
  }

  NumericMatrix log_Z(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  log_Z = interp_from_grid_lin_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  log_lambda = interp_from_grid_lin_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many columns as we have nodes and
  // as many columns as we have items

  NumericMatrix r(n_nodes, m);
  NumericMatrix f(n_nodes, m);
  NumericMatrix h(n_nodes, m);
  for (int k=0;k<n_nodes;k++) {
    for (int j=0;j<m;j++) {
      r(k, j) = 0;
      f(k, j) = 0;
      h(k, j) = 0;
    }
  }

  // cmp log density
  for (int j=0;j<m;j++){
    NumericMatrix item_pp(N, n_nodes);
    NumericMatrix item_log_dens(N, n_nodes);
    NumericMatrix item_pp_num(N, n_nodes);
    NumericVector item_marg_prob(N);
    // compute densities for all persons and nodes (for the item we're currently on)
    // and compute (log) marginal probability for all persons (for the item we're currently on)
    for(int i=0;i<N;i++){
      item_marg_prob[i] = 0;
      for (int k=0;k<n_nodes;k++){
        item_log_dens(i,k) = data(i,j)*log_lambda(k,j) - log_Z(k,j) - disps[j]*lgamma(data(i,j)+1);
        item_pp_num(i,k) = exp(item_log_dens(i,k) + log(weights[k]));
        item_marg_prob[i] = item_marg_prob[i] + item_pp_num(i,k);
      }
      // then compute pp for all persons and nodes (for the item we're currently on)
      // because pp are person and node specific, we need to again loop over nodes, but
      // we needed to complete first loop over nodes first because we need the marg probs
      for (int k=0;k<n_nodes;k++){
        item_pp(i,k) = item_pp_num(i,k) / item_marg_prob[i];
        r(k, j) = r(k, j) + data(i,j) * item_pp(i,k);
        f(k, j) = f(k, j) + item_pp(i,k);
        h(k, j) = h(k, j) + logFactorial(data(i,j)) * item_pp(i,k);
      }
    }
  }
  // pp and density are node, item and person specific
  // marginal prob is only item and person specific (summed over nodes)

  NumericMatrix out(3*n_nodes, m);
  // fill up output vector
  for(int k=0;k<n_nodes;k++){
    for(int j=0;j<m;j++) {
      out(k, j ) = r(k, j);
      out(k + n_nodes, j) = f(k, j);
      out(k + 2*n_nodes, j) = h(k, j);
    }
  }

  return(out);
}

// [[Rcpp::export]]
double marg_ll_cpp (NumericMatrix data,
                    NumericVector alphas,
                    NumericVector deltas,
                    NumericVector disps,
                    NumericVector nodes,
                    NumericVector weights,
                    NumericVector grid_mus,
                    NumericVector grid_nus,
                    NumericVector grid_logZ_long,
                    NumericVector grid_log_lambda_long,
                    double max_mu,
                    double min_mu) {

  int N = data.nrow();
  int M = data.ncol();
  int K = nodes.size();

  NumericMatrix mu(K, M);
  NumericMatrix mu_interp(K, M);
  NumericMatrix disp_interp(K, M);
  for(int i=0;i<M;i++){
    // loop over items (columns)
    for(int j=0;j<K;j++) {
      // loop over nodes (rows)
      mu(j,i) = exp(alphas[i] * nodes[j] + deltas[i]);
      mu_interp(j,i) = mu(j,i);
      if (mu(j,i) > max_mu) { mu_interp(j,i) = max_mu; }
      if (mu(j,i) < min_mu) { mu_interp(j,i) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(j,i) = disps[i];
    }
  }

  NumericMatrix log_Z(K, M);
  NumericMatrix log_lambda(K, M);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many columns as we have nodes and
  // as many columns as we have items

  double log_marg_prob = 0;

  for(int i=0;i<N;i++){
    double integral = 0;
    for(int k=0;k<K;k++) {
      // qudrature over nodes
      double f = 0;
      for(int j=0;j<M;j++) {
        f = f + data(i,j)*log_lambda(k,j) - log_Z(k,j) - disps[j]*lgamma(data(i,j)+1);
      }
      integral = integral + exp(f + log(weights[k]));
    }
    log_marg_prob = log_marg_prob + log(integral);
  }

  return(log_marg_prob);
}

// [[Rcpp::export]]
double marg_ll_cmp_with_pcov_cpp (NumericMatrix data,
                                  NumericVector alphas,
                                  NumericVector deltas,
                                  NumericVector disps,
                                  NumericVector betas,
                                  NumericMatrix p_cov_data,
                                  NumericVector nodes,
                                  NumericVector weights,
                                  NumericVector grid_mus,
                                  NumericVector grid_nus,
                                  NumericVector grid_logZ_long,
                                  NumericVector grid_log_lambda_long,
                                  double max_mu,
                                  double min_mu) {
  
  int N = data.nrow();
  int M = data.ncol();
  int K = nodes.size();
  int P = betas.size();
  
  // for person covariates, we need mus (and lambdas and Zs) which are person
  // as well as node and item specific
  // so we set up KxM matrices (nodesxitems) which we row-bind under each other 
  // for all N persons
  NumericMatrix mu(K*N, M);
  NumericMatrix mu_interp(K*N, M);
  NumericMatrix disp_interp(K*N, M);
  for (int i=0; i<N; i++) {
    // we are computing node-item specific mus for each person
    for(int j=0; j<M; j++){
      // loop over items (columns)
      for(int k=0; k<K; k++) {
        // loop over nodes (rows)
        double log_mu = alphas[j] * nodes[k] + deltas[j];
        for(int p=0; p<P; p++) {
          // add all the (weighted) covariate values for all covariates for the item j
          // (for the specific person i we are currently looking at)
          log_mu += betas[p] * alphas[j] * p_cov_data(i,p);
        }
        mu(k+i*K,j) = exp(log_mu);
        mu_interp(k+i*K,j) = mu(k+i*K,j);
        if (mu(k+i*K,j) > max_mu) { mu_interp(k+i*K,j) = max_mu; }
        if (mu(k+i*K,j) < min_mu) { mu_interp(k+i*K,j) = min_mu; }
        // we need to set maximum for mu to max_mu so that the interpolation will
        // work, max_mu is the maximum mu value in our grid for interpolation
        disp_interp(k+i*K,j) = disps[j];
      }
    }  // end loop over items
  } // end loop over N
  
  NumericMatrix log_lambda(K*N, M);
  NumericMatrix log_Z(K*N, M);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes*persons and
  // as many columns as we have (same as mu_interp and nu_interp matrices)
  
  double log_marg_prob = 0;
  
  for(int i=0;i<N;i++){
    double integral = 0;
    for(int k=0;k<K;k++) {
      // qudrature over nodes
      double f = 0;
      for(int j=0;j<M;j++) {
        f = f + data(i,j)*log_lambda(k+i*K,j) - log_Z(k+i*K,j) - disps[j]*lgamma(data(i,j)+1);
      }
      integral = integral + exp(f + log(weights[k]));
    }
    log_marg_prob = log_marg_prob + log(integral);
  }
  
  return(log_marg_prob);
}

// [[Rcpp::export]]
double marg_ll_cmp_with_pcov_cat_cpp (NumericMatrix data,
                                  NumericVector alphas,
                                  NumericVector deltas,
                                  NumericVector disps,
                                  NumericVector betas,
                                  NumericMatrix p_cov_data,
                                  NumericMatrix resp_pattern, 
                                  NumericVector nodes,
                                  NumericVector weights,
                                  NumericVector grid_mus,
                                  NumericVector grid_nus,
                                  NumericVector grid_logZ_long,
                                  NumericVector grid_log_lambda_long,
                                  double max_mu,
                                  double min_mu) {
  // assume that p_cov is a matrix of dummy coded categorical predictors
  // resp_pattern is a matrix of the same no. of cols than p_cov
  // and as many rows as we have distinct possible response patterns
  
  int N = data.nrow();
  int M = data.ncol();
  int K = nodes.size();
  int P = betas.size(); 
  int n_resp_patterns = resp_pattern.nrow();
  
  // for person covariates, we need mus (and lambdas and Zs) for each node and
  // and then also for each response pattern
  // so first compute that
  NumericMatrix mu(K*n_resp_patterns, M);
  NumericMatrix mu_interp(K*n_resp_patterns, M);
  NumericMatrix disp_interp(K*n_resp_patterns, M);
  for (int l=0; l<n_resp_patterns; l++) {
    for(int j=0; j<M; j++){
      // loop over items (columns)
      for(int k=0; k<K; k++) {
        // loop over nodes (rows)
        double log_mu = alphas[j] * nodes[k] + deltas[j];
        for(int p=0; p<P; p++) {
          // this works because only includes columns for none-reference categories
          // for all covs in ref categories, resp_pattern will just always be zero in that row
          log_mu += betas[p] * alphas[j] * resp_pattern(l,p);
        }
        
        mu(k+l*K,j) = exp(log_mu);
        mu_interp(k+l*K,j) = mu(k+l*K,j);
        if (mu(k+l*K,j) > max_mu) { mu_interp(k+l*K,j) = max_mu; }
        if (mu(k+l*K,j) < min_mu) { mu_interp(k+l*K,j) = min_mu; }
        // we need to set maximum for mu to max_mu so that the interpolation will
        // work, max_mu is the maximum mu value in our grid for interpolation
        disp_interp(k+l*K,j) = disps[j];
      }
    }  // end loop over items
  }
  
  NumericMatrix log_lambda(K*n_resp_patterns, M);
  NumericMatrix log_Z(K*n_resp_patterns, M);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  
  double log_marg_prob = 0;
  
  for(int i=0;i<N;i++){
    // check what response pattern person i had
    int l = 0;
    bool pattern_match = false;
    while (!pattern_match && l<n_resp_patterns) {
      // the second condition is just for safety that we don't get an eternal while loop
      // but we should always find a pattern match
      for (int p=0; p<P; p++) {
        pattern_match = p_cov_data(i,p) == resp_pattern(l,p);
      }
      // if the rows are the same, i am going to get out the foor loop with
      // pattern_match = TRUE and have l at the row of the pattern matrix
      // otherwise I am going to increase l by 1 and stay in the while loop to see if
      // the next row in the pattern matrix is a match for i
      if (!pattern_match) { l += 1; }
    }
    
    // we now know that person i has a response pattern like in row l of resp_pattern matrix
    // so their mu (and lambda, etc.) should be at row k+l*K
    double integral = 0;
    for(int k=0;k<K;k++) {
      // qudrature over nodes
      double f = 0;
      for(int j=0;j<M;j++) {
        f = f + data(i,j)*log_lambda(k+l*K,j) - log_Z(k+l*K,j) - disps[j]*lgamma(data(i,j)+1);
      }
      integral = integral + exp(f + log(weights[k]));
    }
    log_marg_prob = log_marg_prob + log(integral);
  }
  
  return(log_marg_prob);
}

// [[Rcpp::export]]
double marg_ll_cmp_with_icov_delta_cpp (NumericMatrix data,
                                  NumericVector alphas,
                                  double delta,
                                  NumericVector disps,
                                  NumericVector betas,
                                  NumericMatrix i_cov_data,
                                  NumericVector nodes,
                                  NumericVector weights,
                                  NumericVector grid_mus,
                                  NumericVector grid_nus,
                                  NumericVector grid_logZ_long,
                                  NumericVector grid_log_lambda_long,
                                  double max_mu,
                                  double min_mu) {
  
  int N = data.nrow();
  int M = data.ncol();
  int K = nodes.size();
  int I = betas.size();
  
  // for item covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(K, M);
  NumericMatrix mu_interp(K, M);
  NumericMatrix disp_interp(K, M);
  for(int j=0;j<M;j++){
    // loop over items (columns)
    for(int k=0;k<K;k++) {
      // loop over nodes (rows)
      double log_mu = alphas[j] * nodes[k] + delta;
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += betas[c] * i_cov_data(j,c); // for item j
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,j) = disps[j];
    }
  }  // end loop over items
  
  NumericMatrix log_lambda(K, M);
  NumericMatrix log_Z(K, M);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  double log_marg_prob = 0;
  
  for(int i=0;i<N;i++){
    double integral = 0;
    for(int k=0;k<K;k++) {
      // qudrature over nodes
      double f = 0;
      for(int j=0;j<M;j++) {
        f = f + data(i,j)*log_lambda(k,j) - log_Z(k,j) - disps[j]*lgamma(data(i,j)+1);
      }
      integral = integral + exp(f + log(weights[k]));
    }
    log_marg_prob = log_marg_prob + log(integral);
  }
  
  return(log_marg_prob);
}

// [[Rcpp::export]]
double marg_ll_cmp_with_icov_alpha_cpp (NumericMatrix data,
                                        double alpha,
                                        NumericVector deltas,
                                        NumericVector disps,
                                        NumericVector betas,
                                        NumericMatrix i_cov_data,
                                        NumericVector nodes,
                                        NumericVector weights,
                                        NumericVector grid_mus,
                                        NumericVector grid_nus,
                                        NumericVector grid_logZ_long,
                                        NumericVector grid_log_lambda_long,
                                        double max_mu,
                                        double min_mu) {
  
  int N = data.nrow();
  int M = data.ncol();
  int K = nodes.size();
  int I = betas.size();
  
  // for item covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(K, M);
  NumericMatrix mu_interp(K, M);
  NumericMatrix disp_interp(K, M);
  for(int j=0;j<M;j++){
    // loop over items (columns)
    for(int k=0;k<K;k++) {
      // loop over nodes (rows)
      double log_mu = alpha * nodes[k] + deltas[j];
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += nodes[k] * betas[c] * i_cov_data(j,c); // for item j
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,j) = disps[j];
    }
  }  // end loop over items
  
  NumericMatrix log_lambda(K, M);
  NumericMatrix log_Z(K, M);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  double log_marg_prob = 0;
  
  for(int i=0;i<N;i++){
    double integral = 0;
    for(int k=0;k<K;k++) {
      // qudrature over nodes
      double f = 0;
      for(int j=0;j<M;j++) {
        f = f + data(i,j)*log_lambda(k,j) - log_Z(k,j) - disps[j]*lgamma(data(i,j)+1);
      }
      integral = integral + exp(f + log(weights[k]));
    }
    log_marg_prob = log_marg_prob + log(integral);
  }
  
  return(log_marg_prob);
}

// [[Rcpp::export]]
double marg_ll_cmp_with_icov_nu_cpp (NumericMatrix data,
                                        NumericVector alphas,
                                        NumericVector deltas,
                                        double disp,
                                        NumericVector betas,
                                        NumericMatrix i_cov_data,
                                        NumericVector nodes,
                                        NumericVector weights,
                                        NumericVector grid_mus,
                                        NumericVector grid_nus,
                                        NumericVector grid_logZ_long,
                                        NumericVector grid_log_lambda_long,
                                        double max_mu,
                                        double min_mu,
                                        double max_nu,
                                        double min_nu) {
  
  int N = data.nrow();
  int M = data.ncol();
  int K = nodes.size();
  int I = betas.size();
  
  // for item covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(K, M);
  NumericMatrix mu_interp(K, M);
  NumericMatrix disp_interp(K, M);
  for(int j=0;j<M;j++){
    // loop over items (columns)
    for(int k=0;k<K;k++) {
      // loop over nodes (rows)
      double log_mu = alphas[j] * nodes[k] + deltas[j];
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      double log_disp = log(disp);
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_disp += betas[c] * i_cov_data(j,c); // for item j
      }
      disp_interp(k,j) = exp(log_disp);
      if (disp_interp(k,j) > max_nu) { disp_interp(k,j) = max_nu; }
      if (disp_interp(k,j) < min_nu) { disp_interp(k,j) = min_nu; }
    }
  }  // end loop over items
  
  NumericMatrix log_lambda(K, M);
  NumericMatrix log_Z(K, M);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  double log_marg_prob = 0;
  
  for(int i=0;i<N;i++){
    double integral = 0;
    for(int k=0;k<K;k++) {
      // qudrature over nodes
      double f = 0;
      for(int j=0;j<M;j++) {
        f = f + data(i,j)*log_lambda(k,j) - log_Z(k,j) - disp_interp(k,j)*lgamma(data(i,j)+1);
      }
      integral = integral + exp(f + log(weights[k]));
    }
    log_marg_prob = log_marg_prob + log(integral);
  }
  
  return(log_marg_prob);
}

// [[Rcpp::export]]
double marg_ll_cmp_with_icov_all_cpp (NumericMatrix data,
                                     double alpha,
                                     double delta,
                                     double disp,
                                     NumericVector betas_alpha,
                                     NumericVector betas_delta,
                                     NumericVector betas_logdisp,
                                     NumericMatrix i_cov_data,
                                     NumericVector nodes,
                                     NumericVector weights,
                                     NumericVector grid_mus,
                                     NumericVector grid_nus,
                                     NumericVector grid_logZ_long,
                                     NumericVector grid_log_lambda_long,
                                     double max_mu,
                                     double min_mu,
                                     double max_nu,
                                     double min_nu) {
  
  int N = data.nrow();
  int M = data.ncol();
  int K = nodes.size();
  int I = betas_alpha.size();
  
  // for item covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(K, M);
  NumericMatrix mu_interp(K, M);
  NumericMatrix disp_interp(K, M);
  for(int j=0;j<M;j++){
    // loop over items (columns)
    for(int k=0;k<K;k++) {
      // loop over nodes (rows)
      double log_mu = alpha * nodes[k] + delta;
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += nodes[k] * betas_alpha[c] * i_cov_data(j,c) + 
          betas_delta[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      double log_disp = log(disp);
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_disp += betas_logdisp[c] * i_cov_data(j,c); // for item j
      }
      disp_interp(k,j) = exp(log_disp);
      if (disp_interp(k,j) > max_nu) { disp_interp(k,j) = max_nu; }
      if (disp_interp(k,j) < min_nu) { disp_interp(k,j) = min_nu; }
    }
  }  // end loop over items
  
  NumericMatrix log_lambda(K, M);
  NumericMatrix log_Z(K, M);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  double log_marg_prob = 0;
  
  for(int i=0;i<N;i++){
    double integral = 0;
    for(int k=0;k<K;k++) {
      // qudrature over nodes
      double f = 0;
      for(int j=0;j<M;j++) {
        f += data(i,j)*log_lambda(k,j) - log_Z(k,j) - disp_interp(k,j)*lgamma(data(i,j)+1);
      }
      integral += exp(f + log(weights[k]));
    }
    log_marg_prob += log(integral);
  }
  
  return(log_marg_prob);
}


// [[Rcpp::export]]
double marg_ll_cmp_with_icov_alpha_nu_cpp (NumericMatrix data,
                                      double alpha,
                                      NumericVector deltas,
                                      double disp,
                                      NumericVector betas_alpha,
                                      NumericVector betas_logdisp,
                                      NumericMatrix i_cov_data,
                                      NumericVector nodes,
                                      NumericVector weights,
                                      NumericVector grid_mus,
                                      NumericVector grid_nus,
                                      NumericVector grid_logZ_long,
                                      NumericVector grid_log_lambda_long,
                                      double max_mu,
                                      double min_mu,
                                      double max_nu,
                                      double min_nu) {
  
  int N = data.nrow();
  int M = data.ncol();
  int K = nodes.size();
  int I = betas_alpha.size();
  
  // for item covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(K, M);
  NumericMatrix mu_interp(K, M);
  NumericMatrix disp_interp(K, M);
  for(int j=0;j<M;j++){
    // loop over items (columns)
    for(int k=0;k<K;k++) {
      // loop over nodes (rows)
      double log_mu = alpha * nodes[k] + deltas[j];
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += nodes[k] * betas_alpha[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      double log_disp = log(disp);
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_disp += betas_logdisp[c] * i_cov_data(j,c); // for item j
      }
      disp_interp(k,j) = exp(log_disp);
      if (disp_interp(k,j) > max_nu) { disp_interp(k,j) = max_nu; }
      if (disp_interp(k,j) < min_nu) { disp_interp(k,j) = min_nu; }
    }
  }  // end loop over items
  
  NumericMatrix log_lambda(K, M);
  NumericMatrix log_Z(K, M);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  double log_marg_prob = 0;
  
  for(int i=0;i<N;i++){
    double integral = 0;
    for(int k=0;k<K;k++) {
      // qudrature over nodes
      double f = 0;
      for(int j=0;j<M;j++) {
        f += data(i,j)*log_lambda(k,j) - log_Z(k,j) - disp_interp(k,j)*lgamma(data(i,j)+1);
      }
      integral += exp(f + log(weights[k]));
    }
    log_marg_prob += log(integral);
  }
  
  return(log_marg_prob);
}

// [[Rcpp::export]]
double marg_ll_cmp_with_icov_delta_nu_cpp (NumericMatrix data,
                                      NumericVector alphas,
                                      double delta,
                                      double disp,
                                      NumericVector betas_delta,
                                      NumericVector betas_logdisp,
                                      NumericMatrix i_cov_data,
                                      NumericVector nodes,
                                      NumericVector weights,
                                      NumericVector grid_mus,
                                      NumericVector grid_nus,
                                      NumericVector grid_logZ_long,
                                      NumericVector grid_log_lambda_long,
                                      double max_mu,
                                      double min_mu,
                                      double max_nu,
                                      double min_nu) {
  
  int N = data.nrow();
  int M = data.ncol();
  int K = nodes.size();
  int I = betas_delta.size();
  
  // for item covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(K, M);
  NumericMatrix mu_interp(K, M);
  NumericMatrix disp_interp(K, M);
  for(int j=0;j<M;j++){
    // loop over items (columns)
    for(int k=0;k<K;k++) {
      // loop over nodes (rows)
      double log_mu = alphas[j] * nodes[k] + delta;
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += betas_delta[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      double log_disp = log(disp);
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_disp += betas_logdisp[c] * i_cov_data(j,c); // for item j
      }
      disp_interp(k,j) = exp(log_disp);
      if (disp_interp(k,j) > max_nu) { disp_interp(k,j) = max_nu; }
      if (disp_interp(k,j) < min_nu) { disp_interp(k,j) = min_nu; }
    }
  }  // end loop over items
  
  NumericMatrix log_lambda(K, M);
  NumericMatrix log_Z(K, M);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  double log_marg_prob = 0;
  
  for(int i=0;i<N;i++){
    double integral = 0;
    for(int k=0;k<K;k++) {
      // qudrature over nodes
      double f = 0;
      for(int j=0;j<M;j++) {
        f += data(i,j)*log_lambda(k,j) - log_Z(k,j) - disp_interp(k,j)*lgamma(data(i,j)+1);
      }
      integral += exp(f + log(weights[k]));
    }
    log_marg_prob += log(integral);
  }
  
  return(log_marg_prob);
}

// [[Rcpp::export]]
double marg_ll_cmp_with_icov_alpha_delta_cpp (NumericMatrix data,
                                      double alpha,
                                      double delta,
                                      NumericVector disps,
                                      NumericVector betas_alpha,
                                      NumericVector betas_delta,
                                      NumericMatrix i_cov_data,
                                      NumericVector nodes,
                                      NumericVector weights,
                                      NumericVector grid_mus,
                                      NumericVector grid_nus,
                                      NumericVector grid_logZ_long,
                                      NumericVector grid_log_lambda_long,
                                      double max_mu,
                                      double min_mu) {
  
  int N = data.nrow();
  int M = data.ncol();
  int K = nodes.size();
  int I = betas_alpha.size();
  
  // for item covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(K, M);
  NumericMatrix mu_interp(K, M);
  NumericMatrix disp_interp(K, M);
  for(int j=0;j<M;j++){
    // loop over items (columns)
    for(int k=0;k<K;k++) {
      // loop over nodes (rows)
      double log_mu = alpha * nodes[k] + delta;
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += nodes[k] * betas_alpha[c] * i_cov_data(j,c) + 
          betas_delta[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,j) = disps[j];
    }
  }  // end loop over items
  
  NumericMatrix log_lambda(K, M);
  NumericMatrix log_Z(K, M);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  double log_marg_prob = 0;
  
  for(int i=0;i<N;i++){
    double integral = 0;
    for(int k=0;k<K;k++) {
      // qudrature over nodes
      double f = 0;
      for(int j=0;j<M;j++) {
        f += data(i,j)*log_lambda(k,j) - log_Z(k,j) - disp_interp(k,j)*lgamma(data(i,j)+1);
      }
      integral += exp(f + log(weights[k]));
    }
    log_marg_prob += log(integral);
  }
  
  return(log_marg_prob);
}

// [[Rcpp::export]]
double marg_ll_cpp_lininterp (NumericMatrix data,
                      NumericVector alphas,
                      NumericVector deltas,
                      NumericVector disps,
                      NumericVector nodes,
                      NumericVector weights,
                      NumericVector grid_mus,
                      NumericVector grid_nus,
                      NumericVector grid_logZ_long,
                      NumericVector grid_log_lambda_long) {

  int N = data.nrow();
  int M = data.ncol();
  int K = nodes.size();

  NumericMatrix mu(K, M);
  NumericMatrix mu_interp(K, M);
  NumericMatrix disp_interp(K, M);
  for(int i=0;i<M;i++){
    // loop over items (columns)
    for(int j=0;j<K;j++) {
      // loop over nodes (rows)
      mu(j,i) = exp(alphas[i] * nodes[j] + deltas[i]);
      mu_interp(j,i) = mu(j,i);
      //if (mu(j,i) > max_mu) { mu_interp(j,i) = max_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(j,i) = disps[i];
    }
  }

  NumericMatrix log_Z(K, M);
  NumericMatrix log_lambda(K, M);
  log_Z = interp_from_grid_lin_m(grid_mus, grid_nus,
                               grid_logZ_long,
                               mu_interp, disp_interp);
  log_lambda = interp_from_grid_lin_m(grid_mus, grid_nus,
                                    grid_log_lambda_long,
                                    mu_interp, disp_interp);
    // V and log_lambda are matrices with as many columns as we have nodes and
    // as many columns as we have items

  double log_marg_prob = 0;

  for(int i=0;i<N;i++){
    double integral = 0;
    for(int k=0;k<K;k++) {
      // qudrature over nodes
      double f = 0;
      for(int j=0;j<M;j++) {
        f = f + data(i,j)*log_lambda(k,j) - log_Z(k,j) - disps[j]*lgamma(data(i,j)+1);
      }
      integral = integral + exp(f + log(weights[k]));
    }
    log_marg_prob = log_marg_prob + log(integral);
  }

  return(log_marg_prob);
}


// [[Rcpp::export]]
NumericVector grad_cmp_newem_cpp(NumericVector alphas,
                                 NumericVector deltas,
                                 NumericVector disps,
                                 NumericMatrix data,
                                 NumericVector exp_abilities,
                                 NumericVector grid_mus,
                                 NumericVector grid_nus,
                                 NumericVector grid_cmp_var_long,
                                 NumericVector grid_log_lambda_long,
                                 NumericVector grid_logZ_long,
                                 double max_mu,
                                 double min_mu) {

    // r needs to be a matrix with one column per item and then the r values
    // for this item in the column
    // analogously for f and h

    int m = alphas.size();
    int n = exp_abilities.size();
    NumericVector grad_alphas(m);
    NumericVector grad_deltas(m);
    NumericVector grad_disps(m);
    NumericVector out(3*m);

    // set up mu's and nu's for interpolation function to be computed all in one

    NumericMatrix mu(n, m);
    NumericMatrix mu_interp(n, m);
    NumericMatrix disp_interp(n, m);
    for(int i=0;i<m;i++){
      // loop over items (columns)
      for(int j=0;j<n;j++) {
        // loop over persons (rows)
        mu(j,i) = exp(alphas[i] * exp_abilities[j] + deltas[i]);
        mu_interp(j,i) = mu(j,i);
        if (mu(j,i) > max_mu) { mu_interp(j,i) = max_mu; }
        if (mu(j,i) < min_mu) { mu_interp(j,i) = min_mu; }
        // we need to set maximum for mu to max_mu so that the interpolation will
        // work, max_mu is the maximum mu value in our grid for interpolation
        disp_interp(j,i) = disps[i];
      }
    }

    NumericMatrix V(n, m);
    NumericMatrix log_lambda(n, m);
    NumericMatrix log_Z(n, m);
    V = interp_from_grid_m(grid_mus, grid_nus,
                           grid_cmp_var_long,
                           mu_interp, disp_interp);
    log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                    grid_log_lambda_long,
                                    mu_interp, disp_interp);
    log_Z = interp_from_grid_m(grid_mus, grid_nus,
                               grid_logZ_long,
                               mu_interp, disp_interp);
    // V and log_lambda are matrices with as many rows as we have persons and
    // as many columns as we have items


    for(int i=0;i<m;i++){
      // over items (columns in my matrices)
      // so that we get one gradient per item
      grad_alphas[i] = 0;
      grad_deltas[i] = 0;
      grad_disps[i] = 0;

      for(int j=0;j<n;j++) {
        // over persons (rows in my matrices)

        // compute A and B for dispersion gradient
        double lambda = exp(log_lambda(j,i));
        double A = computeA(lambda, mu_interp(j,i), disps[i], log_Z(j,i), 10);
        double B = computeB(lambda, mu_interp(j,i), disps[i], log_Z(j,i), 10);

        // compute the gradients (summing over persons)
        grad_alphas[i] = grad_alphas[i] +
          (exp_abilities[j]*mu_interp(j,i) / V(j,i))*(data(j,i) - mu_interp(j,i));
        grad_deltas[i] = grad_deltas[i] +
          (mu_interp(j,i) / V(j,i))*(data(j,i) - mu_interp(j,i));
        grad_disps[i] = grad_disps[i] +
          (disps[i]*(A*(data(j,i) - mu_interp(j,i))/V(j,i) - (logFactorial(data(j,i))-B)));
      }
    }

    // fill up output vector
    for(int i=0;i<m;i++){
      out[i] = grad_alphas[i];
      out[i + m] = grad_deltas[i];
      out[i + 2*m] = grad_disps[i];
    }

    return(out);
  }


// [[Rcpp::export]]
NumericVector grad_cmp_newem_cpp2(NumericVector alphas,
                                 NumericVector deltas,
                                 NumericVector disps,
                                 NumericMatrix data,
                                 NumericMatrix PPs,
                                 NumericVector nodes, 
                                 NumericVector grid_mus,
                                 NumericVector grid_nus,
                                 NumericVector grid_cmp_var_long,
                                 NumericVector grid_log_lambda_long,
                                 NumericVector grid_logZ_long,
                                 double max_mu,
                                 double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = alphas.size();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  NumericVector grad_alphas(m);
  NumericVector grad_deltas(m);
  NumericVector grad_disps(m);
  NumericVector out(3*m);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int i=0;i<m;i++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over persons (rows)
      mu(k,i) = exp(alphas[i] * nodes[k] + deltas[i]);
      mu_interp(k,i) = mu(k,i);
      if (mu(k,i) > max_mu) { mu_interp(k,i) = max_mu; }
      if (mu(k,i) < min_mu) { mu_interp(k,i) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,i) = disps[i];
    }
  }
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have persons and
  // as many columns as we have items
  
  
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_alphas[i] = 0;
    grad_deltas[i] = 0;
    grad_disps[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
    // over persons (rows in my matrices)
    
    // compute A and B for dispersion gradient
      double lambda = exp(log_lambda(k,i));
      double A = computeA(lambda, mu_interp(k,i), disps[i], log_Z(k,i), 10);
      double B = computeB(lambda, mu_interp(k,i), disps[i], log_Z(k,i), 10);
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_alphas[i] = grad_alphas[i] +
          PPs(j,k) * (nodes[k]*mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_deltas[i] = grad_deltas[i] +
          PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_disps[i] = grad_disps[i] +
          PPs(j,k) * (disps[i]*(A*(data(j,i) - mu_interp(k,i))/V(k,i) - (logFactorial(data(j,i))-B)));
      }
    }
  }
  
  // fill up output vector
  for(int i=0;i<m;i++){
    out[i] = grad_alphas[i];
    out[i + m] = grad_deltas[i];
    out[i + 2*m] = grad_disps[i];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_pcov_cpp(NumericVector alphas,
                                    NumericVector deltas,
                                    NumericVector disps,
                                    NumericVector betas,
                                    NumericMatrix data,
                                    NumericMatrix p_cov_data,
                                    NumericMatrix PPs,
                                    NumericVector nodes, 
                                    NumericVector grid_mus,
                                    NumericVector grid_nus,
                                    NumericVector grid_cmp_var_long,
                                    NumericVector grid_log_lambda_long,
                                    NumericVector grid_logZ_long,
                                    double max_mu,
                                    double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = alphas.size();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int P = betas.size();
  NumericVector grad_alphas(m);
  NumericVector grad_deltas(m);
  NumericVector grad_disps(m);
  NumericVector grad_betas(P);
  NumericVector out(3*m + P);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are person
  // as well as node and item specific
  // so we set up KxM matrices (nodesxitems) which we row-bind under each other 
  // for all N persons
  NumericMatrix mu(n_nodes*n, m);
  NumericMatrix mu_interp(n_nodes*n, m);
  NumericMatrix disp_interp(n_nodes*n, m);
  for (int i=0; i<n; i++) {
    // we are computing node-item specific mus for each person
    for(int j=0;j<m;j++){
      // loop over items (columns)
      for(int k=0;k<n_nodes;k++) {
        // loop over nodes (rows)
        double log_mu = alphas[j] * nodes[k] + deltas[j];
        for(int p=0; p<P; p++) {
          // add all the (weighted) covariate values for all covariates for the item j
          // (for the specific person i we are currently looking at)
          log_mu += betas[p] * alphas[j] * p_cov_data(i,p);
        }
        mu(k+i*n_nodes,j) = exp(log_mu);
        mu_interp(k+i*n_nodes,j) = mu(k+i*n_nodes,j);
        if (mu(k+i*n_nodes,j) > max_mu) { mu_interp(k+i*n_nodes,j) = max_mu; }
        if (mu(k+i*n_nodes,j) < min_mu) { mu_interp(k+i*n_nodes,j) = min_mu; }
        // we need to set maximum for mu to max_mu so that the interpolation will
        // work, max_mu is the maximum mu value in our grid for interpolation
        disp_interp(k+i*n_nodes,j) = disps[j];
      }
    }  // end loop over items
  } // end loop over N
  
  NumericMatrix V(n_nodes*n, m);
  NumericMatrix log_lambda(n_nodes*n, m);
  NumericMatrix log_Z(n_nodes*n, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes*persons and
  // as many columns as we have (same as mu_interp and nu_interp matrices)

  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_alphas[i] = 0;
    grad_deltas[i] = 0;
    grad_disps[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      // note the index for mu_interp and V_interp (and log_lambda and log_Z) must be adjusted 
      // as we now have (additionally) person specific mus and Vs
      for(int j=0;j<n;j++) {
        // loop over persons
      
        // compute A and B for dispersion gradient
        double lambda = exp(log_lambda(k+j*n_nodes,i));
        double A = computeA(lambda, mu_interp(k+j*n_nodes,i), disps[i], log_Z(k+j*n_nodes,i), 10);
        double B = computeB(lambda, mu_interp(k+j*n_nodes,i), disps[i], log_Z(k+j*n_nodes,i), 10);
        
        // compute the sum over the weightes covariates for the gradient for alpha
        double sum_over_pcov = 0;
        for (int p=0; p<P; p++) {
          sum_over_pcov += betas[p] * p_cov_data(j,p);
        }
        
        // compute the gradients (summing over persons)
        grad_alphas[i] = grad_alphas[i] +
          PPs(j,k) * (mu_interp(k+j*n_nodes,i)*(nodes[k] + sum_over_pcov) /  V(k+j*n_nodes,i))*
          (data(j,i) -  mu_interp(k+j*n_nodes,i));
        grad_deltas[i] = grad_deltas[i] +
          PPs(j,k) * (mu_interp(k+j*n_nodes,i) / V(k+j*n_nodes,i))*
          (data(j,i) - mu_interp(k+j*n_nodes,i));
        grad_disps[i] = grad_disps[i] +
          PPs(j,k) * (disps[i]*(A*(data(j,i) - mu_interp(k+j*n_nodes,i))/V(k+j*n_nodes,i) - 
          (logFactorial(data(j,i))-B)));
      }
    }
  }
  
  // gradients for person covariate weights
  for (int p=0; p<P; p++) { // over covariates
    grad_betas[p] = 0;
    for (int j=0; j<m; j++) { // over items
      for (int k=0;k<n_nodes;k++) { // over nodes (rows in my matrices)
        for (int i=0; i<n; i++) { // over persons
          grad_betas[p] += PPs(i,k) * (mu_interp(k+i*n_nodes,j) * alphas[j] * p_cov_data(i,p) / V(k+i*n_nodes,j)) *
            (data(i,j) - mu_interp(k+i*n_nodes,j));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over items
  } //end loop over P (person covariates)
  
  // fill up output vector
  for(int i=0;i<m;i++){
    out[i] = grad_alphas[i];
    out[i + m] = grad_deltas[i];
    out[i + 2*m] = grad_disps[i];
  }
  for(int p=0; p<P; p++) {
    out[3*m + p] = grad_betas[p];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_pcov_cat_cpp(NumericVector alphas,
                                     NumericVector deltas,
                                     NumericVector disps,
                                     NumericVector betas,
                                     NumericMatrix data,
                                     NumericMatrix p_cov_data,
                                     NumericMatrix resp_pattern,
                                     NumericMatrix PPs,
                                     NumericVector nodes, 
                                     NumericVector grid_mus,
                                     NumericVector grid_nus,
                                     NumericVector grid_cmp_var_long,
                                     NumericVector grid_log_lambda_long,
                                     NumericVector grid_logZ_long,
                                     double max_mu,
                                     double min_mu) {
  
  // assume that p_cov is a matrix of dummy coded categorical predictors
  // resp_pattern is a matrix of the same no. of cols than p_cov
  // and as many rows as we have distinct possible response patterns
  
  int N = data.nrow();
  int M = data.ncol();
  int K = nodes.size();
  int P = betas.size(); 
  int n_resp_patterns = resp_pattern.nrow();
  NumericVector grad_alphas(M);
  NumericVector grad_deltas(M);
  NumericVector grad_disps(M);
  NumericVector grad_betas(P);
  NumericVector out(3*M + P);
  
  // for person covariates, we need mus (and lambdas and Zs) for each node and
  // and then also for each response pattern
  // so first compute that
  NumericMatrix mu(K*n_resp_patterns, M);
  NumericMatrix mu_interp(K*n_resp_patterns, M);
  NumericMatrix disp_interp(K*n_resp_patterns, M);
  for (int l=0; l<n_resp_patterns; l++) {
    for(int j=0; j<M; j++){
      // loop over items (columns)
      for(int k=0; k<K; k++) {
        // loop over nodes (rows)
        double log_mu = alphas[j] * nodes[k] + deltas[j];
        for(int p=0; p<P; p++) {
          // this works because only includes columns for none-reference categories
          // for all covs in ref categories, resp_pattern will just always be zero in that row
          log_mu += betas[p] * alphas[j] * resp_pattern(l,p);
        }
        
        mu(k+l*K,j) = exp(log_mu);
        mu_interp(k+l*K,j) = mu(k+l*K,j);
        if (mu(k+l*K,j) > max_mu) { mu_interp(k+l*K,j) = max_mu; }
        if (mu(k+l*K,j) < min_mu) { mu_interp(k+l*K,j) = min_mu; }
        // we need to set maximum for mu to max_mu so that the interpolation will
        // work, max_mu is the maximum mu value in our grid for interpolation
        disp_interp(k+l*K,j) = disps[j];
      }
    }  // end loop over items
  }
  
  NumericMatrix log_lambda(K*n_resp_patterns, M);
  NumericMatrix log_Z(K*n_resp_patterns, M);
  NumericMatrix V(K*n_resp_patterns, M);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);

  
  for(int i=0;i<M;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_alphas[i] = 0;
    grad_deltas[i] = 0;
    grad_disps[i] = 0;
    
    for(int k=0; k<K; k++) {
      // over nodes (rows in my matrices)
      
      for(int j=0;j<N;j++) {
        // loop over persons
        
        // check what response pattern person i had
        int l = 0;
        bool pattern_match = false;
        while (!pattern_match && l<n_resp_patterns) {
          // the second condition is just for safety that we dont get an eternal while loop
          // but we should always find a pattern match
          for (int p=0; p<P; p++) {
            pattern_match = p_cov_data(j,p) == resp_pattern(l,p);
          }
          // if the rows are the same, i am going to get out the foor loop with
          // pattern_match = TRUE and have l at the row of the pattern matrix
          // otherwise I am going to increase l by 1 and stay in the while loop to see if
          // the next row in the pattern matrix is a match for i
          if (!pattern_match) { l += 1; }
        }
        
        // we now know that person i has a response pattern like in row l of resp_pattern matrix
        // so their mu (and lambda, etc.) should be at row k+l*K
        
        // compute A and B for dispersion gradient
        double lambda = exp(log_lambda(k+l*K,i));
        double A = computeA(lambda, mu_interp(k+l*K,i), disps[i], log_Z(k+l*K,i), 10);
        double B = computeB(lambda, mu_interp(k+l*K,i), disps[i], log_Z(k+l*K,i), 10);
        
        // compute the sum over the weightes covariates for the gradient for alpha
        double sum_over_pcov = 0;
        for (int p=0; p<P; p++) {
          sum_over_pcov += betas[p] * p_cov_data(j,p);
        }
        
        // compute the gradients (summing over persons)
        grad_alphas[i] += PPs(j,k) * (mu_interp(k+l*K,i)*(nodes[k] + sum_over_pcov) /  V(k+l*K,i))*
          (data(j,i) -  mu_interp(k+l*K,i));
        grad_deltas[i] += PPs(j,k) * (mu_interp(k+l*K,i) / V(k+l*K,i))*
          (data(j,i) - mu_interp(k+l*K,i));
        grad_disps[i] += PPs(j,k) * (disps[i]*(A*(data(j,i) - mu_interp(k+l*K,i))/V(k+l*K,i) - 
          (logFactorial(data(j,i))-B)));
      }
    }
  }
  
  // gradients for person covariate weights
  for (int p=0; p<P; p++) { // over covariates
    grad_betas[p] = 0;
    for (int j=0; j<M; j++) { // over items
      for (int k=0; k<K; k++) { // over nodes (rows in my matrices)
        for (int i=0; i<N; i++) { // over persons
          // check what response pattern person i had
          int l = 0;
          bool pattern_match = false;
          while (!pattern_match && l<n_resp_patterns) {
            // the second condition is just for safety that we dont get an eternal while loop
            // but we should always find a pattern match
            for (int p=0; p<P; p++) {
              pattern_match = p_cov_data(i,p) == resp_pattern(l,p);
            }
            // if the rows are the same, i am going to get out the foor loop with
            // pattern_match = TRUE and have l at the row of the pattern matrix
            // otherwise I am going to increase l by 1 and stay in the while loop to see if
            // the next row in the pattern matrix is a match for i
            if (!pattern_match) { l += 1; }
          }
          
          // we now know that person i has a response pattern like in row l of resp_pattern matrix
          // so their mu (and lambda, etc.) should be at row k+l*K
          
          grad_betas[p] += PPs(i,k) * (mu_interp(k+l*K,j) * alphas[j] * p_cov_data(i,p) / V(k+l*K,j)) *
            (data(i,j) - mu_interp(k+l*K,j));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over items
  } //end loop over P (person covariates)
  
  // fill up output vector
  for(int i=0;i<M;i++){
    out[i] = grad_alphas[i];
    out[i + M] = grad_deltas[i];
    out[i + 2*M] = grad_disps[i];
  }
  for(int p=0; p<P; p++) {
    out[3*M + p] = grad_betas[p];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_icov_delta_cpp(NumericVector alphas,
                                     double delta,
                                     NumericVector disps,
                                     NumericVector betas,
                                     NumericMatrix data,
                                     NumericMatrix i_cov_data,
                                     NumericMatrix PPs,
                                     NumericVector nodes, 
                                     NumericVector grid_mus,
                                     NumericVector grid_nus,
                                     NumericVector grid_cmp_var_long,
                                     NumericVector grid_log_lambda_long,
                                     NumericVector grid_logZ_long,
                                     double max_mu,
                                     double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = alphas.size();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int I = betas.size();
  NumericVector grad_alphas(m);
  double grad_delta;
  NumericVector grad_disps(m);
  NumericVector grad_betas(I);
  NumericVector out(2*m + 1 + I);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alphas[j] * nodes[k] + delta; // not
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += betas[c] * i_cov_data(j,c); // for item j
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,j) = disps[j];
    }
  }  // end loop over items
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  grad_delta = 0;
  // gradients for item parameters
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_alphas[i] = 0;
    grad_disps[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      // compute A and B for dispersion gradient
      double lambda = exp(log_lambda(k,i));
      double A = computeA(lambda, mu_interp(k,i), disps[i], log_Z(k,i), 10);
      double B = computeB(lambda, mu_interp(k,i), disps[i], log_Z(k,i), 10);
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_alphas[i] = grad_alphas[i] +
          PPs(j,k) * (nodes[k]*mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_delta += PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_disps[i] = grad_disps[i] +
          PPs(j,k) * (disps[i]*(A*(data(j,i) - mu_interp(k,i))/V(k,i) - (logFactorial(data(j,i))-B)));
      }
    }
  }
  
  // gradients for item covariate weights
  for (int c=0; c<I; c++) {
    // for each gamma of which we have one for each covariate-item combination
    grad_betas[c] = 0;
    for (int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      for (int i=0; i<n; i++) {
        // over persons
        for (int j=0; j<m; j++) {
          // over items (as the betas are only specific to item covariates, not items)
          grad_betas[c] += PPs(i,k) * (mu_interp(k,j)*i_cov_data(j,c) / V(k,j)) *
            (data(i,j) - mu_interp(k,j));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over P (person covariates)
  } // end loop over items
  
  // fill up output vector
  for(int i=0;i<m;i++){
    out[i] = grad_alphas[i];
    out[i + m + 1] = grad_disps[i];
  }
  out[m] = grad_delta;
  for(int c=0; c<I; c++) {
    out[2*m + 1 + c] = grad_betas[c];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_icov_alpha_cpp(double alpha,
                                           NumericVector deltas,
                                           NumericVector disps,
                                           NumericVector betas,
                                           NumericMatrix data,
                                           NumericMatrix i_cov_data,
                                           NumericMatrix PPs,
                                           NumericVector nodes, 
                                           NumericVector grid_mus,
                                           NumericVector grid_nus,
                                           NumericVector grid_cmp_var_long,
                                           NumericVector grid_log_lambda_long,
                                           NumericVector grid_logZ_long,
                                           double max_mu,
                                           double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = data.ncol();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int I = betas.size();
  double grad_alpha;
  NumericVector grad_deltas(m);
  NumericVector grad_disps(m);
  NumericVector grad_betas(I);
  NumericVector out(2*m + 1 + I);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alpha * nodes[k] + deltas[j]; // not
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += nodes[k] * betas[c] * i_cov_data(j,c); // for item j
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,j) = disps[j];
    }
  }  // end loop over items
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  grad_alpha = 0;
  // gradients for item parameters
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_deltas[i] = 0;
    grad_disps[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      // compute A and B for dispersion gradient
      double lambda = exp(log_lambda(k,i));
      double A = computeA(lambda, mu_interp(k,i), disps[i], log_Z(k,i), 10);
      double B = computeB(lambda, mu_interp(k,i), disps[i], log_Z(k,i), 10);
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_alpha += PPs(j,k) * (nodes[k]*mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_deltas[i] += PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_disps[i] += PPs(j,k) * (disps[i]*(A*(data(j,i) - mu_interp(k,i))/V(k,i) - (logFactorial(data(j,i))-B)));
      }
    }
  }
  
  // gradients for item covariate weights
  for (int c=0; c<I; c++) {
    // for each gamma of which we have one for each covariate-item combination
    grad_betas[c] = 0;
    for (int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      for (int i=0; i<n; i++) {
        // over persons
        for (int j=0; j<m; j++) {
          // over items (as the betas are only specific to item covariates, not items)
          grad_betas[c] += PPs(i,k) * (nodes[k]*mu_interp(k,j)*i_cov_data(j,c) / V(k,j)) *
            (data(i,j) - mu_interp(k,j));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over P (person covariates)
  } // end loop over items
  
  // fill up output vector
  out[0] = grad_alpha;
  for(int i=0;i<m;i++){
    out[i + 1] = grad_deltas[i];
    out[i + m + 1] = grad_disps[i];
  }
  for(int c=0; c<I; c++) {
    out[2*m + 1 + c] = grad_betas[c];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_icov_nu_cpp(NumericVector alphas,
                                           NumericVector deltas,
                                           double disp,
                                           NumericVector betas,
                                           NumericMatrix data,
                                           NumericMatrix i_cov_data,
                                           NumericMatrix PPs,
                                           NumericVector nodes, 
                                           NumericVector grid_mus,
                                           NumericVector grid_nus,
                                           NumericVector grid_cmp_var_long,
                                           NumericVector grid_log_lambda_long,
                                           NumericVector grid_logZ_long,
                                           double max_mu,
                                           double min_mu,
                                           double max_nu,
                                           double min_nu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = data.ncol();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int I = betas.size();
  NumericVector grad_alphas(m);
  NumericVector grad_deltas(m);
  double grad_disp;
  NumericVector grad_betas(I);
  NumericVector out(2*m + 1 + I);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alphas[j] * nodes[k] + deltas[j];
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      double log_disp = log(disp);
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_disp += betas[c] * i_cov_data(j,c); // for item j
      }
      disp_interp(k,j) = exp(log_disp);
      if (disp_interp(k,j) > max_nu) { disp_interp(k,j) = max_nu; }
      if (disp_interp(k,j) < min_nu) { disp_interp(k,j) = min_nu; }
    }
  }  // end loop over items
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  grad_disp = 0;
  NumericMatrix A(n_nodes, m);
  NumericMatrix B(n_nodes, m);
  // gradients for item parameters
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_deltas[i] = 0;
    grad_alphas[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      // compute A and B for dispersion gradient
      double lambda = exp(log_lambda(k,i));
      A(k,i) = computeA(lambda, mu_interp(k,i), disp_interp(k,i), log_Z(k,i), 10);
      B(k,i) = computeB(lambda, mu_interp(k,i), disp_interp(k,i), log_Z(k,i), 10);
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_alphas[i] += PPs(j,k) * (nodes[k]*mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_deltas[i] += PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_disp += PPs(j,k) * (disp_interp(k,i)*(A(k,i)*
          (data(j,i) - mu_interp(k,i))/V(k,i) - (logFactorial(data(j,i))-B(k,i))));
      }
    }
  }
  
  // gradients for item covariate weights
  for (int c=0; c<I; c++) {
    // for each gamma of which we have one for each covariate-item combination
    grad_betas[c] = 0;
    for (int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      for (int i=0; i<n; i++) {
        // over persons
        for (int j=0; j<m; j++) {
          // over items (as the betas are only specific to item covariates, not items)
          grad_betas[c] += PPs(i,k) * i_cov_data(j,c) * disp_interp(k,j)*
            (A(k,j)*(data(i,j) - mu_interp(k,j))/V(k,j) - (logFactorial(data(i,j))-B(k,j)));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over P (person covariates)
  } // end loop over items
  
  // fill up output vector
  for(int i=0;i<m;i++){
    out[i] = grad_alphas[i];
    out[i + m] = grad_deltas[i];
  }
  out[2*m] = grad_disp;
  for(int c=0; c<I; c++) {
    out[2*m + 1 + c] = grad_betas[c];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_icov_all_cpp(double alpha,
                                        double delta,
                                        double disp,
                                        NumericVector betas_alpha,
                                        NumericVector betas_delta,
                                        NumericVector betas_logdisp,
                                        NumericMatrix data,
                                        NumericMatrix i_cov_data,
                                        NumericMatrix PPs,
                                        NumericVector nodes, 
                                        NumericVector grid_mus,
                                        NumericVector grid_nus,
                                        NumericVector grid_cmp_var_long,
                                        NumericVector grid_log_lambda_long,
                                        NumericVector grid_logZ_long,
                                        double max_mu,
                                        double min_mu,
                                        double max_nu,
                                        double min_nu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = data.ncol();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int I = betas_alpha.size();
  double grad_alpha;
  double grad_delta;
  double grad_disp;
  NumericVector grad_betas_alpha(I);
  NumericVector grad_betas_delta(I);
  NumericVector grad_betas_logdisp(I);
  NumericVector out(3 + 3*I);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alpha * nodes[k] + delta;
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += nodes[k] * betas_alpha[c] * i_cov_data(j,c) + 
          betas_delta[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      double log_disp = log(disp);
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_disp += betas_logdisp[c] * i_cov_data(j,c); // for item j
      }
      disp_interp(k,j) = exp(log_disp);
      if (disp_interp(k,j) > max_nu) { disp_interp(k,j) = max_nu; }
      if (disp_interp(k,j) < min_nu) { disp_interp(k,j) = min_nu; }
    }
  }  // end loop over items
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  grad_disp = 0;
  grad_delta = 0;
  grad_alpha = 0;
  NumericMatrix A(n_nodes, m);
  NumericMatrix B(n_nodes, m);
  // gradients for item parameters
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      // compute A and B for dispersion gradient
      double lambda = exp(log_lambda(k,i));
      A(k,i) = computeA(lambda, mu_interp(k,i), disp_interp(k,i), log_Z(k,i), 10);
      B(k,i) = computeB(lambda, mu_interp(k,i), disp_interp(k,i), log_Z(k,i), 10);
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_alpha += PPs(j,k) * (nodes[k]*mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_delta += PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_disp += PPs(j,k) * (disp_interp(k,i)*(A(k,i)*
          (data(j,i) - mu_interp(k,i))/V(k,i) - (logFactorial(data(j,i))-B(k,i))));
      }
    }
  }
  
  // gradients for item covariate weights
  for (int c=0; c<I; c++) {
    // for each gamma of which we have one for each covariate-item combination
    grad_betas_alpha[c] = 0;
    grad_betas_delta[c] = 0;
    grad_betas_logdisp[c] = 0;
    for (int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      for (int i=0; i<n; i++) {
        // over persons
        for (int j=0; j<m; j++) {
          // over items (as the betas are only specific to item covariates, not items)
          grad_betas_alpha[c] += PPs(i,k) * (nodes[k]*mu_interp(k,j)*i_cov_data(j,c) / V(k,j)) *
            (data(i,j) - mu_interp(k,j));
          grad_betas_delta[c] += PPs(i,k) * (mu_interp(k,j)*i_cov_data(j,c) / V(k,j)) *
            (data(i,j) - mu_interp(k,j));
          grad_betas_logdisp[c] += PPs(i,k) * i_cov_data(j,c) * disp_interp(k,j)*
            (A(k,j)*(data(i,j) - mu_interp(k,j))/V(k,j) - (logFactorial(data(i,j))-B(k,j)));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over P (person covariates)
  } // end loop over items
  
  // fill up output vector
  out[0] = grad_alpha;
  out[1] = grad_delta;
  out[2] = grad_disp;
  for(int c=0; c<I; c++) {
    out[3 + c] = grad_betas_alpha[c];
    out[3 + c + I] = grad_betas_delta[c];
    out[3 + c + 2*I] = grad_betas_logdisp[c];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_icov_alpha_nu_cpp(double alpha,
                                         NumericVector deltas,
                                         double disp,
                                         NumericVector betas_alpha,
                                         NumericVector betas_logdisp,
                                         NumericMatrix data,
                                         NumericMatrix i_cov_data,
                                         NumericMatrix PPs,
                                         NumericVector nodes, 
                                         NumericVector grid_mus,
                                         NumericVector grid_nus,
                                         NumericVector grid_cmp_var_long,
                                         NumericVector grid_log_lambda_long,
                                         NumericVector grid_logZ_long,
                                         double max_mu,
                                         double min_mu,
                                         double max_nu,
                                         double min_nu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = data.ncol();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int I = betas_alpha.size();
  double grad_alpha;
  NumericVector grad_deltas(m);
  double grad_disp;
  NumericVector grad_betas_alpha(I);
  NumericVector grad_betas_logdisp(I);
  NumericVector out(2 + m + 2*I);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alpha * nodes[k] + deltas[j];
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += nodes[k] * betas_alpha[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      double log_disp = log(disp);
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_disp += betas_logdisp[c] * i_cov_data(j,c); // for item j
      }
      disp_interp(k,j) = exp(log_disp);
      if (disp_interp(k,j) > max_nu) { disp_interp(k,j) = max_nu; }
      if (disp_interp(k,j) < min_nu) { disp_interp(k,j) = min_nu; }
    }
  }  // end loop over items
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  grad_disp = 0;
  grad_alpha = 0;
  NumericMatrix A(n_nodes, m);
  NumericMatrix B(n_nodes, m);
  // gradients for item parameters
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_deltas[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      // compute A and B for dispersion gradient
      double lambda = exp(log_lambda(k,i));
      A(k,i) = computeA(lambda, mu_interp(k,i), disp_interp(k,i), log_Z(k,i), 10);
      B(k,i) = computeB(lambda, mu_interp(k,i), disp_interp(k,i), log_Z(k,i), 10);
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_alpha += PPs(j,k) * (nodes[k]*mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_deltas[i] += PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_disp += PPs(j,k) * (disp_interp(k,i)*(A(k,i)*
          (data(j,i) - mu_interp(k,i))/V(k,i) - (logFactorial(data(j,i))-B(k,i))));
      }
    }
  }
  
  // gradients for item covariate weights
  for (int c=0; c<I; c++) {
    // for each gamma of which we have one for each covariate-item combination
    grad_betas_alpha[c] = 0;
    grad_betas_logdisp[c] = 0;
    for (int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      for (int i=0; i<n; i++) {
        // over persons
        for (int j=0; j<m; j++) {
          // over items (as the betas are only specific to item covariates, not items)
          grad_betas_alpha[c] += PPs(i,k) * (nodes[k]*mu_interp(k,j)*i_cov_data(j,c) / V(k,j)) *
            (data(i,j) - mu_interp(k,j));
          grad_betas_logdisp[c] += PPs(i,k) * i_cov_data(j,c) * disp_interp(k,j)*
            (A(k,j)*(data(i,j) - mu_interp(k,j))/V(k,j) - (logFactorial(data(i,j))-B(k,j)));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over P (person covariates)
  } // end loop over items
  
  // fill up output vector
  out[0] = grad_alpha;
  for (int i=0; i<m; i++) {
    out[i + 1] = grad_deltas[i];
  }
  out[m + 1] = grad_disp;
  for(int c=0; c<I; c++) {
    out[2 + m + c] = grad_betas_alpha[c];
    out[2 + m + c + I] = grad_betas_logdisp[c];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_icov_delta_nu_cpp(NumericVector alphas,
                                         double delta,
                                         double disp,
                                         NumericVector betas_delta,
                                         NumericVector betas_logdisp,
                                         NumericMatrix data,
                                         NumericMatrix i_cov_data,
                                         NumericMatrix PPs,
                                         NumericVector nodes, 
                                         NumericVector grid_mus,
                                         NumericVector grid_nus,
                                         NumericVector grid_cmp_var_long,
                                         NumericVector grid_log_lambda_long,
                                         NumericVector grid_logZ_long,
                                         double max_mu,
                                         double min_mu,
                                         double max_nu,
                                         double min_nu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = data.ncol();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int I = betas_delta.size();
  NumericVector grad_alphas(m);
  double grad_delta;
  double grad_disp;
  NumericVector grad_betas_delta(I);
  NumericVector grad_betas_logdisp(I);
  NumericVector out(m + 2 + 2*I);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alphas[j] * nodes[k] + delta;
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += betas_delta[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      double log_disp = log(disp);
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_disp += betas_logdisp[c] * i_cov_data(j,c); // for item j
      }
      disp_interp(k,j) = exp(log_disp);
      if (disp_interp(k,j) > max_nu) { disp_interp(k,j) = max_nu; }
      if (disp_interp(k,j) < min_nu) { disp_interp(k,j) = min_nu; }
    }
  }  // end loop over items
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  grad_disp = 0;
  grad_delta = 0;
  NumericMatrix A(n_nodes, m);
  NumericMatrix B(n_nodes, m);
  // gradients for item parameters
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_alphas[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      // compute A and B for dispersion gradient
      double lambda = exp(log_lambda(k,i));
      A(k,i) = computeA(lambda, mu_interp(k,i), disp_interp(k,i), log_Z(k,i), 10);
      B(k,i) = computeB(lambda, mu_interp(k,i), disp_interp(k,i), log_Z(k,i), 10);
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_alphas[i] += PPs(j,k) * (nodes[k]*mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_delta += PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_disp += PPs(j,k) * (disp_interp(k,i)*(A(k,i)*
          (data(j,i) - mu_interp(k,i))/V(k,i) - (logFactorial(data(j,i))-B(k,i))));
      }
    }
  }
  
  // gradients for item covariate weights
  for (int c=0; c<I; c++) {
    // for each gamma of which we have one for each covariate-item combination
    grad_betas_delta[c] = 0;
    grad_betas_logdisp[c] = 0;
    for (int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      for (int i=0; i<n; i++) {
        // over persons
        for (int j=0; j<m; j++) {
          // over items (as the betas are only specific to item covariates, not items)
          grad_betas_delta[c] += PPs(i,k) * (mu_interp(k,j)*i_cov_data(j,c) / V(k,j)) *
            (data(i,j) - mu_interp(k,j));
          grad_betas_logdisp[c] += PPs(i,k) * i_cov_data(j,c) * disp_interp(k,j)*
            (A(k,j)*(data(i,j) - mu_interp(k,j))/V(k,j) - (logFactorial(data(i,j))-B(k,j)));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over P (person covariates)
  } // end loop over items
  
  // fill up output vector
  for (int i=0; i<m; i++) {
    out[i] = grad_alphas[i];
  }
  out[m] = grad_delta;
  out[m + 1] = grad_disp;
  for(int c=0; c<I; c++) {
    out[m + 2 + c] = grad_betas_delta[c];
    out[m + 2 + c + I] = grad_betas_logdisp[c];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_icov_alpha_delta_cpp(double alpha,
                                         double delta,
                                         NumericVector disps,
                                         NumericVector betas_alpha,
                                         NumericVector betas_delta,
                                         NumericMatrix data,
                                         NumericMatrix i_cov_data,
                                         NumericMatrix PPs,
                                         NumericVector nodes, 
                                         NumericVector grid_mus,
                                         NumericVector grid_nus,
                                         NumericVector grid_cmp_var_long,
                                         NumericVector grid_log_lambda_long,
                                         NumericVector grid_logZ_long,
                                         double max_mu,
                                         double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = data.ncol();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int I = betas_alpha.size();
  double grad_alpha;
  double grad_delta;
  NumericVector grad_disps(m);
  NumericVector grad_betas_alpha(I);
  NumericVector grad_betas_delta(I);
  NumericVector out(2 + m + 2*I);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alpha * nodes[k] + delta;
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += nodes[k] * betas_alpha[c] * i_cov_data(j,c) + 
          betas_delta[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,j) = disps[j];
    }
  }  // end loop over items
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  grad_delta = 0;
  grad_alpha = 0;
  NumericMatrix A(n_nodes, m);
  NumericMatrix B(n_nodes, m);
  // gradients for item parameters
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_disps[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      // compute A and B for dispersion gradient
      double lambda = exp(log_lambda(k,i));
      A(k,i) = computeA(lambda, mu_interp(k,i), disp_interp(k,i), log_Z(k,i), 10);
      B(k,i) = computeB(lambda, mu_interp(k,i), disp_interp(k,i), log_Z(k,i), 10);
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_alpha += PPs(j,k) * (nodes[k]*mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_delta += PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_disps[i] += PPs(j,k) * (disp_interp(k,i)*(A(k,i)*
          (data(j,i) - mu_interp(k,i))/V(k,i) - (logFactorial(data(j,i))-B(k,i))));
      }
    }
  }
  
  // gradients for item covariate weights
  for (int c=0; c<I; c++) {
    // for each gamma of which we have one for each covariate-item combination
    grad_betas_alpha[c] = 0;
    grad_betas_delta[c] = 0;
    for (int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      for (int i=0; i<n; i++) {
        // over persons
        for (int j=0; j<m; j++) {
          // over items (as the betas are only specific to item covariates, not items)
          grad_betas_alpha[c] += PPs(i,k) * (nodes[k]*mu_interp(k,j)*i_cov_data(j,c) / V(k,j)) *
            (data(i,j) - mu_interp(k,j));
          grad_betas_delta[c] += PPs(i,k) * (mu_interp(k,j)*i_cov_data(j,c) / V(k,j)) *
            (data(i,j) - mu_interp(k,j));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over P (person covariates)
  } // end loop over items
  
  // fill up output vector
  out[0] = grad_alpha;
  out[1] = grad_delta;
  for (int i=0; i<m; i++) {
    out[i + 2] = grad_disps[i];
  }
  for(int c=0; c<I; c++) {
    out[2 + m + c] = grad_betas_alpha[c];
    out[2 + m + I] = grad_betas_delta[c];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_fixdisps_newem_cpp(NumericVector alphas,
                                          NumericVector deltas,
                                          NumericVector disps,
                                          NumericMatrix data,
                                          NumericMatrix PPs,
                                          NumericVector nodes, 
                                          NumericVector grid_mus,
                                          NumericVector grid_nus,
                                          NumericVector grid_cmp_var_long,
                                          NumericVector grid_log_lambda_long,
                                          NumericVector grid_logZ_long,
                                          double max_mu,
                                          double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = alphas.size();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  NumericVector grad_alphas(m);
  NumericVector grad_deltas(m);
  NumericVector out(2*m);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int i=0;i<m;i++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over persons (rows)
      mu(k,i) = exp(alphas[i] * nodes[k] + deltas[i]);
      mu_interp(k,i) = mu(k,i);
      if (mu(k,i) > max_mu) { mu_interp(k,i) = max_mu; }
      if (mu(k,i) < min_mu) { mu_interp(k,i) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,i) = disps[i];
    }
  }
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have persons and
  // as many columns as we have items
  
  
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_alphas[i] = 0;
    grad_deltas[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over persons (rows in my matrices)
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_alphas[i] = grad_alphas[i] +
          PPs(j,k) * (nodes[k]*mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_deltas[i] = grad_deltas[i] +
          PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
      }
    }
  }
  
  // fill up output vector
  for(int i=0;i<m;i++){
    out[i] = grad_alphas[i];
    out[i + m] = grad_deltas[i];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_pcov_fixdisps_cpp(NumericVector alphas,
                                     NumericVector deltas,
                                     NumericVector disps,
                                     NumericVector betas,
                                     NumericMatrix data,
                                     NumericMatrix p_cov_data,
                                     NumericMatrix PPs,
                                     NumericVector nodes, 
                                     NumericVector grid_mus,
                                     NumericVector grid_nus,
                                     NumericVector grid_cmp_var_long,
                                     NumericVector grid_log_lambda_long,
                                     NumericVector grid_logZ_long,
                                     double max_mu,
                                     double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = alphas.size();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int P = betas.size();
  NumericVector grad_alphas(m);
  NumericVector grad_deltas(m);
  NumericVector grad_betas(P);
  NumericVector out(2*m + P);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are person
  // as well as node and item specific
  // so we set up KxM matrices (nodesxitems) which we row-bind under each other 
  // for all N persons
  NumericMatrix mu(n_nodes*n, m);
  NumericMatrix mu_interp(n_nodes*n, m);
  NumericMatrix disp_interp(n_nodes*n, m);
  for (int i=0; i<n; i++) {
    // we are computing node-item specific mus for each person
    for(int j=0;j<m;j++){
      // loop over items (columns)
      for(int k=0;k<n_nodes;k++) {
        // loop over nodes (rows)
        double log_mu = alphas[j] * nodes[k] + deltas[j];
        for(int p=0; p<P; p++) {
          // add all the (weighted) covariate values for all covariates for the item j
          // (for the specific person i we are currently looking at)
          log_mu += betas[p] * alphas[j] * p_cov_data(i,p);
        }
        mu(k+i*n_nodes,j) = exp(log_mu);
        mu_interp(k+i*n_nodes,j) = mu(k+i*n_nodes,j);
        if (mu(k+i*n_nodes,j) > max_mu) { mu_interp(k+i*n_nodes,j) = max_mu; }
        if (mu(k+i*n_nodes,j) < min_mu) { mu_interp(k+i*n_nodes,j) = min_mu; }
        // we need to set maximum for mu to max_mu so that the interpolation will
        // work, max_mu is the maximum mu value in our grid for interpolation
        disp_interp(k+i*n_nodes,j) = disps[j];
      }
    }  // end loop over items
  } // end loop over N
  
  NumericMatrix V(n_nodes*n, m);
  NumericMatrix log_lambda(n_nodes*n, m);
  NumericMatrix log_Z(n_nodes*n, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes*persons and
  // as many columns as we have (same as mu_interp and nu_interp matrices)
  
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_alphas[i] = 0;
    grad_deltas[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      // note the index for mu_interp and V_interp (and log_lambda and log_Z) must be adjusted 
      // as we now have (additionally) person specific mus and Vs
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the sum over the weightes covariates for the gradient for alpha
        double sum_over_pcov = 0;
        for (int p=0; p<P; p++) {
          sum_over_pcov += betas[p] * p_cov_data(j,p);
        }
        
        // compute the gradients (summing over persons)
        grad_alphas[i] = grad_alphas[i] +
          PPs(j,k) * (mu_interp(k+j*n_nodes,i)*(nodes[k] + sum_over_pcov) /  V(k+j*n_nodes,i))*
          (data(j,i) -  mu_interp(k+j*n_nodes,i));
        grad_deltas[i] = grad_deltas[i] +
          PPs(j,k) * (mu_interp(k+j*n_nodes,i) / V(k+j*n_nodes,i))*(data(j,i) - 
          mu_interp(k+j*n_nodes,i));
      }
    }
  }
  
  // gradients for person covariate weights
  for (int p=0; p<P; p++) { // over covariates
    grad_betas[p] = 0;
    for (int j=0; j<m; j++) { // over items
      for (int k=0;k<n_nodes;k++) { // over nodes (rows in my matrices)
        for (int i=0; i<n; i++) { // over persons
          grad_betas[p] += PPs(i,k) * (mu_interp(k+i*n_nodes,j) * alphas[j] * p_cov_data(i,p) / V(k+i*n_nodes,j)) *
            (data(i,j) - mu_interp(k+i*n_nodes,j));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over items
  } //end loop over P (person covariates)
  
  // fill up output vector
  for(int i=0;i<m;i++){
    out[i] = grad_alphas[i];
    out[i + m] = grad_deltas[i];
  }
  for(int p=0; p<P; p++) {
    out[2*m + p] = grad_betas[p];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_pcov_cat_fixdisps_cpp(NumericVector alphas,
                                              NumericVector deltas,
                                              NumericVector disps,
                                              NumericVector betas,
                                              NumericMatrix data,
                                              NumericMatrix p_cov_data,
                                              NumericMatrix resp_pattern,
                                              NumericMatrix PPs,
                                              NumericVector nodes, 
                                              NumericVector grid_mus,
                                              NumericVector grid_nus,
                                              NumericVector grid_cmp_var_long,
                                              NumericVector grid_log_lambda_long,
                                              NumericVector grid_logZ_long,
                                              double max_mu,
                                              double min_mu) {
  
  // assume that p_cov is a matrix of dummy coded categorical predictors
  // resp_pattern is a matrix of the same no. of cols than p_cov
  // and as many rows as we have distinct possible response patterns
  
  int N = data.nrow();
  int M = data.ncol();
  int K = nodes.size();
  int P = betas.size(); 
  int n_resp_patterns = resp_pattern.nrow();
  NumericVector grad_alphas(M);
  NumericVector grad_deltas(M);
  NumericVector grad_betas(P);
  NumericVector out(2*M + P);
  
  // for person covariates, we need mus (and lambdas and Zs) for each node and
  // and then also for each response pattern
  // so first compute that
  NumericMatrix mu(K*n_resp_patterns, M);
  NumericMatrix mu_interp(K*n_resp_patterns, M);
  NumericMatrix disp_interp(K*n_resp_patterns, M);
  for (int l=0; l<n_resp_patterns; l++) {
    for(int j=0; j<M; j++){
      // loop over items (columns)
      for(int k=0; k<K; k++) {
        // loop over nodes (rows)
        double log_mu = alphas[j] * nodes[k] + deltas[j];
        for(int p=0; p<P; p++) {
          // this works because only includes columns for none-reference categories
          // for all covs in ref categories, resp_pattern will just always be zero in that row
          log_mu += betas[p] * alphas[j] * resp_pattern(l,p);
        }
        
        mu(k+l*K,j) = exp(log_mu);
        mu_interp(k+l*K,j) = mu(k+l*K,j);
        if (mu(k+l*K,j) > max_mu) { mu_interp(k+l*K,j) = max_mu; }
        if (mu(k+l*K,j) < min_mu) { mu_interp(k+l*K,j) = min_mu; }
        // we need to set maximum for mu to max_mu so that the interpolation will
        // work, max_mu is the maximum mu value in our grid for interpolation
        disp_interp(k+l*K,j) = disps[j];
      }
    }  // end loop over items
  }
  
  NumericMatrix log_lambda(K*n_resp_patterns, M);
  NumericMatrix log_Z(K*n_resp_patterns, M);
  NumericMatrix V(K*n_resp_patterns, M);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  
  for(int i=0;i<M;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_alphas[i] = 0;
    grad_deltas[i] = 0;
    
    for(int k=0; k<K; k++) {
      // over nodes (rows in my matrices)
      
      for(int j=0; j<N; j++) {
        // loop over persons
        
        // check what response pattern person i had
        int l = 0;
        bool pattern_match = false;
        while (!pattern_match && l<n_resp_patterns) {
          // the second condition is just for safety that we dont get an eternal while loop
          // but we should always find a pattern match
          for (int p=0; p<P; p++) {
            pattern_match = p_cov_data(j,p) == resp_pattern(l,p);
          }
          // if the rows are the same, i am going to get out the foor loop with
          // pattern_match = TRUE and have l at the row of the pattern matrix
          // otherwise I am going to increase l by 1 and stay in the while loop to see if
          // the next row in the pattern matrix is a match for i
          if (!pattern_match) { l += 1; }
        }
        
        // we now know that person i has a response pattern like in row l of resp_pattern matrix
        // so their mu (and lambda, etc.) should be at row k+l*K
        
        // compute the sum over the weightes covariates for the gradient for alpha
        double sum_over_pcov = 0;
        for (int p=0; p<P; p++) {
          sum_over_pcov += betas[p] * p_cov_data(j,p);
        }
        
        // compute the gradients (summing over persons)
        grad_alphas[i] += PPs(j,k) * (mu_interp(k+l*K,i)*(nodes[k] + sum_over_pcov) /  V(k+l*K,i))*
          (data(j,i) -  mu_interp(k+l*K,i));
        grad_deltas[i] += PPs(j,k) * (mu_interp(k+l*K,i) / V(k+l*K,i))*(data(j,i) - 
          mu_interp(k+l*K,i));
      }
    }
  }
  
  // gradients for person covariate weights
  for (int p=0; p<P; p++) { // over covariates
    grad_betas[p] = 0;
    for (int j=0; j<M; j++) { // over items
      for (int k=0; k<K; k++) { // over nodes (rows in my matrices)
        for (int i=0; i<N; i++) { // over persons
          // check what response pattern person i had
          int l = 0;
          bool pattern_match = false;
          while (!pattern_match && l<n_resp_patterns) {
            // the second condition is just for safety that we dont get an eternal while loop
            // but we should always find a pattern match
            for (int p=0; p<P; p++) {
              pattern_match = p_cov_data(i,p) == resp_pattern(l,p);
            }
            // if the rows are the same, i am going to get out the foor loop with
            // pattern_match = TRUE and have l at the row of the pattern matrix
            // otherwise I am going to increase l by 1 and stay in the while loop to see if
            // the next row in the pattern matrix is a match for i
            if (!pattern_match) { l += 1; }
          }
          
          // we now know that person i has a response pattern like in row l of resp_pattern matrix
          // so their mu (and lambda, etc.) should be at row k+l*K
          
          grad_betas[p] += PPs(i,k) * (mu_interp(k+l*K,j) * alphas[j] * p_cov_data(i,p) / V(k+l*K,j)) *
            (data(i,j) - mu_interp(k+l*K,j));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over items
  } //end loop over P (person covariates)
  
  // fill up output vector
  for(int i=0;i<M;i++){
    out[i] = grad_alphas[i];
    out[i + M] = grad_deltas[i];
  }
  for(int p=0; p<P; p++) {
    out[2*M + p] = grad_betas[p];
  }
  
  return(out);
}


// [[Rcpp::export]]
NumericVector grad_cmp_with_icov_delta_fixdisps_cpp(NumericVector alphas,
                                     double delta,
                                     NumericVector disps,
                                     NumericVector betas,
                                     NumericMatrix data,
                                     NumericMatrix i_cov_data,
                                     NumericMatrix PPs,
                                     NumericVector nodes, 
                                     NumericVector grid_mus,
                                     NumericVector grid_nus,
                                     NumericVector grid_cmp_var_long,
                                     NumericVector grid_log_lambda_long,
                                     NumericVector grid_logZ_long,
                                     double max_mu,
                                     double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = alphas.size();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int I = betas.size();
  NumericVector grad_alphas(m);
  double grad_delta;
  NumericVector grad_betas(I);
  NumericVector out(m + 1 + I);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for item covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alphas[j] * nodes[k] + delta;
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += betas[c] * i_cov_data(j,c); // for item j
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,j) = disps[j];
    }
  }  // end loop over items
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  grad_delta = 0;
  // gradients for item parameters
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_alphas[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_alphas[i] = grad_alphas[i] +
          PPs(j,k) * (nodes[k]*mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_delta += PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
      }
    }
  }
  
  // gradients for item covariate weights
  for (int c=0; c<I; c++) {
    // for each gamma of which we have one for each covariate-item combination
    grad_betas[c] = 0;
    for (int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      for (int i=0; i<n; i++) {
        // over persons
        for (int j=0; j<m; j++) {
          // over items (as the betas are only specific to item covariates, not items)
          grad_betas[c] += PPs(i,k) * (mu_interp(k,j)*i_cov_data(j,c) / V(k,j)) *
            (data(i,j) - mu_interp(k,j));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over P (person covariates)
  } // end loop over items
  
  // fill up output vector
  for(int i=0;i<m;i++){
    out[i] = grad_alphas[i];
  }
  out[m] = grad_delta;
  for(int c=0; c<I; c++) {
    out[m + 1 + c] = grad_betas[c];
  }
  
  return(out);
}


// [[Rcpp::export]]
NumericVector grad_cmp_with_icov_alpha_fixdisps_cpp(double alpha,
                                                    NumericVector deltas,
                                                    NumericVector disps,
                                                    NumericVector betas,
                                                    NumericMatrix data,
                                                    NumericMatrix i_cov_data,
                                                    NumericMatrix PPs,
                                                    NumericVector nodes, 
                                                    NumericVector grid_mus,
                                                    NumericVector grid_nus,
                                                    NumericVector grid_cmp_var_long,
                                                    NumericVector grid_log_lambda_long,
                                                    NumericVector grid_logZ_long,
                                                    double max_mu,
                                                    double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = data.ncol();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int I = betas.size();
  double grad_alpha;
  NumericVector grad_deltas(m);
  NumericVector grad_betas(I);
  NumericVector out(m + 1 + I);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for item covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alpha * nodes[k] + deltas[j];
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += nodes[k] * betas[c] * i_cov_data(j,c); // for item j
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,j) = disps[j];
    }
  }  // end loop over items
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  grad_alpha = 0;
  // gradients for item parameters
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_deltas[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_alpha += PPs(j,k) * (nodes[k]*mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_deltas[i] += PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
      }
    }
  }
  
  // gradients for item covariate weights
  for (int c=0; c<I; c++) {
    // for each gamma of which we have one for each covariate-item combination
    grad_betas[c] = 0;
    for (int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      for (int i=0; i<n; i++) {
        // over persons
        for (int j=0; j<m; j++) {
          // over items (as the betas are only specific to item covariates, not items)
          grad_betas[c] += PPs(i,k) * (nodes[k]*mu_interp(k,j)*i_cov_data(j,c) / V(k,j)) *
            (data(i,j) - mu_interp(k,j));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over P (person covariates)
  } // end loop over items
  
  // fill up output vector
  out[0] = grad_alpha;
  for(int i=0;i<m;i++){
    out[i+1] = grad_deltas[i];
  }
  for(int c=0; c<I; c++) {
    out[m + 1 + c] = grad_betas[c];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_icov_alpha_delta_fixdisps_cpp(double alpha,
                                                           double delta,
                                                           NumericVector disps,
                                                           NumericVector betas_alpha,
                                                           NumericVector betas_delta,
                                                           NumericMatrix data,
                                                           NumericMatrix i_cov_data,
                                                           NumericMatrix PPs,
                                                           NumericVector nodes, 
                                                           NumericVector grid_mus,
                                                           NumericVector grid_nus,
                                                           NumericVector grid_cmp_var_long,
                                                           NumericVector grid_log_lambda_long,
                                                           NumericVector grid_logZ_long,
                                                           double max_mu,
                                                           double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = data.ncol();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int I = betas_alpha.size();
  double grad_alpha;
  double grad_delta;
  NumericVector grad_betas_alpha(I);
  NumericVector grad_betas_delta(I);
  NumericVector out(2 + 2*I);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alpha * nodes[k] + delta;
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += nodes[k] * betas_alpha[c] * i_cov_data(j,c) + 
          betas_delta[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,j) = disps[j];
    }
  }  // end loop over items
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  grad_delta = 0;
  grad_alpha = 0;
  // gradients for item parameters
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_alpha += PPs(j,k) * (nodes[k]*mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_delta += PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
      }
    }
  }
  
  // gradients for item covariate weights
  for (int c=0; c<I; c++) {
    // for each gamma of which we have one for each covariate-item combination
    grad_betas_alpha[c] = 0;
    grad_betas_delta[c] = 0;
    for (int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      for (int i=0; i<n; i++) {
        // over persons
        for (int j=0; j<m; j++) {
          // over items (as the betas are only specific to item covariates, not items)
          grad_betas_alpha[c] += PPs(i,k) * (nodes[k]*mu_interp(k,j)*i_cov_data(j,c) / V(k,j)) *
            (data(i,j) - mu_interp(k,j));
          grad_betas_delta[c] += PPs(i,k) * (mu_interp(k,j)*i_cov_data(j,c) / V(k,j)) *
            (data(i,j) - mu_interp(k,j));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over P (person covariates)
  } // end loop over items
  
  // fill up output vector
  out[0] = grad_alpha;
  out[1] = grad_delta;
  for(int c=0; c<I; c++) {
    out[2 + c] = grad_betas_alpha[c];
    out[2 + c + I] = grad_betas_delta[c];
  }
  
  return(out);
}


// [[Rcpp::export]]
NumericVector grad_cmp_fixalphas_newem_cpp(NumericVector alphas,
                                  NumericVector deltas,
                                  NumericVector disps,
                                  NumericMatrix data,
                                  NumericMatrix PPs,
                                  NumericVector nodes, 
                                  NumericVector grid_mus,
                                  NumericVector grid_nus,
                                  NumericVector grid_cmp_var_long,
                                  NumericVector grid_log_lambda_long,
                                  NumericVector grid_logZ_long,
                                  double max_mu,
                                  double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = alphas.size();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  NumericVector grad_deltas(m);
  NumericVector grad_disps(m);
  NumericVector out(2*m);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int i=0;i<m;i++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over persons (rows)
      mu(k,i) = exp(alphas[i] * nodes[k] + deltas[i]);
      mu_interp(k,i) = mu(k,i);
      if (mu(k,i) > max_mu) { mu_interp(k,i) = max_mu; }
      if (mu(k,i) < min_mu) { mu_interp(k,i) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,i) = disps[i];
    }
  }
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have persons and
  // as many columns as we have items
  
  
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_deltas[i] = 0;
    grad_disps[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over persons (rows in my matrices)
      
      // compute A and B for dispersion gradient
      double lambda = exp(log_lambda(k,i));
      double A = computeA(lambda, mu_interp(k,i), disps[i], log_Z(k,i), 10);
      double B = computeB(lambda, mu_interp(k,i), disps[i], log_Z(k,i), 10);
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_deltas[i] = grad_deltas[i] +
          PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_disps[i] = grad_disps[i] +
          PPs(j,k) * (disps[i]*(A*(data(j,i) - mu_interp(k,i))/V(k,i) - (logFactorial(data(j,i))-B)));
      }
    }
  }
  
  // fill up output vector
  for(int i=0;i<m;i++){
    out[i] = grad_deltas[i];
    out[i + m] = grad_disps[i];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_pcov_fixalphas_cpp(NumericVector alphas,
                                     NumericVector deltas,
                                     NumericVector disps,
                                     NumericVector betas,
                                     NumericMatrix data,
                                     NumericMatrix p_cov_data,
                                     NumericMatrix PPs,
                                     NumericVector nodes, 
                                     NumericVector grid_mus,
                                     NumericVector grid_nus,
                                     NumericVector grid_cmp_var_long,
                                     NumericVector grid_log_lambda_long,
                                     NumericVector grid_logZ_long,
                                     double max_mu,
                                     double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = alphas.size();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int P = betas.size();
  NumericVector grad_deltas(m);
  NumericVector grad_disps(m);
  NumericVector grad_betas(P);
  NumericVector out(2*m + P);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are person
  // as well as node and item specific
  // so we set up KxM matrices (nodesxitems) which we row-bind under each other 
  // for all N persons
  NumericMatrix mu(n_nodes*n, m);
  NumericMatrix mu_interp(n_nodes*n, m);
  NumericMatrix disp_interp(n_nodes*n, m);
  for (int i=0; i<n; i++) {
    // we are computing node-item specific mus for each person
    for(int j=0;j<m;j++){
      // loop over items (columns)
      for(int k=0;k<n_nodes;k++) {
        // loop over nodes (rows)
        double log_mu = alphas[j] * nodes[k] + deltas[j];
        for(int p=0; p<P; p++) {
          // add all the (weighted) covariate values for all covariates for the item j
          // (for the specific person i we are currently looking at)
          log_mu += betas[p] * alphas[j] * p_cov_data(i,p);
        }
        mu(k+i*n_nodes,j) = exp(log_mu);
        mu_interp(k+i*n_nodes,j) = mu(k+i*n_nodes,j);
        if (mu(k+i*n_nodes,j) > max_mu) { mu_interp(k+i*n_nodes,j) = max_mu; }
        if (mu(k+i*n_nodes,j) < min_mu) { mu_interp(k+i*n_nodes,j) = min_mu; }
        // we need to set maximum for mu to max_mu so that the interpolation will
        // work, max_mu is the maximum mu value in our grid for interpolation
        disp_interp(k+i*n_nodes,j) = disps[j];
      }
    }  // end loop over items
  } // end loop over N
  
  NumericMatrix V(n_nodes*n, m);
  NumericMatrix log_lambda(n_nodes*n, m);
  NumericMatrix log_Z(n_nodes*n, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes*persons and
  // as many columns as we have (same as mu_interp and nu_interp matrices)
  
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_deltas[i] = 0;
    grad_disps[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      // note the index for mu_interp and V_interp (and log_lambda and log_Z) must be adjusted 
      // as we now have (additionally) person specific mus and Vs
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute A and B for dispersion gradient
        double lambda = exp(log_lambda(k+j*n_nodes,i));
        double A = computeA(lambda, mu_interp(k+j*n_nodes,i), disps[i], log_Z(k+j*n_nodes,i), 10);
        double B = computeB(lambda, mu_interp(k+j*n_nodes,i), disps[i], log_Z(k+j*n_nodes,i), 10);
        
        grad_deltas[i] = grad_deltas[i] +
          PPs(j,k) * (mu_interp(k+j*n_nodes,i) / V(k+j*n_nodes,i))*(data(j,i) - 
          mu_interp(k+j*n_nodes,i));
        grad_disps[i] = grad_disps[i] +
          PPs(j,k) * (disps[i]*(A*(data(j,i) - mu_interp(k+j*n_nodes,i))/V(k+j*n_nodes,i) - 
          (logFactorial(data(j,i))-B)));
      }
    }
  }
  
  // gradients for person covariate weights
  for (int p=0; p<P; p++) { // over covariates
    grad_betas[p] = 0;
    for (int j=0; j<m; j++) { // over items
      for (int k=0;k<n_nodes;k++) { // over nodes (rows in my matrices)
        for (int i=0; i<n; i++) { // over persons
          grad_betas[p] += PPs(i,k) * (mu_interp(k+i*n_nodes,j) * alphas[j] * p_cov_data(i,p) / V(k+i*n_nodes,j)) *
            (data(i,j) - mu_interp(k+i*n_nodes,j));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over items
  } //end loop over P (person covariates)
  
  // fill up output vector
  for(int i=0;i<m;i++){
    out[i] = grad_deltas[i];
    out[i + m] = grad_disps[i];
  }
  for(int p=0; p<P; p++) {
    out[2*m + p] = grad_betas[p];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_pcov_cat_fixalphas_cpp(NumericVector alphas,
                                               NumericVector deltas,
                                               NumericVector disps,
                                               NumericVector betas,
                                               NumericMatrix data,
                                               NumericMatrix p_cov_data,
                                               NumericMatrix resp_pattern,
                                               NumericMatrix PPs,
                                               NumericVector nodes, 
                                               NumericVector grid_mus,
                                               NumericVector grid_nus,
                                               NumericVector grid_cmp_var_long,
                                               NumericVector grid_log_lambda_long,
                                               NumericVector grid_logZ_long,
                                               double max_mu,
                                               double min_mu) {
  
  // assume that p_cov is a matrix of dummy coded categorical predictors
  // resp_pattern is a matrix of the same no. of cols than p_cov
  // and as many rows as we have distinct possible response patterns
  
  int N = data.nrow();
  int M = data.ncol();
  int K = nodes.size();
  int P = betas.size(); 
  int n_resp_patterns = resp_pattern.nrow();
  NumericVector grad_deltas(M);
  NumericVector grad_disps(M);
  NumericVector grad_betas(P);
  NumericVector out(2*M + P);
  
  // for person covariates, we need mus (and lambdas and Zs) for each node and
  // and then also for each response pattern
  // so first compute that
  NumericMatrix mu(K*n_resp_patterns, M);
  NumericMatrix mu_interp(K*n_resp_patterns, M);
  NumericMatrix disp_interp(K*n_resp_patterns, M);
  for (int l=0; l<n_resp_patterns; l++) {
    for(int j=0; j<M; j++){
      // loop over items (columns)
      for(int k=0; k<K; k++) {
        // loop over nodes (rows)
        double log_mu = alphas[j] * nodes[k] + deltas[j];
        for(int p=0; p<P; p++) {
          // this works because only includes columns for none-reference categories
          // for all covs in ref categories, resp_pattern will just always be zero in that row
          log_mu += betas[p] * alphas[j] * resp_pattern(l,p);
        }
        
        mu(k+l*K,j) = exp(log_mu);
        mu_interp(k+l*K,j) = mu(k+l*K,j);
        if (mu(k+l*K,j) > max_mu) { mu_interp(k+l*K,j) = max_mu; }
        if (mu(k+l*K,j) < min_mu) { mu_interp(k+l*K,j) = min_mu; }
        // we need to set maximum for mu to max_mu so that the interpolation will
        // work, max_mu is the maximum mu value in our grid for interpolation
        disp_interp(k+l*K,j) = disps[j];
      }
    }  // end loop over items
  }
  
  NumericMatrix log_lambda(K*n_resp_patterns, M);
  NumericMatrix log_Z(K*n_resp_patterns, M);
  NumericMatrix V(K*n_resp_patterns, M);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  
  for(int i=0; i<M; i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_deltas[i] = 0;
    grad_disps[i] = 0;
    
    for(int k=0; k<K; k++) {
      // over nodes (rows in my matrices)
      
      // note the index for mu_interp and V_interp (and log_lambda and log_Z) must be adjusted 
      // as we now have (additionally) person specific mus and Vs
      for(int j=0;j<N;j++) {
        // loop over persons
        
        // check what response pattern person j had
        int l = 0;
        bool pattern_match = false;
        while (!pattern_match && l<n_resp_patterns) {
          // the second condition is just for safety that we dont get an eternal while loop
          // but we should always find a pattern match
          for (int p=0; p<P; p++) {
            pattern_match = p_cov_data(j,p) == resp_pattern(l,p);
          }
          // if the rows are the same, i am going to get out the foor loop with
          // pattern_match = TRUE and have l at the row of the pattern matrix
          // otherwise I am going to increase l by 1 and stay in the while loop to see if
          // the next row in the pattern matrix is a match for i
          if (!pattern_match) { l += 1; }
        }
        
        // we now know that person i has a response pattern like in row l of resp_pattern matrix
        // so their mu (and lambda, etc.) should be at row k+l*K
        
        // compute A and B for dispersion gradient
        double lambda = exp(log_lambda(k+l*K,i));
        double A = computeA(lambda, mu_interp(k+l*K,i), disps[i], log_Z(k+l*K,i), 10);
        double B = computeB(lambda, mu_interp(k+l*K,i), disps[i], log_Z(k+l*K,i), 10);
        
        grad_deltas[i] += PPs(j,k) * (mu_interp(k+l*K,i) / V(k+l*K,i))*(data(j,i) - 
          mu_interp(k+l*K,i));
        grad_disps[i] += PPs(j,k) * (disps[i]*(A*(data(j,i) - mu_interp(k+l*K,i))/V(k+l*K,i) - 
          (logFactorial(data(j,i))-B)));
      }
    }
  }
  
  // gradients for person covariate weights
  for (int p=0; p<P; p++) { // over covariates
    grad_betas[p] = 0;
    for (int j=0; j<M; j++) { // over items
      for (int k=0;k<K;k++) { // over nodes (rows in my matrices)
        for (int i=0; i<N; i++) { // over persons
          
          // check what response pattern person j had
          int l = 0;
          bool pattern_match = false;
          while (!pattern_match && l<n_resp_patterns) {
            // the second condition is just for safety that we dont get an eternal while loop
            // but we should always find a pattern match
            for (int p=0; p<P; p++) {
              pattern_match = p_cov_data(i,p) == resp_pattern(l,p);
            }
            // if the rows are the same, i am going to get out the foor loop with
            // pattern_match = TRUE and have l at the row of the pattern matrix
            // otherwise I am going to increase l by 1 and stay in the while loop to see if
            // the next row in the pattern matrix is a match for i
            if (!pattern_match) { l += 1; }
          }
          
          grad_betas[p] += PPs(i,k) * (mu_interp(k+l*K,j) * alphas[j] * p_cov_data(i,p) / V(k+l*K,j)) *
            (data(i,j) - mu_interp(k+l*K,j));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over items
  } //end loop over P (person covariates)
  
  // fill up output vector
  for(int i=0;i<M;i++){
    out[i] = grad_deltas[i];
    out[i + M] = grad_disps[i];
  }
  for(int p=0; p<P; p++) {
    out[2*M + p] = grad_betas[p];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_icov_delta_fixalphas_cpp(NumericVector alphas,
                                     double delta,
                                     NumericVector disps,
                                     NumericVector betas,
                                     NumericMatrix data,
                                     NumericMatrix i_cov_data,
                                     NumericMatrix PPs,
                                     NumericVector nodes, 
                                     NumericVector grid_mus,
                                     NumericVector grid_nus,
                                     NumericVector grid_cmp_var_long,
                                     NumericVector grid_log_lambda_long,
                                     NumericVector grid_logZ_long,
                                     double max_mu,
                                     double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = alphas.size();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int I = betas.size();
  double grad_delta;
  NumericVector grad_disps(m);
  NumericVector grad_betas(I);
  NumericVector out(1 + m + I);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alphas[j] * nodes[k] + delta;
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += betas[c] * i_cov_data(j,c); // for item j
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,j) = disps[j];
    }
  }  // end loop over items
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  grad_delta = 0;
  // gradients for item parameters
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_disps[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      // compute A and B for dispersion gradient
      double lambda = exp(log_lambda(k,i));
      double A = computeA(lambda, mu_interp(k,i), disps[i], log_Z(k,i), 10);
      double B = computeB(lambda, mu_interp(k,i), disps[i], log_Z(k,i), 10);
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_delta += PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_disps[i] = grad_disps[i] +
          PPs(j,k) * (disps[i]*(A*(data(j,i) - mu_interp(k,i))/V(k,i) - (logFactorial(data(j,i))-B)));
      }
    }
  }
  
  // gradients for item covariate weights
  for (int c=0; c<I; c++) {
    // for each gamma of which we have one for each covariate-item combination
    grad_betas[c] = 0;
    for (int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      for (int i=0; i<n; i++) {
        // over persons
        for (int j=0; j<m; j++) {
          // over items (as the betas are only specific to item covariates, not items)
          grad_betas[c] += PPs(i,k) * (mu_interp(k,j)*i_cov_data(j,c) / V(k,j)) *
            (data(i,j) - mu_interp(k,j));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over P (person covariates)
  } // end loop over items
  
  // fill up output vector
  out[0] = grad_delta;
  for(int i=0;i<m;i++){
    out[i + 1] = grad_disps[i];
  }
  for(int c=0; c<I; c++) {
    out[m + 1 + c] = grad_betas[c];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_icov_nu_fixalphas_cpp(NumericVector alphas,
                                                     NumericVector deltas,
                                                     double disp,
                                                     NumericVector betas,
                                                     NumericMatrix data,
                                                     NumericMatrix i_cov_data,
                                                     NumericMatrix PPs,
                                                     NumericVector nodes, 
                                                     NumericVector grid_mus,
                                                     NumericVector grid_nus,
                                                     NumericVector grid_cmp_var_long,
                                                     NumericVector grid_log_lambda_long,
                                                     NumericVector grid_logZ_long,
                                                     double max_mu,
                                                     double min_mu,
                                                     double max_nu,
                                                     double min_nu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = data.ncol();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int I = betas.size();
  NumericVector grad_deltas(m);
  double grad_disp;
  NumericVector grad_betas(I);
  NumericVector out(1 + m + I);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alphas[j] * nodes[k] + deltas[j];
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      double log_disp = log(disp);
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_disp += betas[c] * i_cov_data(j,c); // for item j
      }
      disp_interp(k,j) = exp(log_disp);
      if (disp_interp(k,j) > max_nu) { disp_interp(k,j) = max_nu; }
      if (disp_interp(k,j) < min_nu) { disp_interp(k,j) = min_nu; }
    }
  }  // end loop over items
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  grad_disp = 0;
  NumericMatrix A(n_nodes, m);
  NumericMatrix B(n_nodes, m);
  // gradients for item parameters
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_deltas[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      // compute A and B for dispersion gradient
      double lambda = exp(log_lambda(k,i));
      A(k,i) = computeA(lambda, mu_interp(k,i), disp_interp(k,i), log_Z(k,i), 10);
      B(k,i) = computeB(lambda, mu_interp(k,i), disp_interp(k,i), log_Z(k,i), 10);
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_deltas[i] += PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_disp += PPs(j,k) * (disp_interp(k,i)*(A(k,i)*
          (data(j,i) - mu_interp(k,i))/V(k,i) - (logFactorial(data(j,i))-B(k,i))));
      }
    }
  }
  
  // gradients for item covariate weights
  for (int c=0; c<I; c++) {
    // for each gamma of which we have one for each covariate-item combination
    grad_betas[c] = 0;
    for (int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      for (int i=0; i<n; i++) {
        // over persons
        for (int j=0; j<m; j++) {
          // over items (as the betas are only specific to item covariates, not items)
          grad_betas[c] += PPs(i,k) * i_cov_data(j,c) * disp_interp(k,j)*
            (A(k,j)*(data(i,j) - mu_interp(k,j))/V(k,j) - (logFactorial(data(i,j))-B(k,j)));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over P (person covariates)
  } // end loop over items
  
  // fill up output vector
  for(int i=0;i<m;i++){
    out[i] = grad_deltas[i];
  }
  out[m] = grad_disp;
  for(int c=0; c<I; c++) {
    out[m + 1 + c] = grad_betas[c];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_icov_delta_nu_fixalphas_cpp(NumericVector alphas,
                                                         double delta,
                                                         double disp,
                                                         NumericVector betas_delta,
                                                         NumericVector betas_logdisp,
                                                         NumericMatrix data,
                                                         NumericMatrix i_cov_data,
                                                         NumericMatrix PPs,
                                                         NumericVector nodes, 
                                                         NumericVector grid_mus,
                                                         NumericVector grid_nus,
                                                         NumericVector grid_cmp_var_long,
                                                         NumericVector grid_log_lambda_long,
                                                         NumericVector grid_logZ_long,
                                                         double max_mu,
                                                         double min_mu,
                                                         double max_nu,
                                                         double min_nu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = data.ncol();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int I = betas_delta.size();
  double grad_delta;
  double grad_disp;
  NumericVector grad_betas_delta(I);
  NumericVector grad_betas_logdisp(I);
  NumericVector out(2 + 2*I);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alphas[j] * nodes[k] + delta;
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += betas_delta[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      double log_disp = log(disp);
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_disp += betas_logdisp[c] * i_cov_data(j,c); // for item j
      }
      disp_interp(k,j) = exp(log_disp);
      if (disp_interp(k,j) > max_nu) { disp_interp(k,j) = max_nu; }
      if (disp_interp(k,j) < min_nu) { disp_interp(k,j) = min_nu; }
    }
  }  // end loop over items
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  grad_disp = 0;
  grad_delta = 0;
  NumericMatrix A(n_nodes, m);
  NumericMatrix B(n_nodes, m);
  // gradients for item parameters
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      // compute A and B for dispersion gradient
      double lambda = exp(log_lambda(k,i));
      A(k,i) = computeA(lambda, mu_interp(k,i), disp_interp(k,i), log_Z(k,i), 10);
      B(k,i) = computeB(lambda, mu_interp(k,i), disp_interp(k,i), log_Z(k,i), 10);
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_delta += PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_disp += PPs(j,k) * (disp_interp(k,i)*(A(k,i)*
          (data(j,i) - mu_interp(k,i))/V(k,i) - (logFactorial(data(j,i))-B(k,i))));
      }
    }
  }
  
  // gradients for item covariate weights
  for (int c=0; c<I; c++) {
    // for each gamma of which we have one for each covariate-item combination
    grad_betas_delta[c] = 0;
    grad_betas_logdisp[c] = 0;
    for (int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      for (int i=0; i<n; i++) {
        // over persons
        for (int j=0; j<m; j++) {
          // over items (as the betas are only specific to item covariates, not items)
          grad_betas_delta[c] += PPs(i,k) * (mu_interp(k,j)*i_cov_data(j,c) / V(k,j)) *
            (data(i,j) - mu_interp(k,j));
          grad_betas_logdisp[c] += PPs(i,k) * i_cov_data(j,c) * disp_interp(k,j)*
            (A(k,j)*(data(i,j) - mu_interp(k,j))/V(k,j) - (logFactorial(data(i,j))-B(k,j)));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over P (person covariates)
  } // end loop over items
  
  // fill up output vector
  out[0] = grad_delta;
  out[1] = grad_disp;
  for(int c=0; c<I; c++) {
    out[2 + c] = grad_betas_delta[c];
    out[2 + c + I] = grad_betas_logdisp[c];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_samedisps_newem_cpp(NumericVector alphas,
                                  NumericVector deltas,
                                  NumericVector disps,
                                  NumericMatrix data,
                                  NumericMatrix PPs,
                                  NumericVector nodes, 
                                  NumericVector grid_mus,
                                  NumericVector grid_nus,
                                  NumericVector grid_cmp_var_long,
                                  NumericVector grid_log_lambda_long,
                                  NumericVector grid_logZ_long,
                                  double max_mu,
                                  double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = alphas.size();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  NumericVector grad_alphas(m);
  NumericVector grad_deltas(m);
  double grad_disp;
  NumericVector out(2*m + 1);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int i=0;i<m;i++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over persons (rows)
      mu(k,i) = exp(alphas[i] * nodes[k] + deltas[i]);
      mu_interp(k,i) = mu(k,i);
      if (mu(k,i) > max_mu) { mu_interp(k,i) = max_mu; }
      if (mu(k,i) < min_mu) { mu_interp(k,i) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,i) = disps[i];
    }
  }
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have persons and
  // as many columns as we have items
  
  grad_disp = 0;
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_alphas[i] = 0;
    grad_deltas[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over persons (rows in my matrices)
      
      // compute A and B for dispersion gradient
      double lambda = exp(log_lambda(k,i));
      double A = computeA(lambda, mu_interp(k,i), disps[i], log_Z(k,i), 10);
      double B = computeB(lambda, mu_interp(k,i), disps[i], log_Z(k,i), 10);
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_alphas[i] = grad_alphas[i] +
          PPs(j,k) * (nodes[k]*mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_deltas[i] = grad_deltas[i] +
          PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_disp = grad_disp +
          PPs(j,k) * (disps[i]*(A*(data(j,i) - mu_interp(k,i))/V(k,i) - (logFactorial(data(j,i))-B)));
      }
    }
  }
  
  // fill up output vector
  for(int i=0;i<m;i++){
    out[i] = grad_alphas[i];
    out[i + m] = grad_deltas[i];
  }
  out[2*m] = grad_disp;
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_pcov_samedisps_cpp(NumericVector alphas,
                                     NumericVector deltas,
                                     NumericVector disps,
                                     NumericVector betas,
                                     NumericMatrix data,
                                     NumericMatrix p_cov_data,
                                     NumericMatrix PPs,
                                     NumericVector nodes, 
                                     NumericVector grid_mus,
                                     NumericVector grid_nus,
                                     NumericVector grid_cmp_var_long,
                                     NumericVector grid_log_lambda_long,
                                     NumericVector grid_logZ_long,
                                     double max_mu,
                                     double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = alphas.size();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int P = betas.size();
  NumericVector grad_alphas(m);
  NumericVector grad_deltas(m);
  double grad_disp;
  NumericVector grad_betas(P);
  NumericVector out(2*m + 1 + P);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are person
  // as well as node and item specific
  // so we set up KxM matrices (nodesxitems) which we row-bind under each other 
  // for all N persons
  NumericMatrix mu(n_nodes*n, m);
  NumericMatrix mu_interp(n_nodes*n, m);
  NumericMatrix disp_interp(n_nodes*n, m);
  for (int i=0; i<n; i++) {
    // we are computing node-item specific mus for each person
    for(int j=0;j<m;j++){
      // loop over items (columns)
      for(int k=0;k<n_nodes;k++) {
        // loop over nodes (rows)
        double log_mu = alphas[j] * nodes[k] + deltas[j];
        for(int p=0; p<P; p++) {
          // add all the (weighted) covariate values for all covariates for the item j
          // (for the specific person i we are currently looking at)
          log_mu += betas[p] * alphas[j] * p_cov_data(i,p);
        }
        mu(k+i*n_nodes,j) = exp(log_mu);
        mu_interp(k+i*n_nodes,j) = mu(k+i*n_nodes,j);
        if (mu(k+i*n_nodes,j) > max_mu) { mu_interp(k+i*n_nodes,j) = max_mu; }
        if (mu(k+i*n_nodes,j) < min_mu) { mu_interp(k+i*n_nodes,j) = min_mu; }
        // we need to set maximum for mu to max_mu so that the interpolation will
        // work, max_mu is the maximum mu value in our grid for interpolation
        disp_interp(k+i*n_nodes,j) = disps[j];
      }
    }  // end loop over items
  } // end loop over N
  
  NumericMatrix V(n_nodes*n, m);
  NumericMatrix log_lambda(n_nodes*n, m);
  NumericMatrix log_Z(n_nodes*n, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes*persons and
  // as many columns as we have (same as mu_interp and nu_interp matrices)
  
  grad_disp = 0;
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_alphas[i] = 0;
    grad_deltas[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      // note the index for mu_interp and V_interp (and log_lambda and log_Z) must be adjusted 
      // as we now have (additionally) person specific mus and Vs
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute A and B for dispersion gradient
        double lambda = exp(log_lambda(k+j*n_nodes,i));
        double A = computeA(lambda, mu_interp(k+j*n_nodes,i), disps[i], log_Z(k+j*n_nodes,i), 10);
        double B = computeB(lambda, mu_interp(k+j*n_nodes,i), disps[i], log_Z(k+j*n_nodes,i), 10);
        
        // compute the sum over the weightes covariates for the gradient for alpha
        double sum_over_pcov = 0;
        for (int p=0; p<P; p++) {
          sum_over_pcov += betas[p] * p_cov_data(j,p);
        }
        
        // compute the gradients (summing over persons)
        grad_alphas[i] = grad_alphas[i] +
          PPs(j,k) * (mu_interp(k+j*n_nodes,i)*(nodes[k] + sum_over_pcov) /  V(k+j*n_nodes,i))*
          (data(j,i) -  mu_interp(k+j*n_nodes,i));
        grad_deltas[i] = grad_deltas[i] +
          PPs(j,k) * (mu_interp(k+j*n_nodes,i) / V(k+j*n_nodes,i))*(data(j,i) - 
          mu_interp(k+j*n_nodes,i));
        grad_disp = grad_disp +
          PPs(j,k) * (disps[i]*(A*(data(j,i) - mu_interp(k+j*n_nodes,i))/V(k+j*n_nodes,i) - 
          (logFactorial(data(j,i))-B)));
      }
    }
  }
  
  // gradients for person covariate weights
  for (int p=0; p<P; p++) { // over covariates
    grad_betas[p] = 0;
    for (int j=0; j<m; j++) { // over items
      for (int k=0;k<n_nodes;k++) { // over nodes (rows in my matrices)
        for (int i=0; i<n; i++) { // over persons
          grad_betas[p] += PPs(i,k) * (mu_interp(k+i*n_nodes,j) * alphas[j] * p_cov_data(i,p) / V(k+i*n_nodes,j)) *
            (data(i,j) - mu_interp(k+i*n_nodes,j));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over items
  } //end loop over P (person covariates)
  
  // fill up output vector
  for(int i=0;i<m;i++){
    out[i] = grad_alphas[i];
    out[i + m] = grad_deltas[i];
  }
  out[2*m] = grad_disp;
  for(int p=0; p<P; p++) {
    out[2*m + 1 + p] = grad_betas[p];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_pcov_cat_samedisps_cpp(NumericVector alphas,
                                               NumericVector deltas,
                                               NumericVector disps,
                                               NumericVector betas,
                                               NumericMatrix data,
                                               NumericMatrix p_cov_data,
                                               NumericMatrix resp_pattern,
                                               NumericMatrix PPs,
                                               NumericVector nodes, 
                                               NumericVector grid_mus,
                                               NumericVector grid_nus,
                                               NumericVector grid_cmp_var_long,
                                               NumericVector grid_log_lambda_long,
                                               NumericVector grid_logZ_long,
                                               double max_mu,
                                               double min_mu) {
  
  // assume that p_cov is a matrix of dummy coded categorical predictors
  // resp_pattern is a matrix of the same no. of cols than p_cov
  // and as many rows as we have distinct possible response patterns
  
  int N = data.nrow();
  int M = data.ncol();
  int K = nodes.size();
  int P = betas.size(); 
  int n_resp_patterns = resp_pattern.nrow();
  NumericVector grad_alphas(M);
  NumericVector grad_deltas(M);
  double grad_disp;
  NumericVector grad_betas(P);
  NumericVector out(2*M + 1 + P);
  
  // for person covariates, we need mus (and lambdas and Zs) for each node and
  // and then also for each response pattern
  // so first compute that
  NumericMatrix mu(K*n_resp_patterns, M);
  NumericMatrix mu_interp(K*n_resp_patterns, M);
  NumericMatrix disp_interp(K*n_resp_patterns, M);
  for (int l=0; l<n_resp_patterns; l++) {
    for(int j=0; j<M; j++){
      // loop over items (columns)
      for(int k=0; k<K; k++) {
        // loop over nodes (rows)
        double log_mu = alphas[j] * nodes[k] + deltas[j];
        for(int p=0; p<P; p++) {
          // this works because only includes columns for none-reference categories
          // for all covs in ref categories, resp_pattern will just always be zero in that row
          log_mu += betas[p] * alphas[j] * resp_pattern(l,p);
        }
        
        mu(k+l*K,j) = exp(log_mu);
        mu_interp(k+l*K,j) = mu(k+l*K,j);
        if (mu(k+l*K,j) > max_mu) { mu_interp(k+l*K,j) = max_mu; }
        if (mu(k+l*K,j) < min_mu) { mu_interp(k+l*K,j) = min_mu; }
        // we need to set maximum for mu to max_mu so that the interpolation will
        // work, max_mu is the maximum mu value in our grid for interpolation
        disp_interp(k+l*K,j) = disps[j];
      }
    }  // end loop over items
  }
  
  NumericMatrix log_lambda(K*n_resp_patterns, M);
  NumericMatrix log_Z(K*n_resp_patterns, M);
  NumericMatrix V(K*n_resp_patterns, M);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  
  grad_disp = 0;
  for(int i=0;i<M;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_alphas[i] = 0;
    grad_deltas[i] = 0;
    
    for(int k=0;k<K;k++) {
      // over nodes (rows in my matrices)
      
      for(int j=0;j<N;j++) {
        // loop over persons
        
        // check what response pattern person j had
        int l = 0;
        bool pattern_match = false;
        while (!pattern_match && l<n_resp_patterns) {
          // the second condition is just for safety that we dont get an eternal while loop
          // but we should always find a pattern match
          for (int p=0; p<P; p++) {
            pattern_match = p_cov_data(j,p) == resp_pattern(l,p);
          }
          // if the rows are the same, i am going to get out the foor loop with
          // pattern_match = TRUE and have l at the row of the pattern matrix
          // otherwise I am going to increase l by 1 and stay in the while loop to see if
          // the next row in the pattern matrix is a match for i
          if (!pattern_match) { l += 1; }
        }
        
        // we now know that person i has a response pattern like in row l of resp_pattern matrix
        // so their mu (and lambda, etc.) should be at row k+l*K
        
        // compute A and B for dispersion gradient
        double lambda = exp(log_lambda(k+l*K,i));
        double A = computeA(lambda, mu_interp(k+l*K,i), disps[i], log_Z(k+l*K,i), 10);
        double B = computeB(lambda, mu_interp(k+l*K,i), disps[i], log_Z(k+l*K,i), 10);
        
        // compute the sum over the weightes covariates for the gradient for alpha
        double sum_over_pcov = 0;
        for (int p=0; p<P; p++) {
          sum_over_pcov += betas[p] * p_cov_data(j,p);
        }
        
        // compute the gradients (summing over persons)
        grad_alphas[i] += PPs(j,k) * (mu_interp(k+l*K,i)*(nodes[k] + sum_over_pcov) /  V(k+l*K,i))*
          (data(j,i) -  mu_interp(k+l*K,i));
        grad_deltas[i] += PPs(j,k) * (mu_interp(k+l*K,i) / V(k+l*K,i))*(data(j,i) - 
          mu_interp(k+l*K,i));
        grad_disp += PPs(j,k) * (disps[i]*(A*(data(j,i) - mu_interp(k+l*K,i))/V(k+l*K,i) - 
          (logFactorial(data(j,i))-B)));
      }
    }
  }
  
  // gradients for person covariate weights
  for (int p=0; p<P; p++) { // over covariates
    grad_betas[p] = 0;
    for (int j=0; j<M; j++) { // over items
      for (int k=0; k<K; k++) { // over nodes (rows in my matrices)
        for (int i=0; i<N; i++) { // over persons
          // check what response pattern person i had
          int l = 0;
          bool pattern_match = false;
          while (!pattern_match && l<n_resp_patterns) {
            // the second condition is just for safety that we dont get an eternal while loop
            // but we should always find a pattern match
            for (int p=0; p<P; p++) {
              pattern_match = p_cov_data(i,p) == resp_pattern(l,p);
            }
            // if the rows are the same, i am going to get out the foor loop with
            // pattern_match = TRUE and have l at the row of the pattern matrix
            // otherwise I am going to increase l by 1 and stay in the while loop to see if
            // the next row in the pattern matrix is a match for i
            if (!pattern_match) { l += 1; }
          }
          
          // we now know that person i has a response pattern like in row l of resp_pattern matrix
          // so their mu (and lambda, etc.) should be at row k+l*K
          
          grad_betas[p] += PPs(i,k) * (mu_interp(k+l*K,j) * alphas[j] * p_cov_data(i,p) / V(k+l*K,j)) *
            (data(i,j) - mu_interp(k+l*K,j));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over items
  } //end loop over P (person covariates)
  
  // fill up output vector
  for(int i=0;i<M;i++){
    out[i] = grad_alphas[i];
    out[i + M] = grad_deltas[i];
  }
  out[2*M] = grad_disp;
  for(int p=0; p<P; p++) {
    out[2*M + 1 + p] = grad_betas[p];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_icov_delta_samedisps_cpp(NumericVector alphas,
                                     double delta,
                                     NumericVector disps,
                                     NumericVector betas,
                                     NumericMatrix data,
                                     NumericMatrix i_cov_data,
                                     NumericMatrix PPs,
                                     NumericVector nodes, 
                                     NumericVector grid_mus,
                                     NumericVector grid_nus,
                                     NumericVector grid_cmp_var_long,
                                     NumericVector grid_log_lambda_long,
                                     NumericVector grid_logZ_long,
                                     double max_mu,
                                     double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = alphas.size();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int I = betas.size();
  NumericVector grad_alphas(m);
  double grad_delta;
  double grad_disp;
  NumericVector grad_betas(I);
  NumericVector out(m + 2 + I);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alphas[j] * nodes[k] + delta;
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += betas[c] * i_cov_data(j,c); // for item j
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,j) = disps[j];
    }
  }  // end loop over items
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  grad_disp = 0;
  grad_delta = 0;
  // gradients for item parameters
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_alphas[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      // compute A and B for dispersion gradient
      double lambda = exp(log_lambda(k,i));
      double A = computeA(lambda, mu_interp(k,i), disps[i], log_Z(k,i), 10);
      double B = computeB(lambda, mu_interp(k,i), disps[i], log_Z(k,i), 10);
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_alphas[i] = grad_alphas[i] +
          PPs(j,k) * (nodes[k]*mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_delta += PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_disp += PPs(j,k) * (disps[i]*(A*(data(j,i) - mu_interp(k,i))/V(k,i) - (logFactorial(data(j,i))-B)));
      }
    }
  }
  
  // gradients for item covariate weights
  for (int c=0; c<I; c++) {
    // for each gamma of which we have one for each covariate-item combination
    grad_betas[c] = 0;
    for (int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      for (int i=0; i<n; i++) {
        // over persons
        for (int j=0; j<m; j++) {
          // over items (as the betas are only specific to item covariates, not items)
          grad_betas[c] += PPs(i,k) * (mu_interp(k,j)*i_cov_data(j,c) / V(k,j)) *
            (data(i,j) - mu_interp(k,j));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over P (person covariates)
  } // end loop over items
  
  // fill up output vector
  for(int i=0;i<m;i++){
    out[i] = grad_alphas[i];
  }
  out[m] = grad_delta;
  out[m + 1] = grad_disp;
  for(int c=0; c<I; c++) {
    out[m + 2 + c] = grad_betas[c];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_icov_alpha_samedisps_cpp(double alpha,
                                                     NumericVector deltas,
                                                     NumericVector disps,
                                                     NumericVector betas,
                                                     NumericMatrix data,
                                                     NumericMatrix i_cov_data,
                                                     NumericMatrix PPs,
                                                     NumericVector nodes, 
                                                     NumericVector grid_mus,
                                                     NumericVector grid_nus,
                                                     NumericVector grid_cmp_var_long,
                                                     NumericVector grid_log_lambda_long,
                                                     NumericVector grid_logZ_long,
                                                     double max_mu,
                                                     double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = data.ncol();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int I = betas.size();
  double grad_alpha;
  NumericVector grad_deltas(m);
  double grad_disp;
  NumericVector grad_betas(I);
  NumericVector out(m + 2 + I);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alpha * nodes[k] + deltas[j];
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += nodes[k] * betas[c] * i_cov_data(j,c); // for item j
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,j) = disps[j];
    }
  }  // end loop over items
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  grad_disp = 0;
  grad_alpha = 0;
  // gradients for item parameters
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_deltas[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      // compute A and B for dispersion gradient
      double lambda = exp(log_lambda(k,i));
      double A = computeA(lambda, mu_interp(k,i), disps[i], log_Z(k,i), 10);
      double B = computeB(lambda, mu_interp(k,i), disps[i], log_Z(k,i), 10);
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_alpha += PPs(j,k) * (nodes[k]*mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_deltas[j] += PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_disp += PPs(j,k) * (disps[i]*(A*(data(j,i) - mu_interp(k,i))/V(k,i) - (logFactorial(data(j,i))-B)));
      }
    }
  }
  
  // gradients for item covariate weights
  for (int c=0; c<I; c++) {
    // for each gamma of which we have one for each covariate-item combination
    grad_betas[c] = 0;
    for (int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      for (int i=0; i<n; i++) {
        // over persons
        for (int j=0; j<m; j++) {
          // over items (as the betas are only specific to item covariates, not items)
          grad_betas[c] += PPs(i,k) * (nodes[k]*mu_interp(k,j)*i_cov_data(j,c) / V(k,j)) *
            (data(i,j) - mu_interp(k,j));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over P (person covariates)
  } // end loop over items
  
  // fill up output vector
  out[0] = grad_alpha;
  for(int i=0;i<m;i++){
    out[i+1] = grad_deltas[i];
  }
  out[m + 1] = grad_disp;
  for(int c=0; c<I; c++) {
    out[m + 2 + c] = grad_betas[c];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_icov_alpha_delta_samedisps_cpp(double alpha,
                                         double delta,
                                         NumericVector disps,
                                         NumericVector betas_alpha,
                                         NumericVector betas_delta,
                                         NumericMatrix data,
                                         NumericMatrix i_cov_data,
                                         NumericMatrix PPs,
                                         NumericVector nodes, 
                                         NumericVector grid_mus,
                                         NumericVector grid_nus,
                                         NumericVector grid_cmp_var_long,
                                         NumericVector grid_log_lambda_long,
                                         NumericVector grid_logZ_long,
                                         double max_mu,
                                         double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = data.ncol();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int I = betas_alpha.size();
  double grad_alpha;
  double grad_delta;
  double grad_disp;
  NumericVector grad_betas_alpha(I);
  NumericVector grad_betas_delta(I);
  NumericVector out(3 + 2*I);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alpha * nodes[k] + delta;
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += nodes[k] * betas_alpha[c] * i_cov_data(j,c) + 
          betas_delta[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,j) = disps[j];
    }
  }  // end loop over items
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  grad_disp = 0;
  grad_delta = 0;
  grad_alpha = 0;
  NumericMatrix A(n_nodes, m);
  NumericMatrix B(n_nodes, m);
  // gradients for item parameters
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      // compute A and B for dispersion gradient
      double lambda = exp(log_lambda(k,i));
      A(k,i) = computeA(lambda, mu_interp(k,i), disp_interp(k,i), log_Z(k,i), 10);
      B(k,i) = computeB(lambda, mu_interp(k,i), disp_interp(k,i), log_Z(k,i), 10);
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_alpha += PPs(j,k) * (nodes[k]*mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_delta += PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_disp += PPs(j,k) * (disp_interp(k,i)*(A(k,i)*
          (data(j,i) - mu_interp(k,i))/V(k,i) - (logFactorial(data(j,i))-B(k,i))));
      }
    }
  }
  
  // gradients for item covariate weights
  for (int c=0; c<I; c++) {
    // for each gamma of which we have one for each covariate-item combination
    grad_betas_alpha[c] = 0;
    grad_betas_delta[c] = 0;
    for (int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      for (int i=0; i<n; i++) {
        // over persons
        for (int j=0; j<m; j++) {
          // over items (as the betas are only specific to item covariates, not items)
          grad_betas_alpha[c] += PPs(i,k) * (nodes[k]*mu_interp(k,j)*i_cov_data(j,c) / V(k,j)) *
            (data(i,j) - mu_interp(k,j));
          grad_betas_delta[c] += PPs(i,k) * (mu_interp(k,j)*i_cov_data(j,c) / V(k,j)) *
            (data(i,j) - mu_interp(k,j));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over P (person covariates)
  } // end loop over items
  
  // fill up output vector
  out[0] = grad_alpha;
  out[1] = grad_delta;
  out[2] = grad_disp;
  for(int c=0; c<I; c++) {
    out[3 + c] = grad_betas_alpha[c];
    out[3 + c + I] = grad_betas_delta[c];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_samealphas_newem_cpp(NumericVector alphas,
                                  NumericVector deltas,
                                  NumericVector disps,
                                  NumericMatrix data,
                                  NumericMatrix PPs,
                                  NumericVector nodes, 
                                  NumericVector grid_mus,
                                  NumericVector grid_nus,
                                  NumericVector grid_cmp_var_long,
                                  NumericVector grid_log_lambda_long,
                                  NumericVector grid_logZ_long,
                                  double max_mu,
                                  double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = alphas.size();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  double grad_alpha;
  NumericVector grad_deltas(m);
  NumericVector grad_disps(m);
  NumericVector out(2*m + 1);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int i=0;i<m;i++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over persons (rows)
      mu(k,i) = exp(alphas[i] * nodes[k] + deltas[i]);
      mu_interp(k,i) = mu(k,i);
      if (mu(k,i) > max_mu) { mu_interp(k,i) = max_mu; }
      if (mu(k,i) < min_mu) { mu_interp(k,i) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,i) = disps[i];
    }
  }
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have persons and
  // as many columns as we have items
  
  grad_alpha = 0;
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_deltas[i] = 0;
    grad_disps[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over persons (rows in my matrices)
      
      // compute A and B for dispersion gradient
      double lambda = exp(log_lambda(k,i));
      double A = computeA(lambda, mu_interp(k,i), disps[i], log_Z(k,i), 10);
      double B = computeB(lambda, mu_interp(k,i), disps[i], log_Z(k,i), 10);
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_alpha = grad_alpha +
          PPs(j,k) * (nodes[k]*mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_deltas[i] = grad_deltas[i] +
          PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_disps[i] = grad_disps[i] +
          PPs(j,k) * (disps[i]*(A*(data(j,i) - mu_interp(k,i))/V(k,i) - (logFactorial(data(j,i))-B)));
      }
    }
  }
  
  // fill up output vector
  out[0] = grad_alpha;
  for(int i=1;i<m+1;i++){
    out[i] = grad_deltas[i-1];
    out[i + m] = grad_disps[i-1];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_pcov_samealphas_cpp(NumericVector alphas,
                                     NumericVector deltas,
                                     NumericVector disps,
                                     NumericVector betas,
                                     NumericMatrix data,
                                     NumericMatrix p_cov_data,
                                     NumericMatrix PPs,
                                     NumericVector nodes, 
                                     NumericVector grid_mus,
                                     NumericVector grid_nus,
                                     NumericVector grid_cmp_var_long,
                                     NumericVector grid_log_lambda_long,
                                     NumericVector grid_logZ_long,
                                     double max_mu,
                                     double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = alphas.size();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int P = betas.size();
  double grad_alpha;
  NumericVector grad_deltas(m);
  NumericVector grad_disps(m);
  NumericVector grad_betas(P);
  NumericVector out(2*m + 1 + P);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are person
  // as well as node and item specific
  // so we set up KxM matrices (nodesxitems) which we row-bind under each other 
  // for all N persons
  NumericMatrix mu(n_nodes*n, m);
  NumericMatrix mu_interp(n_nodes*n, m);
  NumericMatrix disp_interp(n_nodes*n, m);
  for (int i=0; i<n; i++) {
    // we are computing node-item specific mus for each person
    for(int j=0;j<m;j++){
      // loop over items (columns)
      for(int k=0;k<n_nodes;k++) {
        // loop over nodes (rows)
        double log_mu = alphas[j] * nodes[k] + deltas[j];
        for(int p=0; p<P; p++) {
          // add all the (weighted) covariate values for all covariates for the item j
          // (for the specific person i we are currently looking at)
          log_mu += betas[p] * alphas[j] * p_cov_data(i,p);
        }
        mu(k+i*n_nodes,j) = exp(log_mu);
        mu_interp(k+i*n_nodes,j) = mu(k+i*n_nodes,j);
        if (mu(k+i*n_nodes,j) > max_mu) { mu_interp(k+i*n_nodes,j) = max_mu; }
        if (mu(k+i*n_nodes,j) < min_mu) { mu_interp(k+i*n_nodes,j) = min_mu; }
        // we need to set maximum for mu to max_mu so that the interpolation will
        // work, max_mu is the maximum mu value in our grid for interpolation
        disp_interp(k+i*n_nodes,j) = disps[j];
      }
    }  // end loop over items
  } // end loop over N
  
  NumericMatrix V(n_nodes*n, m);
  NumericMatrix log_lambda(n_nodes*n, m);
  NumericMatrix log_Z(n_nodes*n, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes*persons and
  // as many columns as we have (same as mu_interp and nu_interp matrices)
  
  grad_alpha = 0;
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_deltas[i] = 0;
    grad_disps[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      // note the index for mu_interp and V_interp (and log_lambda and log_Z) must be adjusted 
      // as we now have (additionally) person specific mus and Vs
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute A and B for dispersion gradient
        double lambda = exp(log_lambda(k+j*n_nodes,i));
        double A = computeA(lambda, mu_interp(k+j*n_nodes,i), disps[i], log_Z(k+j*n_nodes,i), 10);
        double B = computeB(lambda, mu_interp(k+j*n_nodes,i), disps[i], log_Z(k+j*n_nodes,i), 10);
        
        // compute the sum over the weightes covariates for the gradient for alpha
        double sum_over_pcov = 0;
        for (int p=0; p<P; p++) {
          sum_over_pcov += betas[p] * p_cov_data(j,p);
        }
        
        // compute the gradients (summing over persons)
        grad_alpha = grad_alpha +
          PPs(j,k) * (nodes[k]*mu_interp(k+j*n_nodes,i) / V(k+j*n_nodes,i))*(data(j,i) - 
          mu_interp(k+j*n_nodes,i));
        grad_deltas[i] = grad_deltas[i] +
          PPs(j,k) * (mu_interp(k+j*n_nodes,i) / V(k+j*n_nodes,i))*(data(j,i) - 
          mu_interp(k+j*n_nodes,i));
        grad_disps[i] = grad_disps[i] +
          PPs(j,k) * (disps[i]*(A*(data(j,i) - mu_interp(k+j*n_nodes,i))/V(k+j*n_nodes,i) - 
          (logFactorial(data(j,i))-B)));
      }
    }
  }
  
  // gradients for person covariate weights
  for (int p=0; p<P; p++) { // over covariates
    grad_betas[p] = 0;
    for (int j=0; j<m; j++) { // over items
      for (int k=0;k<n_nodes;k++) { // over nodes (rows in my matrices)
        for (int i=0; i<n; i++) { // over persons
          grad_betas[p] += PPs(i,k) * (mu_interp(k+i*n_nodes,j) * alphas[j] * p_cov_data(i,p) / V(k+i*n_nodes,j)) *
            (data(i,j) - mu_interp(k+i*n_nodes,j));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over items
  } //end loop over P (person covariates)
  
  // fill up output vector
  out[0] = grad_alpha;
  for(int i=0;i<m;i++){
    out[1 + i] = grad_deltas[i];
    out[1 + i + m] = grad_disps[i];
  }
  for(int p=0; p<P; p++) {
    out[2*m + 1 + p] = grad_betas[p];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_pcov_cat_samealphas_cpp(NumericVector alphas,
                                                NumericVector deltas,
                                                NumericVector disps,
                                                NumericVector betas,
                                                NumericMatrix data,
                                                NumericMatrix p_cov_data,
                                                NumericMatrix resp_pattern,
                                                NumericMatrix PPs,
                                                NumericVector nodes, 
                                                NumericVector grid_mus,
                                                NumericVector grid_nus,
                                                NumericVector grid_cmp_var_long,
                                                NumericVector grid_log_lambda_long,
                                                NumericVector grid_logZ_long,
                                                double max_mu,
                                                double min_mu) {
  
  // assume that p_cov is a matrix of dummy coded categorical predictors
  // resp_pattern is a matrix of the same no. of cols than p_cov
  // and as many rows as we have distinct possible response patterns
  
  int N = data.nrow();
  int M = data.ncol();
  int K = nodes.size();
  int P = betas.size(); 
  int n_resp_patterns = resp_pattern.nrow();
  double grad_alpha;
  NumericVector grad_deltas(M);
  NumericVector grad_disps(M);
  NumericVector grad_betas(P);
  NumericVector out(2*M + 1 + P);
  
  // for person covariates, we need mus (and lambdas and Zs) for each node and
  // and then also for each response pattern
  // so first compute that
  NumericMatrix mu(K*n_resp_patterns, M);
  NumericMatrix mu_interp(K*n_resp_patterns, M);
  NumericMatrix disp_interp(K*n_resp_patterns, M);
  for (int l=0; l<n_resp_patterns; l++) {
    for(int j=0; j<M; j++){
      // loop over items (columns)
      for(int k=0; k<K; k++) {
        // loop over nodes (rows)
        double log_mu = alphas[j] * nodes[k] + deltas[j];
        for(int p=0; p<P; p++) {
          // this works because only includes columns for none-reference categories
          // for all covs in ref categories, resp_pattern will just always be zero in that row
          log_mu += betas[p] * alphas[j] * resp_pattern(l,p);
        }
        
        mu(k+l*K,j) = exp(log_mu);
        mu_interp(k+l*K,j) = mu(k+l*K,j);
        if (mu(k+l*K,j) > max_mu) { mu_interp(k+l*K,j) = max_mu; }
        if (mu(k+l*K,j) < min_mu) { mu_interp(k+l*K,j) = min_mu; }
        // we need to set maximum for mu to max_mu so that the interpolation will
        // work, max_mu is the maximum mu value in our grid for interpolation
        disp_interp(k+l*K,j) = disps[j];
      }
    }  // end loop over items
  }
  
  NumericMatrix log_lambda(K*n_resp_patterns, M);
  NumericMatrix log_Z(K*n_resp_patterns, M);
  NumericMatrix V(K*n_resp_patterns, M);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  
  grad_alpha = 0;
  for(int i=0;i<M;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_deltas[i] = 0;
    grad_disps[i] = 0;
    
    for(int k=0;k<K;k++) {
      // over nodes (rows in my matrices)
      
      // note the index for mu_interp and V_interp (and log_lambda and log_Z) must be adjusted 
      // as we now have (additionally) person specific mus and Vs
      for(int j=0;j<N;j++) {
        // loop over persons
        
        // check what response pattern person i had
        int l = 0;
        bool pattern_match = false;
        while (!pattern_match && l<n_resp_patterns) {
          // the second condition is just for safety that we dont get an eternal while loop
          // but we should always find a pattern match
          for (int p=0; p<P; p++) {
            pattern_match = p_cov_data(j,p) == resp_pattern(l,p);
          }
          // if the rows are the same, i am going to get out the foor loop with
          // pattern_match = TRUE and have l at the row of the pattern matrix
          // otherwise I am going to increase l by 1 and stay in the while loop to see if
          // the next row in the pattern matrix is a match for i
          if (!pattern_match) { l += 1; }
        }
        
        // we now know that person i has a response pattern like in row l of resp_pattern matrix
        // so their mu (and lambda, etc.) should be at row k+l*K
        
        // compute A and B for dispersion gradient
        double lambda = exp(log_lambda(k+l*K,i));
        double A = computeA(lambda, mu_interp(k+l*K,i), disps[i], log_Z(k+l*K,i), 10);
        double B = computeB(lambda, mu_interp(k+l*K,i), disps[i], log_Z(k+l*K,i), 10);
        
        // compute the sum over the weightes covariates for the gradient for alpha
        double sum_over_pcov = 0;
        for (int p=0; p<P; p++) {
          sum_over_pcov += betas[p] * p_cov_data(j,p);
        }
        
        // compute the gradients (summing over persons)
        grad_alpha += PPs(j,k) * (nodes[k]*mu_interp(k+l*K,i) / V(k+l*K,i))*(data(j,i) - 
          mu_interp(k+l*K,i));
        grad_deltas[i] += PPs(j,k) * (mu_interp(k+l*K,i) / V(k+l*K,i))*(data(j,i) - 
          mu_interp(k+l*K,i));
        grad_disps[i] += PPs(j,k) * (disps[i]*(A*(data(j,i) - mu_interp(k+l*K,i))/V(k+l*K,i) - 
          (logFactorial(data(j,i))-B)));
      }
    }
  }
  
  // gradients for person covariate weights
  for (int p=0; p<P; p++) { // over covariates
    grad_betas[p] = 0;
    for (int j=0; j<M; j++) { // over items
      for (int k=0;k<K;k++) { // over nodes (rows in my matrices)
        for (int i=0; i<N; i++) { // over persons
          
          // check what response pattern person i had
          int l = 0;
          bool pattern_match = false;
          while (!pattern_match && l<n_resp_patterns) {
            // the second condition is just for safety that we dont get an eternal while loop
            // but we should always find a pattern match
            for (int p=0; p<P; p++) {
              pattern_match = p_cov_data(i,p) == resp_pattern(l,p);
            }
            // if the rows are the same, i am going to get out the foor loop with
            // pattern_match = TRUE and have l at the row of the pattern matrix
            // otherwise I am going to increase l by 1 and stay in the while loop to see if
            // the next row in the pattern matrix is a match for i
            if (!pattern_match) { l += 1; }
          }
          
          // we now know that person i has a response pattern like in row l of resp_pattern matrix
          // so their mu (and lambda, etc.) should be at row k+l*K
          
          grad_betas[p] += PPs(i,k) * (mu_interp(k+l*K,j) * alphas[j] * p_cov_data(i,p) / V(k+l*K,j)) *
            (data(i,j) - mu_interp(k+l*K,j));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over items
  } //end loop over P (person covariates)
  
  // fill up output vector
  out[0] = grad_alpha;
  for(int i=0;i<M;i++){
    out[1 + i] = grad_deltas[i];
    out[1 + i + M] = grad_disps[i];
  }
  for(int p=0; p<P; p++) {
    out[2*M + 1 + p] = grad_betas[p];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_icov_delta_samealphas_cpp(NumericVector alphas,
                                     double delta,
                                     NumericVector disps,
                                     NumericVector betas,
                                     NumericMatrix data,
                                     NumericMatrix i_cov_data,
                                     NumericMatrix PPs,
                                     NumericVector nodes, 
                                     NumericVector grid_mus,
                                     NumericVector grid_nus,
                                     NumericVector grid_cmp_var_long,
                                     NumericVector grid_log_lambda_long,
                                     NumericVector grid_logZ_long,
                                     double max_mu,
                                     double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = alphas.size();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int I = betas.size();
  double grad_alpha;
  double grad_delta;
  NumericVector grad_disps(m);
  NumericVector grad_betas(I);
  NumericVector out(m + 2 + I);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alphas[j] * nodes[k] + delta;
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += betas[c] * i_cov_data(j,c); // for item j
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,j) = disps[j];
    }
  }  // end loop over items
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  grad_alpha = 0;
  grad_delta = 0;
  // gradients for item parameters
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_disps[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      // compute A and B for dispersion gradient
      double lambda = exp(log_lambda(k,i));
      double A = computeA(lambda, mu_interp(k,i), disps[i], log_Z(k,i), 10);
      double B = computeB(lambda, mu_interp(k,i), disps[i], log_Z(k,i), 10);
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_alpha += PPs(j,k) * (nodes[k]*mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_delta += PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_disps[i] = grad_disps[i] +
          PPs(j,k) * (disps[i]*(A*(data(j,i) - mu_interp(k,i))/V(k,i) - (logFactorial(data(j,i))-B)));
      }
    }
  }
  
  // gradients for item covariate weights
  for (int c=0; c<I; c++) {
    // for each gamma of which we have one for each covariate-item combination
    grad_betas[c] = 0;
    for (int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      for (int i=0; i<n; i++) {
        // over persons
        for (int j=0; j<m; j++) {
          // over items (as the betas are only specific to item covariates, not items)
          grad_betas[c] += PPs(i,k) * (mu_interp(k,j)*i_cov_data(j,c) / V(k,j)) *
            (data(i,j) - mu_interp(k,j));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over P (person covariates)
  } // end loop over items
  
  // fill up output vector
  out[0] = grad_alpha;
  out[1] = grad_delta;
  for(int i=0;i<m;i++){
    out[i + 2] = grad_disps[i];
  }
  for(int c=0; c<I; c++) {
    out[2 + m + c] = grad_betas[c];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_icov_nu_samealphas_cpp(NumericVector alphas,
                                                      NumericVector deltas,
                                                      double disp,
                                                      NumericVector betas,
                                                      NumericMatrix data,
                                                      NumericMatrix i_cov_data,
                                                      NumericMatrix PPs,
                                                      NumericVector nodes, 
                                                      NumericVector grid_mus,
                                                      NumericVector grid_nus,
                                                      NumericVector grid_cmp_var_long,
                                                      NumericVector grid_log_lambda_long,
                                                      NumericVector grid_logZ_long,
                                                      double max_mu,
                                                      double min_mu,
                                                      double max_nu,
                                                      double min_nu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = alphas.size();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int I = betas.size();
  double grad_alpha;
  NumericVector grad_deltas(m);
  double grad_disp;
  NumericVector grad_betas(I);
  NumericVector out(m + 2 + I);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alphas[j] * nodes[k] + deltas[j];
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      double log_disp = log(disp);
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_disp += betas[c] * i_cov_data(j,c); // for item j
      }
      disp_interp(k,j) = exp(log_disp);
      if (disp_interp(k,j) > max_nu) { disp_interp(k,j) = max_nu; }
      if (disp_interp(k,j) < min_nu) { disp_interp(k,j) = min_nu; }
    }
  }  // end loop over items
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  grad_alpha = 0;
  grad_disp = 0;
  NumericMatrix A(n_nodes, m);
  NumericMatrix B(n_nodes, m);
  // gradients for item parameters
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    grad_deltas[i] = 0;
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      // compute A and B for dispersion gradient
      double lambda = exp(log_lambda(k,i));
      A(k,i) = computeA(lambda, mu_interp(k,i), disp_interp(k,i), log_Z(k,i), 10);
      B(k,i) = computeB(lambda, mu_interp(k,i),disp_interp(k,i), log_Z(k,i), 10);
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_alpha += PPs(j,k) * (nodes[k]*mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_deltas[i] += PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_disp += PPs(j,k) * (disp_interp(k,i)*(A(k,i)*
          (data(j,i) - mu_interp(k,i))/V(k,i) - (logFactorial(data(j,i))-B(k,i))));
      }
    }
  }
  
  // gradients for item covariate weights
  for (int c=0; c<I; c++) {
    // for each gamma of which we have one for each covariate-item combination
    grad_betas[c] = 0;
    for (int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      for (int i=0; i<n; i++) {
        // over persons
        for (int j=0; j<m; j++) {
          // over items (as the betas are only specific to item covariates, not items)
          grad_betas[c] += PPs(i,k) * i_cov_data(j,c) * disp_interp(k,j)*
            (A(k,j)*(data(i,j) - mu_interp(k,j))/V(k,j) - (logFactorial(data(i,j))-B(k,j)));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over P (person covariates)
  } // end loop over items
  
  // fill up output vector
  out[0] = grad_alpha;
  for(int i=0;i<m;i++){
    out[i + 1] = grad_deltas[i];
  }
  out[m + 1] = grad_disp;
  for(int c=0; c<I; c++) {
    out[2 + m + c] = grad_betas[c];
  }
  
  return(out);
}

// [[Rcpp::export]]
NumericVector grad_cmp_with_icov_delta_nu_samealphas_cpp(NumericVector alphas,
                                                         double delta,
                                                         double disp,
                                                         NumericVector betas_delta,
                                                         NumericVector betas_logdisp,
                                                         NumericMatrix data,
                                                         NumericMatrix i_cov_data,
                                                         NumericMatrix PPs,
                                                         NumericVector nodes, 
                                                         NumericVector grid_mus,
                                                         NumericVector grid_nus,
                                                         NumericVector grid_cmp_var_long,
                                                         NumericVector grid_log_lambda_long,
                                                         NumericVector grid_logZ_long,
                                                         double max_mu,
                                                         double min_mu,
                                                         double max_nu,
                                                         double min_nu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = data.ncol();
  int n = PPs.nrow();
  int n_nodes = nodes.size();
  int I = betas_delta.size();
  double grad_alpha;
  double grad_delta;
  double grad_disp;
  NumericVector grad_betas_delta(I);
  NumericVector grad_betas_logdisp(I);
  NumericVector out(3 + 2*I);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are node and item specific
  // contrary to the case of person covariates, we don't have to make them person specific
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alphas[j] * nodes[k] + delta;
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += betas_delta[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      double log_disp = log(disp);
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_disp += betas_logdisp[c] * i_cov_data(j,c); // for item j
      }
      disp_interp(k,j) = exp(log_disp);
      if (disp_interp(k,j) > max_nu) { disp_interp(k,j) = max_nu; }
      if (disp_interp(k,j) < min_nu) { disp_interp(k,j) = min_nu; }
    }
  }  // end loop over items
  
  NumericMatrix V(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  NumericMatrix log_Z(n_nodes, m);
  V = interp_from_grid_m(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have
  
  grad_disp = 0;
  grad_delta = 0;
  grad_alpha = 0;
  NumericMatrix A(n_nodes, m);
  NumericMatrix B(n_nodes, m);
  // gradients for item parameters
  for(int i=0;i<m;i++){
    // over items (columns in my matrices)
    // so that we get one gradient per item
    
    for(int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      
      // compute A and B for dispersion gradient
      double lambda = exp(log_lambda(k,i));
      A(k,i) = computeA(lambda, mu_interp(k,i), disp_interp(k,i), log_Z(k,i), 10);
      B(k,i) = computeB(lambda, mu_interp(k,i), disp_interp(k,i), log_Z(k,i), 10);
      
      for(int j=0;j<n;j++) {
        // loop over persons
        
        // compute the gradients (summing over persons)
        grad_alpha += PPs(j,k) * (nodes[k]*mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_delta += PPs(j,k) * (mu_interp(k,i) / V(k,i))*(data(j,i) - mu_interp(k,i));
        grad_disp += PPs(j,k) * (disp_interp(k,i)*(A(k,i)*
          (data(j,i) - mu_interp(k,i))/V(k,i) - (logFactorial(data(j,i))-B(k,i))));
      }
    }
  }
  
  // gradients for item covariate weights
  for (int c=0; c<I; c++) {
    // for each gamma of which we have one for each covariate-item combination
    grad_betas_delta[c] = 0;
    grad_betas_logdisp[c] = 0;
    for (int k=0;k<n_nodes;k++) {
      // over nodes (rows in my matrices)
      for (int i=0; i<n; i++) {
        // over persons
        for (int j=0; j<m; j++) {
          // over items (as the betas are only specific to item covariates, not items)
          grad_betas_delta[c] += PPs(i,k) * (mu_interp(k,j)*i_cov_data(j,c) / V(k,j)) *
            (data(i,j) - mu_interp(k,j));
          grad_betas_logdisp[c] += PPs(i,k) * i_cov_data(j,c) * disp_interp(k,j)*
            (A(k,j)*(data(i,j) - mu_interp(k,j))/V(k,j) - (logFactorial(data(i,j))-B(k,j)));
        } // end loop over m (items)
      } // end loop of n_nodes
    } // end loop over P (person covariates)
  } // end loop over items
  
  // fill up output vector
  out[0] = grad_alpha;
  out[1] = grad_delta;
  out[2] = grad_disp;
  for(int c=0; c<I; c++) {
    out[3 + c] = grad_betas_delta[c];
    out[3 + c + I] = grad_betas_logdisp[c];
  }
  
  return(out);
}

// [[Rcpp::export]]
double grad_ll_cmp_ability_1P_cpp(double ability,
                                  NumericVector alphas,
                                  NumericVector deltas,
                                  NumericVector disps,
                                  NumericVector data,
                                  NumericVector grid_mus,
                                  NumericVector grid_nus,
                                  NumericVector grid_cmp_var_long,
                                  NumericVector grid_log_lambda_long,
                                  double max_mu,
                                  double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = alphas.size();
  double grad;
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  NumericVector mu(m);
  NumericVector mu_interp(m);
  for(int j=0;j<m;j++){
    // loop over items
    mu[j] = exp(alphas[j] * ability + deltas[j]);
    mu_interp[j] = mu[j];
    if (mu[j] > max_mu) { mu_interp[j] = max_mu; }
    if (mu[j] < min_mu) { mu_interp[j] = min_mu; }
    // we need to set maximum for mu to max_mu so that the interpolation will
    // work, max_mu is the maximum mu value in our grid for interpolation
  }
  
  NumericVector V(m);
  NumericVector log_lambda(m);
  V = interp_from_grid_v(grid_mus, grid_nus,
                         grid_cmp_var_long,
                         mu_interp, disps);
  log_lambda = interp_from_grid_v(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disps);
  // V and log_lambda are matrices with as many rows as we have persons and
  // as many columns as we have items
  
  grad = 0;
  NumericVector mu_ttpo_two(m);
  for(int j=0;j<m;j++){
    // over items 
    mu_ttpo_two[j] = mu_interp[j]*mu_interp[j];
    grad = grad +
          ((data[j]*alphas[j]*mu_interp[j])/V[j]) - 
          ((alphas[j]*mu_ttpo_two[j])/exp(log_lambda[j]));
  }
  
  return(grad);
}
  


// [[Rcpp::export]]
double ell_cmp_newem_cpp (NumericVector alphas,
                          NumericVector deltas,
                          NumericVector disps,
                          NumericMatrix data,
                          NumericVector exp_abilities,
                          NumericVector grid_mus,
                          NumericVector grid_nus,
                          NumericVector grid_cmp_var_long,
                          NumericVector grid_log_lambda_long,
                          NumericVector grid_logZ_long,
                          double max_mu,
                          double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int m = alphas.size();
  int n = exp_abilities.size();
  // NumericVector grad_alphas(m);
  // NumericVector grad_deltas(m);
  // NumericVector grad_disps(m);
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  NumericMatrix mu(n, m);
  NumericMatrix mu_interp(n, m);
  NumericMatrix disp_interp(n, m);
  for(int i=0;i<m;i++){
    // loop over items (columns)
    for(int j=0;j<n;j++) {
      // loop over persons (rows)
      mu(j,i) = exp(alphas[i] * exp_abilities[j] + deltas[i]);
      mu_interp(j,i) = mu(j,i);
      if (mu(j,i) > max_mu) { mu_interp(j,i) = max_mu; }
      if (mu(j,i) < min_mu) { mu_interp(j,i) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(j,i) = disps[i];
    }
  }
  
  NumericMatrix log_lambda(n, m);
  NumericMatrix log_Z(n, m);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have persons and
  // as many columns as we have items
  
  double out = 0;
  for(int i=0;i<m;i++){
    double sum_for_item = 0;
    for(int j=0;j<n;j++) {
      sum_for_item = sum_for_item + 
        data(j,i)*log_lambda(j,i) - log_Z(j,i) - disps[i]*logFactorial(data(j,i));
//(r(j,i)*log_lambda(j,i) - disps[i]*h(j,i) - f(j,i)*log_Z(j,i));
    }
    out = out + sum_for_item;
  }
  
  return(out);
}

// [[Rcpp::export]]
double ell_cmp_with_pcov_cpp (NumericVector alphas,
                          NumericVector deltas,
                          NumericVector disps,
                          NumericVector betas,
                          NumericMatrix data,
                          NumericMatrix p_cov_data,
                          NumericVector PPs,
                          NumericVector weights,
                          NumericVector nodes,
                          NumericVector grid_mus,
                          NumericVector grid_nus,
                          NumericVector grid_cmp_var_long,
                          NumericVector grid_log_lambda_long,
                          NumericVector grid_logZ_long,
                          double max_mu,
                          double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int K = nodes.size();
  int m = alphas.size();
  int N = data.nrow();
  int P = betas.size();
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for person covariates, we need mus (and lambdas and Zs) which are person
  // as well as node and item specific
  // so then I extend my nu and mu matrices for interpolation accordingly
  // so that i can interpolate lambdas and Zs person-node-item specifically
  // but still only work with matrices so that i can use interp_from_grid_m
  // here i just chain KxM matrices (like I had for no covariates) below each other
  // (so rbind basically), for all N person so that for the first K rows,
  // we have the KxM matrix for person 1, for rows K+1 - K+K we have the
  // KxM matrix for person 1, etc.
  NumericMatrix mu(K*N, m);
  NumericMatrix mu_interp(K*N, m);
  NumericMatrix disp_interp(K*N, m);
  for (int i=0; i<N; i++) {
    // we are computing node-item specific mus for each person
    for(int j=0;j<m;j++){
      // loop over items (columns)
      for(int k=0;k<K;k++) {
        // loop over nodes (rows)
        double log_mu = alphas[j] * nodes[k] + deltas[j];
        // my nodes are here my epsilon
        for(int p=0; p<P; p++) {
          // add all the (weighted) covariate values for that specific item j
          // and person
          log_mu += betas[p] * alphas[j] * p_cov_data(i,p);
        }
        mu(k+i*K,j) = exp(log_mu);
        mu_interp(k+i*K,j) = mu(k+i*K,j);
        if (mu(k+i*K,j) > max_mu) { mu_interp(k+i*K,j) = max_mu; }
        if (mu(k+i*K,j) < min_mu) { mu_interp(k+i*K,j) = min_mu; }
        // we need to set maximum for mu to max_mu so that the interpolation will
        // work, max_mu is the maximum mu value in our grid for interpolation
        disp_interp(k+i*K,j) = disps[j];
      }
    }  // end loop over items
  } // end loop over N
  
  NumericMatrix log_Z(K*N, m);
  NumericMatrix log_lambda(K*N, m);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes*persons and
  // as many columns as we have 
  // they have the same structure as mu_interp and nu_interp matrices above 
  // (where the structure is explained in more detail)
  
  double out = 0;
  for (int k=0; k<K; k++) { // nodes
    for(int i=0;i<N;i++) { // persons
      for(int j=0;j<m;j++) { // items
        out += (data(i,j)*log_lambda(k+i*K,j) - log_Z(k+i*K,j) - 
          disps[j]*logFactorial(data(i,j))) * PPs(i,k);
      }
    }
  } // end loops over K nodes
  
  return(out);
}

// [[Rcpp::export]]
double ell_cmp_with_pcov_cat_cpp (NumericVector alphas,
                              NumericVector deltas,
                              NumericVector disps,
                              NumericVector betas,
                              NumericMatrix data,
                              NumericMatrix p_cov_data,
                              NumericMatrix resp_pattern,
                              NumericVector PPs,
                              NumericVector weights,
                              NumericVector nodes,
                              NumericVector grid_mus,
                              NumericVector grid_nus,
                              NumericVector grid_cmp_var_long,
                              NumericVector grid_log_lambda_long,
                              NumericVector grid_logZ_long,
                              double max_mu,
                              double min_mu) {
  
  // assume that p_cov is a matrix of dummy coded categorical predictors
  // resp_pattern is a matrix of the same no. of cols than p_cov
  // and as many rows as we have distinct possible response patterns
  
  int N = data.nrow();
  int M = data.ncol();
  int K = nodes.size();
  int P = betas.size(); 
  int n_resp_patterns = resp_pattern.nrow();
  
  // for person covariates, we need mus (and lambdas and Zs) for each node and
  // and then also for each response pattern
  // so first compute that
  NumericMatrix mu(K*n_resp_patterns, M);
  NumericMatrix mu_interp(K*n_resp_patterns, M);
  NumericMatrix disp_interp(K*n_resp_patterns, M);
  for (int l=0; l<n_resp_patterns; l++) {
    for(int j=0; j<M; j++){
      // loop over items (columns)
      for(int k=0; k<K; k++) {
        // loop over nodes (rows)
        double log_mu = alphas[j] * nodes[k] + deltas[j];
        for(int p=0; p<P; p++) {
          // this works because only includes columns for none-reference categories
          // for all covs in ref categories, resp_pattern will just always be zero in that row
          log_mu += betas[p] * alphas[j] * resp_pattern(l,p);
        }
        
        mu(k+l*K,j) = exp(log_mu);
        mu_interp(k+l*K,j) = mu(k+l*K,j);
        if (mu(k+l*K,j) > max_mu) { mu_interp(k+l*K,j) = max_mu; }
        if (mu(k+l*K,j) < min_mu) { mu_interp(k+l*K,j) = min_mu; }
        // we need to set maximum for mu to max_mu so that the interpolation will
        // work, max_mu is the maximum mu value in our grid for interpolation
        disp_interp(k+l*K,j) = disps[j];
      }
    }  // end loop over items
  }
  
  NumericMatrix log_lambda(K*n_resp_patterns, M);
  NumericMatrix log_Z(K*n_resp_patterns, M);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  
  double out = 0;
  for (int k=0; k<K; k++) { // nodes
    for(int i=0;i<N;i++) { // persons
      // check what response pattern person i had
      int l = 0;
      bool pattern_match = false;
      while (!pattern_match && l<n_resp_patterns) {
        // the second condition is just for safety that we dont get an eternal while loop
        // but we should always find a pattern match
        for (int p=0; p<P; p++) {
          pattern_match = p_cov_data(i,p) == resp_pattern(l,p);
        }
        // if the rows are the same, i am going to get out the foor loop with
        // pattern_match = TRUE and have l at the row of the pattern matrix
        // otherwise I am going to increase l by 1 and stay in the while loop to see if
        // the next row in the pattern matrix is a match for i
        if (!pattern_match) { l += 1; }
      }
      
      // we now know that person i has a response pattern like in row l of resp_pattern matrix
      // so their mu (and lambda, etc.) should be at row k+l*K
      
      for(int j=0;j<M;j++) { // items
        out += (data(i,j)*log_lambda(k+l*K,j) - log_Z(k+l*K,j) - 
          disps[j]*logFactorial(data(i,j))) * PPs(i,k);
      }
    }
  } // end loops over K nodes
  
  return(out);
}

// [[Rcpp::export]]
double ell_cmp_with_icov_delta_cpp (NumericVector alphas,
                              double delta,
                              NumericVector disps,
                              NumericVector betas,
                              NumericMatrix data,
                              NumericMatrix i_cov_data,
                              NumericVector PPs,
                              NumericVector weights,
                              NumericVector nodes,
                              NumericVector grid_mus,
                              NumericVector grid_nus,
                              NumericVector grid_cmp_var_long,
                              NumericVector grid_log_lambda_long,
                              NumericVector grid_logZ_long,
                              double max_mu,
                              double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int K = nodes.size();
  int m = alphas.size();
  int N = data.nrow();
  int I = betas.size();
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for item covariates, we can do normal mu_interp and nu_interp, so just
  // of dumension KxM, they don't need to be person specific
  NumericMatrix mu(K, m);
  NumericMatrix mu_interp(K, m);
  NumericMatrix disp_interp(K, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<K;k++) {
      // loop over nodes (rows)
      double log_mu = alphas[j] * nodes[k] + delta;
      // my nodes are here my epsilon
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values
        log_mu += betas[c] * i_cov_data(j,c);
        // i_cov_data has a row for each item and a column for each item covariate
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,j) = disps[j];
    }
  }  // end loop over items
  
  NumericMatrix log_Z(K, m);
  NumericMatrix log_lambda(K, m);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes*persons and
  // as many columns as we have 
  // they have the same structure as mu_interp and nu_interp matrices above 
  // (where the structure is explained in more detail)
  
  double out = 0;
  for (int k=0; k<K; k++) { // nodes
    for(int i=0;i<N;i++) { // persons
      for(int j=0;j<m;j++) { // items
        out += (data(i,j)*log_lambda(k,j) - log_Z(k,j) - 
          disps[j]*logFactorial(data(i,j))) * PPs(i,k);
      }
    }
  } // end loops over K nodes
  
  return(out);
}

// [[Rcpp::export]]
double ell_cmp_with_icov_alpha_cpp (double alpha,
                                    NumericVector deltas,
                                    NumericVector disps,
                                    NumericVector betas,
                                    NumericMatrix data,
                                    NumericMatrix i_cov_data,
                                    NumericVector PPs,
                                    NumericVector weights,
                                    NumericVector nodes,
                                    NumericVector grid_mus,
                                    NumericVector grid_nus,
                                    NumericVector grid_cmp_var_long,
                                    NumericVector grid_log_lambda_long,
                                    NumericVector grid_logZ_long,
                                    double max_mu,
                                    double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int K = nodes.size();
  int m = data.ncol();
  int N = data.nrow();
  int I = betas.size();
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for item covariates, we can do normal mu_interp and nu_interp, so just
  // of dumension KxM, they don't need to be person specific
  NumericMatrix mu(K, m);
  NumericMatrix mu_interp(K, m);
  NumericMatrix disp_interp(K, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<K;k++) {
      // loop over nodes (rows)
      double log_mu = alpha * nodes[k] + deltas[j];
      // my nodes are here my epsilon
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values
        log_mu += nodes[k] * betas[c] * i_cov_data(j,c);
        // i_cov_data has a row for each item and a column for each item covariate
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,j) = disps[j];
    }
  }  // end loop over items
  
  NumericMatrix log_Z(K, m);
  NumericMatrix log_lambda(K, m);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes*persons and
  // as many columns as we have 
  // they have the same structure as mu_interp and nu_interp matrices above 
  // (where the structure is explained in more detail)
  
  double out = 0;
  for (int k=0; k<K; k++) { // nodes
    for(int i=0;i<N;i++) { // persons
      for(int j=0;j<m;j++) { // items
        out += (data(i,j)*log_lambda(k,j) - log_Z(k,j) - 
          disps[j]*logFactorial(data(i,j))) * PPs(i,k);
      }
    }
  } // end loops over K nodes
  
  return(out);
}

// [[Rcpp::export]]
double ell_cmp_with_icov_nu_cpp (NumericVector alphas,
                                    NumericVector deltas,
                                    double disp,
                                    NumericVector betas,
                                    NumericMatrix data,
                                    NumericMatrix i_cov_data,
                                    NumericVector PPs,
                                    NumericVector weights,
                                    NumericVector nodes,
                                    NumericVector grid_mus,
                                    NumericVector grid_nus,
                                    NumericVector grid_cmp_var_long,
                                    NumericVector grid_log_lambda_long,
                                    NumericVector grid_logZ_long,
                                    double max_mu,
                                    double min_mu,
                                    double max_nu,
                                    double min_nu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int K = nodes.size();
  int m = data.ncol();
  int N = data.nrow();
  int I = betas.size();
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for item covariates, we can do normal mu_interp and nu_interp, so just
  // of dumension KxM, they don't need to be person specific
  NumericMatrix mu(K, m);
  NumericMatrix mu_interp(K, m);
  NumericMatrix disp_interp(K, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<K;k++) {
      // loop over nodes (rows)
      double log_mu = alphas[j] * nodes[k] + deltas[j];
      // my nodes are here my epsilon
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      double log_disp = log(disp);
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_disp += betas[c] * i_cov_data(j,c); // for item j
      }
      disp_interp(k,j) = exp(log_disp);
      if (disp_interp(k,j) > max_nu) { disp_interp(k,j) = max_nu; }
      if (disp_interp(k,j) < min_nu) { disp_interp(k,j) = min_nu; }
    }
  }  // end loop over items
  
  NumericMatrix log_Z(K, m);
  NumericMatrix log_lambda(K, m);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes*persons and
  // as many columns as we have 
  // they have the same structure as mu_interp and nu_interp matrices above 
  // (where the structure is explained in more detail)
  
  double out = 0;
  for (int k=0; k<K; k++) { // nodes
    for(int i=0;i<N;i++) { // persons
      for(int j=0;j<m;j++) { // items
        out += (data(i,j)*log_lambda(k,j) - log_Z(k,j) - 
          disp_interp(k,j)*logFactorial(data(i,j))) * PPs(i,k);
      }
    }
  } // end loops over K nodes
  
  return(out);
}

// [[Rcpp::export]]
double ell_cmp_with_icov_all_cpp (double alpha,
                                 double delta,
                                 double disp,
                                 NumericVector betas_alpha,
                                 NumericVector betas_delta,
                                 NumericVector betas_logdisp,
                                 NumericMatrix data,
                                 NumericMatrix i_cov_data,
                                 NumericVector PPs,
                                 NumericVector weights,
                                 NumericVector nodes,
                                 NumericVector grid_mus,
                                 NumericVector grid_nus,
                                 NumericVector grid_cmp_var_long,
                                 NumericVector grid_log_lambda_long,
                                 NumericVector grid_logZ_long,
                                 double max_mu,
                                 double min_mu,
                                 double max_nu,
                                 double min_nu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int K = nodes.size();
  int m = data.ncol();
  int N = data.nrow();
  int I = betas_alpha.size();
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for item covariates, we can do normal mu_interp and nu_interp, so just
  // of dumension KxM, they don't need to be person specific
  NumericMatrix mu(K, m);
  NumericMatrix mu_interp(K, m);
  NumericMatrix disp_interp(K, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<K;k++) {
      // loop over nodes (rows)
      double log_mu = alpha * nodes[k] + delta;
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += nodes[k] * betas_alpha[c] * i_cov_data(j,c) + 
          betas_delta[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      double log_disp = log(disp);
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_disp += betas_logdisp[c] * i_cov_data(j,c); // for item j
      }
      disp_interp(k,j) = exp(log_disp);
      if (disp_interp(k,j) > max_nu) { disp_interp(k,j) = max_nu; }
      if (disp_interp(k,j) < min_nu) { disp_interp(k,j) = min_nu; }
    }
  }  // end loop over items
  
  NumericMatrix log_Z(K, m);
  NumericMatrix log_lambda(K, m);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes*persons and
  // as many columns as we have 
  // they have the same structure as mu_interp and nu_interp matrices above 
  // (where the structure is explained in more detail)
  
  double out = 0;
  for (int k=0; k<K; k++) { // nodes
    for(int i=0;i<N;i++) { // persons
      for(int j=0;j<m;j++) { // items
        out += (data(i,j)*log_lambda(k,j) - log_Z(k,j) - 
          disp_interp(k,j)*logFactorial(data(i,j))) * PPs(i,k);
      }
    }
  } // end loops over K nodes
  
  return(out);
}

// [[Rcpp::export]]
double ell_cmp_with_icov_alpha_nu_cpp (double alpha,
                                  NumericVector deltas,
                                  double disp,
                                  NumericVector betas_alpha,
                                  NumericVector betas_logdisp,
                                  NumericMatrix data,
                                  NumericMatrix i_cov_data,
                                  NumericVector PPs,
                                  NumericVector weights,
                                  NumericVector nodes,
                                  NumericVector grid_mus,
                                  NumericVector grid_nus,
                                  NumericVector grid_cmp_var_long,
                                  NumericVector grid_log_lambda_long,
                                  NumericVector grid_logZ_long,
                                  double max_mu,
                                  double min_mu,
                                  double max_nu,
                                  double min_nu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int K = nodes.size();
  int m = data.ncol();
  int N = data.nrow();
  int I = betas_alpha.size();
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for item covariates, we can do normal mu_interp and nu_interp, so just
  // of dumension KxM, they don't need to be person specific
  NumericMatrix mu(K, m);
  NumericMatrix mu_interp(K, m);
  NumericMatrix disp_interp(K, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<K;k++) {
      // loop over nodes (rows)
      double log_mu = alpha * nodes[k] + deltas[j];
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += nodes[k] * betas_alpha[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      double log_disp = log(disp);
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_disp += betas_logdisp[c] * i_cov_data(j,c); // for item j
      }
      disp_interp(k,j) = exp(log_disp);
      if (disp_interp(k,j) > max_nu) { disp_interp(k,j) = max_nu; }
      if (disp_interp(k,j) < min_nu) { disp_interp(k,j) = min_nu; }
    }
  }  // end loop over items
  
  NumericMatrix log_Z(K, m);
  NumericMatrix log_lambda(K, m);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes*persons and
  // as many columns as we have 
  // they have the same structure as mu_interp and nu_interp matrices above 
  // (where the structure is explained in more detail)
  
  double out = 0;
  for (int k=0; k<K; k++) { // nodes
    for(int i=0;i<N;i++) { // persons
      for(int j=0;j<m;j++) { // items
        out += (data(i,j)*log_lambda(k,j) - log_Z(k,j) - 
          disp_interp(k,j)*logFactorial(data(i,j))) * PPs(i,k);
      }
    }
  } // end loops over K nodes
  
  return(out);
}

// [[Rcpp::export]]
double ell_cmp_with_icov_delta_nu_cpp (NumericVector alphas,
                                  double delta,
                                  double disp,
                                  NumericVector betas_delta,
                                  NumericVector betas_logdisp,
                                  NumericMatrix data,
                                  NumericMatrix i_cov_data,
                                  NumericVector PPs,
                                  NumericVector weights,
                                  NumericVector nodes,
                                  NumericVector grid_mus,
                                  NumericVector grid_nus,
                                  NumericVector grid_cmp_var_long,
                                  NumericVector grid_log_lambda_long,
                                  NumericVector grid_logZ_long,
                                  double max_mu,
                                  double min_mu,
                                  double max_nu,
                                  double min_nu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int K = nodes.size();
  int m = data.ncol();
  int N = data.nrow();
  int I = betas_delta.size();
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for item covariates, we can do normal mu_interp and nu_interp, so just
  // of dumension KxM, they don't need to be person specific
  NumericMatrix mu(K, m);
  NumericMatrix mu_interp(K, m);
  NumericMatrix disp_interp(K, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<K;k++) {
      // loop over nodes (rows)
      double log_mu = alphas[j] * nodes[k] + delta;
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += betas_delta[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      double log_disp = log(disp);
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_disp += betas_logdisp[c] * i_cov_data(j,c); // for item j
      }
      disp_interp(k,j) = exp(log_disp);
      if (disp_interp(k,j) > max_nu) { disp_interp(k,j) = max_nu; }
      if (disp_interp(k,j) < min_nu) { disp_interp(k,j) = min_nu; }
    }
  }  // end loop over items
  
  NumericMatrix log_Z(K, m);
  NumericMatrix log_lambda(K, m);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes*persons and
  // as many columns as we have 
  // they have the same structure as mu_interp and nu_interp matrices above 
  // (where the structure is explained in more detail)
  
  double out = 0;
  for (int k=0; k<K; k++) { // nodes
    for(int i=0;i<N;i++) { // persons
      for(int j=0;j<m;j++) { // items
        out += (data(i,j)*log_lambda(k,j) - log_Z(k,j) - 
          disp_interp(k,j)*logFactorial(data(i,j))) * PPs(i,k);
      }
    }
  } // end loops over K nodes
  
  return(out);
}

// [[Rcpp::export]]
double ell_cmp_with_icov_alpha_delta_cpp (double alpha,
                                  double delta,
                                  NumericVector disps,
                                  NumericVector betas_alpha,
                                  NumericVector betas_delta,
                                  NumericMatrix data,
                                  NumericMatrix i_cov_data,
                                  NumericVector PPs,
                                  NumericVector weights,
                                  NumericVector nodes,
                                  NumericVector grid_mus,
                                  NumericVector grid_nus,
                                  NumericVector grid_cmp_var_long,
                                  NumericVector grid_log_lambda_long,
                                  NumericVector grid_logZ_long,
                                  double max_mu,
                                  double min_mu) {
  
  // r needs to be a matrix with one column per item and then the r values
  // for this item in the column
  // analogously for f and h
  
  int K = nodes.size();
  int m = data.ncol();
  int N = data.nrow();
  int I = betas_alpha.size();
  
  // set up mu's and nu's for interpolation function to be computed all in one
  
  // for item covariates, we can do normal mu_interp and nu_interp, so just
  // of dumension KxM, they don't need to be person specific
  NumericMatrix mu(K, m);
  NumericMatrix mu_interp(K, m);
  NumericMatrix disp_interp(K, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<K;k++) {
      // loop over nodes (rows)
      double log_mu = alpha * nodes[k] + delta;
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += nodes[k] * betas_alpha[c] * i_cov_data(j,c) + 
          betas_delta[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,j) = disps[j];
    }
  }  // end loop over items
  
  NumericMatrix log_Z(K, m);
  NumericMatrix log_lambda(K, m);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes*persons and
  // as many columns as we have 
  // they have the same structure as mu_interp and nu_interp matrices above 
  // (where the structure is explained in more detail)
  
  double out = 0;
  for (int k=0; k<K; k++) { // nodes
    for(int i=0;i<N;i++) { // persons
      for(int j=0;j<m;j++) { // items
        out += (data(i,j)*log_lambda(k,j) - log_Z(k,j) - 
          disp_interp(k,j)*logFactorial(data(i,j))) * PPs(i,k);
      }
    }
  } // end loops over K nodes
  
  return(out);
}

// [[Rcpp::export]]
NumericVector e_values_newem_cpp (NumericMatrix data,
                                  NumericVector alphas,
                                  NumericVector deltas,
                                  NumericVector disps,
                                  NumericVector nodes,
                                  NumericVector weights,
                                  NumericVector grid_mus,
                                  NumericVector grid_nus,
                                  NumericVector grid_logZ_long,
                                  NumericVector grid_log_lambda_long,
                                  double max_mu,
                                  double min_mu) {

  int m = alphas.size();
  int n_nodes = nodes.size();
  int N = data.nrow();

  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int i=0;i<m;i++){
    // loop over items (columns)
    for(int j=0;j<n_nodes;j++) {
      // loop over nodes (rows)
      mu(j,i) = exp(alphas[i] * nodes[j] + deltas[i]);
      mu_interp(j,i) = mu(j,i);
      if (mu(j,i) > max_mu) { mu_interp(j,i) = max_mu; }
      if (mu(j,i) < min_mu) { mu_interp(j,i) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(j,i) = disps[i];
    }
  }

  NumericMatrix log_Z(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many columns as we have nodes and
  // as many columns as we have items

  NumericVector exp_abilities(N);
  NumericVector marg_prob(N);

  for(int i=0;i<N;i++){
    // the exp_abilities vector has as many elements as persons
    // and is the weighted mean of nodes (weighted with posterior probabilities)

    // first, we need to compute the posterior probabilities

    // to this end, we need the marginal probabilities (for the denominator)
    // these are person specific (and summed over nodes and probabilities
    // are products over items because they are probs for response vectors,
    // i.e., responses of one person to all items)
    marg_prob(i) = 0;
    NumericVector log_resp_vector_prob(n_nodes);
    for (int k=0;k<n_nodes;k++){
      log_resp_vector_prob(k) = 0;
      for (int j=0;j<m;j++) {
        log_resp_vector_prob(k) += data(i,j)*log_lambda(k,j) -
          log_Z(k,j) - disps[j]*lgamma(data(i,j)+1);
      }
      marg_prob(i) += exp(log_resp_vector_prob(k) + log(weights[k]));
    }

    // compute the numerators and then the posterior probs
    // which are person and node specific (because the numerators are node specific)
    // and prep for the computation of the post. prob weighted mean
    NumericVector post_prob_i(n_nodes);
    double sum_across_nodes_i = 0; // numerator for weighted mean
    double sum_across_post_probs_i = 0; // denominator for weighted mean
    for (int k=0;k<n_nodes;k++){
      post_prob_i(k) = (exp(log_resp_vector_prob(k) + log(weights[k]))) / marg_prob(i);
      sum_across_nodes_i += post_prob_i(k) * nodes[k];
      sum_across_post_probs_i += post_prob_i(k);
    }

    // compute the posterior prob weighted mean across nodes for person i
    exp_abilities(i) = sum_across_nodes_i / sum_across_post_probs_i;
  }

  return(exp_abilities);
}


// [[Rcpp::export]]
NumericMatrix e_values_newem_cpp2(NumericMatrix data,
                                  NumericVector alphas,
                                  NumericVector deltas,
                                  NumericVector disps,
                                  NumericVector nodes,
                                  NumericVector weights,
                                  NumericVector grid_mus,
                                  NumericVector grid_nus,
                                  NumericVector grid_logZ_long,
                                  NumericVector grid_log_lambda_long,
                                  double max_mu,
                                  double min_mu) {
  
  int m = alphas.size();
  int n_nodes = nodes.size();
  int N = data.nrow();
  
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int i=0;i<m;i++){
    // loop over items (columns)
    for(int j=0;j<n_nodes;j++) {
      // loop over nodes (rows)
      mu(j,i) = exp(alphas[i] * nodes[j] + deltas[i]);
      mu_interp(j,i) = mu(j,i);
      if (mu(j,i) > max_mu) { mu_interp(j,i) = max_mu; }
      if (mu(j,i) < min_mu) { mu_interp(j,i) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(j,i) = disps[i];
    }
  }
  
  NumericMatrix log_Z(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many columns as we have nodes and
  // as many columns as we have items
  
  // NumericVector exp_abilities(N);
  NumericVector marg_prob(N);
  NumericMatrix PPs(N, n_nodes);
  
  for(int i=0;i<N;i++){
    // compute the marginal probability for each person 
    // (which we need for the denominator of the posterior probabilities)
    marg_prob(i) = 0;
    NumericVector log_resp_vector_prob(n_nodes);
    for (int k=0;k<n_nodes;k++){
      log_resp_vector_prob(k) = 0;
      for (int j=0;j<m;j++) {
        log_resp_vector_prob(k) += data(i,j)*log_lambda(k,j) -
          log_Z(k,j) - disps[j]*lgamma(data(i,j)+1);
      }
      marg_prob(i) += exp(log_resp_vector_prob(k) + log(weights[k]));
    }
    
    // compute the numerators and then the posterior probs
    // which are person and node specific (because the numerators are node specific)
    for (int k=0;k<n_nodes;k++){
      PPs(i, k) = (exp(log_resp_vector_prob(k) + log(weights[k]))) / marg_prob(i);
      // sum_across_nodes_i += post_prob_i(k) * nodes[k];
      // sum_across_post_probs_i += post_prob_i(k);
    }
    
    // compute the posterior prob weighted mean across nodes for person i
    // exp_abilities(i) = sum_across_nodes_i / sum_across_post_probs_i;
  }
  return(PPs);
}

// [[Rcpp::export]]
NumericMatrix estep_cmp_with_icov_delta_cpp(NumericMatrix data,
                                  NumericVector alphas,
                                  double delta,
                                  NumericVector disps,
                                  NumericVector betas,
                                  NumericMatrix i_cov_data,
                                  NumericVector nodes,
                                  NumericVector weights,
                                  NumericVector grid_mus,
                                  NumericVector grid_nus,
                                  NumericVector grid_logZ_long,
                                  NumericVector grid_log_lambda_long,
                                  double max_mu,
                                  double min_mu) {
  
  int m = alphas.size();
  int n_nodes = nodes.size();
  int N = data.nrow();
  int I = betas.size();
  
  // for item covariates we don't need person specificness (as we do for the person covariates)
  // so our mu_interp and disp_interp can just be of the dimension KxM
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alphas[j] * nodes[k] + delta; // deltas is a scalar for item covariates
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        // (for the specific item j we are currently looking at)
        log_mu += betas[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,j) = disps[j];
    }
  }  // end loop over items
  
  NumericMatrix log_Z(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have 
  
  NumericVector marg_prob(N);
  NumericMatrix PPs(N, n_nodes);
  
  for(int i=0;i<N;i++){
    // compute the marginal probability for each person 
    // (which we need for the denominator of the posterior probabilities)
    marg_prob(i) = 0;
    NumericVector log_resp_vector_prob(n_nodes); // created anew for each person i
    for (int k=0;k<n_nodes;k++){
      log_resp_vector_prob(k) = 0;
      for (int j=0;j<m;j++) {
        // when we access nodes here note that we need to access the nodes for person i
        // here because our lambda and Z values are not only node and item specific but
        // also additionally person specific
        log_resp_vector_prob(k) += data(i,j)*log_lambda(k,j) -
          log_Z(k,j) - disps[j]*lgamma(data(i,j)+1);
      }
      marg_prob(i) += exp(log_resp_vector_prob(k) + log(weights[k]));
    }
    
    // compute the numerators and then the posterior probs
    // which are person and node specific (because the numerators are node specific)
    for (int k=0;k<n_nodes;k++){
      PPs(i, k) = (exp(log_resp_vector_prob(k) + log(weights[k]))) / marg_prob(i);
    }
  }
  return(PPs);
}

// [[Rcpp::export]]
NumericMatrix estep_cmp_with_icov_alpha_cpp(NumericMatrix data,
                                            double alpha,
                                            NumericVector deltas,
                                            NumericVector disps,
                                            NumericVector betas,
                                            NumericMatrix i_cov_data,
                                            NumericVector nodes,
                                            NumericVector weights,
                                            NumericVector grid_mus,
                                            NumericVector grid_nus,
                                            NumericVector grid_logZ_long,
                                            NumericVector grid_log_lambda_long,
                                            double max_mu,
                                            double min_mu) {
  
  int m = data.ncol();
  int n_nodes = nodes.size();
  int N = data.nrow();
  int I = betas.size();
  
  // for item covariates we don't need person specificness (as we do for the person covariates)
  // so our mu_interp and disp_interp can just be of the dimension KxM
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alpha * nodes[k] + deltas[j]; 
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        // (for the specific item j we are currently looking at)
        log_mu += nodes[k] * betas[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,j) = disps[j];
    }
  }  // end loop over items
  
  NumericMatrix log_Z(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have 
  
  NumericVector marg_prob(N);
  NumericMatrix PPs(N, n_nodes);
  
  for(int i=0;i<N;i++){
    // compute the marginal probability for each person 
    // (which we need for the denominator of the posterior probabilities)
    marg_prob(i) = 0;
    NumericVector log_resp_vector_prob(n_nodes); // created anew for each person i
    for (int k=0;k<n_nodes;k++){
      log_resp_vector_prob(k) = 0;
      for (int j=0;j<m;j++) {
        // when we access nodes here note that we need to access the nodes for person i
        // here because our lambda and Z values are not only node and item specific but
        // also additionally person specific
        log_resp_vector_prob(k) += data(i,j)*log_lambda(k,j) -
          log_Z(k,j) - disps[j]*lgamma(data(i,j)+1);
      }
      marg_prob(i) += exp(log_resp_vector_prob(k) + log(weights[k]));
    }
    
    // compute the numerators and then the posterior probs
    // which are person and node specific (because the numerators are node specific)
    for (int k=0;k<n_nodes;k++){
      PPs(i, k) = (exp(log_resp_vector_prob(k) + log(weights[k]))) / marg_prob(i);
    }
  }
  return(PPs);
}

// [[Rcpp::export]]
NumericMatrix estep_cmp_with_icov_nu_cpp(NumericMatrix data,
                                            NumericVector alphas,
                                            NumericVector deltas,
                                            double disp,
                                            NumericVector betas,
                                            NumericMatrix i_cov_data,
                                            NumericVector nodes,
                                            NumericVector weights,
                                            NumericVector grid_mus,
                                            NumericVector grid_nus,
                                            NumericVector grid_logZ_long,
                                            NumericVector grid_log_lambda_long,
                                            double max_mu,
                                            double min_mu,
                                            double max_nu,
                                            double min_nu) {
  
  int m = data.ncol();
  int n_nodes = nodes.size();
  int N = data.nrow();
  int I = betas.size();
  
  // for item covariates we don't need person specificness (as we do for the person covariates)
  // so our mu_interp and disp_interp can just be of the dimension KxM
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alphas[j] * nodes[k] + deltas[j];
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      double log_disp = log(disp);
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_disp += betas[c] * i_cov_data(j,c); // for item j
      }
      disp_interp(k,j) = exp(log_disp);
      if (disp_interp(k,j) > max_nu) { disp_interp(k,j) = max_nu; }
      if (disp_interp(k,j) < min_nu) { disp_interp(k,j) = min_nu; }
    }
  }  // end loop over items
  
  NumericMatrix log_Z(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have 
  
  NumericVector marg_prob(N);
  NumericMatrix PPs(N, n_nodes);
  
  for(int i=0;i<N;i++){
    // compute the marginal probability for each person 
    // (which we need for the denominator of the posterior probabilities)
    marg_prob(i) = 0;
    NumericVector log_resp_vector_prob(n_nodes); // created anew for each person i
    for (int k=0;k<n_nodes;k++){
      log_resp_vector_prob(k) = 0;
      for (int j=0;j<m;j++) {
        // when we access nodes here note that we need to access the nodes for person i
        // here because our lambda and Z values are not only node and item specific but
        // also additionally person specific
        log_resp_vector_prob(k) += data(i,j)*log_lambda(k,j) -
          log_Z(k,j) - disp_interp(k,j)*lgamma(data(i,j)+1);
      }
      marg_prob(i) += exp(log_resp_vector_prob(k) + log(weights[k]));
    }
    
    // compute the numerators and then the posterior probs
    // which are person and node specific (because the numerators are node specific)
    for (int k=0;k<n_nodes;k++){
      PPs(i, k) = (exp(log_resp_vector_prob(k) + log(weights[k]))) / marg_prob(i);
    }
  }
  return(PPs);
}

// [[Rcpp::export]]
NumericMatrix estep_cmp_with_icov_all_cpp(NumericMatrix data,
                                         double alpha,
                                         double delta,
                                         double disp,
                                         NumericVector betas_alpha,
                                         NumericVector betas_delta,
                                         NumericVector betas_logdisp,
                                         NumericMatrix i_cov_data,
                                         NumericVector nodes,
                                         NumericVector weights,
                                         NumericVector grid_mus,
                                         NumericVector grid_nus,
                                         NumericVector grid_logZ_long,
                                         NumericVector grid_log_lambda_long,
                                         double max_mu,
                                         double min_mu,
                                         double max_nu,
                                         double min_nu) {
  
  int m = data.ncol();
  int n_nodes = nodes.size();
  int N = data.nrow();
  int I = betas_alpha.size();
  
  // for item covariates we don't need person specificness (as we do for the person covariates)
  // so our mu_interp and disp_interp can just be of the dimension KxM
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alpha * nodes[k] + delta;
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += nodes[k] * betas_alpha[c] * i_cov_data(j,c) + 
          betas_delta[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      double log_disp = log(disp);
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_disp += betas_logdisp[c] * i_cov_data(j,c); // for item j
      }
      disp_interp(k,j) = exp(log_disp);
      if (disp_interp(k,j) > max_nu) { disp_interp(k,j) = max_nu; }
      if (disp_interp(k,j) < min_nu) { disp_interp(k,j) = min_nu; }
    }
  }  // end loop over items
  
  NumericMatrix log_Z(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have 
  
  NumericVector marg_prob(N);
  NumericMatrix PPs(N, n_nodes);
  
  for(int i=0;i<N;i++){
    // compute the marginal probability for each person 
    // (which we need for the denominator of the posterior probabilities)
    marg_prob(i) = 0;
    NumericVector log_resp_vector_prob(n_nodes); // created anew for each person i
    for (int k=0;k<n_nodes;k++){
      log_resp_vector_prob(k) = 0;
      for (int j=0;j<m;j++) {
        // when we access nodes here note that we need to access the nodes for person i
        // here because our lambda and Z values are not only node and item specific but
        // also additionally person specific
        log_resp_vector_prob(k) += data(i,j)*log_lambda(k,j) -
          log_Z(k,j) - disp_interp(k,j)*lgamma(data(i,j)+1);
      }
      marg_prob(i) += exp(log_resp_vector_prob(k) + log(weights[k]));
    }
    
    // compute the numerators and then the posterior probs
    // which are person and node specific (because the numerators are node specific)
    for (int k=0;k<n_nodes;k++){
      PPs(i, k) = (exp(log_resp_vector_prob(k) + log(weights[k]))) / marg_prob(i);
    }
  }
  return(PPs);
}

// [[Rcpp::export]]
NumericMatrix estep_cmp_with_icov_alpha_nu_cpp(NumericMatrix data,
                                          double alpha,
                                          NumericVector deltas,
                                          double disp,
                                          NumericVector betas_alpha,
                                          NumericVector betas_logdisp,
                                          NumericMatrix i_cov_data,
                                          NumericVector nodes,
                                          NumericVector weights,
                                          NumericVector grid_mus,
                                          NumericVector grid_nus,
                                          NumericVector grid_logZ_long,
                                          NumericVector grid_log_lambda_long,
                                          double max_mu,
                                          double min_mu,
                                          double max_nu,
                                          double min_nu) {
  
  int m = data.ncol();
  int n_nodes = nodes.size();
  int N = data.nrow();
  int I = betas_alpha.size();
  
  // for item covariates we don't need person specificness (as we do for the person covariates)
  // so our mu_interp and disp_interp can just be of the dimension KxM
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alpha * nodes[k] + deltas[j];
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += nodes[k] * betas_alpha[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      double log_disp = log(disp);
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_disp += betas_logdisp[c] * i_cov_data(j,c); // for item j
      }
      disp_interp(k,j) = exp(log_disp);
      if (disp_interp(k,j) > max_nu) { disp_interp(k,j) = max_nu; }
      if (disp_interp(k,j) < min_nu) { disp_interp(k,j) = min_nu; }
    }
  }  // end loop over items
  
  NumericMatrix log_Z(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have 
  
  NumericVector marg_prob(N);
  NumericMatrix PPs(N, n_nodes);
  
  for(int i=0;i<N;i++){
    // compute the marginal probability for each person 
    // (which we need for the denominator of the posterior probabilities)
    marg_prob(i) = 0;
    NumericVector log_resp_vector_prob(n_nodes); // created anew for each person i
    for (int k=0;k<n_nodes;k++){
      log_resp_vector_prob(k) = 0;
      for (int j=0;j<m;j++) {
        // when we access nodes here note that we need to access the nodes for person i
        // here because our lambda and Z values are not only node and item specific but
        // also additionally person specific
        log_resp_vector_prob(k) += data(i,j)*log_lambda(k,j) -
          log_Z(k,j) - disp_interp(k,j)*lgamma(data(i,j)+1);
      }
      marg_prob(i) += exp(log_resp_vector_prob(k) + log(weights[k]));
    }
    
    // compute the numerators and then the posterior probs
    // which are person and node specific (because the numerators are node specific)
    for (int k=0;k<n_nodes;k++){
      PPs(i, k) = (exp(log_resp_vector_prob(k) + log(weights[k]))) / marg_prob(i);
    }
  }
  return(PPs);
}

// [[Rcpp::export]]
NumericMatrix estep_cmp_with_icov_delta_nu_cpp(NumericMatrix data,
                                          NumericVector alphas,
                                          double delta,
                                          double disp,
                                          NumericVector betas_delta,
                                          NumericVector betas_logdisp,
                                          NumericMatrix i_cov_data,
                                          NumericVector nodes,
                                          NumericVector weights,
                                          NumericVector grid_mus,
                                          NumericVector grid_nus,
                                          NumericVector grid_logZ_long,
                                          NumericVector grid_log_lambda_long,
                                          double max_mu,
                                          double min_mu,
                                          double max_nu,
                                          double min_nu) {
  
  int m = data.ncol();
  int n_nodes = nodes.size();
  int N = data.nrow();
  int I = betas_delta.size();
  
  // for item covariates we don't need person specificness (as we do for the person covariates)
  // so our mu_interp and disp_interp can just be of the dimension KxM
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alphas[j] * nodes[k] + delta;
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += betas_delta[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      double log_disp = log(disp);
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_disp += betas_logdisp[c] * i_cov_data(j,c); // for item j
      }
      disp_interp(k,j) = exp(log_disp);
      if (disp_interp(k,j) > max_nu) { disp_interp(k,j) = max_nu; }
      if (disp_interp(k,j) < min_nu) { disp_interp(k,j) = min_nu; }
    }
  }  // end loop over items
  
  NumericMatrix log_Z(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have 
  
  NumericVector marg_prob(N);
  NumericMatrix PPs(N, n_nodes);
  
  for(int i=0;i<N;i++){
    // compute the marginal probability for each person 
    // (which we need for the denominator of the posterior probabilities)
    marg_prob(i) = 0;
    NumericVector log_resp_vector_prob(n_nodes); // created anew for each person i
    for (int k=0;k<n_nodes;k++){
      log_resp_vector_prob(k) = 0;
      for (int j=0;j<m;j++) {
        // when we access nodes here note that we need to access the nodes for person i
        // here because our lambda and Z values are not only node and item specific but
        // also additionally person specific
        log_resp_vector_prob(k) += data(i,j)*log_lambda(k,j) -
          log_Z(k,j) - disp_interp(k,j)*lgamma(data(i,j)+1);
      }
      marg_prob(i) += exp(log_resp_vector_prob(k) + log(weights[k]));
    }
    
    // compute the numerators and then the posterior probs
    // which are person and node specific (because the numerators are node specific)
    for (int k=0;k<n_nodes;k++){
      PPs(i, k) = (exp(log_resp_vector_prob(k) + log(weights[k]))) / marg_prob(i);
    }
  }
  return(PPs);
}

// [[Rcpp::export]]
NumericMatrix estep_cmp_with_icov_alpha_delta_cpp(NumericMatrix data,
                                          double alpha,
                                          double delta,
                                          NumericVector disps,
                                          NumericVector betas_alpha,
                                          NumericVector betas_delta,
                                          NumericMatrix i_cov_data,
                                          NumericVector nodes,
                                          NumericVector weights,
                                          NumericVector grid_mus,
                                          NumericVector grid_nus,
                                          NumericVector grid_logZ_long,
                                          NumericVector grid_log_lambda_long,
                                          double max_mu,
                                          double min_mu,
                                          double max_nu,
                                          double min_nu) {
  
  int m = data.ncol();
  int n_nodes = nodes.size();
  int N = data.nrow();
  int I = betas_alpha.size();
  
  // for item covariates we don't need person specificness (as we do for the person covariates)
  // so our mu_interp and disp_interp can just be of the dimension KxM
  NumericMatrix mu(n_nodes, m);
  NumericMatrix mu_interp(n_nodes, m);
  NumericMatrix disp_interp(n_nodes, m);
  for(int j=0;j<m;j++){
    // loop over items (columns)
    for(int k=0;k<n_nodes;k++) {
      // loop over nodes (rows)
      double log_mu = alpha * nodes[k] + delta;
      for(int c=0; c<I; c++) {
        // add all the (weighted) covariate values for all covariates
        log_mu += nodes[k] * betas_alpha[c] * i_cov_data(j,c) + 
          betas_delta[c] * i_cov_data(j,c);
      }
      mu(k,j) = exp(log_mu);
      mu_interp(k,j) = mu(k,j);
      if (mu(k,j) > max_mu) { mu_interp(k,j) = max_mu; }
      if (mu(k,j) < min_mu) { mu_interp(k,j) = min_mu; }
      // we need to set maximum for mu to max_mu so that the interpolation will
      // work, max_mu is the maximum mu value in our grid for interpolation
      disp_interp(k,j) = disps[j];
    }
  }  // end loop over items
  
  NumericMatrix log_Z(n_nodes, m);
  NumericMatrix log_lambda(n_nodes, m);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes and
  // as many columns as we have 
  
  NumericVector marg_prob(N);
  NumericMatrix PPs(N, n_nodes);
  
  for(int i=0;i<N;i++){
    // compute the marginal probability for each person 
    // (which we need for the denominator of the posterior probabilities)
    marg_prob(i) = 0;
    NumericVector log_resp_vector_prob(n_nodes); // created anew for each person i
    for (int k=0;k<n_nodes;k++){
      log_resp_vector_prob(k) = 0;
      for (int j=0;j<m;j++) {
        // when we access nodes here note that we need to access the nodes for person i
        // here because our lambda and Z values are not only node and item specific but
        // also additionally person specific
        log_resp_vector_prob(k) += data(i,j)*log_lambda(k,j) -
          log_Z(k,j) - disp_interp(k,j)*lgamma(data(i,j)+1);
      }
      marg_prob(i) += exp(log_resp_vector_prob(k) + log(weights[k]));
    }
    
    // compute the numerators and then the posterior probs
    // which are person and node specific (because the numerators are node specific)
    for (int k=0;k<n_nodes;k++){
      PPs(i, k) = (exp(log_resp_vector_prob(k) + log(weights[k]))) / marg_prob(i);
    }
  }
  return(PPs);
}

// [[Rcpp::export]]
NumericMatrix estep_cmp_with_pcov_cpp(NumericMatrix data,
                                     NumericVector alphas,
                                     NumericVector deltas,
                                     NumericVector disps,
                                     NumericVector betas,
                                     NumericMatrix p_cov_data,
                                     NumericVector nodes,
                                     NumericVector weights,
                                     NumericVector grid_mus,
                                     NumericVector grid_nus,
                                     NumericVector grid_logZ_long,
                                     NumericVector grid_log_lambda_long,
                                     double max_mu,
                                     double min_mu) {
  
  int m = alphas.size();
  int n_nodes = nodes.size();
  int N = data.nrow();
  int P = betas.size();
  
  // for person covariates, we need mus (and lambdas and Zs) which are person
  // as well as node and item specific
  // so then I extend my nu and mu matrices for interpolation accordingly
  // so that i can interpolate lambdas and Zs person-node-item specifically
  // but still only work with matrices so that i can use interp_from_grid_m
  // here i just chain KxM matrices (like I had for no covariates) below each other
  // (so rbind basically), for all N person so that for the first K rows,
  // we have the KxM matrix for person 1, for rows K+1 - K+K we have the
  // KxM matrix for person 1, etc.
  NumericMatrix mu(n_nodes*N, m);
  NumericMatrix mu_interp(n_nodes*N, m);
  NumericMatrix disp_interp(n_nodes*N, m);
  for (int i=0; i<N; i++) {
    // we are computing node-item specific mus for each person
    for(int j=0;j<m;j++){
      // loop over items (columns)
      for(int k=0;k<n_nodes;k++) {
        // loop over nodes (rows)
        double log_mu = alphas[j] * nodes[k] + deltas[j];
        // my nodes are here my epsilon
        for(int p=0; p<P; p++) {
          // add all the (weighted) covariate values for that specific item j and person i
          log_mu += betas[p] * alphas[j] * p_cov_data(i,p);
        }
        mu(k+i*n_nodes,j) = exp(log_mu);
        mu_interp(k+i*n_nodes,j) = mu(k+i*n_nodes,j);
        if (mu(k+i*n_nodes,j) > max_mu) { mu_interp(k+i*n_nodes,j) = max_mu; }
        if (mu(k+i*n_nodes,j) < min_mu) { mu_interp(k+i*n_nodes,j) = min_mu; }
        // we need to set maximum for mu to max_mu so that the interpolation will
        // work, max_mu is the maximum mu value in our grid for interpolation
        disp_interp(k+i*n_nodes,j) = disps[j];
      }
    }  // end loop over items
  } // end loop over N
  
  NumericMatrix log_Z(n_nodes*N, m);
  NumericMatrix log_lambda(n_nodes*N, m);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  // V and log_lambda are matrices with as many rows as we have nodes*persons and
  // as many columns as we have 
  // they have the same structure as mu_interp and nu_interp matrices above 
  // (where the structure is explained in more detail)
  
  NumericVector marg_prob(N);
  NumericMatrix PPs(N, n_nodes);
  
  for(int i=0;i<N;i++){
    // compute the marginal probability for each person 
    // (which we need for the denominator of the posterior probabilities)
    marg_prob(i) = 0;
    NumericVector log_resp_vector_prob(n_nodes); // created anew for each person i
    for (int k=0;k<n_nodes;k++){
      log_resp_vector_prob(k) = 0;
      for (int j=0;j<m;j++) {
        // when we access nodes here note that we need to access the nodes for person i
        // here because our lambda and Z values are not only node and item specific but
        // also additionally person specific
        log_resp_vector_prob(k) += data(i,j)*log_lambda(k+i*n_nodes,j) -
          log_Z(k+i*n_nodes,j) - disps[j]*lgamma(data(i,j)+1);
      }
      marg_prob(i) += exp(log_resp_vector_prob(k) + log(weights[k]));
    }
    
    // compute the numerators and then the posterior probs
    // which are person and node specific (because the numerators are node specific)
    for (int k=0;k<n_nodes;k++){
      PPs(i, k) = (exp(log_resp_vector_prob(k) + log(weights[k]))) / marg_prob(i);
    }
  }
  return(PPs);
}


// [[Rcpp::export]]
NumericMatrix estep_cmp_with_pcov_cat_cpp(NumericMatrix data,
                                      NumericVector alphas,
                                      NumericVector deltas,
                                      NumericVector disps,
                                      NumericVector betas,
                                      NumericMatrix p_cov_data,
                                      NumericMatrix resp_pattern,
                                      NumericVector nodes,
                                      NumericVector weights,
                                      NumericVector grid_mus,
                                      NumericVector grid_nus,
                                      NumericVector grid_logZ_long,
                                      NumericVector grid_log_lambda_long,
                                      double max_mu,
                                      double min_mu) {
  
  // assume that p_cov is a matrix of dummy coded categorical predictors
  // resp_pattern is a matrix of the same no. of cols than p_cov
  // and as many rows as we have distinct possible response patterns
  
  int N = data.nrow();
  int M = data.ncol();
  int K = nodes.size();
  int P = betas.size(); 
  int n_resp_patterns = resp_pattern.nrow();
  
  // for person covariates, we need mus (and lambdas and Zs) for each node and
  // and then also for each response pattern
  // so first compute that
  NumericMatrix mu(K*n_resp_patterns, M);
  NumericMatrix mu_interp(K*n_resp_patterns, M);
  NumericMatrix disp_interp(K*n_resp_patterns, M);
  for (int l=0; l<n_resp_patterns; l++) {
    for(int j=0; j<M; j++){
      // loop over items (columns)
      for(int k=0; k<K; k++) {
        // loop over nodes (rows)
        double log_mu = alphas[j] * nodes[k] + deltas[j];
        for(int p=0; p<P; p++) {
          // this works because only includes columns for none-reference categories
          // for all covs in ref categories, resp_pattern will just always be zero in that row
          log_mu += betas[p] * alphas[j] * resp_pattern(l,p);
        }
        
        mu(k+l*K,j) = exp(log_mu);
        mu_interp(k+l*K,j) = mu(k+l*K,j);
        if (mu(k+l*K,j) > max_mu) { mu_interp(k+l*K,j) = max_mu; }
        if (mu(k+l*K,j) < min_mu) { mu_interp(k+l*K,j) = min_mu; }
        // we need to set maximum for mu to max_mu so that the interpolation will
        // work, max_mu is the maximum mu value in our grid for interpolation
        disp_interp(k+l*K,j) = disps[j];
      }
    }  // end loop over items
  }
  
  NumericMatrix log_lambda(K*n_resp_patterns, M);
  NumericMatrix log_Z(K*n_resp_patterns, M);
  log_lambda = interp_from_grid_m(grid_mus, grid_nus,
                                  grid_log_lambda_long,
                                  mu_interp, disp_interp);
  log_Z = interp_from_grid_m(grid_mus, grid_nus,
                             grid_logZ_long,
                             mu_interp, disp_interp);
  
  NumericVector marg_prob(N);
  NumericMatrix PPs(N, K);
  
  for(int i=0;i<N;i++){
    // compute the marginal probability for each person 
    // (which we need for the denominator of the posterior probabilities)
    
    // check what response pattern person i had
    int l = 0;
    bool pattern_match = false;
    while (!pattern_match && l<n_resp_patterns) {
      // the second condition is just for safety that we dont get an eternal while loop
      // but we should always find a pattern match
      for (int p=0; p<P; p++) {
        pattern_match = p_cov_data(i,p) == resp_pattern(l,p);
      }
      // if the rows are the same, i am going to get out the foor loop with
      // pattern_match = TRUE and have l at the row of the pattern matrix
      // otherwise I am going to increase l by 1 and stay in the while loop to see if
      // the next row in the pattern matrix is a match for i
      if (!pattern_match) { l += 1; }
    }
    
    // we now know that person i has a response pattern like in row l of resp_pattern matrix
    // so their mu (and lambda, etc.) should be at row k+l*K
    
    marg_prob(i) = 0;
    NumericVector log_resp_vector_prob(K); // created anew for each person i
    for (int k=0;k<K;k++){
      log_resp_vector_prob(k) = 0;
      for (int j=0;j<M;j++) {
        // when we access nodes here note that we need to access the nodes for person i
        // here because our lambda and Z values are not only node and item specific but
        // also additionally person specific
        log_resp_vector_prob(k) += data(i,j)*log_lambda(k+l*K,j) -
          log_Z(k+l*K,j) - disps[j]*lgamma(data(i,j)+1);
      }
      marg_prob(i) += exp(log_resp_vector_prob(k) + log(weights[k]));
    }
    
    // compute the numerators and then the posterior probs
    // which are person and node specific (because the numerators are node specific)
    for (int k=0;k<K;k++){
      PPs(i, k) = (exp(log_resp_vector_prob(k) + log(weights[k]))) / marg_prob(i);
    }
  }
  return(PPs);
}
