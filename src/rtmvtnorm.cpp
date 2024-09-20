#include "RcppArmadillo.h"


/*** 
 * The code in this file was copied from the source code of the package 
 * 'tmvtnsim' v1.0.3 by Kaifeng Lu on 2024-09-20 distributed under the GPL-3
 * license and subsequently modified by Tomasz Woźniak.
*/


double norm_rej(const double a, const double b) {
  double x;
  do {
    x = R::rnorm(0.0, 1.0);  
  } while (x < a || x > b);
  return x;
}


double unif_rej(const double a, const double b) {
  double x, u, r;
  do {
    x = R::runif(a, b);
    u = R::runif(0.0, 1.0);
    if (a <= 0.0 && b >= 0.0) 
      r = exp(-x*x/2.0);
    else if (a > 0.0)
      r = exp(-(x*x-a*a)/2.0);
    else 
      r = exp(-(x*x-b*b)/2.0);
  } while (u > r);
  return x;
}


double halfnorm_rej(const double a, const double b) {
  double x;
  do {
    x = fabs(R::rnorm(0.0, 1.0));  
  } while (x < a || x > b);
  return x;
}


double exp_rej(const double a, const double b) {
  double lambda = (a+sqrt(a*a+4.0))/2.0;
  double x, u, r;
  do {
    x = a+R::rweibull(1, 1.0/lambda);
    u = R::runif(0.0, 1.0);
    r = exp(-pow(x-lambda,2)/2.0);
  } while (u > r || x > b);
  return x;
}


arma::vec rtnormcpp(const arma::vec& mean, 
                    const double sd, 
                    const arma::vec& lower, 
                    const arma::vec& upper) {
  const unsigned int n=mean.n_elem, n1=lower.n_elem, n2=upper.n_elem;
  
  // check dimensions
  if (!(n==n1 || n1==1)) {
    Rcpp::stop("lower must be a row vector or have the same number of rows as mean");
  }
  
  if (!(n1==n2)) {
    Rcpp::stop("lower and upper must have the same number of rows");
  }
  
  // check boundary conditions
  if (any(lower >= upper)) {
    Rcpp::stop("lower bound must be smaller than upper bound");
  }
  
  static const double pi = 3.141592653589793238462643383279;
  
  auto imp_case1 = [](double a, double b) { 
    double w;
    if (a < 0.0) w = norm_rej(a, b);
    else if (a < 0.25696) w = halfnorm_rej(a, b);
    else w = exp_rej(a, b);
    return w;
  };
  
  auto imp_case2 = [](double a, double b) { 
    double w;
    if (b <= a + sqrt(2*pi)) w = unif_rej(a, b);
    else w = norm_rej(a, b);
    return w;
  };
  
  auto imp_case3 = [](double a, double b) { 
    double w, lambda;
    if (a <= 0.25696) {
      if (b <= a+sqrt(pi/2)*exp(a*a/2)) w = unif_rej(a,b);
      else w = halfnorm_rej(a,b);
    } else {
      lambda = (a+sqrt(a*a+4.0))/2.0;      
      if (b <= a+1/lambda*exp((a*a-a*sqrt(a*a+4))/4+0.5)) w = unif_rej(a,b);
      else w = exp_rej(a,b);
    }
    return w;
  };
  
  double a, b, w;
  arma::vec x(n);
  for (arma::uword i=0; i<n; i++) {
    if (n1==n) {
      a = (lower(i) - mean(i))/sd;
      b = (upper(i) - mean(i))/sd;
    } else {
      a = (lower(0) - mean(i))/sd;
      b = (upper(0) - mean(i))/sd;
    }
    
    if (std::isinf(a) || std::isinf(b)) {
      if (std::isinf(b)) w = imp_case1(a,b);
      else w = -imp_case1(-b,-a); // case 4
    } else {
      if (a<0.0 && b>0.0) w = imp_case2(a,b);
      else if (a>=0.0) w = imp_case3(a,b);
      else w = -imp_case3(-b,-a); // case 5
    }
    
    x(i) = mean(i)+sd*w;
  }
  return x;
}



arma::mat rtmvnormcpp(const arma::mat& mean, 
                      const arma::mat& sigma, 
                      const arma::mat& blc,
                      const arma::mat& lower, 
                      const arma::mat& upper,
                      const arma::mat& init,
                      const arma::uword burn = 10) {
  const unsigned int n=mean.n_rows, p=mean.n_cols;
  arma::mat x(n,p); // output samples
  
  // draw from truncated univariate normal
  if (p==1) {
    if (blc(0,0) > 0.0) { 
      x.col(0) = rtnormcpp(mean.col(0), sqrt(sigma(0,0)), 
            lower.col(0)/blc(0,0), 
            upper.col(0)/blc(0,0));
    } else if (blc(0,0) < 0.0) {
      x.col(0) = rtnormcpp(mean.col(0), sqrt(sigma(0,0)), 
            upper.col(0)/blc(0,0), 
            lower.col(0)/blc(0,0));
    } else {
      arma::vec z = Rcpp::rnorm(n);
      x.col(0) = z*sqrt(sigma(0,0)) + mean.col(0);   
    }
    return x;
  }
  
  // check boundary conditions
  if (sum(any(lower >= upper)) > 0) {
    Rcpp::stop("lower bound must be smaller than upper bound");
  }

  // check initial values, and generate automatically if needed
  const unsigned int n1=lower.n_rows, m=blc.n_rows;
  arma::mat blct = blc.t();
  arma::mat initc = init*blct;
  arma::mat initx(n1,p);
  if (sum(all(initc >= lower + 1e-8 && initc <= upper - 1e-8)) < m) {
    arma::mat estimate(n1,m);
    for (arma::uword i=0; i<n1; i++) {
      for (arma::uword j=0; j<m; j++) {
        if (std::isinf(lower(i,j)) && std::isinf(upper(i,j))) {
          estimate(i,j) = 0;  
        } else if (std::isinf(lower(i,j))) {
          estimate(i,j) = upper(i,j) - 1e-8;  
        } else if (std::isinf(upper(i,j))) {
          estimate(i,j) = lower(i,j) + 1e-8;  
        } else {
          estimate(i,j) = 0.5*(lower(i,j) + upper(i,j));  
        }
      }
    }
    initx = estimate * arma::pinv(blct);
  }
  
  
  // check whether a matrix has identical rows
  auto f = [](const arma::mat& y) {
    const unsigned int n=y.n_rows, p=y.n_cols; 
    bool identical=1;
    if (n>1) {
      for (arma::uword i=1; i<n; i++) {
        for (arma::uword j=0; j<p; j++) {
          if (y(i,j) != y(i-1,j)) {
            identical = 0;
            break;
          }
        }
        if (identical==0) break;
      }
    }
    return identical;
  };
  
  
  // transform to the problem with identity covariance matrix
  arma::mat cholt = arma::trans(arma::chol(sigma));
  arma::mat R = blc*cholt;
  arma::uvec js(p);
  std::iota(js.begin(), js.end(), 0);
  arma::vec mu(1); mu.fill(0); 
  
  // Gibbs step for sampling truncated multivariate normal
  auto g = [p, m, R, js, mu](arma::vec a, arma::vec b, arma::vec& z) {
    for (arma::uword j=0; j<p; j++) {
      // set up linear inequality constraints for z(j)
      arma::vec rj = R.col(j);
      arma::uvec j2 = arma::find(js != j);
      arma::mat Rj = R.cols(j2);
      arma::vec zj = z(j2);
      arma::vec atemp = a - Rj*zj;
      arma::vec btemp = b - Rj*zj;
      
      // determine lower and upper bounds for z(j)
      double lowerj=R_NegInf, upperj=R_PosInf;
      for (arma::uword k=0; k<m; k++) {
        if (rj(k) != 0) {
          double ak=atemp(k)/rj(k), bk=btemp(k)/rj(k);
          if (rj(k) > 0) {
            if (ak > lowerj) lowerj = ak;
            if (bk < upperj) upperj = bk;
          } else {
            if (bk > lowerj) lowerj = bk;
            if (ak < upperj) upperj = ak;
          }
        }
      }
      
      // generate z(j) for truncated univariate normal
      arma::vec lowerj1(1); lowerj1(0) = lowerj;
      arma::vec upperj1(1); upperj1(0) = upperj;
      z(j) = rtnormcpp(mu, 1, lowerj1, upperj1)(0);
    }
  };
  
  
  arma::vec mean1(p), lower1(m), upper1(m), init1(p);
  arma::vec a(m), b(m), z(p);
  
  // obtain burn + n samples in case of identical means, lower and upper bounds 
  if (f(mean) && f(lower) && f(upper)) {
    mean1 = arma::trans(mean.row(0));
    lower1 = arma::trans(lower.row(0));
    upper1 = arma::trans(upper.row(0));
    init1 = arma::trans(initx.row(0));
    
    a = lower1 - blc*mean1;
    b = upper1 - blc*mean1;
    z = solve(cholt, init1-mean1);
    
    for (arma::uword i=0; i<burn+n; i++) {
      g(a, b, z);
      if (i>=burn) x.row(i-burn) = arma::trans(cholt*z + mean1);
    }
    return x;
  }
  
  // obtain (burn+1)*n samples for non-identical means, lower, or upper bounds
  if (n1==1) {
    lower1 = arma::trans(lower.row(0));
    upper1 = arma::trans(upper.row(0));
    init1 = arma::trans(initx.row(0));
    for (arma::uword i=0; i<n; i++) {
      mean1 = arma::trans(mean.row(i));
      
      a = lower1 - blc*mean1;
      b = upper1 - blc*mean1;
      z = solve(cholt, init1-mean1);
      
      for (arma::uword i2=0; i2<burn+1; i2++) {
        g(a, b, z);
      }
      x.row(i) = arma::trans(cholt*z + mean1);
    }
  } else { // n1==n
    for (arma::uword i=0; i<n; i++) {
      mean1 = arma::trans(mean.row(i));
      lower1 = arma::trans(lower.row(i));
      upper1 = arma::trans(upper.row(i));
      init1 = arma::trans(initx.row(i));

      a = lower1 - blc*mean1;
      b = upper1 - blc*mean1;
      z = solve(cholt, init1-mean1);
      
      for (arma::uword i2=0; i2<burn+1; i2++) {
        g(a, b, z);
      }
      x.row(i) = arma::trans(cholt*z + mean1);
    }
  }
  
  return x;
}

