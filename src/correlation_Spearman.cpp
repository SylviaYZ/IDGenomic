#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// https://stackoverflow.com/questions/30822729/create-ranking-for-vector-of-double
// https://ideone.com/Y3k79K

arma::vec rankAVG(arma::vec v){
  std::vector<std::size_t> w(v.size());
  std::iota(begin(w), end(w), 0);
  std::sort(begin(w), end(w),
            [&v](std::size_t i, std::size_t j) { return v[i] < v[j]; });

  std::vector<double> r(w.size());
  for (std::size_t n, i = 0; i < w.size(); i += n)
  {
    n = 1;
    while (i + n < w.size() && v[w[i]] == v[w[i+n]]) ++n;
    for (std::size_t k = 0; k < n; ++k)
    {
      r[w[i+k]] = i + (n + 1) / 2.0; // average rank of n tied values
    }
  }
  return r;
}

// [[Rcpp::export]]
arma::mat corr_Spearman(arma::mat x, arma::mat y) {
  int xcol = x.n_cols;
  int ycol = y.n_cols;
  //int nRows = x.n_rows;
  //arma::mat rho(xcol, ycol);

  arma::mat rankX(x.n_rows, xcol);
  arma::mat rankY(y.n_rows, ycol);
  for(int i = 0; i < xcol; ++i){
    //arma::vec sampleX = x.col(i);
    rankX.col(i) = rankAVG(x.col(i));
  }
  for(int j = 0; j < ycol; ++j){
    //arma::vec sampleY = y.col(j);
    rankY.col(j) = rankAVG(y.col(j));
    }

  return arma::cor(rankX, rankY);
}

