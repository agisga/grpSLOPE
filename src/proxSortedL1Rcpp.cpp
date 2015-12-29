// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

//' Fast prox for the Sorted L1 norm
//' 
//' A stack-based algorithm for the prox for the Sorted L1 norm,
//' which solves the problem in O(n) flops.
//'
//' @export
// [[Rcpp::export]]
Eigen::VectorXd proxSortedL1Rcpp(const Eigen::Map<Eigen::VectorXd>& y, const Eigen::Map<Eigen::VectorXd>& lambda)
{
  int n(y.size());
  Eigen::VectorXd x(y-lambda);
  Eigen::VectorXi block_start(n);
  Eigen::VectorXi block_end(n);
  int block_num(0);
  int block_length(1);
  double block_mean(0.0);
  bool decreasing = false;
  for(int k=0; k<n; k++)
  {
    block_num = 0;
    block_start = Eigen::VectorXi::Constant(n, -999);
    block_end = Eigen::VectorXi::Constant(n, -999);
    for (int i=0; i<n; i++)
    {
      block_start(block_num) = i;
      block_end(block_num) = i;
      if ((i>0) && (x(i-1) <= x(i)))
      {
        if(x(i-1) != x(i)) { decreasing = false; }
        block_end(block_num-1) = i;
        block_num--;
      }
      block_num++;
    }

    if(decreasing==true) break;
    decreasing = true;

    for (int j=0; j<block_num; j++)
    {
      block_mean = 0.0;
      block_length = block_end(j) - block_start(j) + 1;
      block_mean = x.segment(block_start(j), block_length).mean();
      x.segment(block_start(j), block_length) = Eigen::VectorXd::Constant(block_length, block_mean);
    }
  }
  for (int k=0; k<n; k++) { if (x(k) < 0.0){ x(k) = 0.0; } }
  return(x);
}
