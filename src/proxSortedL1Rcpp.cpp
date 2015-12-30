#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd proxSortedL1Rcpp(const Eigen::Map<Eigen::VectorXd>& y, const Eigen::Map<Eigen::VectorXd>& lambda)
{
    int n(y.size());
    Eigen::VectorXd x(n);
    Eigen::VectorXi block_start(n);
    Eigen::VectorXi block_end(n);
    Eigen::VectorXd block_sum(n);
    Eigen::VectorXd block_mean(n);
    int block_num(0);
    int block_length(1);
    block_num = 0;
    for(int k=0; k<n; k++)
    {
        block_start(block_num) = k;
        block_end(block_num) = k;
        block_sum(block_num) = y(k) - lambda(k);
        block_mean(block_num) = block_sum(block_num);
        while((block_num > 0) && (block_mean(block_num-1) <= block_mean(block_num)))
        {
            block_end(block_num-1) = block_end(block_num);
            block_sum(block_num-1) += block_sum(block_num);
            block_mean(block_num-1) = block_sum(block_num-1) / (block_end(block_num-1) - block_start(block_num-1) + 1);
            block_num--;
        }
        block_num++;
    }
    for(int j=0; j<block_num; j++)
    {
        block_length = block_end(j) - block_start(j) + 1;
        if(block_mean(j) < 0.0){ block_mean(j) = 0.0; }
        x.segment(block_start(j), block_length) = Eigen::VectorXd::Constant(block_length, block_mean(j));
    }
    return(x);
}
