#include <iostream>
#include <cmath>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand

using namespace Eigen;
using namespace Rcpp;

double MC_correction(const MatrixXd& X, const VectorXd& lambda,
        int s, int number_of_drawings = 5000)
{
    MatrixXd Xs(X.rows(),s);
    VectorXd Xi(X.rows());
    VectorXd lambda_s(s);
    lambda_s = lambda.head(s);
    VectorXd tmp_vec(VectorXd::Zero(s));
    MatrixXd LinSys(s,s);
    VectorXd RHS(s);
    double v(0.0);
    double correction(0.0);

    std::srand ( unsigned ( std::time(0) ) );

    std::vector<int> indices;
    //std::iota(indices.begin(), indices.end(), 0);
    for (int k=0; k<X.cols(); k++)
    {
        indices.push_back(k);
    }

    for (int i=0; i<number_of_drawings; i++)
    {
        //Randomly select s indices (columns of matrix X)
        std::random_shuffle(indices.begin(), indices.end());
        //Fill Xs with the selected columns
        for (int j=0; j<s; j++)
        {
            Xs.col(j) = X.col(indices[j]);
        }
        //Xi is a column of X not in Xs
        Xi = X.col(indices[s]);
        //Update the sum. Need to compute
        //v^2 = (Xi.transpose() * Xs * (Xs.transpose() * Xs).inverse() * lambda_s)^2.
        //Then Monte Carlo averages over the v's to compute the correction term. 
        LinSys = Xs.transpose() * Xs;
        RHS = Xs.transpose() * Xi;
        tmp_vec = LinSys.ldlt().solve(RHS);
        v = tmp_vec.dot(lambda_s);
        correction += pow(v, 2.0) / ((double)(number_of_drawings));
    }

    return correction;
}

//This function identifies strictly increasing sequences in the given vector
//and replaces them with their average. 
//Algorithm based on Algorithm 4 in Bogdan, van den Berg, Sabatti, Su, 
//Candes: "SLOPE -- Adaptive variable selection via convex optimization".
VectorXd MakeNonIncreasing(const VectorXd& cy)
{
    int n(cy.size());
    VectorXd y(cy);
    VectorXd x(VectorXd::Constant(n,-999));
    VectorXi block_start(VectorXi::Constant(n,-999));
    VectorXi block_end(VectorXi::Constant(n,-999));
    VectorXd block_sum(VectorXd::Constant(n,-999));
    VectorXd block_mean(VectorXd::Constant(n,-999));
    int block_num(0);
    int block_length(1);
    block_num = 0;

    //Make the tail of y constant
    VectorXd::Index minIndex;
    double minValue = y.minCoeff(&minIndex);
    y.tail(n-minIndex) = VectorXd::Constant(n-minIndex, minValue);
    //Average strictly increasing subsequences in y
    for(int k=0; k<n; k++)
    {
        block_start(block_num) = k;
        block_end(block_num) = k;
        block_sum(block_num) = y(k);
        block_mean(block_num) = block_sum(block_num);
        while((block_num > 0) && (block_mean(block_num-1)<=block_mean(block_num)))
        {
            block_end(block_num-1) = block_end(block_num);
            block_sum(block_num-1) += block_sum(block_num);
            block_mean(block_num-1) = block_sum(block_num-1) / (block_end(block_num-1) - block_start(block_num-1) + 1);
            block_num--;
        }
        block_num++;
    }
    //Save the result in x
    for(int j=0; j<block_num; j++)
    {
        block_length = block_end(j) - block_start(j) + 1;
        x.segment(block_start(j), block_length) = VectorXd::Constant(block_length, block_mean(j));
    }

    return(x);
}

// Via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do.
// [[Rcpp::depends(RcppEigen)]]
//
//' Monte Carlo based adjustment of the tuning parameters
//'
//' \code{lambdaMC} adjusts the SLOPE regularizing sequence for correlations in 
//'    the data via a Monte Carlo approach which assumes normality of the error terms.
//'
//' @param lambda_BH The regualrizing sequence as used in Theorem 1.1 in Bogdan et. al. (2015)
//' @param X The model matrix
//' @param lambda_length The corrections of the entries of \code{lambda_BH} will be 
//'    computed up to the index given by \code{lambda_length} only.
//' @param number_of_drawings The number of iterations in the Monte Carlo procedure
//'
//' @references M. Bogdan, E. van den Berg, C. Sabatti, W. Su, E. Candes (2015), \emph{SLOPE - Adaptive variable selection via convex optimization}, \url{http://arxiv.org/abs/1407.3824}
//'
//' @export
// [[Rcpp::export]]
Eigen::VectorXd lambdaMC(const Eigen::Map<Eigen::VectorXd>& lambda_BH, 
        const Eigen::Map<Eigen::MatrixXd>& X, int lambda_length, 
        int number_of_drawings=5000)
{
    int p(lambda_length);
    Eigen::VectorXd lambda_MC(Eigen::VectorXd::Constant(p,-999.999));

    lambda_MC(0) = lambda_BH(0);
    for(int i=1; i<p; i++)
    {
        lambda_MC(i) = lambda_BH(i) * pow(1.0 + MC_correction(X, lambda_MC, i, 
                    number_of_drawings), 0.5);
    }

    Eigen::VectorXd lambda_MC_nonincreasing(p);
    lambda_MC_nonincreasing = MakeNonIncreasing(lambda_MC);

    return lambda_MC_nonincreasing;
}
