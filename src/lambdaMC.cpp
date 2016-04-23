///////////////////////////////////////////////////////////////////////////////
//
//    grpSLOPE: Group SLOPE (Group Sorted L1 Penalized Estimation)
//    Copyright (C) 2016 Alexej Gossmann
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>    // std::random_shuffle, sort
#include <utility>      // std::pair
#include <vector>       // std::vector

using namespace Eigen;
using namespace Rcpp;

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
// (source: http://gallery.rcpp.org/articles/stl-random-shuffle )
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

typedef std::pair<double,int> value_index_pair;

bool comparator(const value_index_pair& l, const value_index_pair& r)
{
    return l.first > r.first; 
}

void order(VectorXd& x, VectorXi& ind_order)
{
    int n(x.size());
    std::vector<value_index_pair> xvec(n);
    value_index_pair p;
    for(int i=0; i<n; i++)
    {
        p = std::make_pair(x(i),i);
        xvec[i] = p;
    }
    sort(xvec.begin(), xvec.end(), comparator);
    for(int j=0; j<n; j++)
    {
        x(j) = xvec[j].first;
        ind_order(j) = xvec[j].second;
    }
}

double correctViaMCGaussian(const MatrixXd& X, const VectorXd& lambda,
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

    std::vector<int> indices;
    //std::iota(indices.begin(), indices.end(), 0);
    for (int k=0; k<X.cols(); k++)
    {
        indices.push_back(k);
    }

    for (int i=0; i<number_of_drawings; i++)
    {
        //Randomly select s indices (columns of matrix X)
        std::random_shuffle(indices.begin(), indices.end(), randWrapper);
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
VectorXd makeNonIncreasing(const VectorXd& cy)
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

//' Monte Carlo based adjustment of the SLOPE tuning parameters
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
//' @references M. Bogdan, E. van den Berg, C. Sabatti, W. Su, E. Candes (2015), \emph{SLOPE -- Adaptive variable selection via convex optimization}, \url{http://arxiv.org/abs/1407.3824}
//'
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
        lambda_MC(i) = lambda_BH(i) * pow(1.0 + correctViaMCGaussian(X, lambda_MC, i, 
                    number_of_drawings), 0.5);
    }

    Eigen::VectorXd lambda_MC_nonincreasing(p);
    lambda_MC_nonincreasing = makeNonIncreasing(lambda_MC);

    return lambda_MC_nonincreasing;
}


//' Monte Carlo based Group SLOPE tuning parameter correction
//'
//' \code{lambdaChiMCAdjustment} approximates the variance of (2.25) in Brzyski et. al. (2015)
//'    via Monte Carlo, in order to adjust the lambda sequence for correlations in the  data.
//'
//' The adjustment is computed for the (s+1)st coefficient of lambda, assuming 
//' that the first s coefficients are known. It is required that rank(X) is greater than the 
//' sum of the number of elements in any s groups. 
//'
//' @param y The response vector 
//' @param X The model matrix
//' @param group_id A list obtained from \code{\link{getGroupID}}
//' @param lambda A vector containing the first s entries of lambda
//' @param w A vector of weights per group
//' @param number_of_drawings The number of iterations in the Monte Carlo procedure
//'
//' @references D. Brzyski, W. Su, M. Bogdan (2015), \emph{Group SLOPE -- adaptive selection of groups of predictors}, \url{http://arxiv.org/abs/1511.09078}
//' @references \url{http://www.alexejgossmann.com/grpSLOPE/Lambda/}
//'
// [[Rcpp::export]]
double lambdaChiMCAdjustment(const Eigen::Map<Eigen::VectorXd>& y,
        const Eigen::Map<Eigen::MatrixXd>& X, const Rcpp::List group_id,
        const Eigen::Map<Eigen::VectorXd>& lambda, const Eigen::Map<Eigen::VectorXd>& w,
        int number_of_drawings=5000)
{
    // number of groups
    int p(group_id.size());
    // Array of vectors of indices of group membership
    NumericVector *groups = new NumericVector[p];
    for (int i=0; i < p; i++)
    {
        groups[i] = group_id[i];
        // adjust for 0-based indexing
        for (int j=0; j < groups[i].length(); j++) { groups[i][j] = groups[i][j] - 1; }
    }

    // number of known enries of lambda
    int s(lambda.size());

    std::vector<int> indices;
    for (int k=0; k<p; k++)
    {
        indices.push_back(k);
    }


    // Monte Carlo loop
    int total_summands = 0;
    double MC_sum = 0.0;
    double wi;
    for (int i=0; i<number_of_drawings; i++)
    {
        // Randomly select s indices (s groups of columns of matrix X)
        std::random_shuffle(indices.begin(), indices.end(), randWrapper);

        // Create matrix Xs filled with the selected columns
        int ncol_Xs = 0;
        for (int j=0; j<s; j++)
        {
            ncol_Xs += groups[indices[j]].length();
        }

        Eigen::MatrixXd Xs(X.rows(),ncol_Xs);

        int colcount = 0;
        for (int j=0; j<s; j++)
        {
            for (int k=0; k < groups[indices[j]].length(); k++)
            {
                Xs.col(colcount) = X.col( groups[indices[j]][k] );
                colcount++;
            }
        }

        //Xi is a group of columns of X not in Xs
        int ncol_Xi = groups[indices[s]].length();
        Eigen::MatrixXd Xi(X.rows(), ncol_Xi);
        for (int j=0; j < ncol_Xi; j++)
        {
           Xi.col(j) = X.col( groups[indices[s]][j] ); 
        }

        wi = w(indices[s]);
       
        // Compute least squares solution on Xs
        Eigen::VectorXd beta(ncol_Xs);
        beta = (Xs.transpose() * Xs).ldlt().solve(Xs.transpose() * y);

        Eigen::VectorXd beta_norms(s);
        int start_ind = 0;
        int group_length = 0;
        for (int j=0; j<s; j++)
        {
            group_length = groups[indices[j]].length();
            beta_norms(j) = beta.segment(start_ind, group_length).norm();
            start_ind += group_length;
        }

        // Sort the beta_norms in decreasing order, and save permutation
        VectorXi id_order(s);
        order(beta_norms, id_order);
       
        // Compute H
        Eigen::VectorXd H(ncol_Xs);
        start_ind = 0;
        int beta_start_ind = 0;
        int rank_j;
        group_length = 0;
        for (int j=0; j<s; j++)
        {
            rank_j = id_order(j);
            group_length = groups[indices[rank_j]].length();
            beta_start_ind = 0;
            for (int k=0; k<rank_j; k++) { beta_start_ind += groups[indices[k]].length(); }
            H.segment(start_ind, group_length) = lambda(j) / beta_norms(j) * beta.segment(beta_start_ind, group_length);
            start_ind += group_length;
        }

        // Update the sum.
        // Need to compute:
        // v1 = Xi.transpose() * Xs * (Xs.transpose() * Xs).inverse() * wi * H
        // and
        // v2 = Xi.transpose() * ( Identity - Xs * (Xs.transpose() * Xs).inverse() Xs.transpose() ) * z,
        // where z is a standard normal vector.
        // Then the sum is updated with (v1 + v2).squaredNorm()
        Eigen::MatrixXd LinSys = Xs.transpose() * Xs;
        Eigen::MatrixXd RHS = Xs.transpose() * Xi;
        Eigen::MatrixXd tmp_mat = LinSys.ldlt().solve(RHS);
        Eigen::VectorXd v1 = wi * tmp_mat.transpose() * H;

        int nXs = Xs.rows();
        Eigen::MatrixXd tmp_mat2 = LinSys.ldlt().solve(Xs.transpose());
        Eigen::MatrixXd tmp_mat3 = Eigen::MatrixXd::Identity(nXs, nXs) - Xs * tmp_mat2;
        Eigen::VectorXd zvec(nXs);
        for (int j=0; j<nXs; j++) { zvec(j) = R::rnorm(0.0,1.0); }
        Eigen::VectorXd v2 = Xi.transpose() * tmp_mat3 * zvec;
        
        MC_sum += (v1 + v2).squaredNorm();
        total_summands += v1.size();
    }

    delete[] groups;

    double MC_correction = sqrt(MC_sum / ((double)(total_summands)));
    return MC_correction;
}
