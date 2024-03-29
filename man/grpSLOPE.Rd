% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grpslope.R
\name{grpSLOPE}
\alias{grpSLOPE}
\title{Group SLOPE (Group Sorted L-One Penalized Estimation)}
\usage{
grpSLOPE(
  X,
  y,
  group,
  fdr,
  lambda = "corrected",
  sigma = NULL,
  verbose = FALSE,
  orthogonalize = NULL,
  normalize = TRUE,
  max.iter = 10000,
  dual.gap.tol = 1e-06,
  infeas.tol = 1e-06,
  x.init = NULL,
  ...
)
}
\arguments{
\item{X}{The model matrix}

\item{y}{The response variable}

\item{group}{A vector describing the grouping structure. It should 
contain a group id for each predictor variable.}

\item{fdr}{Target group false discovery rate (gFDR)}

\item{lambda}{Method used to obtain the regularizing sequence lambda. Possible
values are "max", "mean", and "corrected" (default).
See \code{\link{lambdaGroupSLOPE}} for detail. Alternatively, any
non-increasing sequence of the correct length can be passed.}

\item{sigma}{Noise level. If ommited, estimated from the data, using Procedure 2 in Brzyski et. al. (2016).}

\item{verbose}{A \code{logical} specifying whether to print output or not}

\item{orthogonalize}{Whether to orthogonalize the model matrix within each group.
Do not set manually unless you are certain that your data is appropriately pre-processed.}

\item{normalize}{Whether to center the input data and re-scale the columns
of the design matrix to have unit norms. Do not disable this unless you
are certain that your data are appropriately pre-processed.}

\item{max.iter}{See \code{\link{proximalGradientSolverGroupSLOPE}}.}

\item{dual.gap.tol}{See \code{\link{proximalGradientSolverGroupSLOPE}}.}

\item{infeas.tol}{See \code{\link{proximalGradientSolverGroupSLOPE}}.}

\item{x.init}{See \code{\link{proximalGradientSolverGroupSLOPE}}.}

\item{...}{Options passed to \code{\link{prox_sorted_L1}}}
}
\value{
A list with members:
  \describe{
    \item{beta}{Solution vector. See Details.}
    \item{c}{Solution vector of the transformed model. See Details.}
    \item{group.norms}{Overall effect of each group. See Details.}
    \item{selected}{Names of selected groups (i.e., groups of predictors with at least one non-zero coefficient estimate)}
    \item{optimal}{Convergence status}
    \item{iter}{Iterations of the proximal gradient method}
    \item{lambda}{Regularizing sequence}
    \item{lambda.method}{Method used to construct the regularizing sequence}
    \item{sigma}{(Estimated) noise level}
    \item{group}{The provided grouping structure (corresponding to \code{beta})}
    \item{group.c}{Grouping structure of the transformed model (corresponding to \code{c})}
    \item{original.scale}{A list containing the estimated intercept and regression coefficients on the original scale. See Details.}
  }
}
\description{
Performs selection of significant groups of predictors and estimation of the
corresponding coefficients using the Group SLOPE method (see Brzyski et. al., 2016).
}
\details{
Multiple methods are available to generate the regularizing sequence \code{lambda},
see \code{\link{lambdaGroupSLOPE}} for detail.
The model matrix is transformed by orthogonalization within each group (see Section 2.1
in Brzyski et. al., 2016), and penalization is imposed on \eqn{\| X_{I_i} \beta_{I_i} \|}.
When \code{orthogonalize = TRUE}, due to within group orthogonalization,
the solution vector \code{beta} cannot be computed, if a group submatrix does not have full
column rank (e.g., if there are more predictors in a selected group than there are observations).
In that case only the solution vector \code{c} of the transformed (orthogonalized) model is returned.
Additionally, in any case the vector \code{group.norms} is returned with its \eqn{i}th entry
being \eqn{\| X_{I_i} \beta_{I_i} \|}, i.e., the overall effect of each group.
Note that all of these results are returned on the scale of the normalized versions of \code{X} and \code{y}.
However, \code{original.scale} contains the regression coefficients transformed to correspond to 
the original (unaltered) \code{X} and \code{y}. In that case, an estimate for the intercept term is also
returned with the other coefficients in \code{original.scale} (while on the normalized scale the estimate
of the intercept is always equal to zero, and is not explicitly provided in the \code{grpSLOPE} output).
}
\examples{
# generate some data
set.seed(1)
A   <- matrix(rnorm(100^2), 100, 100)
grp <- rep(rep(1:20), each=5)
b   <- c(runif(20), rep(0, 80))
# (i.e., groups 1, 2, 3, 4, are truly significant)
y   <- A \%*\% b + rnorm(10) 
fdr <- 0.1 # target false discovery rate
# fit a Group SLOPE model
result <- grpSLOPE(X=A, y=y, group=grp, fdr=fdr)
result$selected
# [1] "1"  "2"  "3"  "4"  "14"
result$sigma
# [1] 0.7968632
head(result$group.norms)
#         1         2         3         4         5         6 
#  2.905449  5.516103  8.964201 10.253792  0.000000  0.000000 

}
\references{
D. Brzyski, A. Gossmann, W. Su, and M. Bogdan (2016) \emph{Group SLOPE -- adaptive selection of groups of predictors}, \url{https://arxiv.org/abs/1610.04960}

D. Brzyski, A. Gossmann, W. Su, and M. Bogdan (2019) \emph{Group SLOPE -- adaptive selection of groups of predictors}. Journal of the American Statistical Association 114 (525): 419–33.
}
