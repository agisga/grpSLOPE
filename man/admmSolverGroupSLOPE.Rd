% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimization.R
\name{admmSolverGroupSLOPE}
\alias{admmSolverGroupSLOPE}
\title{Alternating direction method of multipliers}
\usage{
admmSolverGroupSLOPE(
  y,
  A,
  group,
  wt,
  lambda,
  rho = NULL,
  max.iter = 10000,
  verbose = FALSE,
  absolute.tol = 1e-04,
  relative.tol = 1e-04,
  z.init = NULL,
  u.init = NULL,
  ...
)
}
\arguments{
\item{y}{the response vector}

\item{A}{the model matrix}

\item{group}{A vector describing the grouping structure. It should
contain a group id for each predictor variable.}

\item{wt}{A vector of weights (per coefficient)}

\item{lambda}{A decreasing sequence of regularization parameters \eqn{\lambda}}

\item{rho}{Penalty parameter in the augmented Lagrangian (see Boyd et al., 2011)}

\item{max.iter}{Maximal number of iterations to carry out}

\item{verbose}{A \code{logical} specifying whether to print output or not}

\item{absolute.tol}{The absolute tolerance used in the stopping criteria for the primal and dual feasibility conditions (see Boyd et al., 2011, Sec. 3.3.1)}

\item{relative.tol}{The relative tolerance used in the stopping criteria for the primal and dual feasibility conditions (see Boyd et al., 2011, Sec. 3.3.1)}

\item{z.init}{An optional initial value for the iterative algorithm}

\item{u.init}{An optional initial value for the iterative algorithm}

\item{...}{Options passed to \code{\link{prox_sorted_L1}}}
}
\value{
A list with the entries:
  \describe{
    \item{x}{Solution (n-by-1 matrix)}
    \item{status}{Convergence status: 1 if optimal, 2 if iteration limit reached}
    \item{iter}{Number of iterations of the ADMM method}
  }
}
\description{
Compute the coefficient estimates for the Group SLOPE problem.
}
\details{
\code{admmSolverGroupSLOPE} computes the coefficient estimates
for the Group SLOPE model. The employed optimization algorithm is
the alternating direction method of multipliers (ADMM).
}
\examples{
set.seed(1)
A   <- matrix(runif(100, 0, 1), 10, 10)
grp <- c(0, 0, 1, 1, 2, 2, 2, 2, 2, 3)
wt  <- c(2, 2, 2, 2, 5, 5, 5, 5, 5, 1)
x   <- c(0, 0, 5, 1, 0, 0, 0, 1, 0, 3)
y   <- A \%*\% x
lam <- 0.1 * (10:7)
result <- admmSolverGroupSLOPE(y = y, A = A, group = grp, wt = wt,
                               lambda=lam, rho = 1, verbose = FALSE)
result$x
#           [,1]
#  [1,] 0.000000
#  [2,] 0.000000
#  [3,] 3.856002
#  [4,] 2.080742
#  [5,] 0.000000
#  [6,] 0.000000
#  [7,] 0.000000
#  [8,] 0.000000
#  [9,] 0.000000
# [10,] 3.512829

}
\references{
S. Boyd, N. Parikh, E. Chu, B. Peleato, and J. Eckstein (2011) \emph{Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers.} Foundations and Trends in Machine Learning 3 (1).
}
