#' @useDynLib grpSLOPE
#' @importFrom Rcpp sourceCpp
NULL
#> NULL

#' Fast prox for the Sorted L1 norm
#' 
#' A fast stack-based algorithm for the prox for the Sorted L1 norm.
#'
#' See Algorithm 4 in Bogdan et. al. (2015).
#'
#' @param y A vector
#' @param lambda A vector whose entries should form a nonincreasing sequence
#'
#' @references M. Bogdan, E. van den Berg, C. Sabatti, W. Su, E. Candes (2015), \emph{SLOPE - Adaptive variable selection via convex optimization}, \url{http://arxiv.org/abs/1407.3824}
#'
#' @export
proxSortedL1 <- function (y, lambda) 
{
  if (is.complex(y)) {
    sgn = complex(modulus=1, argument=Arg(y))
    y = Mod(y)
  } else {
    sgn = sign(y)
    y = abs(y)
  }
  y.sorted = sort(y, decreasing=TRUE, index.return=TRUE)
  result <- proxSortedL1Rcpp(as.double(y.sorted$x), as.double(lambda))
  result[y.sorted$ix] <- result
  result <- result * sgn
  return(result)
}

#' Prox for group SLOPE
#'
#' Evaluate the proximal mapping for the group SLOPE problem.
#'
#' \code{proxGroupSortedL1} evaluates the proximal mapping of the group SLOPE problem
#' by reducing it to the prox for the (regular) SLOPE and then applying the fast prox
#' algorithm for the Sorted L1 norm. The argument \code{method} specifies which 
#' implementation of the Sorted L1 norm prox should be used. 
#' Possible values are \code{"rcpp"} (default), \code{"c"}, and \code{"isotone"}.
#' The default option \code{"rcpp"} uses the internal implementation in
#' \code{\link{proxSortedL1}}. The alternative options \code{"c"} and
#' \code{"isotone"} call the function \code{\link[SLOPE]{prox_sorted_L1}} from the
#' package \code{SLOPE} (see there for detail on these two options).
#'
#' @param y The response vector
#' @param group A vector or an object of class \code{groupID} (e.g. as produced by 
#'   \code{\link{getGroups}}), which is describing the grouping structure. If it is
#'   a vector, then it should contain a group id for each predictor variable.
#' @param lambda A decreasing sequence of regularization parameters \eqn{\lambda}
#' @param method Specifies which implementation of the Sorted L1 norm prox should be used. 
#'   Possible values are \code{"rcpp"} (default), \code{"c"}, and \code{"isotone"}. See detail.
#'
#' @examples
#' grp <- c(0,0,0,1,1,0,2,1,0,2)
#' proxGroupSortedL1(y = 1:10, group = grp, lambda = 10:1)
#'
#' @references M. Bogdan, E. van den Berg, C. Sabatti, W. Su, E. Candes (2015), \emph{SLOPE - Adaptive variable selection via convex optimization}, \url{http://arxiv.org/abs/1407.3824}
#'
#' @export
proxGroupSortedL1 <- function(y, group, lambda, method = "rcpp") {
  if (inherits(group, "groupID")) {
    n.group <- length(group)
    group.id <- group
  } else {
    n.group <- length(unique(group))
    group.id <- getGroupID(group)
  }

  # compute Euclidean norms for groups in y
  group.norm <- rep(NA, n.group)
  for (i in 1:n.group){
    selected <- group.id[[i]]
    group.norm[i] <- norm(as.matrix(y[selected]), "f")
  }

  # get Euclidean norms of the solution vector
  if (method == "rcpp") {
    prox.norm <- proxSortedL1(group.norm, lambda)
  } else {
    prox.norm <- SLOPE::prox_sorted_L1(group.norm, lambda, method)
  }

  # compute the solution
  prox.solution <- rep(NA, length(y))
  for (i in 1:n.group){
    selected <- group.id[[i]]
    prox.solution[selected] <- prox.norm[i] / group.norm[i] * y[selected]
  }

  return(prox.solution)
}

#' Proximal gradient method for Group SLOPE
#'
#' Compute the coefficient estimates for the Group SLOPE problem.
#'
#' \code{proximalGradientSolverGroupSLOPE} computes the coefficient estimates
#' for the Group SLOPE model. The employed optimization algorithm is FISTA with 
#' backtracking Lipschitz search.
#'
#' @param y the response vector
#' @param A the model matrix
#' @param group A vector describing the grouping structure. It should 
#'   contain a group id for each predictor variable.
#' @param wt A vector of weights
#' @param lambda A decreasing sequence of regularization parameters \eqn{\lambda}
#' @param max.iter Maximal number of iterations to carry out
#' @param verbose A \code{logical} specifying whether to print output or not
#' @param dual.gap.tol The tolerance used in the stopping criteria for the duality gap
#' @param infeas.tol The tolerance used in the stopping criteria for the infeasibility
#' @param x.init An optional initial value for the iterative algorithm
#' @param method Specifies which implementation of the Sorted L1 norm prox should be used. 
#'   See \code{\link{proxGroupSortedL1}} for detail.
#'
#' @return A list with members:
#'   \describe{
#'     \item{x}{Solution (n-by-1 matrix)}
#'     \item{status}{Convergence status: 1 if optimal, 2 if iteration limit reached}
#'     \item{L}{Approximation of the Lipschitz constant (step size)}
#'     \item{iter}{Iterations of the proximal gradient method}
#'     \item{L.iter}{Total number of iterations spent in Lischitz search}
#'   }
#'
#' @examples
#' set.seed(1)
#' A   <- matrix(runif(100, 0, 1), 10, 10)
#' grp <- c(0, 0, 1, 1, 2, 2, 2, 2, 2, 3)
#' wt  <- c(0.5, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.2, 0.2, 1)
#' x   <- c(0, 0, 5, 1, 0, 0, 0, 1, 0, 3)
#' y   <- A %*% x
#' lam <- 0.1 * (10:7)
#' result <- proximalGradientSolverGroupSLOPE(y=y, A=A, group=grp, wt=wt, lambda=lam, verbose=FALSE)
#' result$x
#'
#' @references A. Gossmann, S. Cao, Y.-P. Wang (2015), \emph{Identification of Significant Genetic Variants via SLOPE, and Its Extension to Group SLOPE}, \url{http://dx.doi.org/10.1145/2808719.2808743}
#' @references D. Brzyski, W. Su, M. Bogdan (2015), \emph{Group SLOPE — adaptive selection of groups of predictors}, \url{http://arxiv.org/abs/1511.09078}
#'
#' @export
proximalGradientSolverGroupSLOPE <- function(y, A, group, wt, lambda, max.iter=1e4,
                                             verbose=TRUE, dual.gap.tol=1e-6, 
                                             infeas.tol=1e-6, x.init=vector(), method="rcpp")
{
  # This function is based on the source code available from
  # http://statweb.stanford.edu/~candes/SortedL1/software.html
  # under the GNU GPL-3 licence.
  #
  # Original implementation: Copyright 2013, M. Bogdan, E. van den Berg, W. Su, and E.J. Candes
  # Modifications: Copyright 2015, Alexej Gossmann
  
  # Initialize ---------------------------------------------------------------

  # Prepare grouping information
  group.id <- getGroupID(group)
  n.group  <- length(group.id)
  
  # Ensure that lambda is non-increasing and of right length
  n.lambda <- length(lambda)
  if ((n.lambda > 1) && any(lambda[2:n.lambda] > lambda[1:n.lambda-1])) {
    stop("Lambda must be non-increasing.")
  }
  if (n.lambda != n.group) {
    stop("Lambda must have exactly as many entries as there are groups.")
  }

  # Adjust matrix for prior weights
  Dinv <- diag(wt)
  A    <- A %*% Dinv
  p    <- ncol(A)

  # Auxilliary function to record output information
  recordResult <- function(b, status, L, iter, L.iter) {
    result           <- list()
    result$x         <- Dinv %*% b
    result$status    <- status
    result$L         <- L
    result$iter      <- iter
    result$L.iter    <- L.iter
    return(result)
  }

  # Run regular SLOPE if all groups are singletons
  if (length(group.id) == p) {
    if (length(x.init) == 0) x.init <- matrix(0,p,1)
    sol <- SLOPE::SLOPE_solver(A=A, b=y, lambda=lambda, initial=x.init,
                               max_iter=max.iter, tol_infeas=infeas.tol,
                               tol_rel_gap=dual.gap.tol)
    status <- 2
    if (sol$optimal) status <- 1
    result <- recordResult(b=matrix(sol$x, c(p,1)), status=status, L=sol$lipschitz,
                           iter=sol$iter, L.iter=NULL)
    return(result)
  }
  
  # Get initial lower bound on the Lipschitz constant
  rand.mat <- matrix(rnorm(p),c(p,1))
  rand.mat <- rand.mat / norm(rand.mat, "f")
  rand.mat <- t(A) %*% (A %*% rand.mat)
  L <- norm(rand.mat, "f")
  
  # Set constants
  STATUS_RUNNING    <- 0
  STATUS_OPTIMAL    <- 1
  STATUS_ITERATIONS <- 2
  STATUS_MSG        <- c('Optimal','Iteration limit reached')
  
  # Initialize parameters and iterates
  if (length(x.init) == 0) x.init <- matrix(0,p,1)
  tt       <- 1
  eta      <- 2
  lambda   <- as.matrix(lambda)
  y        <- as.matrix(y)
  x        <- x.init
  b        <- x
  Ax       <- A %*% x
  f        <- Inf
  f.prev   <- Inf
  iter     <- 0
  L.iter   <- 0
  status   <- STATUS_RUNNING
  duality.gap <- Inf
  infeasibility <- Inf
  
  if (verbose == TRUE) {
    printf <- function(...) invisible(cat(sprintf(...)))
    printf('%5s  %9s  %9s  %9s\n', 'Iter', '||r||_2', 'Dual. gap', 'Infeas.')
  }
  
  # Main loop ---------------------------------------------------------
  while (TRUE)
  {    
    # Compute the gradient
    r <- (A %*% b) - y
    g <- crossprod(A, r)
    
    # Stopping criteria --------------------------------------------

    # Compute duality gap (eq. (1.7) in Appendix I in Brzyski et. al. (2015))
    b.norms <- rep(NA, n.group)
    for (i in 1:n.group) {
      selected <- group.id[[i]]
      b.norms[i] <- norm(as.matrix(b[selected]), "f")
    }
    b.norms.sorted <- sort(b.norms, decreasing=TRUE)
    duality.gap <- crossprod(b, g) + crossprod(lambda, b.norms.sorted)

    # Compute the infeasibility
    # (derivation of this infeasibility criterion: http://www.alexejgossmann.com/grpSLOPE/Infeasibility/)
    g.norms <- rep(NA, n.group)
    for (i in 1:n.group) {
      selected <- group.id[[i]]
      g.norms[i] <- norm(as.matrix(g[selected]), "f")
    }
    g.norms.sorted  <- sort(g.norms, decreasing=TRUE)
    infeasibility <- max(max(cumsum(g.norms.sorted-lambda)),0)

    # Format string
    if (verbose == TRUE) str <- sprintf('   %9.2e  %9.2e', duality.gap, infeasibility)
     
    # Check stopping criteria
    if ((duality.gap < dual.gap.tol) && (infeasibility < infeas.tol)) {
      status <- STATUS_OPTIMAL
    }
    
    if (verbose == TRUE) {
      printf('%5d  %9.2e%s\n', iter, f, str)
    }

    if ((status == 0) && (iter >= max.iter)) {
      status <- STATUS_ITERATIONS
    }

    if (status != 0) {
      if (verbose == TRUE) {
        printf('Exiting with status %d -- %s\n', status, STATUS_MSG[[status]])
      }
      break
    }
    # (stopping criteria ends) ------------------
    
    # Increment iteration count
    iter <- iter + 1
    
    # Keep copies of previous values
    f.prev  <- as.double(crossprod(r)) / 2
    Ax.prev <- Ax
    x.prev  <- x
    b.prev  <- b
    tt.prev <- tt
    
    # Lipschitz search
    while (TRUE) {
      x  <- proxGroupSortedL1(y = b - (1/L)*g, group = group.id, lambda = lambda/L, method=method)
      d  <- x - b
      Ax <- A %*% x
      r  <- Ax-y
      f  <- as.double(crossprod(r))/2
      qq <- f.prev + as.double(crossprod(d,g)) + (L/2)*as.double(crossprod(d))
                              
      L.iter <- L.iter + 1
                              
      if (qq >= f*(1-1e-12))
        break
      else
        L <- L * eta
    }
    
    # Update
    tt <- (1 + sqrt(1 + 4*tt^2)) / 2
    b  <- x + ((tt.prev - 1) / tt) * (x - x.prev)
  } # while (TRUE)
  
  
  result <- recordResult(b=b, status=status, L=L, iter=iter, L.iter=L.iter)
  return(result)
}


#' Regularizing sequence for Group SLOPE
#'
#' Generate the regularizing sequence \code{lambda} for the Group SLOPE
#' problem according to one of multiple methods.
#'
#' Multiple methods are available to generate the regularizing sequence \code{lambda}:
#' \itemize{
#'   \item "BH" -- method of Theorem 1.1 in Bogdan et. al. (2015)
#'   \item "gaussian" -- method of Section 3.2.2 in Bogdan et. al. (2015)
#'   \item "gaussianMC" -- method introduced in Gossmann et. al. (2015)
#'   \item "chiOrthoMax" -- lambdas as in Theorem 2.5 in Brzyski et. al. (2015)
#'   \item "chiOrthoMean" -- lambdas of equation (2.14) in Brzyski et. al. (2015)
#'   \item "chiEqual" -- Procedure 1 in Brzyski et. al. (2015)
#'   \item "chiMean" -- Procedure 2 in Brzyski et. al. (2015)
#'   \item "chiMC" -- (experimental) A Monte Carlo lambda selection method based on equation (2.25)
#'            in Brzyski et. al. (2015). Requires that rank(\code{A}) is greater than
#'            the sum of the number of elements in any \code{n.MC} groups. 
#' }
#'
#' @param fdr Target false discovery rate
#' @param n.groups Number of groups
#' @param group A vector describing the grouping structure. It should 
#'    contain a group id for each predictor variable.
#' @param A The model matrix
#' @param y The response variable
#' @param wt A named vector of weights, one weight per group of predictors (named according to names as in vector \code{group})
#' @param n.obs Number of observations (i.e., number of rows in \code{A})
#' @param method Possible values are "BH", "gaussian", "gaussianMC",
#'    "chiOrthoMax", "chiOrthoMean",  "chiEqual", "chiMean", "chiMC". See under Details.
#' @param n.MC The corrections of the entries of lambda will be 
#'    computed up to the index given by \code{n.MC} only. \code{n.MC} should be
#'    less than or equal to \code{n.group}.
#' @param MC.reps The number of repetitions of the Monte Carlo procedure
#'
#' @examples
#' fdr     <- 0.1
#' n.obs   <- 700
#' n.group <- 90
#' group   <- vector()
#' for (i in 1:30) {
#'   tmp <- rep((i-1)*3+c(1,2,3), c(5,10,20))
#'   group <- c(group, tmp)
#' }
#' wt <- rep(c(sqrt(5), sqrt(10), sqrt(20)), 30)
#' names(wt) <- names(getGroupID(group))
#'
#' lambda.BH <- lambdaGroupSLOPE(fdr=fdr, n.group=n.group, method="BH")
#' lambda.G <- lambdaGroupSLOPE(fdr=fdr, n.group=n.group, n.obs=n.obs, method="gaussian")
#' lambda.max <- lambdaGroupSLOPE(fdr=fdr, group=group, wt=wt, method="chiOrthoMax") 
#' lambda.mean <- lambdaGroupSLOPE(fdr=fdr, group=group, wt=wt, method="chiOrthoMean") 
#' lambda.chi <- lambdaGroupSLOPE(fdr=fdr, n.obs=n.obs, group=group, wt=wt, method="chiMean")
#'
#' @references M. Bogdan, E. van den Berg, C. Sabatti, W. Su, E. Candes (2015), \emph{SLOPE - Adaptive variable selection via convex optimization}, Annals of Applied Statistics
#' @references A. Gossmann, S. Cao, Y.-P. Wang (2015), \emph{Identification of Significant Genetic Variants via SLOPE, and Its Extension to Group SLOPE}, \url{http://dx.doi.org/10.1145/2808719.2808743}
#' @references D. Brzyski, W. Su, M. Bogdan (2015), \emph{Group SLOPE — adaptive selection of groups of predictors}, \url{http://arxiv.org/abs/1511.09078}
#'
#' @export
lambdaGroupSLOPE <- function(fdr=0.1, n.group=NULL, group=NULL,
                             A=NULL, y=NULL, wt=NULL, n.obs=NULL, method,
                             n.MC=n.group, MC.reps=5000)
{
  # Prepare grouping information
  if (!is.null(group)) {
    group.id    <- getGroupID(group)
    n.group     <- length(group.id)
    group.sizes <- sapply(group.id, FUN=length)
  }

  if (is.null(n.group)) {
    stop("Either n.group or group needs to be passed as function argument.")
  }

  if (method %in% c("BH", "gaussian", "gaussianMC")) {

    if (method=="BH") {
      return(lambdaBH(fdr=fdr, n.group=n.group))

    } else if (method=="gaussian") {
      if (is.null(A) && is.null(n.obs)) {
        stop("Either A or n.obs needs to be passed as an argument when method is 'gaussian'.")
      }
      if (is.null(n.obs)) n.obs <- nrow(A)

      return(lambdaGaussian(fdr=fdr, n.group=n.group, n.obs=n.obs))

    } else if (method=="gaussianMC") {
      if (is.null(A) || is.null(group)) {
        stop("A and group need to be passed as arguments when method is 'gaussianMC'.")
      }

      return(lambdaGaussianMC(fdr=fdr, n.group=n.group, group.id=group.id,
                              A=A, n.MC=n.MC, MC.reps=MC.reps))

    }
  } else if (method %in% c("chiOrthoMax", "chiOrthoMean", "chiEqual", "chiMean", "chiMC")) {
    if (is.null(group) || is.null(wt)) {
      stop("Arguments group and wt need to be provided when method is one of 'chiOrthoMax', 'chiOrthoMean', 'chiEqual', 'chiMean', 'chiMC'.")
    }

    if (method=="chiOrthoMax" || method=="chiOrthoMean") {
      return(lambdaChiOrtho(fdr=fdr, n.group=n.group, group.sizes=group.sizes, wt=wt, method=method))

    } else if (method=="chiEqual") {
      if (is.null(A) && is.null(n.obs)) {
        stop("Either A or n.obs needs to be passed as an argument when method is 'chiEqual'")
      }
      if (is.null(n.obs)) n.obs <- nrow(A)

      # Equal group sizes and weights
      if ( (length(unique(group.sizes))!=1) || (length(unique(wt))!=1) ) {
        stop("Method 'chiEqual' requires equal group sizes and equal weights.")
      }
      m <- unique(group.sizes)
      w <- unique(wt)

      return(lambdaChiEqual(fdr=fdr, n.obs=n.obs, n.group=n.group, m=m, w=w))

    } else if (method=="chiMean") {
      if (is.null(A) && is.null(n.obs)) {
        stop("Either A or n.obs needs to be passed as an argument when method is 'chiMean'")
      }
      if (is.null(n.obs)) n.obs <- nrow(A)

      return(lambdaChiMean(fdr=fdr, n.obs=n.obs, n.group=n.group, group.sizes=group.sizes, wt=wt))

    } else if (method=="chiMC") {
      if (is.null(y) || is.null(A) || is.null(group) || is.null(wt)) {
        stop("A, y, group and wt need to be passed as arguments when method is 'chiMC'.")
      }

      return(lambdaChiMC(fdr=fdr, X=A, y=y, group.id=group.id, wt=wt, n.MC=n.MC, MC.reps=MC.reps))

    }
  } else {
    stop(paste(method, "is not a valid method."))
  }
}
