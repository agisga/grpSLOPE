#' @useDynLib grpSLOPE
#' @importFrom Rcpp sourceCpp
NULL
#> NULL

#' Prox for group SLOPE
#'
#' Evaluate the proximal mapping for the group SLOPE problem.
#'
#' \code{proxGroupSortedL1} evaluates the proximal mapping of the group SLOPE problem
#' by reducing it to the prox for the (regular) SLOPE and then applying the function
#' \code{\link[SLOPE]{prox_sorted_L1}}.
#'
#' @param y The response vector
#' @param group A vector or an object of class \code{groupID} (e.g. as produced by 
#'   \code{\link{getGroups}}), which is describeing the grouping structure. If it is
#'   a vector, then it should contain a group id for each predictor variable.
#' @param lambda A decreasing sequence of regularization parameters \eqn{\lambda}.
#' @param method Specifies which implementation of the Sorted L1 norm prox should be used. 
#'   See \code{\link[SLOPE]{prox_sorted_L1}}.
#'
#' @examples
#' grp <- c(0,0,0,1,1,0,2,1,0,2)
#' proxGroupSortedL1(y = 1:10, group = grp, lambda = 10:1)
#'
#' @export
proxGroupSortedL1 <- function(y, group, lambda, method = "c") {
  # TODO:  this function should be probably be rewritten from scratch in C++, in order
  # to remove the dependency on the package SLOPE
  # TODO: write my own prox function that can be selected as one of the options for "method"

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
  prox.norm <- SLOPE::prox_sorted_L1(group.norm, lambda, method)

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
#' @param tolerance The tolerance used in the stopping criteria
#' @param x.init An optional initial value for the iterative algorithm
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
#' @export
proximalGradientSolverGroupSLOPE <- function(y, A, group, wt, lambda, max.iter=1e4,
                                             verbose=TRUE, tolerance=1e-6, x.init=vector())
{
  # TODO: rewrite this function from scratch, where the main loop should be written in C++
  # TODO: after writing my own prox function, add an argument "prox"
  # TODO: check the stopping criteria  # TODO: check the stopping criteria  # TODO: check the stopping criteria

  # This is based on the source code available from
  # http://statweb.stanford.edu/~candes/SortedL1/software.html
  # under the GNU GPL-3 licence.
  #
  # Original implementation: Copyright 2013, M. Bogdan, E. van den Berg, W. Su, and E.J. Candes
  # Modifications: Copyright 2015, Alexej Gossmann
  
  # Initialize ---------------------------------------------------------------
  
  # Ensure that lambda is non-increasing
  n.lambda <- length(lambda)
  if ((n.lambda > 1) && any(lambda[2:n.lambda] > lambda[1:n.lambda-1])) {
    stop("Lambda must be non-increasing.")
  }

  # Adjust matrix for prior weights
  Dinv <- diag(wt)
  A    <- A %*% Dinv
  p    <- ncol(A)

  # Prepare grouping information
  group.id <- getGroupID(group)

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
    # TODO: maybe have to change the "prox" argument later
    if (length(x.init) == 0) x.init <- matrix(0,p,1)
    sol <- SLOPE::SLOPE_solver(A=A, b=y, lambda=lambda, initial=x.init,
                               prox=SLOPE::prox_sorted_L1, max_iter=max.iter,
                               tol_infeas=tolerance, tol_rel_gap=tolerance)
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
  f.prev   <- Inf
  iter     <- 0
  L.iter   <- 0
  status   <- STATUS_RUNNING
  relgap   <- Inf
  
  if (verbose == TRUE) {
    printf <- function(...) invisible(cat(sprintf(...)))
    printf('%5s  %9s  %9s\n','Iter','||r||_2','Rel. gap')
  }
  
  # Main loop ---------------------------------------------------------
  while (TRUE)
  {    
    # Compute the gradient
    r <- (A %*% b) - y
    g <- t(A) %*% r
    
    # Increment iteration count
    iter <- iter + 1
    
    # Stopping criteria
    if ((status == 0) && (iter >= max.iter)) {
      status <- STATUS_ITERATIONS
    }

    if (status != 0) {
      if (verbose == TRUE) {
        printf('Exiting with status %d -- %s\n', status, STATUS_MSG[[status]])
      }
      break
    }
    
    # Keep copies of previous values
    f.prev  <- as.double(crossprod(r)) / 2
    Ax.prev <- Ax
    x.prev  <- x
    b.prev  <- b
    tt.prev <- tt
    
    # Lipschitz search
    while (TRUE) {
      x  <- proxGroupSortedL1(y = b - (1/L)*g, group = group.id, lambda = lambda/L)
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

    # Compute 'relative gap'
    if (sqrt(sum(b.prev^2))==0) {
      relgap <- sqrt(sum((b-b.prev)^2))
    } else {
      relgap <- sqrt(sum((b-b.prev)^2)) / sqrt(sum(b.prev^2))
    }

    # Format string
    if (verbose == TRUE) str <- sprintf('   %9.2e', relgap)
     
    # Check relative gap
    if (relgap < tolerance) status <- STATUS_OPTIMAL
    
    if (verbose == TRUE) {
      printf('%5d  %9.2e%s\n', iter,f,str)
    }
  } # while (TRUE)
  
  
  result <- recordResult(b=b, status=status, L=L, iter=iter, L.iter=L.iter)
  return(result)
}
