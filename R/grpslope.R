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

#' Group SLOPE algorithm
#'
#' Compute the coefficient estimates of the Group SLOPE model  
#'
#' @export
grpSLOPE <- function(X, Dinv, b, group, lambda, opts=list())
{
  # TODO: rewrite this function from scratch, where the main loop should be written in C++
  # TODO: check whether all groups have length 1, then use SLOPE
  # TODO: compute groupID somewhere in here

  # Copyright 2013, M. Bogdan, E. van den Berg, W. Su, and E.J. Candes
  # Copyright 2015, Alexej Gossmann

  # -------------------------------------------------------------
  # Start times
  # -------------------------------------------------------------
  t0 <- proc.time()[3]
  
  # -------------------------------------------------------------
  # Define function for retrieving option fields with defaults
  # -------------------------------------------------------------
  getDefaultField <- function(opts,name,default)
  {  
    if (!is.null(opts[[name]])) {
      return(opts[[name]])
    } else {
      return(default)
    }
  }
  
  # -------------------------------------------------------------
  # Parse parameters
  # -------------------------------------------------------------
  iterations <- getDefaultField(opts,"iterations", 10000)
  verbosity  <- getDefaultField(opts,"verbosity" , 1)
  gradIter   <- getDefaultField(opts,"gradIter"  , 20)
  tolerance  <- getDefaultField(opts,"tolerance" , 1e-6)
  xInit      <- getDefaultField(opts,"xInit"     , vector())
  
  # -------------------------------------------------------------
  # Ensure that lambda is non-increasing
  # -------------------------------------------------------------
  n.lambda <- length(lambda)
  if ((n.lambda > 1) && any(lambda[2:n.lambda] > lambda[1:n.lambda-1])) {
    stop("Lambda must be non-increasing.")
  }
  
  # -------------------------------------------------------------
  # Initialize
  # -------------------------------------------------------------

  A <- X %*% Dinv
  n <- ncol(A)
  
  # Get initial lower bound on the Lipschitz constant
  x <- matrix(rnorm(n),c(n,1))
  x <- x / norm(x, "f")
  x <- t(A) %*% (A %*% x)
  L <- norm(x, "f")
  
  # Set constants
  STATUS_RUNNING    <- 0
  STATUS_OPTIMAL    <- 1
  STATUS_ITERATIONS <- 2
  STATUS_MSG        <- c('Optimal','Iteration limit reached')
  
  # Initialize parameters and iterates
  if (length(xInit) == 0) xInit <- matrix(0,n,1)
  tt      <- 1
  eta     <- 2
  lambda  <- as.matrix(lambda)
  b       <- as.matrix(b)
  x       <- xInit
  y       <- x
  Ax      <- A %*% x
  fPrev   <- Inf
  iter    <- 0
  status  <- STATUS_RUNNING
  Aprods  <- 2
  ATprods <- 1
  relgap  <- Inf
  
  if (verbosity > 0) {
    printf <- function(...) invisible(cat(sprintf(...)))
    printf('%5s  %9s  %9s\n','Iter','||r||_2','Rel. gap')
  }
  
  # -------------------------------------------------------------
  # Main loop
  # -------------------------------------------------------------
  while (TRUE)
  {    
    # Compute the gradient at f(y) (Includes first iterations)
    if ((iter %% gradIter) == 0) {
      r <- (A %*% y) - b
      g <- t(A) %*% r
      f <- as.double(crossprod(r)) / 2
    } else {
      r <- (Ax + ((ttPrev - 1) / tt) * (Ax - AxPrev)) - b
      g <- t(A) %*% r
      f <- as.double(crossprod(r)) / 2
    }
    
    # Increment iteration count
    iter <- iter + 1
    
    # Stopping criteria
    if ((status == 0) && (iter >= iterations)) {
      status <- STATUS_ITERATIONS
    }

    if (status != 0) {
      if (verbosity > 0) {
        printf('Exiting with status %d -- %s\n', status, STATUS_MSG[[status]])
      }
      break
    }
    
    # Keep copies of previous values
    AxPrev <- Ax
    xPrev  <- x
    yPrev  <- y
    fPrev  <- f
    ttPrev <- tt
    
    # Lipschitz search
    while (TRUE) {
      # Compute prox mapping
      x <- proxGroupSortedL1(y = y - (1/L)*g, group = group, lambda = lambda/L)
      d <- x - y
    
      Ax <- A %*% x
      r  <- Ax-b
      f  <- as.double(crossprod(r))/2
      q  <- fPrev + as.double(crossprod(d,g)) + (L/2)*as.double(crossprod(d))
                              
      Aprods <- Aprods + 1
                              
      if (q >= f*(1-1e-12))
        break
      else
        L <- L * eta
    }
    
    # Update
    tt <- (1 + sqrt(1 + 4*tt^2)) / 2
    y  <- x + ((ttPrev - 1) / tt) * (x - xPrev)

    # Compute 'relative gap'
    if (sqrt(sum(yPrev^2))==0) {
      relgap <- sqrt(sum((y-yPrev)^2))
    } else {
      relgap <- sqrt(sum((y-yPrev)^2)) / sqrt(sum(yPrev^2))
    }

    # Format string
    if (verbosity > 0) str <- sprintf('   %9.2e', relgap)
     
    # Check relative gap
    if (relgap < tolerance) status <- STATUS_OPTIMAL
    
    if ( verbosity > 0 ) {
      printf('%5d  %9.2e%s\n', iter,f,str)
    }
  } # While (TRUE)
  
  
  # Information structure
  solution           <- list()
  solution$runtime   <- proc.time()[3] - t0
  solution$Aprods    <- Aprods + ceiling(iter / gradIter)
  solution$ATprods   <- ATprods + iter
  solution$status    <- status
  solution$x         <- Dinv %*% y
  solution$L         <- L
  solution$iter      <- iter

  # Define an S3 class
  class(solution) <- "grpSLOPE"
  
  return(solution)
}
