#' Prox for group SLOPE
#'
#' Evaluate the proximal mapping for the group SLOPE problem.
#'
#' \code{proxGroupSortedL1} evaluates the proximal mapping of the group SLOPE problem
#' by reducing it to the prox for the (regular) SLOPE and then applying the function
#' \code{\link[SLOPE]{prox_sorted_L1}}.
#'
#' @param y The response vector
#' @param group A vector describing the grouping structure. It should contain a group id
#'   for each predictor variable.
#' @param lambda A decreasing sequence of regularization parameters \eqn{\lambda}.
#' @param method Specifies which implementation of the Sorted L1 norm prox should be used. 
#'   See \code{\link[SLOPE]{prox_sorted_L1}}.
#'
#' @examples
#' grp <- c(0,0,0,1,1,0,2,1,0,2)
#' proxGroupSortedL1(y = 1:10, group = grp, lambda = 10:1)
#'
proxGroupSortedL1 <- function(y, group, lambda, method = "c") {
  group.unique <- unique(group)
  n.group <- length(group.unique)
  group.id <- list()
  for (i in 1:n.group){
    group.id[[i]] <- which(group==group.unique[i])
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
