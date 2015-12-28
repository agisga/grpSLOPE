#' Get a groupID object
#'
#' Intended for internal use only.
#'
#' @param group A vector describeing the grouping structure. It should 
#'    contain a group id for each predictor variable.
#'
#' @return An object of class groupID, which is a list, whose members are 
#'    vectors of indices corresponding to each group. The names of
#'    the list members are the corresponding group names.
#'
getGroupID <- function(group) {
  group.unique <- unique(group)
  n.group <- length(group.unique)
  group.id <- list()
  for (i in 1:n.group){
    id <- as.character(group.unique[i])
    group.id[[id]] <- which(group==group.unique[i])
  }
  class(group.id) <- "groupID"
  return(group.id)
}
