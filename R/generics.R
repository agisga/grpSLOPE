###############################################################################
#
#    grpSLOPE: Group SLOPE (Group Sorted L1 Penalized Estimation)
#    Copyright (C) 2016 Alexej Gossmann
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#' Extract model coefficients
#' 
#' Extract the regression coefficients from a \code{grpSLOPE} object, either on the
#' scale of the normalized design matrix (i.e., columns centered and scaled to unit norm),
#' or on the original scale.
#'
#' If \code{scaled} is set to \code{TRUE}, then the coefficients are returned for the 
#' normalized version of the design matrix, which is the scale on which they were computed. 
#' If \code{scaled} is set to \code{FALSE}, then the coefficients are transformed to
#' correspond to the original (unaltered) design matrix.
#' In case that \code{scaled = FALSE}, an estimate for the intercept term is returned with
#' the other coefficients. In case that \code{scaled = TRUE}, the estimate of the intercept 
#' is always equal to zero, and is not explicitly provided.
#'
#' @param object A \code{grpSLOPE} object 
#' @param scaled Should the coefficients be returned for the normalized version of the design matrix?
#' @param ... Potentially further arguments passed to and from methods
#'
#' @examples
#' set.seed(1)
#' A   <- matrix(rnorm(100^2), 100, 100)
#' grp <- rep(rep(1:20), each=5)
#' b   <- c(rep(1, 20), rep(0, 80))
#' y   <- A %*% b + rnorm(10) 
#' result <- grpSLOPE(X=A, y=y, group=grp, fdr=0.1)
#' head(coef(result))
#' #        X1        X2        X3        X4        X5        X6 
#' #  7.942177  7.979269  8.667013  8.514861 10.026664  8.963364 
#' head(coef(result, scaled = FALSE))
#' # (Intercept)          X1          X2          X3          X4          X5 
#' #  -0.4418113   0.8886878   0.8372108   0.8422089   0.8629597   0.8615827 
#' 
#' @export
coef.grpSLOPE <- function(object, scaled = TRUE, ...) {
  if(is.null(object$beta)) {
    stop("beta is set to NULL. Maybe one of the group submatrices did not have full column rank? See documentation to grpSLOPE().")
  }

  if(scaled) {
    coefs <- object$beta
    names(coefs) <- paste0("X", 1:length(coefs))
  } else {
    coefs <- c(object$original.scale$intercept, object$original.scale$beta)
    names(coefs) <- c("(Intercept)", paste0("X", 1:length(object$original.scale$beta)))
  }

  return(coefs)
}

#' Extract (estimated) noise level
#' 
#' Extract the noise level of the \code{grpSLOPE} model.
#'
#' @param object A \code{grpSLOPE} object 
#' @param ... Potentially further arguments passed to and from methods
#'
#' @examples
#' set.seed(1)
#' A   <- matrix(rnorm(100^2), 100, 100)
#' grp <- rep(rep(1:20), each = 5)
#' b   <- c(rep(1, 20), rep(0, 80))
#' y   <- A %*% b + rnorm(10) 
#' # model with unknown noise level
#' result <- grpSLOPE(X = A, y = y, group = grp, fdr = 0.1)
#' sigma(result)
#' # [1] 0.6505558
#' # model with known noise level
#' result <- grpSLOPE(X = A, y = y, group = grp, fdr = 0.1, sigma = 1)
#' sigma(result)
#' # [1] 1
#' 
#' @export
sigma.grpSLOPE <- function(object, ...) {
  return(object$sigma)
}

#predict.grpSLOPE <- function(object, newdata, debiased = TRUE) {
#}
