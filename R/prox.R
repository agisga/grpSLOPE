# SLOPE: Sorted L1 Penalized Estimation (SLOPE)
# Copyright (C) 2015 Malgorzata Bogdan, Ewout van den Berg, Chiara Sabatti,
# Weijie Su, Emmanuel Candes, Evan Patterson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Prox for sorted L1 norm
#'
#' Compute the prox for the sorted L1 norm. That is, given a vector \eqn{x}
#' and a decreasing vector \eqn{\lambda}, compute the unique value of \eqn{y}
#' minimizing
#' \deqn{\frac{1}{2} \Vert x - y \Vert_2^2 +
#'       \sum_{i=1}^n \lambda_i |x|_{(i)}.}
#'
#' @param x input vector
#' @param lambda vector of \eqn{\lambda}'s, sorted in decreasing order
#'
#' @details At present, two methods for computing the sorted L1 prox are
#' supported. By default, we use a fast custom C implementation. Since SLOPE
#' can be viewed as an isotonic regression problem, the prox can also be
#' computed using the \code{isotone} package. This option is provided
#' primarily for testing.
#'
#' @author Malgorzata Bogdan, Ewout van den Berg, Chiara Sabatti, Weijie Su,
#'  Emmanuel Candes, Evan Patterson
#'
prox_sorted_L1 <- function(x, lambda) {
  # Normalize input
  if (is.complex(x)) {
    sign = complex(argument = Arg(x))
    x = Mod(x)
  } else {
    sign = sign(x)
    x = abs(x)
  }

  # Sort input
  s = sort(x, decreasing=TRUE, index.return=TRUE)

  # Compute prox
  result = prox_sorted_L1_C(s$x, lambda)

  # Restore original order and sign
  result[s$ix] <- result
  result * sign
}
