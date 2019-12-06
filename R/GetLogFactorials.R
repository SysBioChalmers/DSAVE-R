#' GetLogFactorials
#'
#' GetLogFactorials
#'
#'   Calculates the logarithm of x! for all x <= n and places those in a vector.
#'   To get around that you cannot index a vector with 0, you should always
#'   index this vector with the value + 1.
#'   The purpose of this function is that x! cannot be calculated for
#'   x > 170, while that is not a problem for log(x!), but you cannot use
#'   the built-in factorial function in matlab and do log afterwards!
#'
#' @param n numeric, the max value needed
#' @author Johan Gustafsson, Juan Inda, <inda@@chalmers.se>
#' @return vector
#'
GetLogFactorials <- function(n){
  stopifnot(is.numeric(n) & n>0)
  return(c(0,cumsum(log(1:n))))
}

