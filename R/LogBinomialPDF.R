#' LogBinomialPDF
#'
#' LogBinomialPDF
#'
#' Calculates the probability density function for a vector of observations
#'   for a binomial distribution.
#'   The idea with this function is that it is not easy to calculate the PDF of
#'   a binomial if n is large, since we cannot practically calculate n! of a
#'   large number. The logarithm of the pdf is fine to calculate however, since
#'   log(n!) = log(n) + log(n-1) + ... + log(1). So we cannot use the
#'   built-in function for the pdf and afterwards log it, but need to implement
#'   the function ourselves instead.
#'
#' @param obs numeric, observed values
#' @param prob numeric, probabilities
#' @param logFacs (optional) Precalculated vertical vector of log factorials,
#' created using the function GetLogFactorials(n),
#' where n is at least the sum of all observations.
#' @author Johan Gustafsson, Juan Inda, <inda@@chalmers.se>
#' @return vector
#'
#' \dontrun{LogBinomialPDF(obs, prob)}
LogBinomialPDF <- function(obs, prob, logFacs=NULL){

  n = sum(obs)

  if(is.null(logFacs)){
    logFacs = GetLogFactorials(n)
  }

  #formula:
  #xi = obs i;
  #p = prob;
  # PDF = n!/(xi! * (n - xi)!) * p^xi * (1-p)^(n-xi)
  # =>
  # log(PDF) = log(n!) - log(xi!) - log((n-xi)!) + log(p^xi) + log((1-p)^(n-xi))
  #              a         b            c              d           e

  a = logFacs[n+1]
  b = logFacs[obs+1]
  c = logFacs[n-obs+1]
  d = obs * log(prob);
  e = (n-obs) * log(1-prob);
  lls = a - b - c + d + e;

  return(lls)
}

