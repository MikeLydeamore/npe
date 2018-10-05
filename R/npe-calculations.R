#' Calculate value of n-phase exponential curve
#'
#' Calculates the value of an n-phase exponential curve with initial sizes, \eqn{\lambda} and
#' rate parameters, \eqn{k} at time \eqn{t}. \cr \cr
#' Note importantly that the rate parameters, \eqn{k}, are \strong{not} negated.
#'
#' @param lambda Vector specifying the initial size of each group
#' @param k Vector specifying the decay rates of each group
#' @param t Time(s) at which to evaluate curve
#'
#' @return A vector containing the value of the function,
#' \eqn{\sum_i \lambda_i exp(k_i \times t)}, at each time \eqn{t}.
#'
#' @examples
#' lambda <- c(300, 200)
#' k <- c(-0.1, -0.5)
#' t <- seq(1:10, by=0.1)
#' calculateNPE(lambda = lambda, k = k, t = t)
calculateNPE <- function(lambda, k, t)
{
  if (length(lambda) != length(k))
    stop("Length of lambda must equal length of k")
  l <- lapply(t, function(j)
    {
      s <- sapply(1:length(lambda), function(i) {
      lambda[i] * exp(k[i] * j)
      })
      return (sum(s))
  })

    ret <- do.call(c, l)

  return (ret)
}

#' Calculate residual some of squares for n-phase exponential curve
#'
#' Calculates the residual sum of squares for an n-phase exponential curve with
#' initial sizes, \eqn{\lambda} and rate parameters, \eqn{k}, for a given set of data. \cr\cr
#' Note importantly that the rate paramters, \eqn{k}, are \strong{not} negated.
#'
#' @param lambda Vector specifying the initial size of each group.
#' @param k Vector specifying the decay rates of each group.
#' @param data Data frame containing the observations.
#' @param time_variable Quoted string specifying the name of the time variable in \code{data}.
#' @param value_variable Quoted string specifying the name of the value variable in \code{data}.
#'
#' @return The residual sum of squared error.
calculateRSS <- function(lambda, k, data, time_variable = "t", value_variable = "value")
{
  if (any(lambda < 0))
    return (Inf)

  if (length(lambda) != length(k))
    stop("Length of lambda must equal length of k")

  errors <- sapply(1:nrow(data), function(i) {
    model <- npe::calculateNPE(lambda, k, data[i, time_variable])

    error <- data[i, value_variable] - model

    return (error^2)
  })
  return (sum(errors))
}



#' Calculate residual some of squares for n-phase exponential curve
#'
#' Calculates the residual sum of squares for an n-phase exponential curve with
#' initial sizes, \eqn{\lambda} and rate parameters, \eqn{k}, for a given set of data. \cr\cr
#' Note importantly that the rate paramters, \eqn{k}, are \strong{not} negated.
#'
#' @param lambda Vector specifying the initial size of each group.
#' @param k Vector specifying the decay rates of each group.
#' @param data Data frame containing the observations.
#' @param time_variable Quoted string specifying the name of the time variable in \code{data}.
#' @param value_variable Quoted string specifying the name of the value variable in \code{data}.
#'
#' @return The residual sum of squared error.
calculateRSSNormalised <- function(lambda, k, data, time_variable = "t", value_variable = "value")
{
  if (length(k) > 1)
    lambda <- c(lambda, 1-sum(lambda))

  if (any(lambda < 0))
    return (Inf)

  if (any(lambda > 1))
    return (Inf)

  errors <- sapply(1:nrow(data), function(i) {
    model <- npe::calculateNPE(lambda, k, data[i, time_variable])

    error <- data[i, value_variable] - model

    return (error^2)
  })
  return (sum(errors))
}
