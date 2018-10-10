#' Fit an n-phase exponential distribution
#'
#' Fits an n-phase distribution
#'
#' Under the hood, this function is optimising the residual sum of squares (see \code{\link{calculateRSS}})
#' using \code{\link[stats]{nlm}}. The only control parameter passed to the optimiser is \code{typsize},
#' which is the same as the initial guess.
#'
#' If not supplied, the initial guess for the initial sizes, \eqn{\lambda}, is 200/\code{num_phases}.
#'
#' If not supplied, the initial guess for the rate parameters, \eqn{k}, is \eqn{10^{-2:-2-(num_phases+1)}}.
#' That is, \eqn{k_0 = 0.01}, and each consecutive \eqn{k} decreases by an order of magnitude
#' from there. These guesses seem to work generally well, but may need tweaking for specific problems.
#'
#' @param num_phases Number of phases for exponential distributions.
#' @param data Data frame containing observations at some time points.
#' @param time_variable Quoted string specifying the name of the time variable in \code{data}.
#' @param value_variable Quoted string specifying the name of the value variable in \code{data}.
#' @param lambda_0 \emph{Optional} Initial guess for the size of \eqn{\lambda}.
#' @param k_0 \emph{Optional} Initial guess for the size of \eqn{k}.
#'
#' @return A list of class "nphasefit", which contains the following elements:
#' \itemize{
#'   \item \code{estimates}: Data frame with two columns: \code{lambda} and \code{k}, containing the optimal estimates of \eqn{\lambda} and \eqn{k} respectively.
#'   \item \code{error}: The residual sum of squares at the optimum.
#'   \item \code{AIC}: The AIC at the optimum.
#' }
fitnphaseexponential <- function(num_phases = 1, data, time_variable = "t", value_variable = "dna", lambda_0, k_0)
{
  #Initial guesses:
  if (missing(lambda_0))
    lambdas <- rep(200, times = num_phases)/num_phases
  else
    lambdas <- lambda_0

  if (missing(k_0))
    rates <- -1 * 10^(-2:(-2 - num_phases + 1))#rep(-0.01, times = num_phases)
  else
    rates <- k_0

  minimising_function <- function(p, data, time_variable, value_variable)
  {
    len <- length(p)/2
    lambdas <- p[1:len]
    rates <- p[(len+1):length(p)]

    npe::calculateRSS(lambda = lambdas, k = rates, data = data, time_variable = time_variable, value_variable = value_variable)
  }

  typsize <- c(lambdas, rates)

  opt <- stats::nlm(f = minimising_function, p = c(lambdas, rates), data = data, time_variable = time_variable, value_variable = value_variable, typsize = typsize)

  res <- list("estimates"=data.frame("lambdas"=opt$estimate[1:num_phases], "k"=opt$estimate[(num_phases+1):length(opt$estimate)]),
              "error"=opt$minimum,
              "AIC"=2*num_phases + nrow(data)*log(opt$minimum))
  class(res) <- "nphasefit"
  return (res)
}



fitnphaseexponentialnormalised <- function(num_phases = 1, data, time_variable = "t", value_variable = "dna", k_0)
{
  #Initial guesses:
  if (missing(k_0))
    rates <- -1 * 10^(-1:(-1 - num_phases + 1))#rep(-0.01, times = num_phases)
  else
    rates <- k_0

  lambdas <- seq(0.3, by = 0.3, length.out=num_phases-1)

  minimising_function <- function(p, data, time_variable, value_variable)
  {
    len <- floor(length(p)/2)
    lambdas <- p[1:len]
    rates <- p[(len+1):length(p)]
    print(lambdas)
    print(rates)

    npe::calculateRSSNormalised(lambda = lambdas, k = rates, data = data, time_variable = time_variable, value_variable = value_variable)
  }

  minimising_function_one <- function(p, data, time_variable, value_variable)
  {
    len <- 1
    rates <- p
    npe::calculateRSSNormalised(lambda = 1, k = rates, data = data, time_variable = time_variable, value_variable = value_variable)
  }

  if (num_phases == 1)
  {
    typsize <- rates

    opt <- stats::nlm(f = minimising_function_one, p = rates, data = data, time_variable = time_variable, value_variable = value_variable, typsize = typsize)

    res <- list("estimates"=data.frame("k"=opt$estimate),
                "error"=opt$minimum,
                "AIC"=2*num_phases + nrow(data)*log(opt$minimum))
  }
  else
  {
    typsize <- c(lambdas, rates)

    opt <- stats::nlm(f = minimising_function, p = c(lambdas, rates), data = data, time_variable = time_variable, value_variable = value_variable, typsize = typsize)
    print(opt)
    res <- list("estimates"=data.frame("lambdas"=c(opt$estimate[1:num_phases-1], 1-sum(opt$estimate[1:num_phases-1])), "k"=opt$estimate[(num_phases):length(opt$estimate)]),
                "error"=opt$minimum,
                "AIC"=2*num_phases + nrow(data)*log(opt$minimum))
  }


  class(res) <- "nphasefit"
  return (res)
}

fitnphaseexponentialnormalisedallinone <- function(num_phases = 1, data, time_variable = "t", value_variable = "dna", group_variable = "treatment_day", k_0, max_iter = 1e6)
{
  #Initial guesses:
  if (missing(k_0))
    rates <- -1 * 10^(-1:(-1 - num_phases + 1))
  else
    rates <- k_0

  num_groups <- length(unique(data[,group_variable]))
  lambdas <- c()
  for (i in 1:num_groups)
  {
    lambdas_all <- runif(n = num_phases, min = 0, max=1)
    lambdas_all <- lambdas_all/sum(lambdas_all)
    lambdas <- c(lambdas, lambdas_all[1:(length(lambdas_all)-1)])
  }

  minimising_function <- function(p, data, time_variable, value_variable, group_variable, num_phases)
  {
    len <- floor(length(p)/2)
    rates <- p[1:num_phases]
    lambdas <- p[(num_phases+1):length(p)]

    npe::calculateRSSAllInOneNormalised(lambda = lambdas, k = rates, data = data, time_variable = time_variable, value_variable = value_variable, group_variable)
  }
  typsize <- c(rates, lambdas)
  opt <- optimx::optimx(fn = minimising_function, par = c(rates, lambdas), data = data,
                    time_variable = time_variable, value_variable = value_variable,
                    group_variable = group_variable, num_phases = num_phases, control=list(all.methods=T))

  best <- which.min(opt$value)
  k <- as.numeric(opt[best, 1:num_phases])
  lambdas <- as.numeric(opt[best, (num_phases+1):(num_phases + length(lambdas))])

  estimates <- list("k"=k, "lambdas"=lambdas)
  res <- list("estimates"=estimates,
              "error"=opt$minimum,
                "AIC"=2*(length(opt$par)) + nrow(data)*log(opt$value),
              "lambda_0"=lambdas,
              "k_0"=rates)


  class(res) <- "nphasefitallinone"
  return (res)
}

