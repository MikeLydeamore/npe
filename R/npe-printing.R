print.nphasefit <- function(x)
{
  phases <- length(x$estimates$lambdas)
  cat("Exponential with",phases, "phases fitted.\n")
  cat("Initial size:",sum(x$estimates$lambdas), "\n")
  cat("Relative proportions (3dp):\n")
  percentages <- round(x$estimates$lambda / sum(x$estimates$lambdas), digits = 3)
  for (i in 1:phases)
  {
    cat(i, ": ", percentages[i],"\n", sep = "")
  }
  cat("Decay rates (3dp):","\n")
  rates <- -1 * x$estimates$k
  inverses <- round(1/rates, digits = 3)

  for (i in 1:phases)
  {
    cat(i, ": ", rates[i]," (i.e. 1/", inverses[i],")\n", sep="")
  }

  cat("Sum of squared errors:",x$error,"\n")
  cat("AIC:",x$AIC,"\n")
}
