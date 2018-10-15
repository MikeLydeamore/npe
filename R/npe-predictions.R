predict.nphasefitallinone <- function(fit, newdata, time_variable="t", group_variable="treatment_day")
{
  groups <- unique(newdata[, group_variable])
  num_phases <- length(fit$estimates$k)

  l <- lapply(groups, function(g) {
    whichgroup <- which(g == groups)
    if (num_phases == 1)
      y <- exp(fit$estimates$k[whichgroup] * newdata[,time_variable])

    else
    {
      y <- 0
      for (i in 1:(num_phases-1))
      {
        magicindex <- (i-1)+1 + (num_phases-1)*(whichgroup-1) #Starting to think this would be easier if i just spat out all the lambdas
        y <- y + fit$estimates$lambda[magicindex] * exp(fit$estimates$k[whichgroup] * newdata[,time_variable])
      }

      startindex <- 1 + (num_phases-1)*(whichgroup-1)
      endindex <- (num_phases-1 - 1) + 1 + (num_phases-1)*(whichgroup-1)
      y <- y + (1 - sum(fit$estimates$k[startindex:endindex])) * exp(fit$estimates$k[whichgroup] * newdata[,time_variable])
    }

    d <- data.frame("t"=newdata[,time_variable], y = y, group_variable = g)
    colnames(d)[1] <- time_variable
    colnames(d)[3] <- group_variable

    return (d)
  })

  return (do.call(rbind, l))
}
