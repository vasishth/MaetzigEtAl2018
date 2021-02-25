#' @title Calculate differences between OR and SR for several simulation runs
#'
#' @description At the moment specialised to the analysis of experiment 1 of
#' Grodner & Gibson (2004). Calculates the difference between the mean reading
#' times of object relative (OR) and subject relative (SR) sentences. Also, at
#' the moment, takes the average of several simulation runs of the same parameter
#' value combination, if the indices vector is longer than one.
#'
#' @param results_list list(data.frame) a list of data.frames with mean and SE
#' reading times
#' @param indices list(numeric[]) a list of indices vectors; the first list item
#' of the output of \code{\link{get_reps_same_params}}
#' @param param_combinations character[] a vector containing unique parameter
#' value combinations of a search-param-space-em simulation run
#' @return effects \code{numeric[]} a vector of differences between mean OR and
#' SR reading times
#'
#' @export

calculate_effect_ORSR <- function(results_list, indices, param_combinations) {
  effects <- rep(0.0, length(indices))
  tmp <- rep(0, 4)

  for (i in 1:length(indices)) {
    for (j in 1:length(indices[[i]])) {
      k <- indices[[i]][j]
      tmp[j] <- results_list[[k]][1,]$`mean(TFT)` - results_list[[k]][2,]$`mean(TFT)`
    }
    effects[i] <- mean(tmp)
  }
  return(effects)
}
