#' @title Get indices for iterations with same parameter value combinations
#'
#' @description For the analysis of multiple sets of \code{paramsearch} experiments
#' in ACT-R (e.g. for model checking), this function returns indices of all the
#' iterations that share common parameter value combinations, so that they can be
#' easily plotted or analysed otherwise.
#'
#' @param results_list list(data.frame), containing means and standard
#' errors by region and condition, in one dataframe per simulation set
#' @return \code{indices_param_combinations} list(list(integer[]), character[]) containing [[1]]: n vectors of indices grouping
#' iterations with identical parameter value combinations; and [[2]]: a character vector containing the unique
#' parameter value combinations in the \code{results_list}
#'
#' @export

get_reps_same_params <- function(results_list) {
  names_split <- strsplit(names(results_list), "-data-")
  param_combinations <- c(length(results_list))  ## no better way than for-loop for now, so preallocating the vector

  for (i in 1:length(names_split)) {
    param_combinations[i] <- names_split[[i]][2]
  }

  param_combinations <- unique(param_combinations)
  indices <- lapply(param_combinations, grep, names(results_list))
  indices_param_combinations <- list(indices, param_combinations)

  return(indices_param_combinations)
}
