#' @title Wrapper for ACTRPRedictions2016 preprocessing functions
#'
#' @description This function performs all preprocessing steps that are necessary for
#' visualising data simulated via the ACT-R sentence processing model, for the
#' first experiment in Grodner & Gibson (2005). This will be generalised in the future
#' to fit more experiments.
#'
#' @param path string, the directory where the \code{paramsearch.txt} and fixation files are located
#' @return \code{results_list} list(list(data.frame, data.frame)), containing means and standard
#' errors by region and condition, in one dataframe per simulation set
#'
#' @export

preprocess_data <- function(path=".") {
  p <- get_parameters(path)
  fixation_list <- get_data(path, p)
  data_list <- extract_main_emb_V(calculate_reading_measures(fixation_list))
  results_list <- mean_se_by_region(data_list)

  return(results_list)
}
