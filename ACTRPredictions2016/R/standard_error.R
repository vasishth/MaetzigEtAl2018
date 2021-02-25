#' @title Standard error of the mean
#'
#' @description Calculate the standard error of the mean.
#'
#' @param values numeric vector
#' @return standard_error_mean (numeric) The standard error of the mean of the
#' input vector.
#'
#' @export

standard_error <- function(values) {
  standard_error_mean <- sqrt(var(values) / length(values))
  return(standard_error_mean)
}
