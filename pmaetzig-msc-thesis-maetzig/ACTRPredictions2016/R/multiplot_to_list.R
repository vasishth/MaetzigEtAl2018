#' @title Create multiple plots from a \code{results_list} and return them as a list
#'
#' @description This function creates multiple plots from data of a \code{results_list}.
#' Which datasets from the \code{results_list} are plotted is determined by a vector of
#' indices. At the moment, the only possible indexing method is by integers; support for
#' other indexing methods may be added in the future.
#'
#' @param results_list list(data.frame), containing means and standard
#' errors by region and condition, in one dataframe per simulation set
#' @param indices numeric vector, indices for subsetting the \code{results_list}
#' @return \code{plot_list}, a list of ggplot objects
#'
#' @export

multiplot_to_list <- function(results_list, indices) {
  plot_list <- vector(mode="list", length=length(indices))

  for (i in 1:length(plot_list)) {
    plot_list[[i]] <- plot_by_region(data=results_list[[indices[i]]],
                                         maintitle=strsplit(names(results_list)[indices[i]], "region_")[[1]][2])
  }

  return(plot_list)
}
