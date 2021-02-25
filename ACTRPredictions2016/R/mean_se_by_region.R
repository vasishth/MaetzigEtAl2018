#' @title Calculate mean and standard error for a reading measure
#'
#' @description Use \code{dplyr} to calculate the mean and the standard error
#' for a given reading measure. The data in the given \code{data_list} will be
#' grouped by region of interest and condition automatically via
#' \code{dplyr::group_by} (this should be changed and generalised in a future
#' version).
#'
#' @param data_list list(data.frame), containing dataframes of reading measures
#' @param measure string, the name of a reading measure (must be a column name
#' in a dataframe within \code{data_list}
#'
#' @return \code{results_list}, list(list(data.frame, data.frame)) containing means and standard
#' errors by region and condition, in one dataframe per simulation set
#'
#' @details This function exploits that ACT-R output is already sorted and
#' grouped by region of interest and condition. In the future, this should be
#' adapted to more general cases where input of ordered and grouped data cannot
#' be ascertained.
#'
#' @export

mean_se_by_region <- function(data_list, measure="TFT") {
  name_mean       <- paste("mean", "(", measure, ")", sep="")
  # name_SE         <- paste("SE", "(", measure, ")", sep="")

  results_list    <- vector(mode="list", length=length(data_list))

  for (i in 1:length(data_list)) {
    data <- dplyr::group_by(data_list[[i]][[1]], roi, condition)
    tmp_df <- dplyr::summarise_(data, name_mean)
    tmp_df <- data.frame(tmp_df, dplyr::summarise(data, standard_error(TFT))[,3], check.names=FALSE)

#    tmp_list <- vector(mode="list", length=2)
#    tmp_df_name <- sub(pattern="data", replacement="results", names(data_list[[i]])[1])
#    names(tmp_list) <- c(tmp_df_name, "params")
#    tmp_list[[1]] <- tmp_df
#    tmp_list[[2]] <- data_list[[i]][[2]]  ## parameter dfs are just carried along, always

    complete_df <- merge(tmp_df, data_list[[i]][[2]])
    results_list[[i]] <- complete_df
  }

  return(results_list)
}
