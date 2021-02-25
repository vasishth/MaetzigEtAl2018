#' @title Reading data files for ACT-R paramsearch
#'
#' @description Read all fixation files from an ACT-R paramsearch and create a list of
#' data.frames from them.
#'
#' @param path string of the path of the fixation files directory
#' @param paramsearch data.frame containing paramsearch metainformation
#' @return all_fixations \code{list(list(data.frame, data.frame))} with fixations of the ACT-R paramsearch simulations.
#'
#' @details This is currently only fit to output of the ACT-R function
#' \code{run-param-search-em} and the output it generates for experiment 1 of
#' Grodner & Gibson (2005). A more flexible framework will be added later.
#'
#' @seealso \code{\link{get_parameters}} for creating a paramsearch data.frame
#'
#' @export

get_data <- function(path=".", paramsearch) {
  column_names  <- c("exp", "iter", "condition", "roi", "word", "dur")
  cond_levels   <- c("OR", "SR")

  # create list of fixation files
  file_list <- list.files(path=path,
                          pattern="[[:digit:]{1}]-fixations[[:punct:]{1}]txt")

  # pre-allocate the final list
  all_fixations <- vector(mode="list", length=length(file_list))

  for (i in 1:length(file_list)) {
    df_name <- paste(paramsearch$set[i],
                     paste("fix", paramsearch$params[i], sep="-"), sep="-")
    df <- read.table(paste(path, file_list[i], sep="/"))
    colnames(df) <- column_names
    levels(df$condition) <- cond_levels

    current_parameters <- paramsearch[i, 6:ncol(paramsearch)]

    tmp_list <- vector(mode="list", length=2)
    names(tmp_list) <- c(df_name, "params")
    tmp_list[[1]] <- df
    tmp_list[[2]] <- current_parameters

    all_fixations[[i]] <- tmp_list #the respective pair of data and params
  }

  return(all_fixations)
}
