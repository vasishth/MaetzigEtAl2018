#' @title Calculate reading measures for ACT-R simulations
#'
#' @description Use \code{\link{em2}} to calculate reading measures for the
#' fixation output of ACT-R simulations.
#'
#' @param fixation_list list(data.frame) with fixation output of ACT-R simulations
#' @return \code{data_list} list(list(data.frame, data.frame)), containing
#' data.frames with reading measures for every region.
#'
#' @seealso \code{\link{get_data}} for creating a \code{fixation_list} from the
#' working directory.
#'
#' @export

calculate_reading_measures <- function(fixation_list) {
  data_list <- vector(mode="list", length=length(fixation_list))

  for (i in 1:length(fixation_list)) {
    data <- em2::em2(fixation_list[[i]][[1]]$roi, fixation_list[[i]][[1]]$dur,
                          fixation_list[[i]][[1]][, c("exp", "iter", "condition")])
    data_name <- sub(pattern="fix", replacement="data", names(fixation_list[[i]])[1])

    tmp_list <- vector(mode="list", length=2)
    names(tmp_list) <- c(data_name, "params")
    tmp_list[[1]] <- data
    tmp_list[[2]] <- fixation_list[[i]][[2]]  ## parameter dfs are just carried along, always

    data_list[[i]] <- tmp_list
  }

  return(data_list)
}
