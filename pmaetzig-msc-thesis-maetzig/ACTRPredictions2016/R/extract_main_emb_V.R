#' @title Extract data for main and embedded verb regions
#'
#' @description Use \code{dplyr::filter} to subset the data of a
#' \code{data_list} to include only the lines for main verb and embedded
#' verb regions. Currently only works for experiment 1 of Grodner & Gibson
#' (2005) (\code{gg-exp1} in ACT-R).
#'
#' @param data_list list(data.frame) containing fixation output of an ACT-R
#' simulation.
#' @return \code{data_list} \code{list(list(data.frame, data.frame))}, containing
#' data.frames with reading measures for the subset regions
#'
#' @seealso \code{\link{get_data}} for creating a \code{data_list} from given
#' fixation files.
#'
#' @export

extract_main_emb_V <- function(data_list) {
  ### From a list of data.frames containing reading measures of ACT-R output,
  ### extracts only the data for the main verb and embedded verb regions,
  ### and then renames the regions in the data.frames
  ###
  ### Args: data_list: list(data.frame)
  ### Value: data_list: list(data.frame)

  for (i in 1:length(data_list)) {
    data_list[[i]][[1]] <- dplyr::filter(data_list[[i]][[1]], (condition == "SR" & roi %in% c(4, 10)) | (condition == "OR" & roi %in% c(6, 10)))
    data_list[[i]][[1]]$roi <- ifelse(data_list[[i]][[1]]$condition == "SR",
                                 ifelse(data_list[[i]][[1]]$roi == 4, "EmbV", "MV"),
                                 ifelse(data_list[[i]][[1]]$roi == 6, "EmbV", "MV"))
    data_list[[i]][[1]]$roi <- as.factor(data_list[[i]][[1]]$roi)
  }

  return(data_list)
}
