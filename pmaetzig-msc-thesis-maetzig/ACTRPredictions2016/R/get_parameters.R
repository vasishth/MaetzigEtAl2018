#' @title Get parameters from paramsearch.txt
#'
#' @description Returns a dataframe containing the information in paramsearch.txt, i.e., the
#' running number of the simulation set, the name of the experiment, the name of
#' the fixations file and the parameter values of the simulation. If \code{paramsearch.txt}
#' is in the current working directory, this function can be called without any
#' arguments.
#'
#' @param path string, the path to where \code{paramsearch.txt} is located; with no "/" in the end
#' @param file string, the name of the file (mostly \code{paramsearch.txt})
#' @return paramsearch data.frame containing paramsearch metainformation.
#'
#' @export

get_parameters <- function(path=".", file="paramsearch.txt") {
  param_file                    <- paste(path, file, sep="/")
  paramsearch                   <- read.table(param_file, header=FALSE)
  colnames(paramsearch)         <- c("experiment", "set", "file", "params")
  paramsearch$file              <- as.character(paramsearch$file)
  paramsearch$simulation_index  <- strsplit(paramsearch$file[1], "/")[[1]][4]
  paramsearch$params            <- tolower(as.character(paramsearch$params))

  # First loop: split one arbitrary entry of paramsearch$params (because each
  # row contains all varied parameters) and create variables in the paramsearch
  # data.frame according to these parameters.
  params <- substr(paramsearch$params[1], 2, nchar(paramsearch$params[1]) - 1)
  params <- strsplit(params, " ")[[1]]
  params <- params[c(TRUE, FALSE)]

  for (i in params) {
    paramsearch[i] <- rep(NA, dim(paramsearch)[1])
  }

  # Second loop: for each row in paramsearch, split the current entry of
  # paramsearch$params, and assign to the corresponding variables (e.g.,
  # paramsearch$ga) the value of this parameter in the current simulation set.
  for (i in 1:dim(paramsearch)[1]) {
    tmp <- substr(paramsearch$params[i], 2, nchar(paramsearch$params[i]) - 1)  # remove parentheses
    tmp <- strsplit(tmp, " ")[[1]]
    paramsearch$params[i] <- paste(tmp, collapse="-")
    params <- tmp[c(TRUE, FALSE)]

    for (j in params) {
      paramsearch[j][i,] <- as.numeric(tmp[which(tmp %in% j)+1])
    }
  }

  return(paramsearch)
}
