get_min_params <- function(data, subject) {
  ref <- which(data[c(subject)]==min(data[c(subject)]))
  params <- data[ref, c("GA", "DAT", "ANS")]
  return(params)
}