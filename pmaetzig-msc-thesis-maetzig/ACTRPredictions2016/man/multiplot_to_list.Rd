% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiplot_to_list.R
\name{multiplot_to_list}
\alias{multiplot_to_list}
\title{Create multiple plots from a \code{results_list} and return them as a list}
\usage{
multiplot_to_list(results_list, indices)
}
\arguments{
\item{results_list}{list(data.frame), containing means and standard
errors by region and condition, in one dataframe per simulation set}

\item{indices}{numeric vector, indices for subsetting the \code{results_list}}
}
\value{
\code{plot_list}, a list of ggplot objects
}
\description{
This function creates multiple plots from data of a \code{results_list}.
Which datasets from the \code{results_list} are plotted is determined by a vector of
indices. At the moment, the only possible indexing method is by integers; support for
other indexing methods may be added in the future.
}

