% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mean_se_by_region.R
\name{mean_se_by_region}
\alias{mean_se_by_region}
\title{Calculate mean and standard error for a reading measure}
\usage{
mean_se_by_region(data_list, measure = "TFT")
}
\arguments{
\item{data_list}{list(data.frame), containing dataframes of reading measures}

\item{measure}{string, the name of a reading measure (must be a column name
in a dataframe within \code{data_list}}
}
\value{
\code{results_list}, list(list(data.frame, data.frame)) containing means and standard
errors by region and condition, in one dataframe per simulation set
}
\description{
Use \code{dplyr} to calculate the mean and the standard error
for a given reading measure. The data in the given \code{data_list} will be
grouped by region of interest and condition automatically via
\code{dplyr::group_by} (this should be changed and generalised in a future
version).
}
\details{
This function exploits that ACT-R output is already sorted and
grouped by region of interest and condition. In the future, this should be
adapted to more general cases where input of ordered and grouped data cannot
be ascertained.
}

