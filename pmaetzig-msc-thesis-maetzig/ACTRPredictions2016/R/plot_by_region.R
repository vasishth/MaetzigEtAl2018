#' @title Plot total reading time by region of interest and condition
#'
#' @description Use \code{ggplot2} to plot the mean total reading time for a
#' region of interest and experimental condition, including error bars.
#'
#' @param data list(data.frame), a \code{results_list} as returned by \code{\link{mean_se_by_region}}
#' @param maintitle string, the plot main title
#' @param label_x string, the label of the x axis
#' @param label_y string, the label of the y axis
#' @param limit_y 2x1 numeric vector, the limits of the y-axis (x-axis is not that important because it's
#' automatically assigned the regions, but y represents reading times)
#' @param d numeric, the offset for \code{position_dodge}
#' @param legendpos, 2x1 numeric vector, the position of the legend
#'
#' @export

plot_by_region <- function(data,
                           maintitle,
                           label_x="Region",
                           label_y="Total Reading Time (ms)",
                           limit_y=c(200, 500),
                           d=0.2,
                           legendpos=c(0, 0.2)) {

  ggplot2::ggplot(data, aes(x=data$roi, y=data$`mean(TFT)`, shape=condition, group=condition)) +
    geom_line(aes(linetype=condition), size=1, position=position_dodge(d)) +
    geom_point(position=position_dodge(d),size=4) +
    scale_linetype_manual(values=c("solid","dashed", "solid","dashed")) +
    theme_bw() +
    geom_errorbar(aes(ymin=data$`mean(TFT)` - 2*data$`standard_error(TFT)`,
                      ymax=data$`mean(TFT)` + 2*data$`standard_error(TFT)`),
                  width=.1, size=1,
                  position=position_dodge(d)) +
    xlab(label_x) +
    ylab(label_y) +
    ylim(limit_y) +
    ggtitle(maintitle)
    #theme(legend.justification=c(0, 0), legend.position=legendpos)
}
