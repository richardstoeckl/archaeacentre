#' Plot a microbial growth curve
#'
#' @description
#' This functions plots a microbial growth curve with the option to add a confidence interval.
#' Under the hood, it is a simple ggplot2 wrapper, with defaults set depending on the type of plot.
#' As such, it can be expanded upon using the usual ggplot2 syntax.
#'
#' @param data A data frame containing the data to plot.
#' @param x The x-axis variable. Usually the time in hours after starting the experiment.
#' @param y The y-axis variable. Usually the concentration of microbes per mL as a raw number
#' @param group The grouping variable. Usually the name of the strain or condition.
#' @param color The color variable. Usually the name of the strain or condition.
#' @param type The type of plot to create. Currently only 'robert' is supported. This sets the default look for the plot.
#' @inherit ggplot2::ggplot seealso
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon scale_y_log10 annotation_logticks theme scale_x_continuous
#' @importFrom ggpubr theme_pubr
#' @importFrom dplyr mutate group_by ungroup
#' @export
plotGrowthCurve <- function(data, x, y, group, color, type = "robert") {
    dataForPlot <- data %>%
        dplyr::group_by(group) %>%
        dplyr::mutate(min = min(y), max = max(y)) %>%
        dplyr::ungroup()

    if (type == "robert") {
        plot <- ggplot(dataForPlot, aes(x = x, y = y, group = group, color = color)) +
            geom_line() +
            geom_point() +
            geom_ribbon(aes(ymin = min, ymax = max, fill = color), alpha = 0.5) +
            scale_y_log10(labels = label_log(digits = 2)) +
            annotation_logticks(sides = "l") +
            theme_pubr() +
            theme(panel.grid.minor = element_blank()) +
            scale_x_continuous(expand = c(0, 1))
    } else {
        print("Unknown plot type. Currently only 'robert' is supported.")
    }
    return(plot)
}
