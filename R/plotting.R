#' Example Dataset for microbial Growth Curves
#'
#' @name growthData
#' @docType data
#' @keywords data
NULL


#' Prepare Data for growth curve plotting
#' @keywords internal
#' @param data A data frame containing the data to plot.
#' @param x The x-axis variable. Usually the time in hours after starting the experiment.
#' @param y The y-axis variable. Usually the concentration of microbes per mL as a raw number
#' @param grouping The grouping variable. Usually the name of the strain or condition.
#' @param color The color variable. Usually the name of the strain or condition.
#' @importFrom dplyr mutate group_by ungroup across all_of
.prepareDataForGrowthCurve <- function(data, x, y, grouping, color) {
    dataForPlot <- data %>%
        dplyr::group_by(across(all_of(grouping))) %>%
        dplyr::mutate(XVAR = as.double({{ x }}), YVAR = as.double({{ y }}), COLORVAR = {{ color }}) %>%
        dplyr::mutate(MINVAR = min(YVAR), MAXVAR = max(YVAR), MEANVAR = mean(YVAR))
    return(dataForPlot)
}

#' Plot a microbial growth curve
#'
#' @description
#' This functions plots a microbial growth curve with the option to add a confidence interval.
#' Under the hood, it is a simple ggplot2 wrapper, with defaults set depending on the type of plot.
#' As such, it can be expanded upon using the usual ggplot2 syntax.
#'
#' @param data A data frame containing the data to plot.
#' @param x The x-axis variable. Usually the time in hours after starting the experiment.
#' @param y The y-axis variable. Usually the concentration of microbes per mL as a `double`
#' @param grouping A character vector containing one or more column names for the grouping variables. Usually "c("timepoint", "organism")".
#' @param color The color variable. Usually the strain or condition.
#' @param type The type of plot to create. Currently only 'robert' is supported. This sets the default look for the plot.
#' @inherit ggplot2::ggplot seealso
#' @examples
#' \dontrun{
#' # load example data from this package
#' growthData <- archaeacentre::growthData
#' # plot the growth curve
#' plotGrowthCurve(growthData, timepoint, concentration, grouping = c("timepoint", "organism"), organism, type = "robert")
#' }
#' @import ggplot2
#' @importFrom ggpubr theme_pubr
#' @importFrom scales label_log
#' @export
plotGrowthCurve <- function(data, x, y, grouping, color, type = "robert") {
    dataForPlot <- .prepareDataForGrowthCurve({{ data }}, {{ x }}, {{ y }}, {{ grouping }}, {{ color }})

    if (type == "robert") {
        plot <- ggplot(dataForPlot, aes(x = XVAR, y = MEANVAR, color = COLORVAR, group = COLORVAR)) +
            geom_line(show.legend = FALSE) +
            geom_point() +
            geom_ribbon(aes(ymin = MINVAR, ymax = MAXVAR, fill = COLORVAR), alpha = 0.5, show.legend = FALSE) +
            scale_y_log10(labels = scales::label_log(digits = 2)) +
            annotation_logticks(sides = "l") +
            ggpubr::theme_pubr() +
            theme(panel.grid.minor = element_blank()) +
            scale_x_continuous(expand = c(0, 1))
    } else {
        print("Unknown plot type. Currently only 'robert' is supported.")
    }
    return(plot)
}
