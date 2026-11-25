#' Colour, point shape and line type of the benchmarked methods
#'
#' A prolfqua method (`prolfqua_<facade>`) takes the style of its facade from `facade_styles.tsv` of prolfqua: the
#' colour of its family in the factory figure of the prolfqua README, and a point shape and line type per method within
#' the family. MSstats, proDA and msqrob2 are dark grey with open shapes. The vignettes style every plot with these
#' scales, so a method looks the same in all of them.
#'
#' @param methods method names as in the benchmark tables and the score history; a name without a style is an error.
#' @return A list of ggplot2 scales for the `colour`, `shape` and `linetype` aesthetics, with the methods in the order
#'   of the factory figure and the external tools last.
#' @export
#' @examples
#' points <- data.frame(x = 1:3, y = 1:3, method = c("prolfqua_lm", "prolfqua_limma", "MSstats"))
#' ggplot2::ggplot(points, ggplot2::aes(x, y, colour = method, shape = method)) +
#'   ggplot2::geom_point(size = 3) +
#'   method_scales(points$method)
method_scales <- function(methods) {
  facades <- utils::read.delim(
    system.file("extdata", "facade_styles.tsv", package = "prolfqua", mustWork = TRUE),
    stringsAsFactors = FALSE
  )
  facades$method <- paste0("prolfqua_", facades$method)
  external <- data.frame(
    method = c("MSstats", "proDA", "msqrob2", "msqrob2_protein", "msqrob2_hurdle"),
    colour = "#4d4d4d",
    shape = c(1, 2, 0, 5, 6),
    linetype = c("solid", "longdash", "dashed", "dotdash", "dotted")
  )
  styles <- rbind(facades[c("method", "colour", "shape", "linetype")], external)
  unknown <- setdiff(methods, styles$method)
  if (length(unknown) > 0) {
    stop("No style for method(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  styles <- styles[styles$method %in% methods, ]
  by_method <- function(column) stats::setNames(styles[[column]], styles$method)
  list(
    ggplot2::scale_colour_manual(values = by_method("colour"), breaks = styles$method),
    ggplot2::scale_shape_manual(values = by_method("shape"), breaks = styles$method),
    ggplot2::scale_linetype_manual(values = by_method("linetype"), breaks = styles$method)
  )
}
