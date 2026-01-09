#' Create Length Bin Structure
#'
#' Generates length bin midpoints, lower bounds, and upper bounds for size
#' transition matrix construction. Extends the first and last bins to capture
#' extreme values.
#'
#' @param start Numeric, starting length for bins (lower bound of first bin midpoint).
#' @param stop Numeric, ending length for bins (upper bound of last bin midpoint).
#' @param gap Numeric, width of each length bin. Default is 2 mm.
#'
#' @return A named list with three numeric vectors:
#' \describe{
#'   \item{lbin}{Length bin midpoints}
#'   \item{lbinL}{Length bin lower bounds (first bin extended by 20 mm)}
#'   \item{lbinU}{Length bin upper bounds (last bin extended by 20 mm)}
#' }
#'
#' @details
#' The first bin's lower bound is extended by 20 mm below the calculated value,
#' and the last bin's upper bound is extended by 20 mm above, to ensure all
#' observed lengths fall within the bin structure.
#'
#' @examples
#' # Create 2mm bins from 20-200mm carapace length
#' bins <- MakeLbin(start = 20, stop = 200, gap = 2)
#' head(bins$lbin)  # Midpoints: 19, 21, 23, ...
#'
#' # Create 5mm bins
#' bins <- MakeLbin(start = 50, stop = 150, gap = 5)
#'
#' @export
MakeLbin <- function(start, stop, gap = 2) {
  lbinU <- seq(start, stop, gap)
  lbinL <- lbinU - gap
  lbin  <- rowMeans(cbind(lbinU, lbinL))
  lbinL[1] <- lbinL[1] - 20
  lbinU[length(lbinU)] <- lbinU[length(lbinU)] + 20
  return(list(lbin = lbin, lbinL = lbinL, lbinU = lbinU))
}
