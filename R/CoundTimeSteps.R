#' Count Time Steps Between Release and Recapture
#'
#' Calculates the total number of time steps at liberty for a tag-recapture
#' observation, accounting for fractional time step allocations and spanning
#' multiple years.
#'
#' @param x Numeric vector of length 6 containing release and recapture information:
#'   element 1 = release year, element 2 = recapture year, element 3 = release month,
#'   element 4 = release half-month, element 5 = recapture month, element 6 = recapture half-month
#'
#' @return Numeric value representing total time steps at liberty, rounded to
#'   the nearest integer.
#'
#' @details
#' This function expects a `times` object in the calling environment (created
#' by `Maketimes()`), which maps month/half-month combinations to time steps
#' with fractional allocations.
#'
#' The calculation handles three cases:
#' \itemize{
#'   \item **Same year** (nyrs = 0): Sum fractional time steps from release to recapture
#'   \item **Multi-year**: Sum from release to year-end, plus complete intermediate
#'     years, plus start of final year to recapture
#' }
#'
#' @examples
#' \dontrun{
#' # After creating times object with Maketimes()
#' # Calculate time steps for tag released Nov 2020, recaptured Mar 2022
#' obs <- c(2020, 2022, 11, 0, 3, 1)  # Release Nov (1st half), Recap Mar (2nd half)
#' tsteps <- cntTsteps(obs)
#' }
#'
#' @export
cntTsteps <- function(x) {
  nyrs <- x[2] - x[1]
  spos <- match(paste(x[3], x[4]), times$id)
  epos <- match(paste(x[5], x[6]), times$id)
  if(nyrs == 0) Tfrac <- sum(times$tfrac[spos:epos])
  if(nyrs > 0) {
    Tfrac <- sum(times$tfrac[spos:nrow(times)]) +
      ((nyrs - 1) * sum(times$tfrac[1:nrow(times)])) +
      sum(times$tfrac[1:epos])
  }
  Tfrac <- round(Tfrac, 0)
  return(Tfrac)
}
