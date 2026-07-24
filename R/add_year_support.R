#' Flag Years with Insufficient Recapture Support for Year-Effect Estimation
#'
#' Computes, for each modelled year, how many "animal-years" of recapture
#' data actually inform that year (i.e. how many tagged animals had at least
#' one moult opportunity while that year was current), and flags which years
#' meet a minimum support threshold. Must be run on \code{datain}
#' \strong{before} calling \code{Makepin(TemporalGrowth = TRUE)}, since
#' \code{Makepin} sizes the \code{Sraw} parameter vector based on the number
#' of supported years.
#'
#' @param datain List. Must contain \code{nyears}, \code{relyr}, \code{relts},
#'   \code{tsteps}, \code{ntsteps}.
#' @param min_support Integer. Minimum number of animal-years of support
#'   required for a year to be freely estimated. Years below this threshold
#'   have their growth-scaling effect fixed at exactly 0 (average growth)
#'   inside \code{\link{growmod}}, rather than being left to float on the
#'   sum-to-zero ridge penalty's degenerate equilibrium — an unsupported year
#'   does not shrink toward 0 under that penalty alone, it collapses onto
#'   whatever value satisfies the constraint given the other years, which can
#'   be arbitrarily large if it happens to land on the position TMB derives
#'   automatically (previously always the last chronological year). Default
#'   is 1 (any support at all counts as supported).
#'
#' @return \code{datain}, with two additional elements:
#' \describe{
#'   \item{yr_support}{Integer vector, length \code{nyears}. Count of
#'     animal-years of recapture support per modelled year.}
#'   \item{yr_supported}{Logical vector, length \code{nyears}. \code{TRUE}
#'     where \code{yr_support >= min_support}.}
#' }
#'
#' @details
#' Uses the same absolute-time-index arithmetic as \code{\link{growmod}}'s
#' likelihood loop (\code{relyr}/\code{relts} -> \code{abs0} -> \code{abst} ->
#' \code{yearvec}), so "supported" here means exactly what \code{growmod(..., TemporalGrowth = TRUE)}
#' would use to decide which years an animal's liberty period passes through
#' — this function is not a separate approximation, it is the same logic run
#' once, ahead of time, on plain (non-AD) data.
#'
#' @examples
#' \dontrun{
#' datain <- add_year_support(datain, min_support = 3)
#' table(datain$yr_supported)          # how many years pass the threshold
#' which(!datain$yr_supported)         # which years will be fixed at S = 0
#' pin <- Makepin(TemporalGrowth = TRUE)   # sizes Sraw from datain$yr_supported
#' }
#'
#' @seealso \code{\link{Makepin}}, \code{\link{growmod}}
#'
#' @export
add_year_support <- function(datain, min_support = 1) {
  # relyr must be a 1-based index into 1..nyears, NOT a raw calendar year.
  # A raw calendar year (e.g. 2002) silently produces nonsense here: every
  # animal gets pmin()-clamped into the last modelled year regardless of its
  # actual release year, giving one wildly "supported" year and everything
  # else showing zero support -- with no error to catch it. Guard against
  # that shape of mistake explicitly rather than letting it fail silently.
  if (any(datain$relyr < 1) || any(datain$relyr > datain$nyears)) {
    stop("datain$relyr has values outside 1..", datain$nyears,
         " (range: ", min(datain$relyr), "-", max(datain$relyr), "). ",
         "growmod/add_year_support expect relyr as a 1-based index into ",
         "the modelled years, not a raw calendar year. Convert first, e.g.:\n",
         "  yr0 <- min(datain$relyr)\n",
         "  datain$relyr <- datain$relyr - yr0 + 1")
  }

  abs0 <- (datain$relyr - 1) * datain$ntsteps + datain$relts
  yr_support <- integer(datain$nyears)

  for (r in seq_along(datain$tsteps)) {
    if (datain$tsteps[r] > 0) {
      abst <- abs0[r] + (0:(datain$tsteps[r] - 1))
      yrs  <- ((abst - 1) %/% datain$ntsteps) + 1
      yrs  <- pmin(yrs, datain$nyears)
      yr_support[unique(yrs)] <- yr_support[unique(yrs)] + 1
    }
  }

  datain$yr_support   <- yr_support
  datain$yr_supported <- yr_support >= min_support

  datain
}
