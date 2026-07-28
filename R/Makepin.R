#' Create Initial Parameter List for Tag-Recapture Growth Model
#'
#' Generates a list of initial parameter values (pin file) for fitting
#' tag-recapture growth models using RTMB/TMB. See previous version's
#' roxygen for avgrowth/LsigError/LsigGrow/LMerrorRelsigma/LMerrorRecsigma/
#' TemporalGrowth -- unchanged.
#'
#' @param avgrowth Numeric, default 2. Not currently implemented; retained
#'   for backward compatibility.
#' @param LsigError Numeric, default log(2). Log SD of baseline measurement
#'   error.
#' @param LsigGrow Numeric, default log(2). Log CV of growth increments.
#' @param LMerrorRelsigma Numeric, default 0. Log SD of individual release
#'   measurement-error random effects.
#' @param LMerrorRecsigma Numeric, default 0. Log SD of individual
#'   recapture measurement-error random effects.
#' @param TemporalGrowth Logical, default FALSE. If TRUE, adds \code{Sraw}
#'   for year-specific growth scaling (requires
#'   \code{datain <- add_year_support(datain)} first).
#'
#' @return A named list of initial parameter values. Key elements:
#' \describe{
#'   \item{growth_vecpar}{Length \code{nlbin * ntsteps}, as before.}
#'   \item{Pmoult_par}{Numeric MATRIX, \code{ntsteps x 2} (intercept, slope
#'     columns) -- one logistic P(moult)-by-size curve per season:
#'     \code{Pmoult(ns, fm) = plogis(Pmoult_par[ns, 1] + Pmoult_par[ns, 2] *
#'     lbin[fm])}, further reparameterised inside \code{growmod} to
#'     asymptote at \code{mpy_floor[ns]} rather than 0 when
#'     \code{datain$mpy > 0} -- see \code{\link{growmod}}. Every row is
#'     initialised to the same starting values, \code{c(0.756, -0.0120)}.
#'     \code{Makemap} determines which rows are estimated independently
#'     versus fixed/shared (see \code{\link{Makemap}}).}
#'   \item{mpy_split_par}{Only present when \code{length(datain$goodts) > 1}.
#'     Numeric vector of length \code{length(datain$goodts) - 1}, initialised
#'     to all zeros (an exactly even split of \code{datain$mpy} across
#'     seasons at the start of optimisation). Estimated by \code{growmod}
#'     via a softmax transform -- see \code{\link{growmod}} for the
#'     minimum-moults-per-year floor this feeds into.}
#'   \item{Sraw}{Only present when \code{TemporalGrowth = TRUE}, as before.}
#' }
#'
#' @seealso \code{\link{growmod}}, \code{\link{Makemap}}
#' @export
Makepin <- function(avgrowth = 2,
                    LsigError = log(2),
                    LsigGrow = log(2),
                    LMerrorRelsigma = 0,
                    LMerrorRecsigma = 0,
                    TemporalGrowth = FALSE) {
  nlbin   <- length(bins$lbin)
  ntsteps <- datain$ntsteps

  pin <- list(
    growth_vecpar = rep(c(rep(log(0.2), nlbin - 1), log(0.01)), ntsteps),
    LsigError = LsigError,
    LsigGrow = LsigGrow,
    MerrorRel = rep(0, nrow(tdat)),
    LMerrorRelsigma = LMerrorRelsigma,
    MerrorRec = rep(0, nrow(tdat)),
    LMerrorRecsigma = LMerrorRecsigma,
    # ntsteps x 2 matrix -- growmod indexes Pmoult_par[ns, 1]/[ns, 2], so
    # this must be an actual matrix, not a length-2 vector. Every row
    # starts identical; Makemap decides which rows are freed independently.
    Pmoult_par = matrix(rep(c(0.756, -0.0120), each = ntsteps),
                         nrow = ntsteps, ncol = 2)
  )

  # Estimated split of datain$mpy across goodts seasons (softmax
  # construction in growmod). Only present when there's more than one
  # goodts season -- must be present at trace time whenever that's true,
  # regardless of what mpy itself is currently set to (mpy can be changed
  # between fits via datain$mpy <- ... without re-tracing; the parameter
  # vector's shape can't). All zeros = an exactly even split at the start
  # of optimisation, unbiased since you don't know in advance which season
  # should carry more of the floor.
  n_goodts <- length(datain$goodts)
  if (n_goodts > 1) {
    pin$mpy_split_par <- rep(0, n_goodts - 1)
  }

  if (TemporalGrowth) {
    if (is.null(datain$yr_supported)) {
      stop("datain$yr_supported not found. Run datain <- add_year_support(datain) ",
           "before calling Makepin(TemporalGrowth = TRUE) -- Sraw is sized from the ",
           "number of years with adequate recapture support, not raw nyears.")
    }
    n_supported <- sum(datain$yr_supported)
    if (n_supported < 2) {
      stop("Fewer than 2 supported years in datain$yr_supported (", n_supported,
           ") -- not enough information to estimate any year effects. Consider ",
           "lowering min_support in add_year_support(), or fitting with TemporalGrowth = FALSE.")
    }
    pin$Sraw <- rep(0, n_supported - 1)
  }

  attr(pin, "TemporalGrowth") <- TemporalGrowth

  return(pin)
}
