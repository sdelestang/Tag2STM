#' Create Initial Parameter List for Tag-Recapture Growth Model
#'
#' Generates a list of initial parameter values (pin file) for fitting
#' tag-recapture growth models using RTMB/TMB. See previous version's
#' roxygen for LsigError/LsigGrow/LMerrorRelsigma/LMerrorRecsigma/
#' TemporalGrowth -- unchanged.
#'
#' @param avgrowth Numeric, default 2. Starting value (mm) for the maximum
#'   growth increment (\code{Amax}) in the \code{Growth_par} double
#'   logistic used by \code{\link{growmod}} -- used directly as the
#'   \code{log(Amax)} starting value for every row of \code{Growth_par}.
#'   (Previously accepted but not implemented, retained only for backward
#'   compatibility; it now does something.)
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
#' @details
#' \strong{Growth_par starting values.} The four shape parameters of the
#' double logistic (\code{P1}, \code{P2}, \code{P3}, \code{P5} -- see
#' \code{\link{growmod}} for the curve itself) are seeded from
#' \code{datain$lbin} rather than hardcoded, so they scale sensibly
#' whatever species/size range is in use:
#' \itemize{
#'   \item \code{P2} starts at \code{median(datain$lbin)} -- an
#'     inflection point near the middle of the observed size range is a
#'     reasonable, uninformative starting guess.
#'   \item \code{P1} starts at \code{diff(range(datain$lbin)) / 20} (a
#'     fairly steep initial logistic scale) and \code{P3} at
#'     \code{diff(range(datain$lbin)) / 5} (shallower). \code{P1 < P3} at
#'     the starting values is only a starting-value convenience, NOT a
#'     constraint enforced anywhere -- \code{growmod} does not care which
#'     of the two ends up larger after fitting; the data decide.
#'   \item \code{P5} (the swap/blend transition width) starts at
#'     \code{diff(range(datain$lbin)) / 20}, i.e. roughly as sharp as the
#'     steeper of the two growth logistics.
#'   \item \code{Amax} starts at \code{avgrowth} (mm), estimated on the
#'     log scale.
#' }
#' Every row of \code{Growth_par} is initialised identically (same
#' convention as \code{Pmoult_par}); \code{\link{Makemap}} decides which
#' rows are estimated independently versus fixed/shared.
#'
#' @return A named list of initial parameter values. Key elements:
#' \describe{
#'   \item{Growth_par}{Numeric MATRIX, \code{ntsteps x 5} -- one
#'     5-parameter double-logistic growth-at-length curve per season,
#'     columns \code{log(Amax)}, \code{P2}, \code{log(P1)}, \code{log(P3)},
#'     \code{log(P5)} (see \code{\link{growmod}} for the curve itself and
#'     Details above for the starting values). Replaces the former
#'     \code{growth_vecpar} (an \code{nlbin * ntsteps}-length per-bin
#'     random walk). \code{Makemap} determines which rows are estimated
#'     independently versus fixed/shared (see \code{\link{Makemap}}).}
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
MakepinPar <- function(avgrowth = 2,
                    LsigError = NULL,
                    LsigGrow = log(2),
                    LMerrorRelsigma = 0,
                    LMerrorRecsigma = 0,
                    TemporalGrowth = NULL) {
  # Both now default to whatever Makedata() stored in datain, so they are
  # set once rather than passed separately to Makedata and Makepin (where
  # they could disagree). An explicit argument still wins, so existing
  # calls like Makepin(TemporalGrowth = TRUE) behave as before.
  if (is.null(TemporalGrowth)) TemporalGrowth <- isTRUE(datain$TemporalGrowth)
  if (is.null(LsigError)) {
    LsigError <- if (!is.null(datain$LsigError_init)) datain$LsigError_init else log(2)
  }

  nlbin   <- length(datain$lbin)
  ntsteps <- datain$ntsteps
  # Prefer datain$nobs (set by Makedata) over nrow(tdat): the random effects
  # must be sized to the SAME rows datain was built from, and reading the
  # global tdat silently breaks if it has been filtered since.
  nobs <- if (!is.null(datain$nobs)) datain$nobs else nrow(tdat)

  # Growth_par starting values -- see Details above. Derived from
  # datain$lbin so they scale with whatever species/size range is in use,
  # rather than hardcoding fixed mm values that would only suit one
  # species. avgrowth now feeds Amax directly (previously accepted but
  # unused).
  lbin_range <- diff(range(datain$lbin))
  growth_start <- c(
    log(avgrowth),           # log(Amax)
    median(datain$lbin),     # P2
    log(lbin_range / 20),    # log(P1)  -- steeper logistic scale
    log(lbin_range / 5),     # log(P3)  -- shallower logistic scale
    log(lbin_range / 20)     # log(P5)  -- swap/blend transition width
  )

  pin <- list(
    # ntsteps x 5 matrix -- growmod indexes Growth_par[ns, 1:5], so this
    # must be an actual matrix, not a length-5 vector. Every row starts
    # identical; Makemap decides which rows are freed independently.
    Growth_par = matrix(rep(growth_start, each = ntsteps),
                        nrow = ntsteps, ncol = 5),
    LsigError = LsigError,
    LsigGrow = LsigGrow,
    MerrorRel = rep(0, nobs),
    LMerrorRelsigma = LMerrorRelsigma,
    MerrorRec = rep(0, nobs),
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

  # Tagging-induced moult suppression as a recovery function of total time
  # at liberty (see growmod). Present only when datain$suppress is TRUE,
  # and must be present at trace time whenever it is -- like mpy_split_par,
  # the parameter vector's shape is fixed when MakeADFun() traces.
  #
  #   r = r0 + (1 - r0) * plogis((liberty_days - lib50) / slope)
  #
  # r0_par    : logit(r0). Starts at logit(0.4) ~= -0.405, i.e. animals
  #             recaptured immediately show ~40% of the expected moult
  #             probability -- roughly the deep-sea crab ratio (2.4 mm at
  #             one opportunity against ~5.5 mm per opportunity later).
  #             For a program with no handling effect the fit should push
  #             r0 toward 1 (r0_par -> +Inf).
  # lib50_par : log(days) at half recovery. Starts at log(300).
  # slope_par : log(days) transition width. Starts at log(60), i.e. most
  #             of the recovery spans roughly +/- 4 months around lib50.
  #             This is the least well-informed of the three -- liberty is
  #             coarsely distributed -- so fix it if its SE is large.
  if (isTRUE(datain$suppress)) {
    pin$r0_par    <- qlogis(0.4)
    pin$lib50_par <- log(300)
    pin$slope_par <- log(60)
  }

  if (TemporalGrowth) {
    if (isTRUE(datain$period_mode)) {
      # One effect per period, last anchored by the sum-to-zero constraint.
      if (is.null(datain$nperiods)) {
        stop("datain$period_mode is TRUE but datain$nperiods is missing -- ",
             "build datain with Makedata(..., period_col = ...).")
      }
      if (datain$nperiods < 2) {
        stop("Only ", datain$nperiods, " period(s) -- nothing to estimate. ",
             "Use period_mode = FALSE, or check the period column.")
      }
      pin$Sraw <- rep(0, datain$nperiods - 1)
    } else {
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
  }

  attr(pin, "TemporalGrowth") <- TemporalGrowth

  return(pin)
}
