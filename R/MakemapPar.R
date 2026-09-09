## Change from previous version (1):
##
## Pmoult_par was excluded from the growth_vecpar-style `turnon` mapping
## sweep and left entirely unmapped, which was fine when it was a length-2
## vector (both elements always estimated, nothing to fix). Now that
## Makepin builds it as an ntsteps x 2 matrix (see Makepin.R), leaving it
## completely unmapped would estimate a separate logistic for every season
## INCLUDING seasons not in goodts, which growmod's Pmoult_fn never even
## evaluates -- unidentifiable parameters wasting degrees of freedom and
## likely causing convergence trouble, same failure mode growth_vecpar's
## turnon mapping already exists to prevent.
##
## New default: goodts rows are freed independently (mirrors growth_vecpar
## exactly), non-goodts rows fixed at their Makepin initial values.
##
## New `Pmoult_shared` argument (default FALSE): set TRUE to instead map
## ALL goodts rows to the SAME two factor levels -- i.e. reproduces the old
## single-shared-curve behaviour exactly. Useful as a fallback if a
## particular goodts season turns out to have too little identifying
## information even with the broader internal identification mixture (see
## growmod's "seasonal selection caveat", now less restrictive than before
## since identification no longer requires single-opportunity liberty
## specifically -- but the caveat isn't eliminated, just relaxed).

## Change from previous version (2):
##
## growth_vecpar (and with it, the entire `pnames`/`turnon` generic
## catch-all block) is gone. That block existed for exactly one purpose:
## growth_vecpar was the sole parameter left over after excluding every
## sig/error-named parameter and the explicitly-listed Pmoult_par /
## mpy_split_par / r0_par / lib50_par / slope_par / Sraw -- so
## `pnames` always resolved to just `"growth_vecpar"`, and the elaborate
## exclusion filter was really just a roundabout way of finding it.
##
## Growth_par (see Makepin.R) is now built and mapped exactly like
## Pmoult_par -- an ntsteps x 5 matrix, own dedicated mapping block below
## -- rather than falling through the generic sweep, so:
##   (a) it is added to the exclusion list alongside Pmoult_par etc., and
##   (b) the `pnames`/`turnon`/`lapply` block is removed outright, since
##       with Growth_par excluded there is nothing left for it to catch --
##       keeping dead code around that computes an unused `turnon` vector
##       would only mislead the next reader into thinking some other
##       parameter still relies on it.
## `map` is now built as a plain empty list and populated by each
## parameter's own block, same as MerrorRel/LsigError/Sraw already were.
##
## New `Growth_shared` argument (default FALSE), exactly mirroring
## `Pmoult_shared`: FALSE frees each goodts row's 5 parameters
## independently (feasible now that Growth_par is 5 parameters per row
## rather than nlbin), TRUE maps every goodts row to the same 5 factor
## levels (a single shared growth curve across seasons) -- a reasonable
## fallback if a diagnostic run shows a particular season's row poorly
## identified on its own.

## -----------------------------------------------------------------------
## Full function
## -----------------------------------------------------------------------

#' Create Parameter Mapping for Selective Time Step Fitting
#'
#' (See previous version's roxygen for the general \code{re}/
#' \code{estTemporalGrowth} behaviour -- unchanged. This version adds
#' explicit handling for \code{Pmoult_par} and \code{Growth_par}, both
#' \code{ntsteps x n} matrices rather than flat vectors -- see
#' \code{\link{Makepin}}.)
#'
#' @param pin List, from \code{Makepin}.
#' @param re Logical, as before.
#' @param estTemporalGrowth Logical, as before.
#' @param Pmoult_shared Logical (default \code{FALSE}). If \code{FALSE}
#'   (default), each \code{goodts} row of \code{Pmoult_par} is estimated
#'   independently (own intercept + slope per season) -- appropriate now
#'   that identification comes from growmod's internal moment-matched
#'   mixture across every animal's \code{goodts} opportunities, not just
#'   literal single-opportunity records, so per-season identification is
#'   generally more available than it used to be. If \code{TRUE}, all
#'   \code{goodts} rows are mapped to the SAME two factor levels (a single
#'   shared logistic across seasons) -- reproduces the old pre-matrix
#'   behaviour exactly, and is a reasonable fallback if a diagnostic run
#'   shows a particular season's row not moving from its initial value /
#'   producing enormous standard errors, i.e. still under-identified even
#'   with the broader mixture. Non-\code{goodts} rows are always fixed
#'   regardless of this argument, since \code{growmod} never evaluates them.
#' @param Growth_shared Logical (default \code{FALSE}). Same idea as
#'   \code{Pmoult_shared}, applied to the 5-column \code{Growth_par}
#'   double-logistic growth-at-length curve (see \code{\link{growmod}}
#'   and \code{\link{Makepin}}). \code{FALSE} (default) frees each
#'   \code{goodts} row's 5 parameters (\code{Amax}, \code{P1}, \code{P2},
#'   \code{P3}, \code{P5}) independently, which is generally feasible now
#'   that the growth curve is 5 parameters per row rather than
#'   \code{nlbin}. \code{TRUE} maps every \code{goodts} row to the same 5
#'   factor levels, i.e. a single shared growth curve across seasons --
#'   use this as a fallback if a particular season is poorly identified
#'   on its own (e.g. very few recaptures in that \code{goodts} window).
#'   Non-\code{goodts} rows are always fixed regardless of this argument,
#'   since \code{growmod} never evaluates them.
#'
#' @export
MakemapPar <- function(pin, re = FALSE, estTemporalGrowth = TRUE,
                    Pmoult_shared = FALSE, Growth_shared = FALSE,
                    estSuppress = TRUE, estSlope = TRUE, estSigError = FALSE) {

  ## map is built up incrementally, one parameter (or parameter group) at
  ## a time, below -- there is no longer a generic catch-all sweep (see
  ## "Change from previous version (2)" above for why).
  map <- list()

  if (max(grepl('MerrorR', names(pin))) == 1) {
    if (re == FALSE) {
      map$MerrorRel <- rep(factor(NA), length(pin$MerrorRel))
      map$MerrorRec <- rep(factor(NA), length(pin$MerrorRec))
    }
  }

  # LsigError is fixed by default. Freeing it under the OLD loose prior
  # (log(2), sd 0.5 -- roughly +/-65% within one sd) collapsed the growth
  # curve toward flat: with one release and one recapture measurement per
  # animal, that much slack let the model explain genuine growth away as
  # measurement noise.
  #
  # estSigError = TRUE frees it. That is only safe when datain carries a
  # data-derived prior from add_sigError(), whose width is the sampling SE
  # of the estimate (typically ~0.05 on the log scale, i.e. +/-5%) -- far
  # too tight for the collapse mechanism to operate. The check below
  # refuses to free it against the loose fallback prior.
  if (estSigError) {
    psd <- datain$LsigError_prior_sd
    if (is.null(psd) || psd > 0.2) {
      stop("estSigError = TRUE needs a data-derived prior on LsigError. ",
           "Build datain with Makedata(..., LsigError = NULL) so that ",
           "add_sigError() estimates it, or leave estSigError = FALSE. ",
           "Freeing LsigError against a loose prior collapses the growth curve.")
    }
  } else {
    map$LsigError <- factor(NA)
  }

  if (re == FALSE) {
    map$LMerrorRelsigma <- factor(NA)
    map$LMerrorRecsigma <- factor(NA)
  }

  if ('Sraw' %in% names(pin)) {
    if (estTemporalGrowth == FALSE) {
      map$Sraw <- rep(factor(NA), length(pin$Sraw))
    }
  }

  ## --- Pmoult_par mapping -------------------------------------------------
  ## pin$Pmoult_par is ntsteps x 2 (see Makepin). Non-goodts rows: fixed at
  ## initial values (NA), same reasoning as Growth_par below. goodts
  ## rows: independent factor levels by default, or all sharing one pair of
  ## levels if Pmoult_shared = TRUE.
  ntsteps <- datain$ntsteps
  Pmoult_map <- matrix(NA_integer_, nrow = ntsteps, ncol = 2)

  if (Pmoult_shared) {
    # every goodts row -> the same two levels (1 = intercept, 2 = slope)
    Pmoult_map[datain$goodts, 1] <- 1
    Pmoult_map[datain$goodts, 2] <- 2
  } else {
    free_idx <- 0
    for (ns in datain$goodts) {
      for (cc in 1:2) {
        free_idx <- free_idx + 1
        Pmoult_map[ns, cc] <- free_idx
      }
    }
  }
  # as.vector() on a matrix is column-major, matching how MakeADFun flattens
  # a matrix-valued parameter -- must stay column-major for map and pin to
  # correspond to the same elements.
  map$Pmoult_par <- as.factor(as.vector(Pmoult_map))

  ## --- Growth_par mapping (new) -------------------------------------------
  ## pin$Growth_par is ntsteps x 5 (see Makepin): log(Amax), P2, log(P1),
  ## log(P3), log(P5). Same convention as Pmoult_par immediately above:
  ## non-goodts rows fixed at their Makepin initial values (growmod's
  ## growth-at-length block never reads them), goodts rows either freed
  ## independently (5 levels per season) or all sharing one set of 5
  ## levels if Growth_shared = TRUE. Independent-by-default is feasible
  ## here in a way it never was for the old growth_vecpar (nlbin
  ## parameters per row) -- 5 parameters per goodts row is a realistic ask
  ## even for a single-season, moderate-sample dataset.
  ncol_growth <- ncol(pin$Growth_par)
  Growth_map <- matrix(NA_integer_, nrow = ntsteps, ncol = ncol_growth)

  if (Growth_shared) {
    # every goodts row -> the same ncol_growth levels
    Growth_map[datain$goodts, ] <- matrix(1:ncol_growth,
                                          nrow = length(datain$goodts),
                                          ncol = ncol_growth, byrow = TRUE)
  } else {
    free_idx <- 0
    for (ns in datain$goodts) {
      for (cc in 1:ncol_growth) {
        free_idx <- free_idx + 1
        Growth_map[ns, cc] <- free_idx
      }
    }
  }
  # Column-major, same reasoning as Pmoult_par above.
  map$Growth_par <- as.factor(as.vector(Growth_map))

  ## --- mpy_split_par mapping -----------------------------------------------
  ## Only present in pin when length(datain$goodts) > 1 (see Makepin). Free
  ## by default -- the whole point is letting the data decide the split.
  ## Auto-fixed at its initial value (0, i.e. inert/even) whenever
  ## datain$mpy is 0: estimating how to split a floor of zero across
  ## seasons is meaningless and just adds an unidentified free direction.
  if ('mpy_split_par' %in% names(pin)) {
    if (is.null(datain$mpy) || datain$mpy == 0) {
      map$mpy_split_par <- rep(factor(NA), length(pin$mpy_split_par))
    }
    # else: not added to map, so estimated freely (default)
  }

  ## --- Suppression / recovery parameters -----------------------------------
  ## Free by default when present. estSuppress = FALSE fixes all three at
  ## their Makepin starting values, which is the right way to compare
  ## against a no-suppression fit. estSlope = FALSE fixes only slope_par --
  ## liberty is coarsely distributed, so the transition width is usually
  ## the least well-informed of the three; fix it if its SE comes back
  ## large or if the fit wanders to an implausible value.
  for (nm in c('r0_par', 'lib50_par', 'slope_par')) {
    if (nm %in% names(pin) && !estSuppress) map[[nm]] <- factor(NA)
  }
  if ('slope_par' %in% names(pin) && estSuppress && !estSlope) {
    map$slope_par <- factor(NA)
  }

  return(map)
}
