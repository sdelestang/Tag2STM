## Change from previous version:
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

## Inside Makemap(), the `pnames` exclusion of Pmoult_par from the turnon
## sweep stays as before (it needs its own mapping logic, not growth_vecpar's).
## Add this block after the existing Sraw block, at the end of the function,
## before `return(map)`:

## -----------------------------------------------------------------------
## Full function
## -----------------------------------------------------------------------

#' Create Parameter Mapping for Selective Time Step Fitting
#'
#' (See previous version's roxygen for the general \code{re}/
#' \code{estTemporalGrowth} behaviour -- unchanged. This version adds
#' explicit handling for \code{Pmoult_par}, now an \code{ntsteps x 2}
#' matrix rather than a length-2 vector -- see \code{\link{Makepin}}.)
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
#'
#' @export
Makemap <- function(pin, re = FALSE, estTemporalGrowth = TRUE,
                     Pmoult_shared = FALSE, estSuppress = TRUE) {
  pnames <- names(pin)[!(grepl('sig', names(pin), ignore.case = TRUE) |
                           grepl('error', names(pin), ignore.case = TRUE) |
                           names(pin) == 'Sraw' |
                           names(pin) == 'Pmoult_par' |
                           names(pin) == 'mpy_split_par' |
                           names(pin) == 'suppress_par' |
                           names(pin) == 'comp_par')]

  turnon <- 1:(datain$ntsteps * datain$nlbin)
  turnon[!(rep(1:(datain$ntsteps), each = datain$nlbin)) %in% datain$goodts] <- NA
  turnon <- as.factor(turnon)

  func1 <- function(x) x = turnon
  map <- lapply(pnames, func1)
  map <- setNames(map, pnames)

  if (max(grepl('MerrorR', names(pin))) == 1) {
    if (re == FALSE) {
      map$MerrorRel <- rep(factor(NA), length(pin$MerrorRel))
      map$MerrorRec <- rep(factor(NA), length(pin$MerrorRec))
    }
  }

  map$LsigError <- factor(NA)

  if (re == FALSE) {
    map$LMerrorRelsigma <- factor(NA)
    map$LMerrorRecsigma <- factor(NA)
  }

  if ('Sraw' %in% names(pin)) {
    if (estTemporalGrowth == FALSE) {
      map$Sraw <- rep(factor(NA), length(pin$Sraw))
    }
  }

  ## --- Pmoult_par mapping (new) -------------------------------------------
  ## pin$Pmoult_par is ntsteps x 2 (see Makepin). Non-goodts rows: fixed at
  ## initial values (NA), same reasoning as growth_vecpar's turnon. goodts
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

  ## --- mpy_split_par mapping (new) -----------------------------------------
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

  ## --- Suppression parameters (new) ---------------------------------------
  ## Free by default when present -- the whole point is estimating how much
  ## of the first post-release moult is lost. Set estSuppress = FALSE to fix
  ## them at their Makepin starting values, which is the right way to run a
  ## likelihood comparison against the no-suppression case (fit both, compare
  ## Pmoult and the growth curve; the objective is not a true log-likelihood
  ## so treat any difference as indicative rather than as a test statistic).
  if ('suppress_par' %in% names(pin) && !estSuppress) {
    map$suppress_par <- factor(NA)
  }
  if ('comp_par' %in% names(pin) && !estSuppress) {
    map$comp_par <- factor(NA)
  }

  return(map)
}
