#' Create Parameter Mapping for Selective Time Step Fitting
#'
#' Generates a mapping list that controls which parameters are estimated in the
#' tag-recapture growth model. This allows growth parameters to be selectively
#' turned on or off based on which time steps have sufficient data (defined in
#' \code{goodts}), while controlling whether measurement error random effects
#' — and, when present, year-specific growth scaling (\code{Sraw}) — are estimated.
#'
#' @param pin List. A pin file created by \code{\link{Makepin}} containing initial
#'   parameter values and structure. If \code{pin} was created with
#'   \code{TemporalGrowth = TRUE}, it will contain an additional \code{Sraw} element,
#'   which this function detects automatically and maps separately from
#'   \code{growth_vecpar} (see Details).
#' @param re Logical. Should individual-level measurement error random effects
#'   (\code{MerrorRel} and \code{MerrorRec}) be estimated? If \code{FALSE} (default),
#'   these parameters are fixed at their initial values (typically 0). If \code{TRUE},
#'   they are removed from the map and estimated as random effects, AND
#'   \code{LMerrorRelsigma}/\code{LMerrorRecsigma} (the SD they are drawn
#'   from) are also freed to be estimated -- fixing that SD while letting
#'   the individual values vary would cap how much they're allowed to move
#'   at an arbitrary default, defeating the point of \code{re = TRUE}.
#'   \code{LsigError} is NOT coupled to \code{re} and stays fixed regardless
#'   (see Details for why). Default is \code{FALSE}.
#' @param estTemporalGrowth Logical. Only relevant when \code{pin} contains \code{Sraw}
#'   (i.e. was built with \code{Makepin(TemporalGrowth = TRUE)}). If \code{TRUE} (default),
#'   \code{Sraw} is left out of the map and estimated as ordinary fixed effects,
#'   one parameter per element. If \code{FALSE}, all \code{Sraw} elements are fixed
#'   at their initial value (0), which collapses \code{growmod} back to
#'   average-year-only growth — useful for a baseline fit or a likelihood-ratio
#'   comparison against the year-effect model. Ignored (with no error) when
#'   \code{pin} has no \code{Sraw} element.
#'
#' @return A named list suitable for the \code{map} argument in \code{MakeADFun()}.
#'   Each element is a factor vector indicating which parameters should be fixed
#'   (\code{NA}) versus estimated (unique integer factor levels). The list contains:
#' \describe{
#'   \item{growth_vecpar}{Factor vector with \code{NA} for time steps not in
#'     \code{goodts} (fixed at initial values), and sequential integer factor levels
#'     for time steps in \code{goodts} (estimated). This prevents estimation of
#'     growth parameters for time steps without data.}
#'   \item{MerrorRel}{Factor vector that is all \code{NA} (fixed at 0) if \code{re = FALSE}.
#'     If \code{re = TRUE}, this element is not included in the map, allowing these
#'     parameters to be estimated as random effects.}
#'   \item{MerrorRec}{Factor vector that is all \code{NA} (fixed at 0) if \code{re = FALSE}.
#'     If \code{re = TRUE}, this element is not included in the map, allowing these
#'     parameters to be estimated as random effects.}
#'   \item{LsigError}{Single \code{factor(NA)} - FIXED at initial value. Modify map
#'     to estimate if needed.}
#'   \item{LMerrorRelsigma}{Single \code{factor(NA)} - FIXED at initial value. Modify
#'     map to estimate if needed.}
#'   \item{LMerrorRecsigma}{Single \code{factor(NA)} - FIXED at initial value. Modify
#'     map to estimate if needed.}
#'   \item{Sraw}{Only present in the returned map (as \code{factor(NA)} for every
#'     element) when \code{pin} contains \code{Sraw} AND \code{estTemporalGrowth = FALSE}.
#'     When \code{pin} contains \code{Sraw} and \code{estTemporalGrowth = TRUE} (default),
#'     \code{Sraw} is deliberately left out of the map so every element is estimated
#'     independently by \code{MakeADFun()}'s default behaviour. When \code{pin} has
#'     no \code{Sraw} element at all, this list has no \code{Sraw} entry either way.}
#' }
#'
#' @details
#' **TMB Mapping Convention:**
#' In TMB/RTMB:
#' \itemize{
#'   \item \code{factor(NA)} = parameter is FIXED at its initial value (not estimated)
#'   \item \code{factor(c(1, 2, 3))} = these parameters ARE estimated
#'   \item Parameters with the same factor level are constrained to be equal
#'   \item Parameters not in the map are estimated by default
#' }
#'
#' **Purpose:**
#' This function implements a critical feature for fitting growth models with sparse
#' or irregular data. By fixing parameters for time steps not in \code{goodts},
#' the model:
#' \itemize{
#'   \item Reduces the number of parameters to estimate, improving convergence
#'   \item Prevents estimation of growth for seasons/periods without recapture data
#'   \item Allows growth to be held constant during periods of interest
#'   \item Focuses estimation power on time steps with sufficient information
#' }
#'
#' **Time Step Mapping:**
#' The function creates a factor vector for \code{growth_vecpar} where:
#' \itemize{
#'   \item Each time step occupies \code{nlbin} consecutive positions
#'   \item Time steps in \code{goodts} receive sequential integer factor levels (1, 2, 3, ...) = ESTIMATED
#'   \item Time steps NOT in \code{goodts} receive \code{NA} = FIXED at initial values
#'   \item This pattern repeats for all \code{ntsteps} time periods
#' }
#'
#' **Why \code{Sraw} is excluded from the time-step mapping:**
#' The \code{turnon} factor vector built for \code{growth_vecpar} has length
#' \code{ntsteps * nlbin} and encodes which length-bin/season combinations have
#' data. \code{Sraw} is a completely different parameter — length
#' \code{nyears - 1}, indexed by year, not by length bin or season — so applying
#' the same \code{turnon} vector to it would silently produce a length mismatch
#' (or, worse, a coincidental length match with the wrong meaning). \code{Sraw}
#' is therefore identified by name and excluded from the \code{pnames} sweep that
#' applies \code{turnon}, and handled explicitly afterward via \code{estTemporalGrowth}.
#'
#' **Random Effects Control:**
#' When \code{re = FALSE} (default):
#' \itemize{
#'   \item \code{MerrorRel} and \code{MerrorRec} are added to map as \code{factor(NA)}
#'   \item Individual measurement errors are fixed at 0 (no individual variation)
#'   \item Model is faster and uses fewer degrees of freedom
#'   \item Appropriate when measurement error is consistent across individuals
#' }
#'
#' When \code{re = TRUE}:
#' \itemize{
#'   \item \code{MerrorRel} and \code{MerrorRec} are NOT added to the map
#'   \item Each individual gets its own estimated release and recapture measurement error
#'   \item Useful when taggers differ in precision or animals have individual characteristics
#'   \item Requires \code{random = c("MerrorRel", "MerrorRec")} in \code{MakeADFun()}
#' }
#'
#' **Sigma Parameters:**
#' \code{LsigError} is ALWAYS fixed at its initial value, regardless of
#' \code{re} -- see Details above for why (freeing it collapses the growth
#' curve, since a single release/recapture measurement pair per animal can't
#' separate measurement noise from real growth). \code{LMerrorRelsigma} and
#' \code{LMerrorRecsigma} are coupled to \code{re}: fixed when \code{re =
#' FALSE} (consistent with \code{MerrorRel}/\code{MerrorRec} themselves being
#' fixed at 0 then, so their SD has nothing to act on anyway), and left free
#' to estimate automatically when \code{re = TRUE} -- no manual
#' \code{map$LMerrorRelsigma <- NULL} step needed.
#'
#' To estimate \code{LsigError} anyway (not recommended given the collapse
#' risk above, but available for sensitivity analysis), remove it from the
#' map explicitly before calling \code{MakeADFun()}:
#' \preformatted{
#' map$LsigError <- NULL
#' }
#'
#' @section Required Environment Variables:
#' The function requires \code{datain} to exist in the calling environment, containing:
#' \itemize{
#'   \item \code{ntsteps}: Total number of time steps in the model
#'   \item \code{nlbin}: Number of length bins
#'   \item \code{goodts}: Integer vector of time steps to estimate growth for
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming datain with goodts defined and pin created:
#'
#' # Example 1: Basic model, LsigError fixed, no random effects
#' pin <- Makepin()
#' map <- Makemap(pin, re = FALSE)
#' obj <- MakeADFun(growmod, pin, map = map)
#'
#' # Example 2: Estimate LsigError anyway (not recommended, see Details)
#' pin <- Makepin()
#' map <- Makemap(pin, re = FALSE)
#' map$LsigError <- NULL  # Remove from map to estimate
#' obj <- MakeADFun(growmod, pin, map = map)
#'
#' # Example 3: Individual measurement error random effects -- re = TRUE
#' # frees MerrorRel/MerrorRec AND their sigmas automatically, no extra step
#' pin <- Makepin(LMerrorRelsigma = log(0.5), LMerrorRecsigma = log(0.5))
#' map <- Makemap(pin, re = TRUE)
#' obj <- MakeADFun(
#'   growmod,
#'   pin,
#'   map = map,
#'   random = c("MerrorRel", "MerrorRec")
#' )
#'
#' # Example 4: Fit only specific seasons
#' datain$goodts <- c(1, 2, 3)  # Summer months only
#' map <- Makemap(pin, re = FALSE)
#' # Growth for time steps 4-12 fixed at initial values
#'
#' # Example 5: Check the mapping structure
#' pin <- Makepin()
#' map <- Makemap(pin, re = FALSE)
#' table(map$growth_vecpar, useNA = "ifany")  # See which parameters are estimated
#' # NA = fixed, integers = estimated
#'
#' # Example 6: Year-specific growth scaling — estimate Sraw
#' pin <- Makepin(TemporalGrowth = TRUE)
#' map <- Makemap(pin, re = FALSE)  # Sraw estimated by default (estTemporalGrowth = TRUE)
#' obj <- make_growmod_obj(pin = pin, map = map)  # dispatches to growmod(..., TemporalGrowth = TRUE)
#'
#' # Example 7: Baseline comparison — fit with year effects fixed off
#' pin <- Makepin(TemporalGrowth = TRUE)
#' map <- Makemap(pin, re = FALSE, estTemporalGrowth = FALSE)  # Sraw fixed at 0
#' obj_baseline <- make_growmod_obj(pin = pin, map = map)
#' }
#'
#' @seealso
#' \code{\link{Makepin}} for creating the initial parameter list
#' \code{\link{growmod}} for the model function (single function, with a
#'   \code{TemporalGrowth} argument, that uses the \code{Sraw} mapping
#'   described here)
#'
#' @export
Makemap <- function(pin, re = FALSE, estTemporalGrowth = TRUE) {
  # Identify parameter names to map (exclude sigma, error, Sraw, and
  # Pmoult_par terms — each is a different length/meaning to growth_vecpar
  # and either doesn't need mapping (Pmoult_par: always estimated, length 2)
  # or gets its own mapping below (Sraw), not the growth_vecpar time-step mapping)
  pnames <- names(pin)[!(grepl('sig', names(pin), ignore.case = TRUE) |
                           grepl('error', names(pin), ignore.case = TRUE) |
                           names(pin) == 'Sraw' |
                           names(pin) == 'Pmoult_par')]

  # Create factor vector for turning parameters on/off by time step
  turnon <- 1:(datain$ntsteps * datain$nlbin)
  turnon[!(rep(1:(datain$ntsteps), each = datain$nlbin)) %in% datain$goodts] <- NA
  turnon <- as.factor(turnon)

  # Apply mapping to identified parameters
  func1 <- function(x) x = turnon
  map <- lapply(pnames, func1)
  map <- setNames(map, pnames)

  # Handle random effects for individual measurement errors
  if (max(grepl('MerrorR', names(pin))) == 1) {
    if (re == FALSE) {
      # Fix random effects at initial values (typically 0)
      map$MerrorRel <- rep(factor(NA), length(pin$MerrorRel))
      map$MerrorRec <- rep(factor(NA), length(pin$MerrorRec))
    }
    # If re == TRUE, MerrorRel/Rec are NOT added to map, so they'll be estimated
  }

  # Fix LsigError at its initial value UNCONDITIONALLY, regardless of re.
  # Freeing LsigError specifically (even alongside re = TRUE) was found to
  # collapse the growth curve toward flat: with typically one release and one
  # recapture measurement per animal, the model can't distinguish "this
  # animal really grew" from "both measurements were just noisy", so a free
  # LsigError can explain away genuine growth signal as measurement error.
  map$LsigError <- factor(NA)

  # LMerrorRelsigma/LMerrorRecsigma (the SD that MerrorRel/MerrorRec are drawn
  # from), by contrast, ARE coupled to re: fixing them while re = TRUE frees
  # the individual MerrorRel/MerrorRec values is self-defeating -- it caps
  # how much individual variation is allowed to matter at a fixed default
  # (exp(0) = 1mm) regardless of what the data actually supports, which
  # silently neutered re = TRUE in practice (individual REs converged to
  # exactly their initial value, never having had room to move). If
  # re == TRUE, leave these two out of map so they estimate; if re == FALSE,
  # fix them (matches the original unconditional-fix behaviour, and is
  # harmless either way since MerrorRel/MerrorRec are themselves fixed at 0
  # when re == FALSE, so their sigma has nothing to act on regardless).
  if (re == FALSE) {
    map$LMerrorRelsigma <- factor(NA)
    map$LMerrorRecsigma <- factor(NA)
  }
  # If re == TRUE, LMerrorRelsigma/LMerrorRecsigma are NOT added to map, so
  # they estimate automatically -- no manual `map$LMerrorRelsigma <- NULL`
  # needed anymore.

  # Handle year-specific growth scaling (Sraw), only present when pin was
  # built via Makepin(TemporalGrowth = TRUE)
  if ('Sraw' %in% names(pin)) {
    if (estTemporalGrowth == FALSE) {
      # Fix all year effects at their initial value (0) -> growmod
      # collapses back to average-year-only growth
      map$Sraw <- rep(factor(NA), length(pin$Sraw))
    }
    # If estTemporalGrowth == TRUE (default), Sraw is NOT added to map, so each
    # element is estimated independently
  }

  return(map)
}
