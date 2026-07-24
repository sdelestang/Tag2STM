#' Create Initial Parameter List for Tag-Recapture Growth Model
#'
#' Generates a list of initial parameter values (pin file) for fitting tag-recapture
#' growth models using RTMB/TMB. The function initializes growth parameters that define
#' a flexible, data-driven growth curve across length bins and time steps, along with
#' measurement error components.
#'
#' @param avgrowth Numeric. Originally intended to represent expected growth over one
#'   time step, but currently not implemented in the model. Retained for backward
#'   compatibility. Default is 2.
#' @param LsigError Numeric. Log of the standard deviation for measurement error.
#'   This represents the baseline uncertainty in length measurements at both release
#'   and recapture. Applied symmetrically in the likelihood calculation via truncated
#'   normal distributions across length bins. Default is log(2) mm.
#' @param LsigGrow Numeric scalar. Log coefficient of variation on growth
#'   increments, used by \code{growmod} (both \code{TemporalGrowth} branches)
#'   when constructing the size transition matrix (STM):
#'   \code{sd_growth = exp(LsigGrow) * mn_growth} -- i.e. growth spread is
#'   directly proportional to mean growth increment at each size/season,
#'   tapering to near-zero as growth itself tapers toward zero (e.g. near
#'   maximum size). This replaced an earlier saturating-ceiling formula that
#'   could not reproduce the size-dependent spread seen in raw single-moult
#'   recapture data (spread staying flat/near-maximal even for large animals
#'   with tiny mean increments, when the data showed spread shrinking
#'   alongside the mean). Default is \code{log(2)}, i.e. sd_growth starts at
#'   roughly 2x the mean increment, matching what raw single-moult data in
#'   this package's typical use case implied.
#' @param LMerrorRelsigma Numeric. Log of the standard deviation for individual-level
#'   random effects on measurement error at release. When non-zero, each tagged animal
#'   gets its own release measurement error drawn from N(0, exp(LMerrorRelsigma)).
#'   A value of 0 indicates no individual-level variation beyond the baseline
#'   \code{LsigError}. Default is 0.
#' @param LMerrorRecsigma Numeric. Log of the standard deviation for individual-level
#'   random effects on measurement error at recapture. When non-zero, each recapture
#'   gets its own measurement error drawn from N(0, exp(LMerrorRecsigma)).
#'   A value of 0 indicates no individual-level variation. Default is 0.
#' @param TemporalGrowth Logical. If \code{TRUE}, initializes the additional \code{Sraw}
#'   parameter vector needed for year-specific proportional growth scaling
#'   (see \code{\link{growmod}}), and marks the returned \code{pin} so that
#'   downstream functions (e.g. \code{make_growmod_obj}) know to fit
#'   \code{growmod(..., TemporalGrowth = TRUE)}. \code{Sraw} is sized from the number of years with adequate
#'   recapture support, \code{sum(datain$yr_supported) - 1} — NOT from raw
#'   \code{nyears} — so \code{datain <- add_year_support(datain)} must be run
#'   first (see \code{\link{add_year_support}}); \code{Makepin} errors otherwise.
#'   Default is \code{FALSE}, i.e. the standard no-interannual-variation model.
#'
#' @return A named list containing initial parameter values for model fitting:
#' \describe{
#'   \item{growth_vecpar}{Numeric vector of length \code{nlbin * ntsteps}. Log-scale
#'     increment parameters that define growth curves. These are transformed into
#'     cumulative growth by length bin within the model: starting from the largest
#'     bin (initialized to log(0.01), near-zero growth), each smaller bin adds
#'     exp(parameter) to create an increasing growth curve. Most values initialized
#'     to log(0.2) mm to represent typical increment sizes. The random walk structure
#'     allows flexible, data-driven growth patterns.}
#'   \item{LsigError}{Scalar. Log standard deviation for baseline measurement error.}
#'   \item{LsigGrow}{Scalar. Log standard deviation for growth dispersion in STM.}
#'   \item{MerrorRel}{Numeric vector of length \code{nrow(tdat)}. Individual random
#'     effects for release measurement error, one per tagged animal. Initialized to 0
#'     and estimated during model fitting if \code{LMerrorRelsigma > 0}.}
#'   \item{LMerrorRelsigma}{Scalar. Log standard deviation governing the distribution
#'     of \code{MerrorRel}.}
#'   \item{MerrorRec}{Numeric vector of length \code{nrow(tdat)}. Individual random
#'     effects for recapture measurement error, one per recapture event. Initialized
#'     to 0 and estimated during model fitting if \code{LMerrorRecsigma > 0}.}
#'   \item{LMerrorRecsigma}{Scalar. Log standard deviation governing the distribution
#'     of \code{MerrorRec}.}
#'   \item{Sraw}{Only present when \code{TemporalGrowth = TRUE}. Numeric vector of length
#'     \code{sum(datain$yr_supported) - 1} (see \code{\link{add_year_support}}),
#'     initialized to all zeros (average growth, no year effects). Inside
#'     \code{growmod}, these values are placed into the *supported* year
#'     slots of the full \code{nyears}-length year-effect vector \code{S} as
#'     \code{c(Sraw, -sum(Sraw))}, enforcing a sum-to-zero constraint among
#'     supported years only; unsupported years are fixed at exactly \code{S = 0}
#'     (average growth) rather than participating in that constraint.}
#'   \item{Pmoult_par}{ALWAYS present (unconditional on \code{TemporalGrowth}
#'     -- both branches of \code{growmod} use the moult hurdle). Numeric
#'     vector of length 2 (intercept, slope), the logistic P(moult)-by-size
#'     coefficients: \code{Pmoult(fm) = plogis(Pmoult_par[1] +
#'     Pmoult_par[2] * lbin[fm])}. Default \code{c(0.756, -0.0120)} is
#'     derived directly from observed near-zero-growth fractions in
#'     single-timestep recapture data, not an arbitrary guess -- see
#'     \code{\link{growmod}} for the full rationale.}
#' }
#'
#' @details
#' This function requires that \code{lbin}, \code{ntsteps}, and \code{tdat}
#' exist in the calling environment (typically loaded as package data or defined
#' in the global environment before model fitting). When \code{TemporalGrowth = TRUE},
#' \code{tdat} must additionally include \code{relyr} and \code{recyr} columns
#' (release and recapture year for each tagged animal) so that \code{nyears} can
#' be derived.
#'
#' **Growth Model Structure:**
#' The \code{growth_vecpar} parameters use a random walk structure where:
#' \itemize{
#'   \item Parameters are reshaped into a matrix with \code{nlbin} columns (length bins)
#'     and \code{ntsteps} rows (time periods)
#'   \item Within each time step, parameters are converted to cumulative growth starting
#'     from the largest length bin
#'   \item The model estimates \code{nlbin} parameters per time step, allowing growth
#'     to vary flexibly with both size and season
#' }
#'
#' **Measurement Error Structure:**
#' \itemize{
#'   \item \code{LsigError}: Applied to all measurements, represents instrument precision
#'     and rounding error
#'   \item \code{MerrorRel/Rec}: Individual-level deviations, useful for capturing
#'     systematic differences between taggers, measuring instruments, or individual
#'     animal characteristics
#' }
#'
#' **Interannual Growth Variation (\code{TemporalGrowth}):**
#' When enabled, \code{Sraw} adds year-specific proportional scaling
#' (\code{exp(S_t)}) of the mean growth increment (and, since sd_growth is
#' proportional to mean growth, the growth variance too) at every length bin
#' and season. See \code{\link{growmod}} for the full model. The
#' \code{TemporalGrowth} flag itself is not part of the numeric parameter
#' list passed to \code{MakeADFun} — it is stored as an attribute on the
#' returned \code{pin} object (\code{attr(pin, "TemporalGrowth")}) so that
#' downstream functions such as \code{make_growmod_obj} can detect it and
#' pass the matching \code{TemporalGrowth} value through to \code{growmod()},
#' without it being treated as an estimable parameter.
#'
#' **Parameter Scale:**
#' All parameters involving standard deviations are on the log scale to ensure
#' positivity during unconstrained optimization and to improve numerical stability.
#'
#' @section Model Components:
#' The pin file feeds into \code{growmod(pin, TemporalGrowth = ...)}, which:
#' \enumerate{
#'   \item Builds size transition matrices (STMs) from growth parameters
#'   \item Projects each tagged animal forward through time using STMs
#'   \item Compares projected size distributions with observed recapture sizes
#'   \item Maximizes multinomial-like likelihood weighted by number of lobsters
#'     in each release cohort
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming lbin, ntsteps, and tdat are loaded:
#'
#' # Create pin file with default values
#' pin_default <- Makepin()
#'
#' # Create pin file with tighter measurement error (more precise measurements)
#' pin_precise <- Makepin(LsigError = log(1.0))
#'
#' # Create pin file with more diffuse growth (accounts for high variability)
#' pin_variable <- Makepin(LsigGrow = log(0.3))
#'
#' # Create pin file with individual-level measurement error variation
#' # (useful when multiple taggers or instruments are involved)
#' pin_random <- Makepin(
#'   LMerrorRelsigma = log(0.5),  # 0.5mm SD between individuals at release
#'   LMerrorRecsigma = log(0.5)   # 0.5mm SD between individuals at recapture
#' )
#'
#' # Create pin file for the year-specific growth scaling model
#' # (requires tdat$relyr and tdat$recyr to be present)
#' pin_var <- Makepin(TemporalGrowth = TRUE)
#'
#' # Fit the model
#' obj <- MakeADFun(growmod, pin_default, random = c("MerrorRel", "MerrorRec"))
#' opt <- nlminb(obj$par, obj$fn, obj$gr)
#' }
#'
#' @seealso \code{\link{growmod}}, the model function that uses this pin file
#'   -- pass \code{TemporalGrowth = TRUE}/\code{FALSE} to match how \code{pin}
#'   was built here.
#'
#' @export
Makepin <- function(avgrowth = 2,
                    LsigError = log(2),
                    LsigGrow = log(2),
                    LMerrorRelsigma = 0,
                    LMerrorRecsigma = 0,
                    TemporalGrowth = FALSE) {
  nlbin <- length(bins$lbin)
  pin <- list(
    growth_vecpar = rep(c(rep(log(0.2), nlbin - 1), log(0.01)), ntsteps),
    LsigError = LsigError,
    LsigGrow = LsigGrow,
    MerrorRel = rep(0, nrow(tdat)),
    LMerrorRelsigma = LMerrorRelsigma,
    MerrorRec = rep(0, nrow(tdat)),
    LMerrorRecsigma = LMerrorRecsigma,
    # Starting values derived directly from the observed near-zero-growth
    # fraction in single-timestep recaptures (~0.41 moult probability at the
    # small end of the size range, ~0.25 at the large end) via
    # plogis(a + b*lbin) = target_Pmoult solved at two representative sizes
    # (rel ~= 93 and rel ~= 154) -- not an arbitrary guess. Unconditional:
    # both TemporalGrowth branches of growmod() now use the moult hurdle.
    Pmoult_par = c(0.756, -0.0120)
  )

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
