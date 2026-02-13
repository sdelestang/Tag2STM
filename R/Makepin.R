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
#' @param LsigGrow Numeric. Log of the standard deviation for growth variability used
#'   in constructing the size transition matrix (STM). Controls how growth probability
#'   spreads across adjacent length bins - smaller values create sharper transitions,
#'   larger values create more diffuse growth distributions. Default is log(0.15) mm.
#' @param LMerrorRelsigma Numeric. Log of the standard deviation for individual-level
#'   random effects on measurement error at release. When non-zero, each tagged animal
#'   gets its own release measurement error drawn from N(0, exp(LMerrorRelsigma)).
#'   A value of 0 indicates no individual-level variation beyond the baseline
#'   \code{LsigError}. Default is 0.
#' @param LMerrorRecsigma Numeric. Log of the standard deviation for individual-level
#'   random effects on measurement error at recapture. When non-zero, each recapture
#'   gets its own measurement error drawn from N(0, exp(LMerrorRecsigma)).
#'   A value of 0 indicates no individual-level variation. Default is 0.
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
#' }
#'
#' @details
#' This function requires that \code{lbin}, \code{ntsteps}, and \code{tdat}
#' exist in the calling environment (typically loaded as package data or defined
#' in the global environment before model fitting).
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
#' **Parameter Scale:**
#' All parameters involving standard deviations are on the log scale to ensure
#' positivity during unconstrained optimization and to improve numerical stability.
#'
#' @section Model Components:
#' The pin file feeds into \code{growmod()}, which:
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
#' # Fit the model
#' obj <- MakeADFun(growmod, pin_default, random = c("MerrorRel", "MerrorRec"))
#' opt <- nlminb(obj$par, obj$fn, obj$gr)
#' }
#'
#' @seealso \code{\link{growmod}} for the main model function that uses this pin file
#'
#' @export
Makepin <- function(avgrowth = 2,
                    LsigError = log(2),
                    LsigGrow = rep(log(1.8),2),
                    LMerrorRelsigma = 0,
                    LMerrorRecsigma = 0) {
  nlbin <- length(bins$lbin)
  pin <- list(
    growth_vecpar = rep(c(rep(log(0.2), nlbin - 1), log(0.01)), ntsteps),
    LsigError = LsigError,
    LsigGrow = LsigGrow,
    MerrorRel = rep(0, nrow(tdat)),
    LMerrorRelsigma = LMerrorRelsigma,
    MerrorRec = rep(0, nrow(tdat)),
    LMerrorRecsigma = LMerrorRecsigma
  )
  return(pin)
}
