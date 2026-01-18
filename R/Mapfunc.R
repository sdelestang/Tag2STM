#' Create Parameter Mapping for Selective Time Step Fitting
#'
#' Generates a mapping list that controls which parameters are estimated in the
#' tag-recapture growth model. This allows growth parameters to be selectively
#' turned on or off based on which time steps have sufficient data (defined in
#' \code{goodts}), while controlling whether measurement error random effects
#' are estimated.
#'
#' @param pin List. A pin file created by \code{\link{Makepin}} containing initial
#'   parameter values and structure.
#' @param re Logical. Should individual-level measurement error random effects
#'   (\code{MerrorRel} and \code{MerrorRec}) be estimated? If \code{FALSE} (default),
#'   these parameters are fixed at their initial values (typically 0). If \code{TRUE},
#'   they are removed from the map and estimated as random effects. Default is \code{FALSE}.
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
#' By default, \code{LsigError}, \code{LMerrorRelsigma}, and \code{LMerrorRecsigma}
#' are FIXED at their initial values. This is appropriate when:
#' \itemize{
#'   \item You want to use expert-specified measurement error values
#'   \item You're conducting sensitivity analyses with different fixed error levels
#'   \item The model has difficulty estimating these simultaneously with other parameters
#' }
#'
#' To estimate these parameters, remove them from the map before calling \code{MakeADFun()}:
#' \preformatted{
#' map$LsigError <- NULL
#' map$LMerrorRelsigma <- NULL
#' map$LMerrorRecsigma <- NULL
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
#' # Example 1: Basic model with fixed sigma, no random effects
#' pin <- Makepin()
#' map <- Makemap(pin, re = FALSE)
#' obj <- MakeADFun(growmod, pin, map = map)
#'
#' # Example 2: Estimate measurement error sigma
#' pin <- Makepin()
#' map <- Makemap(pin, re = FALSE)
#' map$LsigError <- NULL  # Remove from map to estimate
#' obj <- MakeADFun(growmod, pin, map = map)
#'
#' # Example 3: Model with individual measurement error random effects
#' pin <- Makepin(LMerrorRelsigma = log(0.5), LMerrorRecsigma = log(0.5))
#' map <- Makemap(pin, re = TRUE)
#' map$LMerrorRelsigma <- NULL  # Estimate the sigma for random effects
#' map$LMerrorRecsigma <- NULL
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
#' }
#'
#' @seealso
#' \code{\link{Makepin}} for creating the initial parameter list
#' \code{\link{growmod}} for the model function
#'
#' @export
Makemap <- function(pin, re = FALSE) {
  # Identify parameter names to map (exclude sigma and error terms)
  pnames <- names(pin)[!(grepl('sig', names(pin), ignore.case = TRUE) |
                           grepl('error', names(pin), ignore.case = TRUE))]

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

  # Fix sigma parameters at initial values (remove from map to estimate)
  map$LsigError <- factor(NA)
  map$LMerrorRelsigma <- factor(NA)
  map$LMerrorRecsigma <- factor(NA)

  return(map)
}
