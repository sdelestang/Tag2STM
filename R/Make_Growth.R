#' Create growth objective function with data baked in
#'
#' @param datain List containing model data
#' @return Function that takes parameters and returns negative log-likelihood
#' @export
growth_objective <- function(datain) {
  # Validate datain
  required <- c("lbin", "nlbin", "ntsteps", "goodts", "Tlbin",
                "Rlcl", "Rccl", "relts", "tsteps", "nlob", "lbinU", "lbinL")
  missing <- setdiff(required, names(datain))
  if(length(missing) > 0) {
    stop("Missing required data elements: ", paste(missing, collapse = ", "))
  }

  function(pin) growmod(pin, datain)
}

#' Fit growth model (complete workflow)
#'
#' @param parameters Initial parameters list
#' @param datain Model data list
#' @param random Character vector of random effects or NULL for no random effects (default: NULL)
#' @param map Parameter map list for fixing parameters (default: empty list)
#' @param control Optimization control list passed to nlminb (default: empty list)
#' @param silent Suppress RTMB output? (default: FALSE)
#' @return List containing obj (RTMB object), opt (optimization result),
#'   sdr (sdreport if successful), and converged (logical)
#' @examples
#' \dontrun{
#' # Without random effects
#' fit <- fit_growth(parameters, datain)
#'
#' # With random effects
#' fit <- fit_growth(parameters, datain, random = c("MerrorRel", "MerrorRec"))
#'
#' # With parameter mapping
#' map <- list(sigGrow = factor(NA))  # Fix sigGrow
#' fit <- fit_growth(parameters, datain, map = map)
#' }
#' @export
fit_growth <- function(parameters, datain,
                       random = NULL,  # Default: no random effects
                       map = list(),
                       control = list(),
                       silent = FALSE) {

  # Create objective
  obj <- MakeADFun(
    func = growth_objective(datain),
    parameters = parameters,
    random = random,  # Will be NULL by default
    map = map,
    silent = silent
  )

  # Optimize
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = control)

  # Standard errors
  sdr <- try(sdreport(obj), silent = TRUE)

  return(list(
    obj = obj,
    opt = opt,
    sdr = if(!inherits(sdr, "try-error")) sdr else NULL,
    converged = opt$convergence == 0
  ))
}

#' Create RTMB objective object for growth model
#'
#' Lower-level function that just creates the MakeADFun object without optimizing.
#' Useful when you want more control over the optimization process.
#'
#' @param parameters Initial parameters list
#' @param datain Model data list
#' @param random Character vector of random effects or NULL (default: NULL)
#' @param map Parameter map list (default: empty list)
#' @param silent Suppress RTMB output? (default: FALSE)
#' @return RTMB ADFun object
#' @examples
#' \dontrun{
#' # Create object
#' obj <- make_growth_obj(parameters, datain)
#'
#' # Optimize with custom settings
#' opt <- nlminb(obj$par, obj$fn, obj$gr,
#'               control = list(eval.max = 2000))
#' }
#' @export
make_growth_obj <- function(parameters, datain,
                            random = NULL,  # Default: no random effects
                            map = list(),
                            silent = FALSE) {

  obj <- MakeADFun(
    func = growth_objective(datain),
    parameters = parameters,
    random = random,
    map = map,
    silent = silent
  )

  return(obj)
}
