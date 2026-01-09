#' Create RTMB objective from growmod with proper environment
#'
#' This function copies growmod to the global environment to avoid
#' namespace issues with RTMB's automatic differentiation.
#'
#' @param pin Parameter list
#' @param datain Data list (must be in calling environment)
#' @param map Optional parameter map
#' @param random Optional random effects
#' @export
make_growmod_obj <- function(pin, datain = NULL, map = list(), random = NULL) {

  # If datain provided, put it in global env
  if(!is.null(datain)) {
    assign("datain", datain, envir = .GlobalEnv)
  }

  # Copy growmod to global environment to avoid package namespace issues
  growmod_local <- growmod
  environment(growmod_local) <- .GlobalEnv

  # Create objective
  obj <- MakeADFun(
    func = function(p) growmod_local(p),
    parameters = pin,
    map = map,
    random = random,
    silent = FALSE
  )

  return(obj)
}
