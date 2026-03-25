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
make_growmod_obj <- function(pin, datain = NULL, map = list(), random = NULL, Like = 1) {
  if(!is.null(datain)) {
    assign("datain", datain, envir = .GlobalEnv)
  }
  growmod_local <- growmod
  environment(growmod_local) <- .GlobalEnv

  obj <- MakeADFun(
    func = function(p) growmod_local(p, Like = Like),  # Like captured from enclosing scope
    parameters = pin,
    map = map,
    random = random,
    silent = FALSE
  )
  return(obj)
}
