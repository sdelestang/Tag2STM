#' Create RTMB objective from a growmod-family function, with proper environment
#'
#' This function copies the requested model function (see \code{model}) to
#' the global environment to avoid namespace issues with RTMB's automatic
#' differentiation, and calls it with the \code{TemporalGrowth} value implied
#' by \code{pin} (single merged function per model, not a dispatch between
#' separate \code{TemporalGrowth}/non-\code{TemporalGrowth} model functions).
#'
#' @param pin Parameter list, typically created by \code{\link{Makepin}} or
#'   \code{\link{MakepinPar}} (matching whichever \code{model} you pass). If
#'   \code{pin} was built with \code{Makepin(TemporalGrowth = TRUE)} (detected
#'   via \code{attr(pin, "TemporalGrowth")}), the objective is built by
#'   calling \code{<model>(..., TemporalGrowth = TRUE)}. No separate argument
#'   is needed here — the dispatch follows \code{pin} automatically, so
#'   \code{pin} and \code{datain} must agree on whether year effects are in
#'   play (see Details).
#' @param datain Data list. Optional here: if omitted (\code{NULL}, the
#'   default), the function falls back to whatever \code{datain} object
#'   already exists in \code{.GlobalEnv} — consistent with the original
#'   "must be in calling environment" convention. If supplied, it is assigned
#'   into \code{.GlobalEnv}, overwriting any existing global \code{datain}.
#' @param map Optional parameter map, typically from \code{\link{Makemap}} or
#'   \code{\link{MakemapPar}} (matching \code{model}).
#' @param random Optional random effects.
#' @param Like Integer. Likelihood form passed through to \code{model}
#'   (1 = asymmetric KL, 2 = symmetric KL). Default 1.
#' @param model Character, default \code{'growmod'}. Name of the model
#'   function to build the objective from — e.g. \code{'growmod'} (the
#'   original \code{growth_vecpar} random-walk growth curve) or
#'   \code{'growmodPar'} (the \code{Growth_par} double-logistic growth
#'   curve, see \code{\link{growmodPar}}). Any function with the same
#'   \code{(pin, Like, TemporalGrowth)} interface as \code{growmod} works —
#'   this is a name lookup (\code{get(model, mode = "function")}), not a
#'   fixed switch between two hardcoded options, so a third growth
#'   formulation added later under a new name works here without changes to
#'   this function. Pass whichever model's \code{pin}/\code{map} you built
#'   (\code{Makepin}/\code{Makemap} for \code{'growmod'},
#'   \code{MakepinPar}/\code{MakemapPar} for \code{'growmodPar'}) — this is
#'   what makes fitting both side by side for comparison a matter of
#'   swapping this one argument (plus the matching \code{pin}/\code{map}).
#'
#' @details
#' **Dispatch logic:** \code{pin} carries a \code{TemporalGrowth} attribute
#' set by \code{\link{Makepin}}/\code{\link{MakepinPar}}. When \code{TRUE},
#' this function checks that \code{datain} carries the fields the model's
#' temporal branch needs, and that \code{length(pin$Sraw)} is consistent with
#' them — meant to catch \code{pin} and \code{datain} built from mismatched
#' \code{tdat} versions or settings. These checks are identical regardless of
#' which \code{model} is requested: \code{growmod} and \code{growmodPar} share
#' the same \code{TemporalGrowth}/\code{Sraw} machinery and differ only in how
#' the growth-at-length curve is built internally (\code{growth_vecpar} vs
#' \code{Growth_par}), which this function has no need to know about.
#'
#' Those requirements differ by mode. With \strong{annual} effects
#' (\code{datain$period_mode} absent or \code{FALSE}), the model needs
#' \code{nyears}, \code{relyr} and \code{yr_supported}, and \code{Sraw} has
#' length \code{sum(yr_supported) - 1}. With \strong{period} effects
#' (\code{period_mode = TRUE}, set by \code{\link{Makedata}} when
#' \code{period = TRUE}), it needs \code{nyears}, \code{relyr},
#' \code{year_period} and \code{nperiods} instead, and \code{Sraw} has length
#' \code{nperiods - 1}. \code{yr_supported} is deliberately NOT required in
#' period mode: support is a property of the period, not of the individual
#' year, so \code{\link{Makedata}} does not call
#' \code{\link{add_year_support}} at all in that mode and the field
#' legitimately does not exist.
#'
#' The model function named by \code{model} is looked up
#' (\code{get(model, mode = "function")}) and copied into \code{.GlobalEnv}
#' before \code{MakeADFun()} is called, since RTMB's automatic
#' differentiation requires the function's environment to resolve
#' \code{datain} and \code{pin} names directly rather than through this
#' function's local scope. \code{TemporalGrowth} itself is passed as a plain
#' (non-AD) argument to the model function, captured via closure alongside
#' \code{Like} — it is decided once when \code{MakeADFun} traces the
#' function, not re-evaluated per optimization step, so this is safe in the
#' same way \code{Like} always was.
#'
#' @export
make_growmod_obj <- function(pin, datain = NULL, map = list(), random = NULL,
                              Like = 1, model = 'growmod') {
  if (!is.null(datain)) {
    assign("datain", datain, envir = .GlobalEnv)
  } else if (exists("datain", envir = .GlobalEnv)) {
    # No datain argument supplied — fall back to whatever is already in the
    # global environment, consistent with the original documented usage
    # ("datain must be in calling environment").
    datain <- get("datain", envir = .GlobalEnv)
  }

  if (is.null(attr(pin, "TemporalGrowth"))) {
    warning("attr(pin, 'TemporalGrowth') is missing (not FALSE, genuinely absent) -- ",
            "defaulting to TemporalGrowth = FALSE. This commonly happens if pin was ",
            "reconstructed via mod$env$parList(), which does not preserve custom ",
            "attributes. If you intended to fit with year effects, rebuild pin with ",
            "Makepin(TemporalGrowth = TRUE) / MakepinPar(TemporalGrowth = TRUE) ",
            "(optionally using parList()'s values as new starting values via ",
            "pin$Sraw <- ..., pin$Pmoult_par <- ..., pin$Growth_par <- ..., etc., ",
            "rather than replacing pin wholesale).")
  }
  useVar <- isTRUE(attr(pin, "TemporalGrowth"))

  if (useVar) {
    if (is.null(datain)) {
      stop("pin was built with TemporalGrowth = TRUE, but no datain is ",
           "available (neither supplied nor present in .GlobalEnv).")
    }
    if (is.null(pin$Sraw)) {
      stop("pin has attr(pin, 'TemporalGrowth') == TRUE but no Sraw element. ",
           "Rebuild pin with Makepin(TemporalGrowth = TRUE) or ",
           "MakepinPar(TemporalGrowth = TRUE), matching the model you intend to fit.")
    }

    period_mode <- isTRUE(datain$period_mode)

    # Required fields differ by mode -- see @details. yr_supported is not a
    # period-mode concept and Makedata does not create it there. Identical
    # for every model in the growmod family (see @details).
    need <- if (period_mode) {
      c("nyears", "relyr", "year_period", "nperiods")
    } else {
      c("nyears", "relyr", "yr_supported")
    }
    miss <- need[vapply(need, function(nm) is.null(datain[[nm]]), logical(1))]
    if (length(miss)) {
      stop("pin was built with TemporalGrowth = TRUE and datain is in ",
           if (period_mode) "PERIOD" else "ANNUAL", " mode, but datain is missing: ",
           paste(sQuote(miss), collapse = ", "), ". Rebuild datain with Makedata(",
           if (period_mode) "..., period = TRUE)" else "..., TemporalGrowth = TRUE)",
           ".")
    }

    # Sraw length: nperiods - 1 in period mode, sum(yr_supported) - 1 otherwise.
    if (period_mode) {
      if (datain$nperiods != length(pin$Sraw) + 1) {
        stop("datain$nperiods (", datain$nperiods, ") does not match ",
             "length(pin$Sraw) + 1 (", length(pin$Sraw) + 1, "). pin and datain ",
             "are out of sync -- rebuild pin using this exact datain. ",
             "Note that switching between period = TRUE and FALSE, or changing the ",
             "number of periods, changes the shape of Sraw, so pin, map and the ",
             "model object must all be rebuilt.")
      }
    } else {
      n_supported <- sum(datain$yr_supported)
      if (n_supported != length(pin$Sraw) + 1) {
        stop("sum(datain$yr_supported) (", n_supported, ") does not match length(pin$Sraw) + 1 (",
             length(pin$Sraw) + 1, "). pin and datain appear to be out of sync — ",
             "rebuild pin with TemporalGrowth = TRUE using this exact datain ",
             "(after add_year_support()), or re-run add_year_support() if datain changed.")
      }
    }
  }

  # Look up the requested model function by name -- a plain name lookup
  # rather than a hardcoded growmod/growmodPar switch, so a third growth
  # formulation added later under a new name (following the same
  # (pin, Like, TemporalGrowth) interface) works here with no change to this
  # function. mode = "function" so a same-named non-function object
  # elsewhere on the search path can't be picked up by mistake.
  model_fn <- tryCatch(
    get(model, mode = "function"),
    error = function(e) {
      stop("model = '", model, "' was not found as a function. Expected one of ",
           "the growmod-family functions (e.g. 'growmod', 'growmodPar') to already ",
           "be defined/sourced with the same (pin, Like, TemporalGrowth) interface ",
           "as growmod. Check spelling, and that the relevant script has been sourced.")
    }
  )

  growmod_local <- model_fn
  environment(growmod_local) <- .GlobalEnv

  # Explicitly qualified: if TMB is also loaded (e.g. for tmbprofile()),
  # TMB::MakeADFun can shadow RTMB::MakeADFun on the search path, and the two
  # have incompatible signatures (TMB's is MakeADFun(data, parameters, ...),
  # RTMB's is MakeADFun(func, parameters, ...)) -- an unqualified call would
  # then fail with "argument 'data' is missing" regardless of load order.
  obj <- RTMB::MakeADFun(
    func = function(p) growmod_local(p, Like = Like, TemporalGrowth = useVar),  # Like/TemporalGrowth captured from enclosing scope
    parameters = pin,
    map = map,
    random = random,
    silent = FALSE
  )
  return(obj)
}
