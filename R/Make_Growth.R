#' Create RTMB objective from growmod with proper environment
#'
#' This function copies \code{growmod} to the global environment to avoid
#' namespace issues with RTMB's automatic differentiation, and calls it with
#' the \code{TemporalGrowth} value implied by \code{pin} (single merged
#' function, not a dispatch between two separate model functions).
#'
#' @param pin Parameter list, typically created by \code{\link{Makepin}}. If
#'   \code{pin} was built with \code{Makepin(TemporalGrowth = TRUE)} (detected
#'   via \code{attr(pin, "TemporalGrowth")}), the objective is built by
#'   calling \code{growmod(..., TemporalGrowth = TRUE)}. No separate argument
#'   is needed here — the dispatch follows \code{pin} automatically, so
#'   \code{pin} and \code{datain} must agree on whether year effects are in
#'   play (see Details).
#' @param datain Data list. Optional here: if omitted (\code{NULL}, the
#'   default), the function falls back to whatever \code{datain} object
#'   already exists in \code{.GlobalEnv} — consistent with the original
#'   "must be in calling environment" convention. If supplied, it is assigned
#'   into \code{.GlobalEnv}, overwriting any existing global \code{datain}.
#' @param map Optional parameter map, typically from \code{\link{Makemap}}.
#' @param random Optional random effects.
#' @param Like Integer. Likelihood form passed through to \code{growmod}
#'   (1 = asymmetric KL, 2 = symmetric KL). Default 1.
#'
#' @details
#' **Dispatch logic:** \code{pin} carries a \code{TemporalGrowth} attribute
#' set by \code{\link{Makepin}}. When \code{TRUE}, this function checks that
#' \code{datain} carries the fields \code{growmod}'s temporal branch needs,
#' and that \code{length(pin$Sraw)} is consistent with them — meant to catch
#' \code{pin} and \code{datain} built from mismatched \code{tdat} versions or
#' settings.
#'
#' Those requirements differ by mode. With \strong{annual} effects
#' (\code{datain$period_mode} absent or \code{FALSE}), \code{growmod} needs
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
#' The model function is copied into \code{.GlobalEnv} before
#' \code{MakeADFun()} is called, since RTMB's automatic differentiation
#' requires the function's environment to resolve \code{datain} and
#' \code{pin} names directly rather than through this function's local scope.
#' \code{TemporalGrowth} itself is passed as a plain (non-AD) argument to
#' \code{growmod}, captured via closure alongside \code{Like} — it is decided
#' once when \code{MakeADFun} traces the function, not re-evaluated per
#' optimization step, so this is safe in the same way \code{Like} always was.
#'
#' @export
make_growmod_obj <- function(pin, datain = NULL, map = list(), random = NULL, Like = 1) {
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
            "Makepin(TemporalGrowth = TRUE) (optionally using parList()'s values as ",
            "new starting values via pin$Sraw <- ..., pin$Pmoult_par <- ..., etc., ",
            "rather than replacing pin wholesale).")
  }
  useVar <- isTRUE(attr(pin, "TemporalGrowth"))

  if (useVar) {
    if (is.null(datain)) {
      stop("pin was built with Makepin(TemporalGrowth = TRUE), but no datain is ",
           "available (neither supplied nor present in .GlobalEnv).")
    }
    if (is.null(pin$Sraw)) {
      stop("pin has attr(pin, 'TemporalGrowth') == TRUE but no Sraw element. ",
           "Rebuild pin with Makepin(TemporalGrowth = TRUE).")
    }

    period_mode <- isTRUE(datain$period_mode)

    # Required fields differ by mode -- see @details. yr_supported is not a
    # period-mode concept and Makedata does not create it there.
    need <- if (period_mode) {
      c("nyears", "relyr", "year_period", "nperiods")
    } else {
      c("nyears", "relyr", "yr_supported")
    }
    miss <- need[vapply(need, function(nm) is.null(datain[[nm]]), logical(1))]
    if (length(miss)) {
      stop("pin was built with Makepin(TemporalGrowth = TRUE) and datain is in ",
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
             "are out of sync -- rebuild pin with Makepin() using this exact datain. ",
             "Note that switching between period = TRUE and FALSE, or changing the ",
             "number of periods, changes the shape of Sraw, so pin, map and the ",
             "model object must all be rebuilt.")
      }
    } else {
      n_supported <- sum(datain$yr_supported)
      if (n_supported != length(pin$Sraw) + 1) {
        stop("sum(datain$yr_supported) (", n_supported, ") does not match length(pin$Sraw) + 1 (",
             length(pin$Sraw) + 1, "). pin and datain appear to be out of sync — ",
             "rebuild pin with Makepin(TemporalGrowth = TRUE) using this exact datain ",
             "(after add_year_support()), or re-run add_year_support() if datain changed.")
      }
    }
  }

  growmod_local <- growmod
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
