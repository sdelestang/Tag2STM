#' Build the Model Data List for growmod
#'
#' Assembles \code{datain} from a tag-recapture data frame and length-bin
#' structure, applies all model switches in one place, validates the result,
#' and (when \code{TemporalGrowth = TRUE}) runs \code{\link{add_year_support}}
#' so it cannot be forgotten. Replaces hand-writing the \code{datain} list.
#'
#' @param tdat Data frame of tag-recapture records. Required columns:
#'   \code{tag}, \code{relyr}, \code{recyr}, \code{rccl}, \code{rlcl},
#'   \code{relts}, \code{ntstep}, \code{num}.
#' @param bins List with \code{lbin}, \code{lbinL}, \code{lbinU}.
#' @param ntsteps Integer, number of time steps per year.
#' @param goodts Integer vector, which time steps are growth/moult
#'   opportunities.
#' @param M Numeric, natural mortality (used by \code{plotfit} for the
#'   trajectory x-axis, not by \code{growmod} itself).
#' @param smoother Numeric, smoothness penalty weight on growth. Default 2.
#' @param TemporalGrowth Logical. Stored in \code{datain} so that
#'   \code{\link{Makepin}} and \code{\link{growmod}} pick up the same value
#'   rather than it being passed separately to each -- a mismatch between
#'   them was previously possible. Default \code{FALSE}.
#' @param period Logical. When \code{TRUE} (and \code{TemporalGrowth = TRUE}),
#'   growth deviations are estimated once per PERIOD rather than once per
#'   year -- appropriate when the series has too few recaptures to support
#'   annual effects, or has a gap that annual effects would otherwise fill
#'   with zeros. Default \code{FALSE}, i.e. per-year effects. The period
#'   column can stay in \code{tdat} either way, so this is a one-line switch
#'   between the two.
#' @param period_col Character, name of the record-level period column in
#'   \code{tdat}. Default \code{"period"}. Used only when
#'   \code{period = TRUE}. Record-level values are mapped up to a year-level
#'   vector covering every year, gaps filled by carrying the preceding
#'   period forward, so animals whose liberty spans a gap still get a
#'   defined effect.
#' @param LsigError Numeric or \code{NULL}. Log measurement-error SD. When
#'   \code{NULL} (default) it is ESTIMATED from the data by
#'   \code{\link{add_sigError}} -- reflection on negative increments, which
#'   can only be measurement error since moulting always increases length.
#'   Supply a value only to override that; doing so also widens the
#'   \code{PenSigError} prior to 0.5, since a supplied number is an
#'   assertion rather than an estimate with known sampling error.
#' @param mpy Numeric, minimum moults per year (biological floor on Pmoult).
#'   Default 0 (no floor). See \code{\link{growmod}}.
#' @param n_pmoult1 Integer, number of smallest length bins with Pmoult
#'   hard-fixed at exactly 1. Default 1.
#' @param ident_wt Numeric, weight on the internal identification mixture.
#'   Default 1. Set to 0 for a pure marginal fit.
#' @param use_individual_error Logical, default \code{FALSE}.
#' @param suppress Logical, estimate tagging-induced moult suppression as a
#'   recovery function of time at liberty. Default \code{FALSE}. Requires
#'   liberty in days (see \code{liberty_days}).
#' @param liberty_days Numeric vector or \code{NULL}. Time at liberty in
#'   days, one per record. If \code{NULL} (default) it is taken from the
#'   \code{lib_col} column of \code{tdat} when present. Required when
#'   \code{suppress = TRUE}; \code{tsteps} is far too coarse a substitute,
#'   since with annual timesteps everything from ~6 to ~18 months maps to 1.
#' @param lib_col Character, name of the time-at-liberty (days) column in
#'   \code{tdat}. Default \code{"Ldays"}. Ignored when \code{liberty_days}
#'   is supplied directly.
#' @param sigError_max_liberty Numeric, days. Passed to
#'   \code{\link{add_sigError}}: only records at liberty up to this
#'   contribute to the measurement-error estimate, since negatives at long
#'   liberty are more likely tag-matching errors than measurement noise.
#'   Default 730; reduce for faster-growing species. Ignored when
#'   \code{LsigError} is supplied.
#' @param min_support Integer, passed to \code{\link{add_year_support}} when
#'   \code{TemporalGrowth = TRUE}. Default 3.
#' @param quiet Logical, suppress the diagnostic summary. Default
#'   \code{FALSE}.
#'
#' @return A \code{datain} list ready for \code{\link{Makepin}},
#'   \code{\link{Makemap}} and \code{\link{growmod}}.
#'
#' @details
#' \strong{Year indexing.} \code{relyr} is rebased to a 1-based index off
#' the earliest release year, and \code{nyears} is computed as
#' \code{max(c(relyr, recyr)) - min(relyr) + 1}. Note this is NOT the same
#' as \code{length(unique(c(relyr, recyr)))}, which is what hand-written
#' versions of this list have used: that counts distinct years PRESENT, so
#' if any year within the span has neither a release nor a recapture it
#' returns a value smaller than the largest \code{relyr} index, and
#' \code{growmod}'s range check then fails. The two agree only when every
#' year in the span is populated.
#'
#' \strong{Row filtering happens BEFORE this call.} \code{Makepin} sizes the
#' \code{MerrorRel}/\code{MerrorRec} random effects from the tag data, so
#' \code{datain} and the data frame \code{Makepin} sees must contain the
#' same rows. If you want to exclude records (e.g. dropping single-
#' opportunity animals), filter \code{tdat} first and pass the filtered
#' frame here -- do not filter afterwards.
#'
#' The printed summary includes the distribution of moult opportunities per
#' animal, which is worth reading before fitting: a large single-opportunity
#' group is where tagging-induced suppression shows up (see
#' \code{suppress}), and the count of animals with zero opportunities tells
#' you how many records carry no growth information at all.
#'
#' @seealso \code{\link{growmod}}, \code{\link{Makepin}},
#'   \code{\link{Makemap}}, \code{\link{add_year_support}}
#'
#' @export
Makedata <- function(tdat, bins, ntsteps, goodts, M,
                     smoother = 2,
                     TemporalGrowth = FALSE,
                     period = FALSE,
                     period_col = "period",
                     LsigError = NULL,
                     mpy = 0,
                     n_pmoult1 = 1,
                     ident_wt = 1,
                     use_individual_error = FALSE,
                     suppress = FALSE,
                     liberty_days = NULL,
                     lib_col = "Ldays",
                     sigError_max_liberty = 730,
                     min_support = 3,
                     quiet = FALSE) {

  ## --- Validate inputs ----------------------------------------------------
  need <- c("tag", "relyr", "recyr", "rccl", "rlcl", "relts", "ntstep", "num")
  miss <- setdiff(need, names(tdat))
  if (length(miss)) {
    stop("tdat is missing required column(s): ", paste(miss, collapse = ", "))
  }
  if (!all(c("lbin", "lbinL", "lbinU") %in% names(bins))) {
    stop("bins must contain lbin, lbinL and lbinU.")
  }
  if (length(bins$lbin) != length(bins$lbinL) ||
      length(bins$lbin) != length(bins$lbinU)) {
    stop("bins$lbin, bins$lbinL and bins$lbinU must all be the same length.")
  }
  if (!all(goodts %in% seq_len(ntsteps))) {
    stop("goodts contains values outside 1..ntsteps (", ntsteps, ").")
  }
  if (length(goodts) == 0) {
    stop("goodts is empty -- no growth opportunities, nothing to estimate.")
  }
  if (mpy > 1) {
    stop("mpy = ", mpy, " exceeds 1 -- not a valid probability floor.")
  }
  if (n_pmoult1 > length(bins$lbin)) {
    stop("n_pmoult1 (", n_pmoult1, ") exceeds the number of length bins.")
  }
  if (any(tdat$relts < 1 | tdat$relts > ntsteps, na.rm = TRUE)) {
    stop("tdat$relts has values outside 1..ntsteps (", ntsteps, ").")
  }
  if (any(tdat$ntstep < 0, na.rm = TRUE)) {
    stop("tdat$ntstep has negative values.")
  }
  if (anyNA(tdat$rccl) || anyNA(tdat$rlcl)) {
    warning("tdat has NA release or recapture lengths -- these will produce ",
            "NA likelihood contributions. Filter them out before fitting.")
  }

  rng <- range(c(tdat$rlcl, tdat$rccl), na.rm = TRUE)
  if (rng[1] < min(bins$lbinL) || rng[2] > max(bins$lbinU)) {
    warning("Length data (", rng[1], " to ", rng[2], ") fall outside the bin ",
            "structure (", min(bins$lbinL), " to ", max(bins$lbinU),
            ") -- records outside the bins are silently mis-assigned.")
  }

  ## --- Time at liberty, in days -------------------------------------------
  ## Needed by growmod's recovery function when suppress = TRUE. tsteps is
  ## far too coarse for this: with an annual timestep and the half-step
  ## rounding rule, everything from ~6 to ~18 months maps to tsteps = 1,
  ## which is exactly the resolution the recovery function exists to use.
  if (is.null(liberty_days)) {
    if (lib_col %in% names(tdat)) {
      liberty_days <- as.numeric(tdat[[lib_col]])
    } else if (suppress) {
      stop("suppress = TRUE needs time at liberty in days. Either supply ",
           "liberty_days =, or add a '", lib_col, "' column to tdat.")
    } else {
      liberty_days <- rep(NA_real_, nrow(tdat))
    }
  }
  if (length(liberty_days) != nrow(tdat)) {
    stop("liberty_days has length ", length(liberty_days), " but tdat has ",
         nrow(tdat), " rows.")
  }
  if (suppress) {
    if (anyNA(liberty_days)) {
      stop(sum(is.na(liberty_days)), " records have NA time at liberty. ",
           "These cannot be used with suppress = TRUE -- filter them out ",
           "or supply the missing dates.")
    }
    if (any(liberty_days < 0)) {
      stop("liberty_days has negative values -- check the date column order ",
           "(release first, recapture second).")
    }
  }

  ## --- Year indexing ------------------------------------------------------
  ## See @details: NOT length(unique(...)), which undercounts whenever a year
  ## in the span has neither a release nor a recapture.
  yr0    <- min(tdat$relyr, na.rm = TRUE)
  relyr  <- tdat$relyr - yr0 + 1
  nyears <- max(c(tdat$relyr, tdat$recyr), na.rm = TRUE) - yr0 + 1

  ## --- Assemble -----------------------------------------------------------
  datain <- list(
    tagNum   = tdat$tag,
    relyr    = relyr,
    nyears   = nyears,
    Rccl     = tdat$rccl,
    Rlcl     = tdat$rlcl,
    relts    = tdat$relts,
    tsteps   = as.integer(tdat$ntstep),
    nlob     = tdat$num,
    nlbin    = length(bins$lbin),
    ntsteps  = ntsteps,
    lbinL    = bins$lbinL,
    lbinU    = bins$lbinU,
    lbin     = bins$lbin,
    goodts   = goodts,
    M        = M,
    smoother = smoother,
    # model switches
    TemporalGrowth       = TemporalGrowth,
    mpy                  = mpy,
    n_pmoult1            = n_pmoult1,
    ident_wt             = ident_wt,
    use_individual_error = use_individual_error,
    suppress             = suppress,
    liberty              = liberty_days,
    # bookkeeping
    nobs     = nrow(tdat),
    year0    = yr0
  )
  ## --- Period mode --------------------------------------------------------
  ## A record-level period column is mapped up to a year-level vector
  ## covering EVERY year, including those with no records: an animal
  ## released before a gap and recaptured after it needs S defined
  ## throughout, and gap years have no records of their own to label them.
  ## Gaps are filled by carrying the preceding period forward (and back-
  ## filling any leading gap), which is what "periods are contiguous blocks
  ## of time" implies.
  period_mode <- isTRUE(period)
  if (period_mode) {
    if (!TemporalGrowth) {
      stop("period = TRUE but TemporalGrowth = FALSE -- period effects are a ",
           "form of temporal growth variation, so set TemporalGrowth = TRUE.")
    }
    if (is.null(period_col) || !period_col %in% names(tdat)) {
      stop("period = TRUE needs a period column in tdat. Looked for '",
           period_col, "' -- set period_col to its actual name.")
    }
    prec <- as.integer(as.factor(tdat[[period_col]]))   # 1..nperiods, in order
    nperiods <- max(prec, na.rm = TRUE)

    ## Modal period of the records released in each year
    yp <- rep(NA_integer_, nyears)
    for (y in seq_len(nyears)) {
      pv <- prec[relyr == y & !is.na(prec)]
      if (length(pv)) yp[y] <- as.integer(names(which.max(table(pv))))
    }
    if (all(is.na(yp))) stop("No records with a usable period value.")

    ## Carry forward, then back-fill any leading NAs
    for (y in 2:nyears) if (is.na(yp[y])) yp[y] <- yp[y - 1]
    first <- which(!is.na(yp))[1]
    if (first > 1) yp[1:(first - 1)] <- yp[first]

    if (any(diff(yp) < 0)) {
      warning("Periods are not contiguous in time (year-to-period mapping is ",
              "not monotonic). Check the period column -- overlapping periods ",
              "are almost certainly a coding error.")
    }

    datain$year_period <- yp
    datain$nperiods    <- nperiods
    datain$period_mode <- TRUE

    if (!quiet) {
      yrs <- yr0 + seq_len(nyears) - 1
      cat("Period mode: ", nperiods, " periods\n", sep = "")
      for (p in seq_len(nperiods)) {
        inp <- which(yp == p)
        cat("  period ", p, ": ", min(yrs[inp]), "-", max(yrs[inp]),
            "  (", sum(prec == p, na.rm = TRUE), " records)\n", sep = "")
      }
    }
  } else {
    datain$period_mode <- FALSE
  }

  ## --- Measurement error --------------------------------------------------
  ## Estimated from the data unless supplied. A hand-set value has to be
  ## recalculated for every species and dataset, and getting it wrong
  ## biases the moult/no-moult split -- for deep-sea crab the previously
  ## fixed log(2) was nearly double the empirical value.
  if (is.null(LsigError)) {
    datain <- add_sigError(datain, tdat, quiet = quiet,
                            max_liberty = sigError_max_liberty)
  } else {
    datain$LsigError_init       <- LsigError
    datain$LsigError_prior_mean <- LsigError
    datain$LsigError_prior_sd   <- 0.5   # loose: a supplied value is an assertion
  }

  ## add_year_support is an annual-effects concept: in period mode, support
  ## is a property of the period (reported above), not of the year.
  if (TemporalGrowth && !period_mode) {
    datain <- add_year_support(datain, min_support = min_support)
  }

  ## --- Diagnostic summary -------------------------------------------------
  if (!quiet) {
    K <- vapply(seq_len(nrow(tdat)), function(r) {
      nts <- as.integer(tdat$ntstep[r])
      if (is.na(nts) || nts <= 0) return(0L)
      tsv <- c(tdat$relts[r]:ntsteps, rep(1:ntsteps, 30))[1:nts]
      sum(tsv %in% goodts)
    }, integer(1))

    cat("Makedata: ", nrow(tdat), " records, ", nyears, " years, ",
        ntsteps, " timestep(s)/yr, goodts = ",
        paste(goodts, collapse = ","), "\n", sep = "")
    cat("Moult opportunities per animal:\n")
    print(table(K))
    if (sum(K == 1) > 0) {
      cat("  ", sum(K == 1), " animals with exactly one opportunity",
          if (!suppress) " -- consider suppress = TRUE if handling effects are plausible"
          else " (suppression is being estimated)", "\n", sep = "")
    }
    if (sum(K == 0) > 0) {
      cat("  ", sum(K == 0), " animals with NO opportunity (no growth ",
          "information; useful only for measurement error)\n", sep = "")
    }
    if (TemporalGrowth) {
      cat("Supported years: ", sum(datain$yr_supported), " of ", nyears,
          " (min_support = ", min_support, ")\n", sep = "")
    }
  }

  datain
}
