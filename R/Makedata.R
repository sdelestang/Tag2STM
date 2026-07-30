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
#' @param LsigError Numeric or \code{NULL}. Log measurement-error SD, stored
#'   as \code{LsigError_init} for \code{Makepin} to use as its default.
#'   \code{NULL} (default) leaves \code{Makepin}'s own default in force.
#'   This should come from the data rather than being guessed -- for the
#'   deep-sea crab data the fixed \code{log(2)} was nearly double the
#'   empirical value and biased the moult/no-moult split.
#' @param mpy Numeric, minimum moults per year (biological floor on Pmoult).
#'   Default 0 (no floor). See \code{\link{growmod}}.
#' @param n_pmoult1 Integer, number of smallest length bins with Pmoult
#'   hard-fixed at exactly 1. Default 1.
#' @param ident_wt Numeric, weight on the internal identification mixture.
#'   Default 1. Set to 0 for a pure marginal fit.
#' @param use_individual_error Logical, default \code{FALSE}.
#' @param suppress Logical, estimate tagging-induced moult suppression in
#'   each animal's first post-release moult opportunity. Default
#'   \code{FALSE}.
#' @param suppress_compensate Logical, allow part of a suppressed moult to
#'   be recovered in the second opportunity (delay rather than loss).
#'   Default \code{TRUE}; ignored when \code{suppress = FALSE}.
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
                     LsigError = NULL,
                     mpy = 0,
                     n_pmoult1 = 1,
                     ident_wt = 1,
                     use_individual_error = FALSE,
                     suppress = FALSE,
                     suppress_compensate = TRUE,
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
    suppress_compensate  = suppress_compensate,
    # bookkeeping
    nobs     = nrow(tdat),
    year0    = yr0
  )
  if (!is.null(LsigError)) datain$LsigError_init <- LsigError

  if (TemporalGrowth) datain <- add_year_support(datain, min_support = min_support)

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
