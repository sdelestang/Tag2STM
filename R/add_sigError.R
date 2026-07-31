#' Estimate Measurement Error from Tag-Recapture Increments
#'
#' Estimates the per-measurement standard deviation (\code{sigError}) directly
#' from the data and stores it in \code{datain} as a prior centre and width,
#' so that \code{\link{growmod}}'s \code{PenSigError} and
#' \code{\link{Makepin}}'s \code{LsigError} starting value both derive from
#' the data rather than a hand-set constant. Called automatically by
#' \code{\link{Makedata}} when \code{LsigError = NULL}.
#'
#' @param datain A \code{datain} list (from \code{\link{Makedata}}).
#' @param tdat The tag-recapture data frame \code{datain} was built from.
#' @param round_unit Numeric, the recording resolution in mm for
#'   whole-unit-recorded lengths. Default 1. Used for Sheppard's correction.
#' @param max_liberty Numeric, days. Only records at liberty up to this are
#'   used for the reflection estimate. Default 730. A negative increment is
#'   measurement error at any liberty, but the fraction of records that are
#'   genuinely un-moulted falls as liberty grows, so the share of negatives
#'   that are tag-matching or transcription errors rises -- past roughly one
#'   intermoult period those dominate and inflate the estimate. Should scale
#'   with the species: 730 days suits a slow-growing deep-water crab, but a
#'   lobster moulting twice a year needs a much shorter window. The printed
#'   profile of sigError against cutoff shows where the estimate starts to
#'   climb, which is the contamination becoming visible.
#' @param outlier_k Numeric, default 5. Negative increments beyond this
#'   multiple of the robust (MAD) scale are trimmed before the moment
#'   estimate and reported. Set wide deliberately: the aim is to catch a
#'   gross tag mismatch, not to trim a heavy tail.
#' @param quiet Logical, suppress the diagnostic printout. Default
#'   \code{FALSE}.
#'
#' @return \code{datain} with \code{LsigError_init} (log scale, for
#'   \code{Makepin}), \code{LsigError_prior_mean} and
#'   \code{LsigError_prior_sd} (for \code{growmod}'s \code{PenSigError})
#'   added.
#'
#' @details
#' \strong{Why reflection.} A moult swells the animal with water, so it can
#' only ever increase length. Every NEGATIVE observed increment is therefore
#' pure measurement error, whatever the animal's moult history or time at
#' liberty. Since that error is symmetric about zero, the negative tail is
#' a clean half-sample of it:
#'
#'   \code{sd(difference) = sqrt(mean(neg^2))},
#'   \code{sigError = sd(difference) / sqrt(2)}
#'
#' (the difference of two independent measurements carries error at both
#' ends, hence the \code{sqrt(2)}).
#'
#' This is preferred over the more obvious route of using zero-opportunity
#' records (\code{tsteps == 0}), because those are defined by the
#' half-timestep rounding rule rather than by biology: with an annual
#' timestep an animal at liberty five months is recorded as having had no
#' opportunity, but it may well have moulted. That group is therefore
#' contaminated in its POSITIVE tail. The zero-opportunity estimate is
#' still computed and reported as a cross-check, since agreement between
#' the two is the best evidence the number is trustworthy.
#'
#' \strong{Robustness.} A moment estimator is destroyed by a couple of
#' transcription errors -- in the deep-sea crab data two bad records
#' inflated a raw sd from ~1.4 mm to 3.58 mm. Gross outliers are therefore
#' trimmed at \code{outlier_k} times the robust (MAD) scale and reported,
#' with the moment estimator applied to what remains.
#'
#' The estimator is NOT switched to a MAD-based one when the two disagree.
#' MAD lies below the standard deviation for any peaked, heavy-tailed error
#' distribution without contamination of any kind, so a moment/MAD ratio
#' cannot distinguish the two situations. Note also that reflection is
#' already immune to the failure mode above: those crab outliers were
#' POSITIVE, and a negatives-only sample never sees them. Only a gross
#' negative needs handling, which trimming does.
#'
#' \strong{Sheppard's correction.} Lengths recorded to the nearest whole mm
#' carry an extra \code{round_unit^2 / 12} of variance from rounding alone.
#' The correction is applied in proportion to the fraction of records whose
#' release and recapture lengths are both whole numbers, since recording
#' precision is often mixed within one dataset.
#'
#' \strong{Prior width.} \code{LsigError_prior_sd} is the sampling standard
#' error of the estimate on the log scale, approximately
#' \code{1 / sqrt(2 * n_neg)}, not an arbitrary choice. With a few hundred
#' negative increments this is tight (~0.05), which is what makes it safe
#' to estimate \code{LsigError} inside the model rather than fixing it:
#' the historical reason for fixing it was that a loose prior let the model
#' explain real growth away as measurement noise, and a prior this narrow
#' leaves no room for that.
#'
#' @seealso \code{\link{Makedata}}, \code{\link{growmod}},
#'   \code{\link{Makepin}}
#'
#' @export
add_sigError <- function(datain, tdat, round_unit = 1, quiet = FALSE,
                          outlier_k = 5, max_liberty = 730) {

  inc <- tdat$rccl - tdat$rlcl
  lib <- datain$liberty
  if (is.null(lib)) lib <- rep(NA_real_, length(inc))
  ok  <- !is.na(inc)

  ## --- Liberty profile of the estimate (diagnostic) -----------------------
  ## Printed so the max_liberty choice is visible rather than assumed. If
  ## sigError climbs steadily with the cutoff, that is contamination
  ## entering, not measurement error growing.
  cuts <- c(365, 730, 1095, 1460, Inf)
  prof <- do.call(rbind, lapply(cuts, function(L) {
    sel <- ok & inc < 0 & (is.na(lib) | lib <= L)
    n   <- sum(sel)
    data.frame(max_liberty = L, n_neg = n,
               sigError = if (n >= 5) sqrt(mean(inc[sel]^2)) / sqrt(2) else NA_real_)
  }))

  ## --- Restrict to short-liberty records ----------------------------------
  ## A negative increment is measurement error at any liberty, but the
  ## fraction of records that are genuinely un-moulted falls as liberty
  ## grows, so the share of negatives that are tag-matching or transcription
  ## errors rises. Past roughly one intermoult period those errors dominate
  ## and inflate the estimate. max_liberty should scale with the species'
  ## intermoult period -- 730 days suits a slow-growing deep-water crab;
  ## reduce it substantially for a fast-growing species (a lobster moulting
  ## twice a year does not need, and should not use, a two-year window).
  if (all(is.na(lib))) {
    if (!quiet) {
      cat("add_sigError: no liberty data -- using all records ",
          "(max_liberty ignored).\n", sep = "")
    }
    keep <- ok
  } else {
    keep <- ok & !is.na(lib) & lib <= max_liberty
  }

  neg <- inc[keep & inc < 0]
  if (length(neg) < 10) {
    stop("Only ", length(neg), " negative increments within max_liberty = ",
         max_liberty, " -- too few to estimate measurement error by ",
         "reflection. Raise max_liberty or supply LsigError directly.")
  }
  refl <- c(neg, -neg)

  ## --- Estimate, with gross outliers trimmed ------------------------------
  ## Always the moment estimator, computed after removing genuinely extreme
  ## values. An earlier version switched wholesale to the MAD estimate when
  ## moment/MAD exceeded a ratio, which was wrong: MAD sits below the sd for
  ## any peaked (leptokurtic) error distribution, with no contamination
  ## whatever, so the ratio test cannot tell the two apart. On lobster data
  ## it fired at a ratio of 1.44 and returned 0.205 mm when both the liberty
  ## profile and the zero-opportunity cross-check said 0.303, flagging
  ## ordinary -0.9 and -1.0 mm records as outliers.
  ##
  ## Note reflection is already immune to the failure mode that motivated
  ## the fallback: the crab dataset's damaging outliers were POSITIVE
  ## (+8 and +35 mm), and a negatives-only sample never sees them. Only a
  ## gross NEGATIVE -- a tag mismatch reading -30 mm, say -- needs handling,
  ## and trimming at a wide multiple of the robust scale does that without
  ## penalising a heavy-tailed but legitimate error distribution.
  mad_scale <- mad(refl, center = 0)
  sd_moment_all <- sqrt(mean(neg^2))

  if (mad_scale > 0) {
    gross    <- abs(neg) > outlier_k * mad_scale
    neg_used <- neg[!gross]
  } else {
    gross    <- rep(FALSE, length(neg))
    neg_used <- neg
  }
  if (length(neg_used) < 10) {
    warning("Trimming left fewer than 10 negatives -- using the untrimmed ",
            "sample. Check outlier_k.")
    gross    <- rep(FALSE, length(neg))
    neg_used <- neg
  }
  if (sum(gross) > 0.1 * length(neg)) {
    warning(sum(gross), " of ", length(neg), " negatives flagged as gross ",
            "outliers. That is a lot -- the error distribution may be very ",
            "heavy-tailed rather than contaminated. Inspect before trusting.")
  }

  sd_diff <- sqrt(mean(neg_used^2))

  ## --- Sheppard's correction, weighted by the integer-recorded fraction ---
  both_int <- abs(tdat$rlcl - round(tdat$rlcl)) < 1e-8 &
              abs(tdat$rccl - round(tdat$rccl)) < 1e-8
  frac_int <- mean(both_int[keep], na.rm = TRUE)
  shep     <- frac_int * (round_unit^2) / 12

  var_diff <- sd_diff^2 - shep
  if (var_diff <= 0) {
    warning("Sheppard's correction exceeds the estimated variance -- ",
            "recording resolution may be coarser than the true error. ",
            "Using the uncorrected estimate.")
    var_diff <- sd_diff^2
  }
  sigError <- sqrt(var_diff / 2)

  ## --- Cross-check: zero-opportunity records ------------------------------
  ## Contaminated in the positive tail (see @details), so reported not used.
  ## Uses the SAME estimator as the main route -- reflection on negatives,
  ## moment-based -- so the two numbers are comparable. An earlier version
  ## took the MAD of all zero-opportunity increments, which is not the same
  ## quantity: that distribution has a spike at exactly zero from whole-mm
  ## recording, which deflates the MAD, so the cross-check printed a value
  ## ~30% below the main estimate and looked like a disagreement when there
  ## was none.
  K0 <- if (!is.null(datain$tsteps)) datain$tsteps == 0 else rep(FALSE, nrow(tdat))
  sig_K0 <- NA_real_
  if (sum(K0, na.rm = TRUE) >= 10) {
    inc0 <- (tdat$rccl - tdat$rlcl)[K0]
    inc0 <- inc0[!is.na(inc0)]
    neg0 <- inc0[inc0 < 0]
    sig_K0 <- if (length(neg0) >= 5) sqrt(mean(neg0^2)) / sqrt(2) else NA_real_
  }

  ## --- Prior width = sampling SE on the log scale -------------------------
  n_neg    <- length(neg)
  prior_sd <- 1 / sqrt(2 * n_neg)

  datain$LsigError_init       <- log(sigError)
  datain$LsigError_prior_mean <- log(sigError)
  datain$LsigError_prior_sd   <- prior_sd

  if (!quiet) {
    cat("add_sigError:\n")
    cat("  sigError by liberty cutoff (days):\n")
    for (i in seq_len(nrow(prof))) {
      cat("    <= ", formatC(prof$max_liberty[i], width = 5),
          "  n_neg ", formatC(prof$n_neg[i], width = 4),
          "  sigError ",
          if (is.na(prof$sigError[i])) "   --" else formatC(prof$sigError[i], format = "f", digits = 3),
          "\n", sep = "")
    }
    cat("  using max_liberty        : ", max_liberty, " days\n", sep = "")
    cat("  negative increments used : ", n_neg, " of ", sum(ok), "\n", sep = "")
    cat("  difference sd  moment    : ", round(sd_moment_all, 3), "\n", sep = "")
    cat("  robust scale (MAD)       : ", round(mad_scale, 3),
        "   (reference only, not used)\n", sep = "")
    if (sum(gross) > 0) {
      cat("  trimmed ", sum(gross), " gross outlier(s) beyond ", outlier_k,
          " x robust scale: ",
          paste(round(sort(neg[gross]), 1), collapse = ", "), "\n", sep = "")
      cat("    (check these against the raw data -- likely tag mismatches)\n")
      cat("  difference sd  trimmed   : ", round(sd_diff, 3), "\n", sep = "")
    }
    cat("  whole-mm recorded pairs  : ", round(100 * frac_int), "%",
        "  (Sheppard correction ", signif(shep, 2), ")\n", sep = "")
    cat("  => sigError              : ", round(sigError, 3), " mm",
        "   [LsigError = ", round(log(sigError), 3), "]\n", sep = "")
    cat("  prior sd (log scale)     : ", round(prior_sd, 3), "\n", sep = "")
    if (!is.na(sig_K0)) {
      cat("  cross-check, zero-opportunity records: ", round(sig_K0, 3), " mm\n",
          sep = "")
      if (abs(log(sig_K0) - log(sigError)) > 0.4) {
        cat("    NOTE: the two routes disagree by more than ~50%. Worth\n",
            "    checking recording precision or moult contamination before\n",
            "    trusting either.\n", sep = "")
      }
    }
  }

  datain
}
