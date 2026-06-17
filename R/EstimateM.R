#' Estimate Natural Mortality from Tag-Recapture Data Using Size-Structured Survival
#'
#' Uses the relative recapture rates across release size classes to estimate natural
#' mortality (M). The key insight is that under homogenous, size-independent fishing
#' mortality (F), the ratio of recapture probabilities between size classes depends
#' only on M — F cancels completely. Animals at larger sizes at release have accumulated
#' more M to reach that size and have fewer years remaining, making them
#' disproportionately less likely to be recaptured when M is high.
#'
#' @param obs Data frame with two columns:
#'   \itemize{
#'     \item \code{release_len} - Carapace length at release (mm), no NAs permitted
#'     \item \code{isrecap} - Binary recapture outcome (1 = recaptured, 0 = not)
#'   }
#' @param dout List output from the best fitted Tag2STM model. Must contain
#'   \code{$stm} as a 3D array \code{[n_bins x n_bins x n_tsteps]}.
#' @param bins List returned by \code{MakeLbin()}. Must contain \code{$lbin}
#'   (midpoints), \code{$lbinL} (lower bounds), and \code{$lbinU} (upper bounds).
#'   The annual STM is constructed internally via
#'   \code{ClipSTM(lbinL[1]+20, lbinL[n-1], gap, Annual = TRUE)}.
#' @param M_range Numeric vector of length 2 giving the search range for M (annual
#'   rate). Default \code{c(0.01, 1.0)}.
#' @param plot Logical. If TRUE produces four diagnostic plots. Default TRUE.
#'
#' @return Invisibly returns a list with components:
#'   \itemize{
#'     \item \code{M_est} - MLE of M (annual rate)
#'     \item \code{M_ci} - 95\% confidence interval from likelihood profile
#'     \item \code{nll_profile} - Data frame of M values and delta NLL
#'     \item \code{age_at_release} - Relative age (moult-years) at release per animal
#'     \item \code{n_recap} - Number recaptured
#'     \item \code{n_total} - Total animals
#'   }
#'
#' @details
#' The annual STM is constructed from \code{dout$stm} via \code{ClipSTM(Annual = TRUE)},
#' which compounds only the active growth timesteps (e.g. timesteps 1 and 3 where
#' growth is non-zero) and clips to the size range defined by \code{bins}.
#'
#' Relative age at release for each animal is estimated by propagating a unit cohort
#' forward from the smallest size bin through the annual STM, recording the number of
#' annual steps until the modal bin reaches the observed release size.
#'
#' Recapture probability for animal \code{i} with relative age \code{a_i} is:
#'
#' \deqn{p_i = \frac{\exp(-M \cdot a_i)}{\sum_j \exp(-M \cdot a_j)} \times N_{recap}}
#'
#' F cancels in this relative formulation. M is estimated by minimising the Bernoulli
#' negative log-likelihood and profiled for 95\% CI using a chi-squared threshold of
#' 1.92 (= qchisq(0.95, 1) / 2).
#'
#' @note
#' Assumes F and reporting rate are size-independent. Size-selective gear would
#' bias M estimates upward.
#'
#' @examples
#' \dontrun{
#' bins <- MakeLbin(start = 50, stop = 178, gap = 2)
#' obs2 <- obs[, c("LCl", "isrecap")]
#' colnames(obs2) <- c("release_len", "isrecap")
#'
#' result <- EstimateM(obs = obs2, dout = dout, bins = bins)
#' result$M_est
#' result$M_ci
#' }
#'
#' @export
EstimateM <- function(obs,
                      dout,
                      bins,
                      M_range = c(0.01, 1.0),
                      plot    = TRUE) {

  ## ── 0. Input checks ──────────────────────────────────────────────────────────
  if (!all(c("release_len", "isrecap") %in% names(obs)))
    stop("obs must have columns 'release_len' and 'isrecap'.")
  if (any(is.na(obs$release_len)))
    stop("release_len contains NAs — release size must be recorded for all animals.")
  if (!all(obs$isrecap %in% c(0, 1)))
    stop("isrecap must be binary (0 or 1).")
  if (is.null(dout$stm))
    stop("dout must contain $stm [n_bins x n_bins x n_tsteps].")
  if (!all(c("lbin", "lbinL", "lbinU") %in% names(bins)))
    stop("bins must contain $lbin, $lbinL, $lbinU (from MakeLbin()).")

  n_total <- nrow(obs)
  n_recap <- sum(obs$isrecap)
  message(sprintf("EstimateM: %d animals total, %d recaptured (%.1f%%).",
                  n_total, n_recap, 100 * n_recap / n_total))

  ## ── 1. Build annual STM via ClipSTM ──────────────────────────────────────────
  n_bins <- length(bins$lbin)
  low_lb <- bins$lbinL[1] + 20                        # undo padding on first bin
  up_lb  <- bins$lbinL[length(bins$lbinL) - 1]        # second-last lower bound
  gap    <- bins$lbinU[2] - bins$lbinL[2]             # bin width from non-padded bin

  message(sprintf("Building annual STM via ClipSTM(%g, %g, %g, Annual=TRUE).",
                  low_lb, up_lb, gap))
  STM_annual <- ClipSTM(low_lb, up_lb, gap, Annual = TRUE)

  if (!is.matrix(STM_annual) || nrow(STM_annual) != ncol(STM_annual))
    stop("ClipSTM did not return a square matrix — check bins and dout.")

  stm_bins <- nrow(STM_annual)

  ## ── 2. Age-at-release via forward STM propagation ────────────────────────────
  # Starting bin: smallest bin in the clipped STM
  # lbin from bins may be wider than the clipped STM, so build clipped lbin
  lbin_clip <- bins$lbin[bins$lbin >= low_lb & bins$lbin <= (up_lb + gap)]
  if (length(lbin_clip) != stm_bins) {
    # Fallback: construct evenly spaced bins matching STM dimensions
    lbin_clip <- seq(low_lb + gap / 2, by = gap, length.out = stm_bins)
  }
  min_bin <- 1L  # always start from smallest clipped bin

  # Propagate unit cohort forward; return steps until modal bin >= target bin
  .age_from_stm <- function(target_bin, max_age = 40L) {
    dist <- rep(0, stm_bins)
    dist[min_bin] <- 1
    for (yr in seq_len(max_age)) {
      dist <- as.numeric(STM_annual %*% dist)
      dist <- dist / sum(dist)             # renormalise for numerical stability
      if (which.max(dist) >= target_bin) return(yr)
    }
    return(max_age)
  }

  # Map each animal's release size to closest bin in the clipped lbin
  release_bins <- sapply(obs$release_len,
                         function(l) which.min(abs(lbin_clip - l)))

  message("Computing relative age at release for all animals...")
  age_at_release <- sapply(release_bins, .age_from_stm)
  age_at_release <- pmax(age_at_release, 0.5)   # floor at 0.5 to avoid zero

  ## ── 3. Negative log-likelihood ───────────────────────────────────────────────
  # p_i = exp(-M * a_i) / sum_j(exp(-M * a_j)) * n_recap
  # Bernoulli: ll = sum[ isrecap * log(p_i) + (1-isrecap) * log(1 - p_i) ]
  # F cancels in the relative formulation; only M drives size contrast.
  nll <- function(M) {
    surv  <- exp(-M * age_at_release)
    p_raw <- surv / sum(surv) * n_recap       # scale to expected n recaptures
    p_i   <- pmin(pmax(p_raw, 1e-10), 1 - 1e-10)
    -sum(obs$isrecap * log(p_i) + (1 - obs$isrecap) * log(1 - p_i))
  }

  ## ── 4. MLE ───────────────────────────────────────────────────────────────────
  message("Optimising M...")
  opt     <- optimise(nll, interval = M_range, tol = 1e-6)
  M_est   <- opt$minimum
  nll_min <- opt$objective

  ## ── 5. Likelihood profile and 95% CI ─────────────────────────────────────────
  M_seq       <- seq(M_range[1], M_range[2], length.out = 200)
  nll_seq     <- sapply(M_seq, nll)
  nll_profile <- data.frame(M = M_seq, delta_nll = nll_seq - nll_min)

  # 95% CI: delta NLL <= 1.92  (= qchisq(0.95, 1) / 2)
  ci_idx <- which(nll_profile$delta_nll <= 1.92)
  M_ci <- if (length(ci_idx) >= 2) {
    c(M_seq[min(ci_idx)], M_seq[max(ci_idx)])
  } else {
    c(NA_real_, NA_real_)
  }

  message(sprintf("M estimate: %.4f  (95%% CI: %.4f - %.4f)",
                  M_est, M_ci[1], M_ci[2]))

  ## ── 6. Diagnostic plots ───────────────────────────────────────────────────────
  if (plot) {

    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

    # Plot 1: Likelihood profile
    plot(nll_profile$M, nll_profile$delta_nll,
         type = "l", lwd = 2, col = "steelblue",
         xlab = "M (annual rate)", ylab = "Delta NLL",
         main = "Likelihood Profile for M",
         ylim = c(0, min(max(nll_profile$delta_nll), 10)))
    abline(h = 1.92, lty = 2, col = "red")
    abline(v = M_est, lty = 1, col = "darkgreen", lwd = 1.5)
    if (!is.na(M_ci[1])) {
      abline(v = M_ci[1], lty = 3, col = "orange")
      abline(v = M_ci[2], lty = 3, col = "orange")
    }
    legend("topright",
           legend = c(sprintf("MLE = %.3f", M_est),
                      sprintf("95%% CI: %.3f - %.3f", M_ci[1], M_ci[2]),
                      "Chi-sq threshold (1.92)"),
           lty = c(1, 3, 2),
           col = c("darkgreen", "orange", "red"),
           bty = "n", cex = 0.8)

    # Plot 2: Age at release distribution
    hist(age_at_release, breaks = 20,
         col = "steelblue", border = "white",
         xlab = "Relative age at release (moult-years)",
         ylab = "Frequency",
         main = "Age at Release Distribution")

    # Plot 3: Observed recapture rate by release size decile
    size_breaks   <- quantile(obs$release_len, probs = seq(0, 1, by = 0.1))
    size_class    <- cut(obs$release_len, breaks = size_breaks, include.lowest = TRUE)
    recap_by_size <- tapply(obs$isrecap, size_class, mean, na.rm = TRUE)
    midpoints     <- (size_breaks[-1] + size_breaks[-length(size_breaks)]) / 2
    surv_pred     <- exp(-M_est * age_at_release)
    p_pred        <- surv_pred / sum(surv_pred) * n_recap
    pred_by_size  <- tapply(p_pred, size_class, mean, na.rm = TRUE) * 100

    plot(midpoints, recap_by_size * 100,
         type = "b", pch = 16, col = "steelblue", lwd = 1.5,
         xlab = "Release size (mm CL)", ylab = "Recapture rate (%)",
         main = "Recapture Rate by Release Size",
         ylim = c(0, max(recap_by_size * 100, na.rm = TRUE) * 1.25))
    lines(midpoints, pred_by_size, col = "red", lwd = 1.5, lty = 2)
    legend("topright",
           legend = c("Observed", sprintf("Predicted (M = %.3f)", M_est)),
           col = c("steelblue", "red"), lty = c(1, 2), pch = c(16, NA),
           bty = "n", cex = 0.8)

    # Plot 4: Observed vs predicted by predicted-probability decile
    p_all       <- pmin(pmax(p_pred, 1e-10), 1 - 1e-10)
    pred_breaks <- quantile(p_all, probs = seq(0, 1, by = 0.1))
    pred_class  <- cut(p_all, breaks = pred_breaks, include.lowest = TRUE)
    obs_bin     <- tapply(obs$isrecap, pred_class, mean, na.rm = TRUE)
    pred_bin    <- tapply(p_all,       pred_class, mean, na.rm = TRUE)
    plot(pred_bin * 100, obs_bin * 100,
         pch = 16, col = "steelblue",
         xlab = "Predicted recapture prob. (%)",
         ylab = "Observed recapture rate (%)",
         main = "Observed vs Predicted (decile bins)")
    abline(0, 1, col = "red", lty = 2)
  }

  ## ── 7. Return ─────────────────────────────────────────────────────────────────
  invisible(list(
    M_est          = M_est,
    M_ci           = M_ci,
    nll_profile    = nll_profile,
    age_at_release = age_at_release,
    n_recap        = n_recap,
    n_total        = n_total
  ))
}
