#' Estimate Natural Mortality from Tag-Recapture Data Using Size-Structured Survival
#'
#' Estimates natural mortality (M) from binary tag-recapture outcomes using two
#' complementary signals:
#' \enumerate{
#'   \item \strong{Size at release} — converted to relative age via a fitted STM.
#'     Under constant, size-independent F, the ratio of recapture probabilities
#'     across size classes depends only on M (F cancels analytically).
#'   \item \strong{Time at liberty} — used as a weight on each animal's likelihood
#'     contribution. A crab not seen after 4 years is stronger evidence of low
#'     survival than one not seen after 1 year. Liberty is NOT used to estimate Z
#'     absolutely (too many confounders: tag loss, spatial movement, reporting rate
#'     variation) but purely to up-weight long-unseen animals in the Bernoulli NLL.
#' }
#'
#' @param obs Data frame with three columns:
#'   \itemize{
#'     \item \code{release_len} - Carapace length at release (mm), no NAs
#'     \item \code{isrecap} - Binary recapture outcome (1 = recaptured, 0 = not)
#'     \item \code{Liberty} - Days at liberty. For recaptured animals: recapture
#'       date minus release date. For non-recaptured animals: today's date (or
#'       study end date) minus release date (right-censored).
#'   }
#'   Should contain \strong{all} releases including sub-legal animals so that
#'   the auto min-size cut correctly identifies the gear selectivity ramp-up.
#' @param dout List output from the best fitted Tag2STM model. Must contain
#'   \code{$stm} as a 3D array \code{[n_bins x n_bins x n_tsteps]}.
#' @param bins List returned by \code{MakeLbin()}. Must contain \code{$lbin}
#'   (midpoints), \code{$lbinL} (lower bounds), and \code{$lbinU} (upper bounds).
#' @param min_cut Optional numeric. Starting minimum size cut (mm). If NULL
#'   (default), determined automatically as the first size bin after which
#'   recapture rate declines for at least \code{n_consec} consecutive bins.
#' @param n_consec Integer. Consecutive declining bins for auto min_cut. Default 3.
#' @param step_mm Numeric. Increment (mm) for profiling over min_cut. Default 1.
#' @param n_steps Integer. Maximum number of min_cut steps to profile. Default 20.
#' @param min_recap Integer. Stop profiling if fewer recaptures remain. Default 50.
#' @param lambda Numeric. Decay constant controlling how steeply liberty weight
#'   increases with time at liberty (years). Weight for non-recaptured animal i is
#'   \eqn{1 - \exp(-\lambda \cdot t_i)}. Recaptured animals always have weight 1.
#'   Higher lambda gives more weight to long-unseen animals. Default 0.3.
#' @param M_range Numeric vector of length 2. Search range for M. Default
#'   \code{c(0.01, 0.5)}.
#' @param plot Logical. If TRUE produces diagnostic plots. Default TRUE.
#'
#' @return Invisibly returns a list with components:
#'   \itemize{
#'     \item \code{M_est} - Weighted central M estimate across stable plateau steps
#'     \item \code{M_ci} - 95\% CI for M from best-fitting plateau step
#'     \item \code{cut_profile} - Data frame with min_cut, M_est, M_lwr, M_upr,
#'       n_recap, ci_width, ssr at each step
#'     \item \code{min_cut_auto} - Auto-detected starting min_cut (mm)
#'     \item \code{age_at_release} - Fractional ages from the starting min_cut step
#'     \item \code{nll_profile} - NLL profile from the starting min_cut step
#'   }
#'
#' @details
#' The liberty-weighted Bernoulli NLL for animal \eqn{i} is:
#'
#' \deqn{\ell_i = w_i \left[ y_i \log(p_i) + (1 - y_i) \log(1 - p_i) \right]}
#'
#' where \eqn{y_i} is the recapture indicator, \eqn{p_i} is the size-based
#' relative recapture probability, and the liberty weight is:
#'
#' \deqn{w_i = \begin{cases} 1 & \text{if recaptured} \\ 1 - \exp(-\lambda t_i) & \text{if not recaptured} \end{cases}}
#'
#' A crab not seen after 4 years (\eqn{t = 4}) receives weight
#' \eqn{1 - e^{-0.3 \times 4} \approx 0.70} versus \eqn{1 - e^{-0.3 \times 1} \approx 0.26}
#' for one at liberty only 1 year, correctly reflecting that the longer absence
#' is stronger evidence of death.
#'
#' @note
#' Assumes F and reporting rate are size-independent above \code{min_cut}.
#' Liberty is used for weighting only — Z and F are not estimated. Fisher
#' targeting of larger animals and spatial depth refuge will cause upward bias in M.
#'
#' @examples
#' \dontrun{
#' bins <- MakeLbin(start = 50, stop = 198, gap = 2)
#' EstM <- obs %>%
#'   filter(!is.na(LCl)) %>%
#'   rename(release_len = LCl) %>%
#'   mutate(
#'     Ldate   = as.Date(Ldate),
#'     date    = as.Date(date),
#'     date    = if_else(is.na(date), Sys.Date(), date),
#'     Liberty = as.numeric(date - Ldate)
#'   ) %>%
#'   dplyr::select(release_len, isrecap, Liberty)
#'
#' result <- EstimateM(obs = EstM, dout = dout, bins = bins)
#' result$M_est
#' result$M_ci
#' }
#'
#' @export
EstimateM <- function(obs,
                      dout,
                      bins,
                      min_cut   = NULL,
                      n_consec  = 3L,
                      step_mm   = 1,
                      n_steps   = 20L,
                      min_recap = 50L,
                      lambda    = 0.3,
                      M_range   = c(0.01, 0.5),
                      plot      = TRUE) {

  ## ── 0. Input checks ──────────────────────────────────────────────────────────
  if (!all(c("release_len", "isrecap", "Liberty") %in% names(obs)))
    stop("obs must have columns 'release_len', 'isrecap', and 'Liberty'.")
  if (any(is.na(obs$release_len)))
    stop("release_len contains NAs — release size must be recorded for all animals.")
  if (!all(obs$isrecap %in% c(0, 1)))
    stop("isrecap must be binary (0 or 1).")
  if (any(is.na(obs$Liberty) | obs$Liberty < 0))
    stop("Liberty must be non-negative and non-NA for all animals.")
  if (is.null(dout$stm))
    stop("dout must contain $stm [n_bins x n_bins x n_tsteps].")
  if (!all(c("lbin", "lbinL", "lbinU") %in% names(bins)))
    stop("bins must contain $lbin, $lbinL, $lbinU (from MakeLbin()).")

  ## ── 1. Build annual STM via ClipSTM ──────────────────────────────────────────
  low_lb <- bins$lbinL[2]
  up_lb  <- bins$lbinL[length(bins$lbinL) - 1]
  gap    <- bins$lbinU[2] - bins$lbinL[2]

  message(sprintf("Building annual STM via ClipSTM(%g, %g, %g, Annual=TRUE, return=TRUE).",
                  low_lb, up_lb, gap))
  STM_annual <- ClipSTM(low_lb, up_lb, gap, Annual = TRUE, return = TRUE)

  if (!is.matrix(STM_annual) || nrow(STM_annual) != ncol(STM_annual))
    stop("ClipSTM did not return a square matrix — check bins and dout.")

  stm_bins <- nrow(STM_annual)

  ## ── 2. Build clipped lbin to match STM dimensions ────────────────────────────
  lbin_clip <- bins$lbin[bins$lbin >= low_lb & bins$lbin <= (up_lb + gap)]
  if (length(lbin_clip) != stm_bins)
    lbin_clip <- seq(low_lb + gap / 2, by = gap, length.out = stm_bins)

  min_bin <- 1L

  ## ── 3. Auto min-size cut ─────────────────────────────────────────────────────
  .find_min_cut <- function(release_len, isrecap, n_bins = 10, n_consec = 3L) {
    breaks     <- unique(quantile(release_len, probs = seq(0, 1, length.out = n_bins + 1)))
    if (length(breaks) < 3) return(min(release_len))
    size_class <- cut(release_len, breaks = breaks, include.lowest = TRUE)
    minimums   <- as.numeric(tapply(release_len, size_class, min,  na.rm = TRUE))
    recap_rate <- as.numeric(tapply(isrecap,     size_class, mean, na.rm = TRUE))
    n <- length(recap_rate)
    for (i in 1:(n - n_consec)) {
      if (all(diff(recap_rate[i:(i + n_consec)]) <= 0)) return(minimums[i])
    }
    message("Auto min-size cut: no clear declining run found — using min(release_len).")
    return(min(release_len))
  }

  if (is.null(min_cut)) {
    min_cut <- .find_min_cut(obs$release_len, obs$isrecap, n_consec = n_consec)
    message(sprintf(
      "Auto min-size cut: %g mm (first bin with %d consecutive declining recapture rates).",
      min_cut, n_consec))
  } else {
    message(sprintf("Using supplied min_cut: %g mm.", min_cut))
  }
  min_cut_auto     <- min_cut
  min_cut_auto_rnd <- round(min_cut)

  ## ── 4. Fractional age-at-release via forward STM propagation ─────────────────
  .age_from_stm <- function(target_bin, max_age = 40L) {
    dist       <- rep(0, stm_bins)
    dist[min_bin] <- 1
    prev_modal <- min_bin
    for (yr in seq_len(max_age)) {
      dist      <- as.numeric(STM_annual %*% dist)
      dist      <- dist / sum(dist)
      modal_now <- which.max(dist)
      if (modal_now >= target_bin) {
        frac <- (target_bin - prev_modal) / (modal_now - prev_modal)
        frac <- pmin(pmax(frac, 0), 1)
        return((yr - 1) + frac)
      }
      prev_modal <- modal_now
    }
    return(as.numeric(max_age))
  }

  ## ── 5. Core estimation helper ─────────────────────────────────────────────────
  .estimate_at_cut <- function(cut) {
    obs_i <- obs[obs$release_len >= cut, ]
    n_rec <- sum(obs_i$isrecap)
    if (n_rec < min_recap) return(NULL)

    rel_bins <- sapply(obs_i$release_len, function(l) which.min(abs(lbin_clip - l)))
    ages     <- pmax(sapply(rel_bins, .age_from_stm), 0.5)
    t_lib    <- obs_i$Liberty / 365.25   # days to years

    # Liberty weights: recaptured animals always weight 1;
    # non-recaptured animals weighted by 1 - exp(-lambda * t)
    # so longer-unseen animals contribute more to the likelihood
    lib_wt <- ifelse(obs_i$isrecap == 1, 1, 1 - exp(-lambda * t_lib))
    lib_wt <- pmax(lib_wt, 1e-6)   # floor to avoid zero weights

    n_total <- nrow(obs_i)

    # Liberty-weighted Bernoulli NLL — M only, single parameter
    nll <- function(M) {
      surv  <- exp(-M * ages)
      p_raw <- surv / sum(surv) * n_rec
      p_i   <- pmin(pmax(p_raw, 1e-10), 1 - 1e-10)
      -sum(lib_wt * (obs_i$isrecap * log(p_i) +
                       (1 - obs_i$isrecap) * log(1 - p_i)))
    }

    opt     <- optimise(nll, interval = M_range, tol = 1e-6)
    M_est   <- opt$minimum
    nll_min <- opt$objective

    M_seq   <- seq(M_range[1], M_range[2], length.out = 200)
    nll_seq <- sapply(M_seq, nll)
    ci_idx  <- which(nll_seq - nll_min <= 1.92)
    M_ci <- if (length(ci_idx) >= 2) {
      c(M_seq[min(ci_idx)], M_seq[max(ci_idx)])
    } else c(NA_real_, NA_real_)

    # SSR: observed vs predicted recapture rate by size decile
    surv_pred <- exp(-M_est * ages)
    p_pred    <- surv_pred / sum(surv_pred) * n_rec
    sz_breaks <- unique(quantile(obs_i$release_len, probs = seq(0, 1, by = 0.1)))
    if (length(sz_breaks) < 3)
      sz_breaks <- seq(min(obs_i$release_len), max(obs_i$release_len), length.out = 11)
    sz_class  <- cut(obs_i$release_len, breaks = sz_breaks, include.lowest = TRUE)
    obs_rate  <- tapply(obs_i$isrecap, sz_class, mean,   na.rm = TRUE)
    pred_rate <- tapply(p_pred,        sz_class, mean,   na.rm = TRUE)
    ssr       <- sum((obs_rate - pred_rate / 100)^2, na.rm = TRUE)

    list(M_est       = M_est,
         M_ci        = M_ci,
         nll_profile = data.frame(M = M_seq, delta_nll = nll_seq - nll_min),
         ages        = ages,
         t_lib       = t_lib,
         lib_wt      = lib_wt,
         p_pred      = p_pred,
         sz_breaks   = sz_breaks,
         n_recap     = n_rec,
         n_total     = n_total,
         ssr         = ssr)
  }

  ## ── 6. Profile over min_cut ───────────────────────────────────────────────────
  cuts <- seq(min_cut_auto_rnd, by = step_mm, length.out = n_steps)
  message(sprintf("Profiling M over %d min_cut steps from %g mm in %g mm increments...",
                  n_steps, min_cut_auto_rnd, step_mm))

  cut_profile <- data.frame(min_cut  = cuts,
                            M_est    = NA_real_,
                            M_lwr    = NA_real_,
                            M_upr    = NA_real_,
                            n_recap  = NA_integer_,
                            ci_width = NA_real_,
                            ssr      = NA_real_)

  first_result <- NULL
  for (i in seq_along(cuts)) {
    res_i <- .estimate_at_cut(cuts[i])
    if (is.null(res_i)) {
      message(sprintf("  Step %d (%g mm): fewer than %d recaptures — stopping.",
                      i, cuts[i], min_recap))
      break
    }
    cut_profile$M_est[i]    <- res_i$M_est
    cut_profile$M_lwr[i]    <- res_i$M_ci[1]
    cut_profile$M_upr[i]    <- res_i$M_ci[2]
    cut_profile$n_recap[i]  <- res_i$n_recap
    cut_profile$ci_width[i] <- diff(res_i$M_ci)
    cut_profile$ssr[i]      <- res_i$ssr
    message(sprintf(
      "  Step %d (%g mm): M = %.4f  CI [%.4f, %.4f]  SSR = %.5f  n_recap = %d",
      i, cuts[i], res_i$M_est, res_i$M_ci[1], res_i$M_ci[2],
      res_i$ssr, res_i$n_recap))
    if (i == 1) first_result <- res_i
  }

  cut_profile <- cut_profile[!is.na(cut_profile$M_est), ]

  ## ── 7. Weighted central estimate ─────────────────────────────────────────────
  # Stable plateau: all three criteria
  ci_ok    <- cut_profile$ci_width < quantile(cut_profile$ci_width, 0.75, na.rm = TRUE)
  roll_med <- as.numeric(zoo::rollmedian(cut_profile$M_est, k = 5,
                                         fill = NA, align = "center"))
  stab_ok  <- !is.na(roll_med) &
    abs(cut_profile$M_est - roll_med) / (roll_med + 1e-10) < 0.30
  ssr_ok   <- cut_profile$ssr < quantile(cut_profile$ssr, 0.75, na.rm = TRUE)
  stable_idx <- ci_ok & stab_ok & ssr_ok

  if (sum(stable_idx, na.rm = TRUE) < 3) {
    message("  Plateau: relaxing stability criterion — using CI width + SSR only.")
    stable_idx <- ci_ok & ssr_ok
  }
  if (sum(stable_idx, na.rm = TRUE) < 3) {
    message("  Plateau: insufficient steps — using all steps.")
    stable_idx <- rep(TRUE, nrow(cut_profile))
  }
  message(sprintf("  Plateau: %d of %d steps included.", sum(stable_idx), nrow(cut_profile)))

  ssr_stable <- cut_profile$ssr[stable_idx]
  wts        <- 1 / ssr_stable
  wts        <- wts / sum(wts)
  M_central  <- sum(wts * cut_profile$M_est[stable_idx])

  plateau_df   <- cut_profile[stable_idx, ]
  recap_thresh <- median(plateau_df$n_recap, na.rm = TRUE)
  data_ok      <- plateau_df$n_recap >= recap_thresh
  if (sum(data_ok) < 1) data_ok <- rep(TRUE, nrow(plateau_df))
  best_local   <- which(data_ok)[which.min(plateau_df$ssr[data_ok])]
  best_idx     <- which(stable_idx)[best_local]
  M_ci_central <- c(cut_profile$M_lwr[best_idx], cut_profile$M_upr[best_idx])

  message(sprintf("\nWeighted central M estimate (stable plateau, inv-SSR weighted): %.4f",
                  M_central))
  message(sprintf("CI from best-fitting plateau step (%g mm): %.4f - %.4f",
                  cut_profile$min_cut[best_idx], M_ci_central[1], M_ci_central[2]))

  ## ── 8. Diagnostic plots ───────────────────────────────────────────────────────
  if (plot) {

    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

    # Plot 1: M profile over min_cut
    ylim_M <- range(c(cut_profile$M_lwr, cut_profile$M_upr), na.rm = TRUE)
    plot(cut_profile$min_cut, cut_profile$M_est,
         type = "b", pch = 16, col = "steelblue", lwd = 1.5,
         xlab = "Min size cut (mm)", ylab = "M estimate",
         main = "M Profile over Min Size Cut", ylim = ylim_M)
    polygon(c(cut_profile$min_cut, rev(cut_profile$min_cut)),
            c(cut_profile$M_lwr,   rev(cut_profile$M_upr)),
            col = rgb(0.27, 0.51, 0.71, 0.2), border = NA)
    points(cut_profile$min_cut[stable_idx], cut_profile$M_est[stable_idx],
           pch = 16, col = "darkgreen", cex = 1.2)
    abline(h = M_central, col = "red", lty = 2, lwd = 1.5)
    abline(v = cut_profile$min_cut[best_idx], col = "orange", lty = 3)
    legend("topright",
           legend = c(sprintf("Weighted M = %.3f", M_central),
                      "95% CI ribbon", "Stable plateau", "Best fit step"),
           col    = c("red", rgb(0.27,0.51,0.71,0.4), "darkgreen", "orange"),
           lty    = c(2, 1, NA, 3), pch = c(NA, NA, 16, NA),
           lwd    = c(1.5, 6, NA, 1), bty = "n", cex = 0.8)

    # Plot 2: SSR vs min_cut
    plot(cut_profile$min_cut, cut_profile$ssr,
         type = "b", pch = 16, col = "steelblue", lwd = 1.5,
         xlab = "Min size cut (mm)", ylab = "SSR",
         main = "Model Fit (SSR) vs Min Size Cut")
    points(cut_profile$min_cut[stable_idx], cut_profile$ssr[stable_idx],
           pch = 16, col = "darkgreen", cex = 1.2)
    abline(v = cut_profile$min_cut[best_idx], col = "orange", lty = 3)
    legend("topright", legend = c("Stable plateau", "Best fit step"),
           col = c("darkgreen", "orange"), lty = c(NA, 3), pch = c(16, NA),
           bty = "n", cex = 0.8)

    # Plot 3: Recapture rate by size at starting cut
    obs_start  <- obs[obs$release_len >= min_cut_auto_rnd, ]
    sz_breaks  <- first_result$sz_breaks
    if (length(sz_breaks) < 3)
      sz_breaks <- seq(min(obs_start$release_len), max(obs_start$release_len),
                       length.out = 11)
    size_class    <- cut(obs_start$release_len[seq_along(first_result$ages)],
                         breaks = sz_breaks, include.lowest = TRUE)
    recap_by_size <- tapply(obs_start$isrecap[seq_along(first_result$ages)],
                            size_class, mean, na.rm = TRUE)
    midpoints     <- (sz_breaks[-1] + sz_breaks[-length(sz_breaks)]) / 2
    surv_c        <- exp(-M_central * first_result$ages)
    p_central     <- surv_c / sum(surv_c) * first_result$n_recap
    pred_by_size  <- tapply(p_central, size_class, mean, na.rm = TRUE) * 100

    plot(midpoints, recap_by_size * 100,
         type = "b", pch = 16, col = "steelblue", lwd = 1.5,
         xlab = "Release size (mm CL)", ylab = "Recapture rate (%)",
         main = sprintf("Recapture Rate by Size (cut = %g mm)", min_cut_auto_rnd),
         ylim = c(0, max(recap_by_size * 100, na.rm = TRUE) * 1.25))
    lines(midpoints, pred_by_size, col = "red", lwd = 1.5, lty = 2)
    legend("topright",
           legend = c("Observed", sprintf("Predicted (M = %.3f)", M_central)),
           col = c("steelblue", "red"), lty = c(1, 2), pch = c(16, NA),
           bty = "n", cex = 0.8)

    # Plot 4: NLL profile at starting cut
    plot(first_result$nll_profile$M, first_result$nll_profile$delta_nll,
         type = "l", lwd = 2, col = "steelblue",
         xlab = "M (annual rate)", ylab = "Delta NLL",
         main = sprintf("NLL Profile at Starting Cut (%g mm)", min_cut_auto_rnd),
         ylim = c(0, min(max(first_result$nll_profile$delta_nll), 10)))
    abline(h = 1.92, lty = 2, col = "red")
    abline(v = first_result$M_est, lty = 1, col = "steelblue", lwd = 1.5)
    abline(v = M_central,          lty = 2, col = "darkgreen", lwd = 1.5)

    # Annotate with liberty weighting info
    mtext(sprintf("lambda = %g (liberty weighting)", lambda),
          side = 3, line = 0.2, cex = 0.75, col = "grey40")

    legend("topright",
           legend = c(sprintf("Step MLE = %.3f", first_result$M_est),
                      sprintf("Weighted M = %.3f", M_central),
                      "Chi-sq threshold (1.92)"),
           lty = c(1, 2, 2),
           col = c("steelblue", "darkgreen", "red"),
           bty = "n", cex = 0.8)
  }

  ## ── 9. Return ─────────────────────────────────────────────────────────────────
  invisible(list(
    M_est          = M_central,
    M_ci           = M_ci_central,
    cut_profile    = cut_profile,
    min_cut_auto   = min_cut_auto,
    age_at_release = first_result$ages,
    nll_profile    = first_result$nll_profile
  ))
}
