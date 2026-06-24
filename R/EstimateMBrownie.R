#' Estimate Natural Mortality Using Brownie Dead-Recovery Model
#'
#' Estimates natural mortality (M) from tag-recapture data using a Brownie
#' dead-recovery model. Animals are tracked through time with year-specific
#' fishing mortality derived from relative catch rates. Only animals released
#' at or above the minimum legal length are used, ensuring constant gear
#' selectivity. Likelihood profiling is used to derive confidence intervals
#' for M.
#'
#' @param obs Data frame of tagged animals containing columns:
#'   \itemize{
#'     \item \code{tag} - unique tag identifier
#'     \item \code{LCl} - carapace length at release (mm)
#'     \item \code{relyr} - year of release (numeric)
#'     \item \code{recyr} - year of recapture (numeric, NA if not recaptured)
#'     \item \code{isrecap} - 1 if recaptured, 0 if not
#'   }
#'   Multiple recaptures of the same tag are treated as independent events.
#'
#' @param catch Data frame of annual commercial catch containing columns:
#'   \itemize{
#'     \item \code{Year} - fishing year (numeric)
#'     \item \code{Catch} - total catch in weight (kg or tonnes, consistent units)
#'   }
#'   Used to derive relative annual fishing mortality. Absolute scale does not
#'   matter as only relative year-to-year patterns are used.
#'
#' @param MLL Numeric. Minimum legal length (mm). Only animals with
#'   \code{LCl >= MLL} at release are included in the analysis. Smaller animals
#'   are silently dropped. Default is 100 mm.
#'
#' @param tag_loss Numeric between 0 and 1. Annual tag loss rate, assumed known
#'   from double-tagging experiments. Applied as a fixed correction to survival
#'   each year. Default is 0.05 (5% per year).
#'
#' @param report_rate Numeric between 0 and 1. Probability that a recaptured
#'   animal is reported. Assumed constant across years. Cannot be estimated
#'   separately from F without additional data; fix at 1 unless information
#'   is available. Default is 1.
#'
#' @param plot Logical. If TRUE (default), produces diagnostic plots including
#'   the likelihood profile for M, observed vs predicted recovery rates by year,
#'   and residuals by release year.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{M} - MLE estimate of natural mortality
#'     \item \code{M_lo95} - lower 95\% profile confidence interval for M
#'     \item \code{M_hi95} - upper 95\% profile confidence interval for M
#'     \item \code{F_base} - estimated average fishing mortality
#'     \item \code{F_yr} - named vector of annual fishing mortality estimates
#'     \item \code{Z_yr} - named vector of annual total mortality (M + F_yr)
#'     \item \code{nll} - negative log-likelihood at MLE
#'     \item \code{n_legal} - number of animals included (released >= MLL)
#'     \item \code{n_recap} - number of recapture events used
#'     \item \code{pred_recap} - data frame of predicted vs observed recovery
#'       rates by release year and years at liberty
#'     \item \code{profile} - data frame of M values and profile NLL for plotting
#'   }
#'
#' @details
#' **Model Structure:**
#'
#' For each legal-sized tagged animal \eqn{i} released in year \eqn{r}:
#'
#' Annual total mortality: \eqn{Z_t = M + F_t}
#'
#' where \eqn{F_t = F_{base} \times C_t / \bar{C}} and \eqn{C_t} is catch in year \eqn{t}.
#'
#' Effective annual survival (accounting for tag loss):
#' \deqn{S_t = \exp(-Z_t) \times (1 - \delta)}
#' where \eqn{\delta} is the annual tag loss rate.
#'
#' Probability of first recapture in year \eqn{t} (years at liberty = \eqn{t - r}):
#' \deqn{p_{recap,t} = \frac{F_t}{Z_t}(1 - \exp(-Z_t)) \times \rho \times \prod_{s=r}^{t-1} S_s}
#' where \eqn{\rho} is the reporting rate.
#'
#' Probability of never being recaptured:
#' \deqn{p_{never} = \prod_{t=r}^{T} (1 - p_{recap,t} / S_{cumul,t-1})}
#'
#' Log-likelihood:
#' \deqn{\ell = \sum_{recaptured} \log(p_{recap,t}) + \sum_{not recaptured} \log(p_{never})}
#'
#' **Multiple recaptures:** Each recapture event for a given tag is treated as
#' an independent observation. This is conservative (slightly inflates effective
#' sample size) but avoids the complexity of conditioning on previous recaptures.
#'
#' **Identifiability:** M and F_base are identified from the temporal pattern
#' of recaptures — specifically from the ratio of early to late recoveries.
#' Years with high catch give high F_rel, accelerating the decline in recapture
#' probability, which separates M from F. At least 3-4 years of recapture data
#' and meaningful year-to-year variation in catch are needed for reliable
#' separation.
#'
#' **Profile CI:** The 95\% confidence interval is computed by profiling the
#' likelihood over a grid of M values, with F_base re-optimised at each M.
#' The CI boundary is where the profile NLL exceeds the MLE NLL by
#' \eqn{\chi^2_{1,0.95}/2 = 1.92}.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' result <- EstimateMBrownie(obs, catch, MLL = 100, tag_loss = 0.05)
#' cat("M estimate:", result$M, "(", result$M_lo95, "-", result$M_hi95, ")\n")
#'
#' # Sensitivity to tag loss
#' sens_tl <- lapply(c(0, 0.05, 0.10, 0.15), function(tl) {
#'   r <- EstimateMBrownie(obs, catch, MLL = 100, tag_loss = tl, plot = FALSE)
#'   data.frame(tag_loss = tl, M = r$M, lo = r$M_lo95, hi = r$M_hi95)
#' })
#' do.call(rbind, sens_tl)
#'
#' # Sensitivity to MLL
#' sens_mll <- lapply(c(90, 95, 100, 105, 110), function(mll) {
#'   r <- EstimateMBrownie(obs, catch, MLL = mll, tag_loss = 0.05, plot = FALSE)
#'   data.frame(MLL = mll, n = r$n_legal, M = r$M, lo = r$M_lo95, hi = r$M_hi95)
#' })
#' do.call(rbind, sens_mll)
#'
#' # Sensitivity to reporting rate
#' sens_rr <- lapply(c(0.7, 0.8, 0.9, 1.0), function(rr) {
#'   r <- EstimateMBrownie(obs, catch, MLL = 100, tag_loss = 0.05,
#'                         report_rate = rr, plot = FALSE)
#'   data.frame(report_rate = rr, M = r$M, lo = r$M_lo95, hi = r$M_hi95)
#' })
#' do.call(rbind, sens_rr)
#' }
#'
#' @importFrom RTMB MakeADFun
#' @importFrom stats nlminb optimise
#' @importFrom graphics plot lines points abline legend par mtext
#' @export
EstimateMBrownie <- function(obs,
                              catch,
                              MLL         = 100,
                              tag_loss    = 0.05,
                              report_rate = 1,
                              plot        = TRUE) {

  # ── 0. Input checks ────────────────────────────────────────────────────────
  required_obs <- c("tag", "LCl", "relyr", "recyr", "isrecap")
  missing_obs  <- setdiff(required_obs, names(obs))
  if (length(missing_obs) > 0)
    stop("obs is missing columns: ", paste(missing_obs, collapse = ", "))

  required_catch <- c("Year", "Catch")
  missing_catch  <- setdiff(required_catch, names(catch))
  if (length(missing_catch) > 0)
    stop("catch is missing columns: ", paste(missing_catch, collapse = ", "))

  if (tag_loss < 0 || tag_loss >= 1)
    stop("tag_loss must be in [0, 1)")
  if (report_rate <= 0 || report_rate > 1)
    stop("report_rate must be in (0, 1]")

  # ── 1. Filter to legal-sized animals ───────────────────────────────────────
  obs_legal <- obs[!is.na(obs$LCl) & obs$LCl >= MLL, ]
  n_legal   <- nrow(obs_legal)
  n_recap   <- sum(obs_legal$isrecap, na.rm = TRUE)

  if (n_legal == 0)
    stop("No animals with LCl >= MLL = ", MLL, ". Check MLL or LCl column.")
  if (n_recap == 0)
    stop("No recaptured animals with LCl >= MLL = ", MLL,
         ". Cannot estimate M.")

  message(sprintf("Using %d tagged animals (%d recaptured) with LCl >= %g mm",
                  n_legal, n_recap, MLL))

  # ── 2. Build year index ────────────────────────────────────────────────────
  all_years  <- sort(unique(c(obs_legal$relyr,
                               obs_legal$recyr[!is.na(obs_legal$recyr)],
                               catch$Year)))
  yr_min     <- min(all_years)
  yr_max     <- max(all_years)
  years      <- yr_min:yr_max
  nyears     <- length(years)

  # ── 3. Relative F from catch ───────────────────────────────────────────────
  # Fill missing years with mean catch (conservative — low F in gap years)
  catch_vec       <- rep(mean(catch$Catch), nyears)
  names(catch_vec) <- years
  for (i in seq_len(nrow(catch))) {
    yr <- as.character(catch$Year[i])
    if (yr %in% names(catch_vec))
      catch_vec[yr] <- catch$Catch[i]
  }
  F_rel <- catch_vec / mean(catch_vec)   # relative F; mean = 1

  # ── 4. Build recovery array ────────────────────────────────────────────────
  # For each animal: release year index, recap year index (NA if not recapped)
  obs_legal$rel_idx <- match(obs_legal$relyr,  years)
  obs_legal$rec_idx <- match(obs_legal$recyr,  years)

  # ── 5. Negative log-likelihood function ───────────────────────────────────
  nll_fn <- function(log_M, log_F_base) {

    M      <- exp(log_M)
    F_base <- exp(log_F_base)
    F_yr   <- F_base * F_rel          # annual F vector

    nll <- 0

    for (i in seq_len(nrow(obs_legal))) {

      rel_t  <- obs_legal$rel_idx[i]
      rec_t  <- obs_legal$rec_idx[i]
      recap  <- obs_legal$isrecap[i]

      if (is.na(rel_t)) next

      # Maximum year to track this animal
      max_t <- ifelse(is.na(rec_t), nyears, rec_t)

      S_cumul       <- 1    # cumulative survival to start of year t
      p_never_recap <- 1    # probability of never being recaptured
      found         <- FALSE

      for (t in rel_t:max_t) {

        Ft  <- F_yr[t]
        Zt  <- M + Ft
        # Survival this year including tag loss
        St  <- exp(-Zt) * (1 - tag_loss)

        # Probability of recapture in year t (given alive at start of t)
        # = fishing mortality component * reporting rate
        p_recap_yr <- (Ft / Zt) * (1 - exp(-Zt)) * report_rate

        # Probability of recapture in year t (unconditional)
        p_recap_uncond <- p_recap_yr * S_cumul

        if (recap == 1 && !is.na(rec_t) && t == rec_t) {
          # This is the recapture year — add log probability
          nll  <- nll - log(p_recap_uncond + 1e-15)
          found <- TRUE
          break
        }

        # Update probability of never being recaptured
        # Conditional on being alive at start of t: probability of NOT
        # being recaptured this year = 1 - p_recap_yr
        p_never_recap <- p_never_recap * (1 - p_recap_yr)

        # Update cumulative survival
        S_cumul <- S_cumul * St
      }

      # Non-recaptured contribution
      if (recap == 0 && !found) {
        nll <- nll - log(p_never_recap * S_cumul + 1e-15)
      }
    }

    return(nll)
  }

  # ── 6. Optimise jointly over M and F_base ─────────────────────────────────
  # Grid search for starting values
  grid_M      <- seq(log(0.05), log(0.8), length.out = 12)
  grid_F      <- seq(log(0.05), log(1.0), length.out = 12)
  best_nll    <- Inf
  best_start  <- c(log(0.15), log(0.2))

  for (gm in grid_M) {
    for (gf in grid_F) {
      tryCatch({
        v <- nll_fn(gm, gf)
        if (is.finite(v) && v < best_nll) {
          best_nll   <- v
          best_start <- c(gm, gf)
        }
      }, error = function(e) NULL)
    }
  }

  # Full optimisation from best grid start
  opt <- nlminb(best_start,
                function(p) nll_fn(p[1], p[2]),
                lower = c(log(0.001), log(0.001)),
                upper = c(log(5.0),   log(5.0)),
                control = list(iter.max = 1000, eval.max = 2000))

  M_mle      <- exp(opt$par[1])
  F_base_mle <- exp(opt$par[2])
  F_yr_mle   <- F_base_mle * F_rel
  Z_yr_mle   <- M_mle + F_yr_mle
  nll_mle    <- opt$objective

  message(sprintf("MLE: M = %.4f, F_base = %.4f, NLL = %.4f",
                  M_mle, F_base_mle, nll_mle))

  # ── 7. Likelihood profile for M ────────────────────────────────────────────
  # Profile over a grid of M values, re-optimising F_base at each
  M_grid   <- exp(seq(log(max(0.01, M_mle * 0.1)),
                      log(min(2.0,  M_mle * 10)),
                      length.out = 100))
  prof_nll <- numeric(length(M_grid))

  for (k in seq_along(M_grid)) {
    fixed_logM <- log(M_grid[k])
    opt_k <- tryCatch(
      optimise(function(lF) nll_fn(fixed_logM, lF),
               interval = c(log(0.001), log(5)),
               tol = 1e-6),
      error = function(e) list(objective = Inf)
    )
    prof_nll[k] <- opt_k$objective
  }

  # Chi-squared threshold for 95% CI: delta NLL = 1.92
  chi_thresh <- 1.92
  in_ci      <- (prof_nll - nll_mle) <= chi_thresh

  M_lo95 <- if (any(in_ci)) min(M_grid[in_ci]) else NA_real_
  M_hi95 <- if (any(in_ci)) max(M_grid[in_ci]) else NA_real_

  profile_df <- data.frame(M       = M_grid,
                            nll     = prof_nll,
                            delta   = prof_nll - nll_mle,
                            in_ci95 = in_ci)

  # ── 8. Predicted vs observed recovery rates ────────────────────────────────
  # Compute predicted recovery probability by years-at-liberty, pooled over
  # release years (for diagnostic plot)
  max_lib   <- nyears
  pred_rows <- list()

  for (rel_t in unique(obs_legal$rel_idx)) {

    rel_yr    <- years[rel_t]
    n_rel     <- sum(obs_legal$rel_idx == rel_t, na.rm = TRUE)
    if (n_rel == 0) next

    S_cumul <- 1
    for (lib in 1:(nyears - rel_t + 1)) {
      t   <- rel_t + lib - 1
      if (t > nyears) break

      Ft  <- F_yr_mle[t]
      Zt  <- M_mle + Ft
      St  <- exp(-Zt) * (1 - tag_loss)
      p_r <- (Ft / Zt) * (1 - exp(-Zt)) * report_rate * S_cumul

      # Observed recoveries at this liberty time
      n_obs_recap <- sum(obs_legal$rel_idx == rel_t &
                         obs_legal$rec_idx == t &
                         obs_legal$isrecap == 1,
                         na.rm = TRUE)

      pred_rows[[length(pred_rows) + 1]] <- data.frame(
        rel_yr      = rel_yr,
        years_lib   = lib,
        n_released  = n_rel,
        n_recap_obs = n_obs_recap,
        rate_obs    = n_obs_recap / n_rel,
        rate_pred   = p_r
      )
      S_cumul <- S_cumul * St
    }
  }

  pred_df <- do.call(rbind, pred_rows)

  # ── 9. Diagnostic plots ────────────────────────────────────────────────────
  if (plot) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

    # (a) Likelihood profile
    plot(profile_df$M, profile_df$delta,
         type = "l", lwd = 2, col = "#185FA5",
         xlab = "M (natural mortality)",
         ylab = expression(Delta ~ "NLL"),
         main = "Likelihood Profile for M",
         ylim = c(0, max(min(profile_df$delta[is.finite(profile_df$delta)],
                             15), chi_thresh * 2)))
    abline(h = chi_thresh, lty = 2, col = "red")
    abline(v = M_mle,    lty = 1, col = "black", lwd = 1.5)
    if (!is.na(M_lo95)) abline(v = M_lo95, lty = 3, col = "grey40")
    if (!is.na(M_hi95)) abline(v = M_hi95, lty = 3, col = "grey40")
    legend("topright",
           legend = c(sprintf("MLE = %.3f", M_mle),
                      sprintf("95%% CI: %.3f - %.3f", M_lo95, M_hi95),
                      "Chi-sq threshold"),
           lty = c(1, 3, 2),
           col = c("black", "grey40", "red"),
           bty = "n", cex = 0.85)

    # (b) Observed vs predicted recovery rates by years at liberty
    lib_summary <- aggregate(cbind(n_recap_obs, n_released, rate_pred) ~
                               years_lib,
                             data = pred_df, FUN = sum)
    lib_summary$rate_obs  <- lib_summary$n_recap_obs / lib_summary$n_released
    lib_summary$rate_pred <- aggregate(rate_pred ~ years_lib,
                                       data = pred_df, FUN = mean)$rate_pred

    ylim_b <- range(c(lib_summary$rate_obs, lib_summary$rate_pred), na.rm = TRUE)
    plot(lib_summary$years_lib, lib_summary$rate_obs,
         pch = 16, col = "black",
         xlab = "Years at liberty",
         ylab = "Recovery rate",
         main = "Observed vs Predicted Recovery",
         ylim = c(0, max(ylim_b) * 1.1))
    lines(lib_summary$years_lib, lib_summary$rate_pred,
          col = "#185FA5", lwd = 2)
    legend("topright",
           legend = c("Observed", "Predicted"),
           pch    = c(16, NA),
           lty    = c(NA, 1),
           col    = c("black", "#185FA5"),
           bty    = "n", cex = 0.85)

    # (c) Residuals by release year
    rel_yr_res <- pred_df[pred_df$n_recap_obs > 0 | pred_df$rate_pred > 0.001, ]
    rel_yr_res$resid <- rel_yr_res$rate_obs - rel_yr_res$rate_pred
    if (nrow(rel_yr_res) > 0) {
      plot(rel_yr_res$rel_yr, rel_yr_res$resid,
           pch = 16, col = rgb(0.1, 0.4, 0.7, 0.5),
           xlab = "Release year",
           ylab = "Obs - Pred recovery rate",
           main = "Residuals by Release Year")
      abline(h = 0, lty = 2, col = "red")
    }

    # (d) Annual F and M
    yr_df <- data.frame(year = years,
                        F    = F_yr_mle,
                        M    = M_mle,
                        Z    = Z_yr_mle)
    ylim_d <- c(0, max(yr_df$Z) * 1.2)
    plot(yr_df$year, yr_df$Z,
         type = "l", lwd = 2, col = "black",
         xlab = "Year", ylab = "Mortality rate",
         main = "Annual Mortality Components",
         ylim = ylim_d)
    lines(yr_df$year, yr_df$F, col = "#E05C1A", lwd = 1.5, lty = 2)
    abline(h = M_mle, col = "#185FA5", lwd = 1.5, lty = 3)
    legend("topright",
           legend = c("Z (total)", "F (fishing)",
                      sprintf("M = %.3f", M_mle)),
           lty  = c(1, 2, 3),
           col  = c("black", "#E05C1A", "#185FA5"),
           lwd  = c(2, 1.5, 1.5),
           bty  = "n", cex = 0.85)

    mtext(sprintf(
      "EstimateMBrownie: MLL=%g mm, tag_loss=%.2f, report_rate=%.2f  |  n=%d released, %d recaptured",
      MLL, tag_loss, report_rate, n_legal, n_recap),
      side = 3, outer = TRUE, line = -1.2, cex = 0.75)
  }

  # ── 10. Return results ─────────────────────────────────────────────────────
  invisible(list(
    M           = M_mle,
    M_lo95      = M_lo95,
    M_hi95      = M_hi95,
    F_base      = F_base_mle,
    F_yr        = setNames(F_yr_mle,  years),
    Z_yr        = setNames(Z_yr_mle,  years),
    nll         = nll_mle,
    n_legal     = n_legal,
    n_recap     = n_recap,
    MLL         = MLL,
    tag_loss    = tag_loss,
    report_rate = report_rate,
    pred_recap  = pred_df,
    profile     = profile_df
  ))
}
