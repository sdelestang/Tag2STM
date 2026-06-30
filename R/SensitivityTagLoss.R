#' Tag Loss Sensitivity and Survival-to-Liberty Constraint Analysis
#'
#' Two complementary analyses of tag loss uncertainty. Part 1 runs
#' EstimateMBrownie across a grid of fixed tag loss rates. Part 2 is a
#' model-free survival-to-liberty constraint using a 2D grid over
#' (M, tag_loss) with F fixed at the MLE.
#'
#' @param obs Data frame: tag, LCl, relyr, recyr, isrecap.
#' @param brownie_result Output of EstimateMBrownie. If NULL, fitted internally.
#' @param tag_loss_grid Tag loss grid for Part 1. Default seq(0, 0.40, by=0.025).
#' @param tag_loss_ref Reference tag loss if brownie_result is NULL. Default 0.05.
#' @param MLL Minimum legal length (mm). Default 110.
#' @param report_rate Reporting rate. Default 1.
#' @param sigma_F SD on log(F) random walk. Default 0.5.
#' @param M_init Starting M. Default 0.1.
#' @param max_lib Liberty threshold (yr) for survival constraint. Default 6.
#' @param M_grid M values for 2D grid. Default seq(0.01, 0.60, by=0.01).
#' @param tl_grid Tag loss values for 2D grid. Default seq(0.00, 0.60, by=0.01).
#' @param plot Logical. Default TRUE.
#' @return List: sensitivity, nll_ci95, liberty_bound, constraint, results.
#' @export
SensitivityTagLoss <- function(obs,
                               brownie_result = NULL,
                               tag_loss_grid  = seq(0, 0.40, by = 0.025),
                               tag_loss_ref   = 0.05,
                               MLL            = 110,
                               report_rate    = 1,
                               sigma_F        = 0.5,
                               M_init         = 0.1,
                               max_lib        = 6,
                               M_grid         = seq(0.01, 0.60, by = 0.01),
                               tl_grid        = seq(0.00, 0.60, by = 0.01),
                               plot           = TRUE) {

  if (any(tag_loss_grid < 0) || any(tag_loss_grid >= 1))
    stop("tag_loss_grid values must be in [0, 1)")
  tag_loss_grid <- sort(unique(tag_loss_grid))
  n_grid <- length(tag_loss_grid)

  ## ---- Get F MLE -----------------------------------------------------------
  if (is.null(brownie_result)) {
    message("brownie_result not supplied - fitting at tag_loss_ref = ", tag_loss_ref)
    brownie_result <- EstimateMBrownie(obs, MLL = MLL, tag_loss = tag_loss_ref,
                                       report_rate = report_rate, sigma_F = sigma_F, M_init = M_init,
                                       n_boot = 0, plot = FALSE, profile_M = FALSE)
  }
  F_mle <- brownie_result$meanF
  F_lo  <- if (!is.null(brownie_result$meanF_lo95) && !is.na(brownie_result$meanF_lo95))
    brownie_result$meanF_lo95 else F_mle
  F_hi  <- if (!is.null(brownie_result$meanF_hi95) && !is.na(brownie_result$meanF_hi95))
    brownie_result$meanF_hi95 else F_mle
  if (is.na(F_mle)) stop("Could not extract meanF from brownie_result")
  message(sprintf("Using mean F = %.4f (%.4f - %.4f)", F_mle, F_lo, F_hi))

  ## ---- Observed long-liberty proportion ------------------------------------
  obs_sub <- obs[!is.na(obs$LCl) & obs$LCl >= MLL, ]
  n_rel   <- nrow(obs_sub)
  obs_sub$liberty <- obs_sub$recyr - obs_sub$relyr
  n_long   <- sum(obs_sub$isrecap == 1 & !is.na(obs_sub$liberty) &
                    obs_sub$liberty >= max_lib, na.rm = TRUE)
  obs_prop <- n_long / n_rel
  max_obs_lib <- if (any(obs_sub$isrecap == 1 & !is.na(obs_sub$liberty)))
    max(obs_sub$liberty[obs_sub$isrecap == 1 & !is.na(obs_sub$liberty)], na.rm = TRUE)
  else max_lib
  message(sprintf("Observed: %d tags at liberty >= %d yr out of %d (prop = %.5f)",
                  n_long, max_lib, n_rel, obs_prop))
  p_min_vec    <- c(0.001, 0.01, 0.05)
  liberty_bound <- data.frame(min_prob = p_min_vec, max_liberty_yr = max_obs_lib,
                              tag_loss_upper = 1 - p_min_vec^(1 / max_obs_lib))

  ## ---- Part 1: Model sensitivity -------------------------------------------
  message("\n--- Part 1: Model sensitivity ---")
  results  <- vector("list", n_grid)
  sum_rows <- vector("list", n_grid)
  for (i in seq_len(n_grid)) {
    tl <- tag_loss_grid[i]
    message(sprintf("  [%d/%d] tag_loss = %.3f", i, n_grid, tl))
    res <- tryCatch(
      EstimateMBrownie(obs, MLL = MLL, tag_loss = tl, report_rate = report_rate,
                       sigma_F = sigma_F, M_init = M_init, n_boot = 0, plot = FALSE,
                       profile_M = FALSE),
      error = function(e) { message("    FAILED: ", e$message); NULL })
    results[[i]] <- res
    if (!is.null(res)) {
      sum_rows[[i]] <- data.frame(tag_loss = tl, retention_10yr = (1 - tl)^10,
                                  nll = res$nll, delta_nll = NA_real_, M = res$M, M_lo95 = res$M_lo95,
                                  M_hi95 = res$M_hi95, meanF = res$meanF, meanF_lo95 = res$meanF_lo95,
                                  meanF_hi95 = res$meanF_hi95, meanF_cv = res$meanF_cv, converged = TRUE)
    } else {
      sum_rows[[i]] <- data.frame(tag_loss = tl, retention_10yr = (1 - tl)^10,
                                  nll = NA, delta_nll = NA, M = NA, M_lo95 = NA, M_hi95 = NA,
                                  meanF = NA, meanF_lo95 = NA, meanF_hi95 = NA, meanF_cv = NA,
                                  converged = FALSE)
    }
  }
  summ           <- do.call(rbind, sum_rows)
  nll_min        <- min(summ$nll, na.rm = TRUE)
  summ$delta_nll <- summ$nll - nll_min
  chi_thresh     <- 1.92
  in_ci          <- !is.na(summ$delta_nll) & summ$delta_nll <= chi_thresh
  nll_ci95       <- if (any(in_ci))
    c(lo = min(summ$tag_loss[in_ci]), hi = max(summ$tag_loss[in_ci]))
  else c(lo = NA_real_, hi = NA_real_)
  best_tl        <- summ$tag_loss[which.min(summ$nll)]
  summ$boundary_M <- !is.na(summ$M) & summ$M <= 0.005

  ## ---- Part 2: Survival-to-liberty grid ------------------------------------
  message("\n--- Part 2: Survival-to-liberty ---")
  surv_fn  <- function(Mv, tlv, Fv)
    outer(Mv, tlv, function(M, tl) (1 - tl)^max_lib * exp(-(M + Fv) * max_lib))
  P_mat    <- surv_fn(M_grid, tl_grid, F_mle)
  P_mat_lo <- surv_fn(M_grid, tl_grid, F_hi)
  P_mat_hi <- surv_fn(M_grid, tl_grid, F_lo)
  constraint <- list(M_grid = M_grid, tl_grid = tl_grid, P_mat = P_mat,
                     P_mat_lo = P_mat_lo, P_mat_hi = P_mat_hi, obs_prop = obs_prop,
                     n_rel = n_rel, n_long = n_long, F_mle = F_mle, F_lo = F_lo,
                     F_hi = F_hi, max_lib = max_lib)

  ## ---- Plots ---------------------------------------------------------------
  if (plot) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    col_ci  <- "#185FA5"
    col_F   <- "#E05C1A"
    col_ref <- "grey50"
    col_bad <- adjustcolor("firebrick", alpha.f = 0.12)

    ## === Page 1: 2x2 model sensitivity =======================================
    layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))
    par(mar = c(4, 4, 3, 1), oma = c(0, 0, 2.5, 0))
    ok <- !is.na(summ$delta_nll)

    ## 1a: NLL profile
    ylim_nll <- c(0, min(15, max(summ$delta_nll[ok]) * 1.05))
    plot(summ$tag_loss[ok], summ$delta_nll[ok], type = "l", lwd = 2, col = col_ci,
         xlab = "Annual tag loss rate", ylab = expression(Delta ~ NLL),
         main = "Part 1: Likelihood Profile - Tag Loss", ylim = ylim_nll)
    abline(h = chi_thresh, lty = 2, col = "red")
    abline(v = best_tl, col = "black", lwd = 1.5)
    if (!is.na(nll_ci95["lo"])) abline(v = nll_ci95["lo"], lty = 3, col = "grey40")
    if (!is.na(nll_ci95["hi"])) abline(v = nll_ci95["hi"], lty = 3, col = "grey40")
    if (any(summ$boundary_M & ok)) {
      usr <- par("usr")
      rect(min(summ$tag_loss[summ$boundary_M & ok]), usr[3], usr[2], usr[4],
           col = col_bad, border = NA)
      text(mean(c(min(summ$tag_loss[summ$boundary_M & ok]), usr[2])),
           ylim_nll[2] * 0.85, "M -> 0\n(boundary)", cex = 0.65, col = "firebrick")
    }
    legend("topleft",
           legend = c(sprintf("Profile MLE = %.3f", best_tl),
                      sprintf("95%% CI: %.3f-%.3f", nll_ci95["lo"], nll_ci95["hi"]),
                      "Chi-sq (1.92)", "M -> boundary"),
           lty = c(1, 3, 2, NA), fill = c(NA, NA, NA, col_bad),
           col = c("black", "grey40", "red", NA), border = NA, bty = "n", cex = 0.75)

    ## 1b: M vs tag_loss
    ok_M      <- ok & !is.na(summ$M)
    M_hi_vals <- summ$M_hi95[ok_M & !summ$boundary_M]
    M_top     <- if (any(is.finite(M_hi_vals))) max(M_hi_vals, na.rm = TRUE)
    else max(summ$M[ok_M], na.rm = TRUE)
    plot(summ$tag_loss[ok_M], summ$M[ok_M], type = "n",
         xlab = "Annual tag loss rate", ylab = "M (natural mortality)",
         main = "Part 1: M Estimate vs Tag Loss", ylim = c(0, M_top * 1.15))
    if (any(summ$boundary_M)) {
      usr <- par("usr")
      rect(min(summ$tag_loss[summ$boundary_M]), usr[3], usr[2], usr[4],
           col = col_bad, border = NA)
    }
    ci_ok <- ok_M & !is.na(summ$M_lo95) & !is.na(summ$M_hi95) & !summ$boundary_M
    if (any(ci_ok))
      polygon(c(summ$tag_loss[ci_ok], rev(summ$tag_loss[ci_ok])),
              c(summ$M_hi95[ci_ok], rev(summ$M_lo95[ci_ok])),
              col = adjustcolor(col_ci, 0.18), border = NA)
    lines(summ$tag_loss[ok_M], summ$M[ok_M], lwd = 2, col = col_ci)
    abline(v = best_tl, lty = 2, col = col_ref)
    legend("topright", legend = c("M MLE", "M 95% CI", "M -> boundary"),
           lty = c(1, NA, NA), lwd = c(2, NA, NA),
           fill = c(NA, adjustcolor(col_ci, 0.18), col_bad),
           col = c(col_ci, NA, NA), border = NA, bty = "n", cex = 0.75)

    ## 1c: mean F vs tag_loss
    ok_F      <- ok & !is.na(summ$meanF)
    F_hi_vals <- summ$meanF_hi95[ok_F]
    F_top     <- if (any(is.finite(F_hi_vals))) max(F_hi_vals, na.rm = TRUE)
    else max(summ$meanF[ok_F], na.rm = TRUE)
    plot(summ$tag_loss[ok_F], summ$meanF[ok_F], type = "n",
         xlab = "Annual tag loss rate", ylab = "Mean annual F",
         main = "Part 1: Mean F vs Tag Loss", ylim = c(0, F_top * 1.15))
    ci_ok_F <- ok_F & !is.na(summ$meanF_lo95) & !is.na(summ$meanF_hi95)
    if (any(ci_ok_F))
      polygon(c(summ$tag_loss[ci_ok_F], rev(summ$tag_loss[ci_ok_F])),
              c(summ$meanF_hi95[ci_ok_F], rev(summ$meanF_lo95[ci_ok_F])),
              col = adjustcolor(col_F, 0.18), border = NA)
    lines(summ$tag_loss[ok_F], summ$meanF[ok_F], lwd = 2, col = col_F)
    abline(v = best_tl, lty = 2, col = col_ref)
    legend("topleft", legend = c("Mean F MLE", "Mean F 95% CI (delta)"),
           lty = c(1, NA), lwd = c(2, NA), fill = c(NA, adjustcolor(col_F, 0.18)),
           col = c(col_F, NA), border = NA, bty = "n", cex = 0.75)

    ## 1d: M + tag_loss decomposition
    plot(summ$tag_loss[ok_M], summ$M[ok_M] + summ$tag_loss[ok_M],
         type = "l", lwd = 2, col = "grey30", lty = 2,
         xlab = "Annual tag loss rate", ylab = "Rate",
         main = "Part 1: M and Tag Loss Decomposition",
         ylim = c(0, max(summ$M[ok_M] + summ$tag_loss[ok_M],
                         summ$M[ok_M], na.rm = TRUE) * 1.15))
    lines(summ$tag_loss[ok_M], summ$M[ok_M], lwd = 2, col = col_ci)
    lines(summ$tag_loss[ok_M], summ$tag_loss[ok_M], lwd = 1, col = "grey60", lty = 3)
    abline(v = best_tl, lty = 2, col = col_ref)
    legend("topright", legend = c("M", "M + tag_loss", "tag_loss (1:1 ref)"),
           lty = c(1, 2, 3), lwd = c(2, 2, 1),
           col = c(col_ci, "grey30", "grey60"), bty = "n", cex = 0.75)
    mtext(sprintf(
      "Part 1: Model sensitivity  |  MLL=%g mm  sigma_F=%.2f  |  profile MLE=%.3f (%.3f-%.3f)",
      MLL, sigma_F, best_tl, nll_ci95["lo"], nll_ci95["hi"]),
      side = 3, outer = TRUE, line = 1.0, cex = 0.72)

    ## === Page 2: survival-to-liberty constraint ==============================
    dev.new()
    layout(matrix(c(1, 2), nrow = 1))
    par(mar = c(4.5, 4.5, 3.5, 1), oma = c(0, 0, 3, 0))
    M_ref <- brownie_result$M

    ## 2a: 2D image + contours
    image(tl_grid, M_grid, t(P_mat), col = hcl.colors(64, "YlOrRd", rev = TRUE),
          xlab = "Annual tag loss rate", ylab = "M (natural mortality)",
          main = sprintf("P(tag survives >= %d yr) | F = %.3f", max_lib, F_mle))
    contour(tl_grid, M_grid, t(P_mat), levels = obs_prop, add = TRUE, lwd = 2.5,
            col = "black", labels = sprintf("obs p=%.4f", obs_prop), labcex = 0.75)
    if (!isTRUE(all.equal(F_lo, F_mle))) {
      contour(tl_grid, M_grid, t(P_mat_hi), levels = obs_prop, add = TRUE,
              lwd = 1.5, col = "white", lty = 2,
              labels = sprintf("F=%.3f", F_lo), labcex = 0.65)
      contour(tl_grid, M_grid, t(P_mat_lo), levels = obs_prop, add = TRUE,
              lwd = 1.5, col = "white", lty = 2,
              labels = sprintf("F=%.3f", F_hi), labcex = 0.65)
    }
    abline(h = M_ref,        col = col_ci,   lwd = 1.5, lty = 2)
    abline(v = tag_loss_ref, col = "grey80", lwd = 1.5, lty = 2)
    points(tag_loss_ref, M_ref, pch = 21, bg = "white", col = col_ci, cex = 1.5)
    legend("topright",
           legend = c(sprintf("Obs prop contour (F=%.3f)", F_mle),
                      sprintf("F +/- SE (%.3f/%.3f)", F_lo, F_hi),
                      sprintf("Reference fit (tl=%.2f)", tag_loss_ref)),
           lty = c(1, 2, 2), lwd = c(2.5, 1.5, 1.5),
           col = c("black", "white", col_ci), bty = "n", cex = 0.70,
           bg = adjustcolor("grey20", 0.55), text.col = "white")

    ## 2b: Slice at M_ref
    M_lo_ref <- if (!is.null(brownie_result$M_lo95) && !is.na(brownie_result$M_lo95))
      brownie_result$M_lo95 else M_ref
    M_hi_ref <- if (!is.null(brownie_result$M_hi95) && !is.na(brownie_result$M_hi95))
      brownie_result$M_hi95 else M_ref
    P_mle <- (1 - tl_grid)^max_lib * exp(-(M_ref    + F_mle) * max_lib)
    P_Mlo <- (1 - tl_grid)^max_lib * exp(-(M_lo_ref + F_mle) * max_lib)
    P_Mhi <- (1 - tl_grid)^max_lib * exp(-(M_hi_ref + F_mle) * max_lib)
    P_Flo <- (1 - tl_grid)^max_lib * exp(-(M_ref    + F_lo)  * max_lib)
    P_Fhi <- (1 - tl_grid)^max_lib * exp(-(M_ref    + F_hi)  * max_lib)
    ylim_s <- c(0, max(P_Mlo, P_Fhi, na.rm = TRUE) * 1.15)
    plot(tl_grid, P_mle, type = "l", lwd = 2.5, col = col_ci,
         xlab = "Annual tag loss rate",
         ylab = sprintf("P(tag survives >= %d yr)", max_lib),
         main = sprintf("Survival Constraint Slice  (M=%.3f, F=%.3f)", M_ref, F_mle),
         ylim = ylim_s)
    polygon(c(tl_grid, rev(tl_grid)), c(P_Mlo, rev(P_Mhi)),
            col = adjustcolor(col_ci, 0.15), border = NA)
    polygon(c(tl_grid, rev(tl_grid)), c(P_Flo, rev(P_Fhi)),
            col = adjustcolor(col_F, 0.15), border = NA)
    lines(tl_grid, P_mle, lwd = 2.5, col = col_ci)
    abline(h = obs_prop, lwd = 2, col = "black")
    cross_idx <- which(diff(P_mle < obs_prop) != 0)
    if (length(cross_idx) > 0) {
      tl_cross <- tl_grid[cross_idx[1]]
      abline(v = tl_cross, lwd = 1.5, col = "black", lty = 3)
      usr <- par("usr")
      rect(tl_cross, usr[3], usr[2], usr[4],
           col = adjustcolor("firebrick", 0.08), border = NA)
      text(tl_cross + diff(usr[1:2]) * 0.015, ylim_s[2] * 0.92,
           sprintf("max plausible\ntag_loss ~ %.2f", tl_cross),
           adj = 0, cex = 0.72, col = "firebrick")
    }
    abline(v = tag_loss_ref, lty = 2, col = col_ref)
    legend("topright",
           legend = c(sprintf("P at M=%.3f, F=%.3f", M_ref, F_mle),
                      "M 95% CI band", "F 95% CI band",
                      sprintf("Observed: %d/%d = %.4f", n_long, n_rel, obs_prop),
                      "Max plausible tag_loss", "Implausible region"),
           lty  = c(1, NA, NA, 1, 3, NA), lwd = c(2.5, NA, NA, 2, 1.5, NA),
           fill = c(NA, adjustcolor(col_ci, 0.15), adjustcolor(col_F, 0.15),
                    NA, NA, adjustcolor("firebrick", 0.08)),
           col  = c(col_ci, NA, NA, "black", "black", NA),
           border = NA, bty = "n", cex = 0.72)
    mtext(sprintf(
      "Part 2: Survival-to-liberty  |  threshold=%d yr  |  F=%.3f (%.3f-%.3f)  |  obs=%.4f",
      max_lib, F_mle, F_lo, F_hi, obs_prop),
      side = 3, outer = TRUE, line = 1.0, cex = 0.70)
  }

  invisible(list(sensitivity = summ, nll_ci95 = nll_ci95,
                 liberty_bound = liberty_bound, constraint = constraint,
                 results = results))
}
