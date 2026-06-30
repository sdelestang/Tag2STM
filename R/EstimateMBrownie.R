#' Estimate Natural Mortality Using Brownie Dead-Recovery Model
#'
#' Estimates natural mortality (M) from tag-recapture data using a Brownie
#' dead-recovery model. Follows a matrix-based formulation where rows are
#' release cohorts (by year) and columns are recapture years, with a final
#' 'nsa' (not seen again) column. Year-specific fishing mortality is estimated
#' directly from the recovery pattern — no catch or effort data required.
#' Multiple release cohorts sharing the same F in a given year provide the
#' cross-cohort signal that identifies M. Delta-method uncertainty is computed
#' for all F quantities from the joint Hessian at the MLE.
#'
#' @param obs Data frame of tagged animals containing columns:
#'   \itemize{
#'     \item \code{tag}     - unique tag identifier
#'     \item \code{LCl}     - carapace length at release (mm)
#'     \item \code{relyr}   - year of release (numeric)
#'     \item \code{recyr}   - year of recapture (numeric, NA if not recaptured)
#'     \item \code{isrecap} - 1 if recaptured, 0 if not
#'   }
#'   Multiple recaptures of the same tag are treated as independent events.
#'
#' @param MLL Numeric. Minimum legal length (mm). Only animals with
#'   LCl >= MLL at release are included. Default 110 mm.
#'
#' @param tag_loss Numeric matrix or scalar. Annual tag loss rate(s).
#'   If scalar, applied uniformly. If a matrix with rows = release cohorts
#'   and cols = recapture years, allows tag loss to vary by cohort and time.
#'   Fixed from double-tagging experiments. Default 0.05.
#'
#' @param report_rate Numeric vector or scalar. Reporting rate by recapture
#'   year (phi in Brownie formulation). If scalar, applied uniformly across
#'   years. Default 1.
#'
#' @param sigma_F Numeric. SD of random-walk penalty on log(F) between
#'   consecutive years. Controls smoothness of F trajectory.
#'   Lower = smoother F (more signal attributed to M).
#'   Higher = F varies more freely.
#'   Set to Inf to remove the penalty entirely. Default 0.5.
#'
#' @param M_init Numeric. Starting value for M in optimisation. Default 0.1.
#'
#' @param n_boot Integer. Number of bootstrap resamples for M uncertainty.
#'   Set to 0 to skip bootstrapping and use likelihood profiling only.
#'   Default 0.
#'
#' @param plot Logical. Produce diagnostic plots. Default TRUE.
#'
#' @param profile_M Logical. Compute likelihood profile CI for M? This is the
#'   slowest part of the function (40 inner optimisations). Set to FALSE when
#'   calling inside a sensitivity loop (e.g. \code{\link{SensitivityTagLoss}})
#'   where only the MLE and delta-method F uncertainty are needed. Default TRUE.
#'
#' @return Invisibly returns a list with:
#'   \itemize{
#'     \item \code{M}           - MLE of natural mortality (annual rate)
#'     \item \code{M_lo95}      - lower 95\% profile CI for M
#'     \item \code{M_hi95}      - upper 95\% profile CI for M
#'     \item \code{M_boot}      - bootstrap distribution of M (if n_boot > 0)
#'     \item \code{F_yr}        - named vector of annual F estimates
#'     \item \code{F_lo95}      - named vector of lower 95\% delta-method CI for F
#'     \item \code{F_hi95}      - named vector of upper 95\% delta-method CI for F
#'     \item \code{F_cv}        - named vector of CV for each annual F estimate
#'     \item \code{meanF}       - arithmetic mean of annual F estimates
#'     \item \code{meanF_lo95}  - lower 95\% delta-method CI for mean F
#'     \item \code{meanF_hi95}  - upper 95\% delta-method CI for mean F
#'     \item \code{meanF_cv}    - CV of mean F estimate
#'     \item \code{Z_yr}        - named vector of annual Z = M + F
#'     \item \code{nll}         - NLL at MLE
#'     \item \code{n_legal}     - number of animals included (LCl >= MLL)
#'     \item \code{n_recap}     - number of recapture events
#'     \item \code{obs_mat}     - observed recovery matrix (cohort x year + nsa)
#'     \item \code{pred_mat}    - predicted recovery matrix (expected counts)
#'     \item \code{profile}     - likelihood profile data frame for M
#'   }
#'
#' @details
#' **Likelihood (following Brownie et al. 1985):**
#'
#' For release cohort r recaptured in year t:
#' \deqn{p_{r,t} = \phi_t \cdot \lambda_{r,t} \cdot u(F_t) \cdot
#'   \exp\left(-\sum_{s=r}^{t-1} F_s - (t-r) \cdot M\right)}
#'
#' where \eqn{u(F) = \frac{F}{F+M}(1-e^{-(F+M)})} is the probability of
#' dying from fishing, \eqn{\phi_t} is reporting rate, and
#' \eqn{\lambda_{r,t}} is tag retention probability.
#'
#' Probability of not being seen again:
#' \deqn{p_{r,nsa} = 1 - \sum_{t \geq r} p_{r,t}}
#'
#' Log-likelihood:
#' \deqn{\ell = \sum_{r,t} n_{r,t} \log(p_{r,t}) +
#'              \sum_r n_{r,nsa} \log(p_{r,nsa})}
#'
#' **Identification of M:** In any year t, cohorts r1 < r2 both contribute
#' recaptures. They share F_t but differ in cumulative survival. The ratio
#' of their expected recovery rates is:
#' \deqn{\frac{p_{r1,t}}{p_{r2,t}} = \exp\left(-\sum_{s=r1}^{r2-1}(F_s+M)\right)}
#' which pins down M once F is accounted for.
#'
#' **F uncertainty (delta method):** The Hessian of the joint NLL with respect
#' to \eqn{(\log M, \log F_1, \ldots, \log F_T)} is computed via
#' \code{optimHess()} at the MLE and inverted to give the parameter covariance
#' matrix. Per-year 95\% CIs are \eqn{F_t \exp(\pm 1.96\, \widehat{SE}_{\log F_t})},
#' which are asymmetric and always positive. The CV on the natural scale is
#' \eqn{\sqrt{\exp(\widehat{SE}^2_{\log F_t}) - 1}}. Mean F uncertainty is
#' propagated via the gradient \eqn{g_t = F_t / T} applied to the log-F
#' covariance submatrix, then log-transformed for the CI.
#' Note that because the Hessian includes the \code{sigma_F} random-walk
#' penalty, the F CIs reflect both data information and the smoothness prior.
#'
#' @references
#' Brownie, C., Anderson, D.R., Burnham, K.P. & Robson, D.S. (1985).
#' \emph{Statistical Inference from Band Recovery Data: A Handbook}, 2nd ed.
#' US Fish and Wildlife Service Resource Publication 156.
#'
#' @examples
#' \dontrun{
#' result <- EstimateMBrownie(obs, MLL = 110, tag_loss = 0.05)
#' cat("M:", result$M, "(", result$M_lo95, "-", result$M_hi95, ")\n")
#' cat("Mean F:", result$meanF, "(", result$meanF_lo95, "-", result$meanF_hi95,
#'     ")  CV =", round(result$meanF_cv, 3), "\n")
#'
#' # Access per-year F table
#' data.frame(year  = as.integer(names(result$F_yr)),
#'            F     = round(result$F_yr,   3),
#'            lo95  = round(result$F_lo95, 3),
#'            hi95  = round(result$F_hi95, 3),
#'            cv    = round(result$F_cv,   3))
#'
#' # Sensitivity to sigma_F
#' lapply(c(0.2, 0.5, 1.0, Inf), function(sf) {
#'   r <- EstimateMBrownie(obs, MLL = 110, tag_loss = 0.05, sigma_F = sf,
#'                         plot = FALSE)
#'   data.frame(sigma_F = sf, M = round(r$M, 3),
#'              M_lo = round(r$M_lo95, 3), M_hi = round(r$M_hi95, 3),
#'              meanF = round(r$meanF, 3),
#'              meanF_lo = round(r$meanF_lo95, 3),
#'              meanF_hi = round(r$meanF_hi95, 3))
#' }) |> do.call(what = rbind)
#'
#' # Sensitivity to tag loss
#' lapply(c(0, 0.05, 0.10, 0.20), function(tl) {
#'   r <- EstimateMBrownie(obs, MLL = 110, tag_loss = tl, plot = FALSE)
#'   data.frame(tag_loss = tl, M = round(r$M, 3),
#'              lo = round(r$M_lo95, 3), hi = round(r$M_hi95, 3))
#' }) |> do.call(what = rbind)
#'
#' # Sensitivity to reporting rate
#' lapply(c(0.7, 0.8, 0.9, 1.0), function(rr) {
#'   r <- EstimateMBrownie(obs, MLL = 110, tag_loss = 0.05, report_rate = rr,
#'                         plot = FALSE)
#'   data.frame(report_rate = rr, M = round(r$M, 3),
#'              lo = round(r$M_lo95, 3), hi = round(r$M_hi95, 3))
#' }) |> do.call(what = rbind)
#' }
#'
#' @export
EstimateMBrownie <- 
function (obs, MLL = 110, tag_loss = 0.05, report_rate = 1, sigma_F = 0.5, 
          M_init = 0.1, n_boot = 0, plot = TRUE, profile_M = TRUE) 
{
  required <- c("tag", "LCl", "relyr", "recyr", "isrecap")
  missing <- setdiff(required, names(obs))
  if (length(missing) > 0) 
    stop("obs missing columns: ", paste(missing, collapse = ", "))
  if (any(tag_loss < 0) || any(tag_loss >= 1)) 
    stop("tag_loss must be in [0, 1)")
  if (any(report_rate <= 0) || any(report_rate > 1)) 
    stop("report_rate must be in (0, 1]")
  obs <- obs[!is.na(obs$LCl) & obs$LCl >= MLL, ]
  n_legal <- nrow(obs)
  n_recap <- sum(obs$isrecap, na.rm = TRUE)
  if (n_legal == 0) 
    stop("No animals with LCl >= MLL = ", MLL)
  if (n_recap == 0) 
    stop("No recaptured animals with LCl >= MLL = ", MLL)
  message(sprintf("Using %d animals (%d recaptured) with LCl >= %g mm", 
                  n_legal, n_recap, MLL))
  rel_years <- sort(unique(obs$relyr))
  rec_years <- sort(unique(obs$recyr[!is.na(obs$recyr) & obs$isrecap == 
                                       1]))
  years <- sort(unique(c(rel_years, rec_years)))
  yr_min <- min(years)
  yr_max <- max(years)
  years <- yr_min:yr_max
  nyears <- length(years)
  nrel <- length(rel_years)
  yr_idx <- function(y) match(y, years)
  ry_idx <- function(y) match(y, rel_years)
  obs_mat <- matrix(0L, nrow = nrel, ncol = nyears + 1L, dimnames = list(rel_years, 
                                                                         c(years, "nsa")))
  n_rel_vec <- setNames(integer(nrel), rel_years)
  for (ry in rel_years) {
    sub <- obs[obs$relyr == ry, ]
    ri <- ry_idx(ry)
    n_rel_vec[as.character(ry)] <- nrow(sub)
    recs <- sub[sub$isrecap == 1 & !is.na(sub$recyr), ]
    for (i in seq_len(nrow(recs))) {
      cy <- as.character(recs$recyr[i])
      if (cy %in% colnames(obs_mat)) 
        obs_mat[ri, cy] <- obs_mat[ri, cy] + 1L
    }
    obs_mat[ri, "nsa"] <- nrow(sub) - sum(obs_mat[ri, as.character(years)])
  }
  phi_mat <- matrix(1, nrow = nrel, ncol = nyears, dimnames = list(rel_years, 
                                                                   years))
  lambda_mat <- matrix(1, nrow = nrel, ncol = nyears, dimnames = list(rel_years, 
                                                                      years))
  phi_vec <- if (length(report_rate) == 1) 
    rep(report_rate, nyears)
  else {
    if (length(report_rate) != nyears) 
      stop("report_rate must be scalar or length nyears = ", 
           nyears)
    report_rate
  }
  for (ri in seq_len(nrel)) phi_mat[ri, ] <- phi_vec
  if (is.matrix(tag_loss)) {
    if (!identical(dim(tag_loss), c(nrel, nyears))) 
      stop("tag_loss matrix must be nrel x nyears = ", 
           nrel, " x ", nyears)
    lambda_mat <- tag_loss
  }
  else {
    tl_annual <- tag_loss[1]
    for (ri in seq_len(nrel)) {
      rel_t <- yr_idx(rel_years[ri])
      for (t in seq_len(nyears)) {
        lib <- t - rel_t
        lambda_mat[ri, t] <- if (lib < 0) 
          0
        else (1 - tl_annual)^lib
      }
    }
  }
  u_fm <- function(F, M) {
    val <- (F/(F + M)) * (1 - exp(-(F + M)))
    pmax(val, 1e-10)
  }
  nll_fn <- function(log_M, log_F) {
    M <- exp(log_M)
    F_t <- exp(log_F)
    pred <- matrix(0, nrow = nrel, ncol = nyears, dimnames = list(rel_years, 
                                                                  years))
    for (ri in seq_len(nrel)) {
      rel_t <- yr_idx(rel_years[ri])
      cum_F <- 0
      cum_M <- 0
      for (t in rel_t:nyears) {
        lib <- t - rel_t
        if (lib > 0) {
          cum_F <- sum(F_t[rel_t:(t - 1)])
          cum_M <- M * lib
        }
        else {
          cum_F <- 0
          cum_M <- 0
        }
        pred[ri, t] <- phi_mat[ri, t] * lambda_mat[ri, 
                                                   t] * u_fm(F_t[t], M) * exp(-cum_F - cum_M)
      }
    }
    p_nsa <- pmax(1 - rowSums(pred), 1e-10)
    ll_recap <- sum(obs_mat[, as.character(years)] * log(pmax(pred, 
                                                              1e-10)), na.rm = TRUE)
    ll_nsa <- sum(obs_mat[, "nsa"] * log(p_nsa), na.rm = TRUE)
    nll <- -(ll_recap + ll_nsa)
    if (is.finite(sigma_F) && sigma_F > 0 && nyears > 1) 
      nll <- nll + sum(diff(log_F)^2)/(2 * sigma_F^2)
    return(nll)
  }
  par_init <- c(log(M_init), rep(log(0.1), nyears))
  message(sprintf("Optimising M + %d year-specific F values...", 
                  nyears))
  opt <- nlminb(par_init, function(p) nll_fn(p[1], p[-1]), 
                lower = rep(log(1e-04), 1 + nyears), upper = rep(log(5), 
                                                                 1 + nyears), control = list(iter.max = 2000, eval.max = 5000, 
                                                                                             rel.tol = 1e-09))
  M_mle <- exp(opt$par[1])
  F_yr <- setNames(exp(opt$par[-1]), years)
  Z_yr <- setNames(M_mle + F_yr, years)
  nll_mle <- opt$objective
  message(sprintf("MLE: M = %.4f, mean(F) = %.4f, NLL = %.2f", 
                  M_mle, mean(F_yr), nll_mle))

  ## ---- Delta-method uncertainty for F ----
  message("Computing delta-method SEs for F...")
  H <- tryCatch(
    optimHess(opt$par, function(p) nll_fn(p[1], p[-1])),
    error = function(e) { message("Hessian failed: ", e$message); NULL }
  )
  F_lo95 <- F_hi95 <- rep(NA_real_, nyears)
  meanF <- mean(F_yr)
  meanF_lo95 <- meanF_hi95 <- meanF_cv <- NA_real_
  F_cv <- rep(NA_real_, nyears)

  if (!is.null(H)) {
    V <- tryCatch(solve(H), error = function(e) NULL)
    if (!is.null(V) && all(is.finite(diag(V))) && all(diag(V) > 0)) {
      ## log(F_t) parameters are indices 2:(nyears+1) in the par vector
      logF_idx  <- 2:(nyears + 1L)
      se_logF   <- sqrt(diag(V)[logF_idx])
      F_lo95    <- setNames(F_yr * exp(-1.96 * se_logF), years)
      F_hi95    <- setNames(F_yr * exp( 1.96 * se_logF), years)
      F_cv      <- setNames(sqrt(exp(se_logF^2) - 1), years)   # CV on natural scale

      ## Mean F delta method: grad = (1/nyears) * F_t  (chain rule through exp)
      ## Var(mean_F) = g' V_logF g,  g_t = F_t / nyears
      g         <- F_yr / nyears
      V_logF    <- V[logF_idx, logF_idx]
      var_meanF <- as.numeric(t(g) %*% V_logF %*% g)
      se_meanF  <- sqrt(max(var_meanF, 0))
      ## CI via log transform of mean F to keep positive
      se_log_meanF <- sqrt(log(1 + (se_meanF / meanF)^2))
      meanF_lo95   <- meanF * exp(-1.96 * se_log_meanF)
      meanF_hi95   <- meanF * exp( 1.96 * se_log_meanF)
      meanF_cv     <- se_meanF / meanF
      message(sprintf("Mean F = %.4f  (95%% CI: %.4f - %.4f,  CV = %.3f)",
                      meanF, meanF_lo95, meanF_hi95, meanF_cv))
    } else {
      message("Variance matrix not positive definite — F CIs not computed")
    }
  }
  pred_mat <- matrix(0, nrow = nrel, ncol = nyears + 1, dimnames = list(rel_years, 
                                                                        c(years, "nsa")))
  for (ri in seq_len(nrel)) {
    rel_t <- yr_idx(rel_years[ri])
    for (t in rel_t:nyears) {
      lib <- t - rel_t
      cum_F <- if (lib > 0) 
        sum(F_yr[rel_t:(t - 1)])
      else 0
      cum_M <- M_mle * lib
      pred_mat[ri, t] <- phi_mat[ri, t] * lambda_mat[ri, 
                                                     t] * u_fm(F_yr[t], M_mle) * exp(-cum_F - cum_M)
    }
    pred_mat[ri, "nsa"] <- max(1 - sum(pred_mat[ri, as.character(years)]), 
                               1e-10)
  }
  pred_n <- pred_mat * n_rel_vec
  chi_thresh <- 1.92
  if (profile_M) {
    message("Computing likelihood profile for M...")
    M_grid <- exp(seq(log(max(0.005, M_mle * 0.1)), log(min(3,
                                                            M_mle * 10)), length.out = 40))
    prof_nll <- numeric(length(M_grid))
    mle_k <- which.min(abs(M_grid - M_mle))
    lF_warm <- opt$par[-1]
    for (k in mle_k:length(M_grid)) {
      opt_k <- tryCatch(nlminb(lF_warm, function(lF) nll_fn(log(M_grid[k]),
                                                            lF), lower = rep(log(1e-04), nyears), upper = rep(log(5),
                                                                                                              nyears), control = list(iter.max = 200, eval.max = 400,
                                                                                                                                      rel.tol = 1e-06)), error = function(e) list(objective = Inf,
                                                                                                                                                                                  par = lF_warm))
      prof_nll[k] <- opt_k$objective
      if (is.finite(opt_k$objective))
        lF_warm <- opt_k$par
    }
    lF_warm <- opt$par[-1]
    for (k in mle_k:1) {
      opt_k <- tryCatch(nlminb(lF_warm, function(lF) nll_fn(log(M_grid[k]),
                                                            lF), lower = rep(log(1e-04), nyears), upper = rep(log(5),
                                                                                                              nyears), control = list(iter.max = 200, eval.max = 400,
                                                                                                                                      rel.tol = 1e-06)), error = function(e) list(objective = Inf,
                                                                                                                                                                                  par = lF_warm))
      prof_nll[k] <- opt_k$objective
      if (is.finite(opt_k$objective))
        lF_warm <- opt_k$par
    }
    in_ci <- (prof_nll - nll_mle) <= chi_thresh
    M_lo95 <- if (any(in_ci)) min(M_grid[in_ci]) else NA_real_
    M_hi95 <- if (any(in_ci)) max(M_grid[in_ci]) else NA_real_
    profile_df <- data.frame(M = M_grid, nll = prof_nll,
                             delta = prof_nll - nll_mle, in_ci95 = in_ci)
    message(sprintf("M = %.4f  (95%% CI: %.4f - %.4f)", M_mle, M_lo95, M_hi95))
  } else {
    message(sprintf("M = %.4f  (profile skipped)", M_mle))
    M_lo95 <- M_hi95 <- NA_real_
    profile_df <- data.frame(M = M_mle, nll = nll_mle, delta = 0, in_ci95 = TRUE)
  }
  M_boot <- NULL
  if (n_boot > 0) {
    message(sprintf("Running %d bootstrap resamples...", 
                    n_boot))
    M_boot <- numeric(n_boot)
    obs_orig <- obs
    for (b in seq_len(n_boot)) {
      set.seed(b)
      obs_b <- do.call(rbind, lapply(rel_years, function(ry) {
        sub <- obs_orig[obs_orig$relyr == ry, ]
        sub[sample(nrow(sub), replace = TRUE), ]
      }))
      tryCatch({
        obs_mat_b <- matrix(0L, nrow = nrel, ncol = nyears + 
                              1, dimnames = list(rel_years, c(years, "nsa")))
        for (ry in rel_years) {
          ri <- ry_idx(ry)
          sub <- obs_b[obs_b$relyr == ry, ]
          obs_mat_b[ri, "nsa"] <- nrow(sub)
          recs <- sub[sub$isrecap == 1 & !is.na(sub$recyr), ]
          for (i in seq_len(nrow(recs))) {
            cy <- as.character(recs$recyr[i])
            if (cy %in% colnames(obs_mat_b)) {
              obs_mat_b[ri, cy] <- obs_mat_b[ri, cy] + 1L
              obs_mat_b[ri, "nsa"] <- obs_mat_b[ri, "nsa"] - 1L
            }
          }
        }
        obs_mat_save <- obs_mat
        obs_mat <<- obs_mat_b
        opt_b <- nlminb(opt$par, function(p) nll_fn(p[1], p[-1]),
                        lower = rep(log(1e-04), 1 + nyears), 
                        upper = rep(log(5), 1 + nyears),
                        control = list(iter.max = 500))
        M_boot[b] <- exp(opt_b$par[1])
        obs_mat <<- obs_mat_save
      }, error = function(e) {
        M_boot[b] <<- NA_real_
      })
    }
    M_boot <- M_boot[!is.na(M_boot)]
    message(sprintf("Bootstrap M: median=%.4f, 95%% CI: %.4f-%.4f", 
                    median(M_boot), quantile(M_boot, 0.025), quantile(M_boot, 0.975)))
  }

  if (plot) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)

    ## Layout: slots 1-4 as 2x2 top-left, slot 5 spans full bottom row
    layout(matrix(c(1, 2, 3,
                    4, 5, 5), nrow = 2, byrow = TRUE),
           widths = c(1, 1, 1), heights = c(1, 1))
    par(mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

    ## ---- Panel 1: Likelihood profile for M ----
    if (profile_M) {
      ylim_a <- c(0, min(15, max(profile_df$delta[is.finite(profile_df$delta)])))
      plot(profile_df$M, profile_df$delta, type = "l", lwd = 2,
           col = "#185FA5", xlab = "M (natural mortality)",
           ylab = expression(Delta ~ NLL), main = "Likelihood Profile for M",
           ylim = ylim_a)
      abline(h = chi_thresh, lty = 2, col = "red")
      abline(v = M_mle, col = "black", lwd = 1.5)
      if (!is.na(M_lo95)) abline(v = M_lo95, lty = 3, col = "grey40")
      if (!is.na(M_hi95)) abline(v = M_hi95, lty = 3, col = "grey40")
      legend("topright", legend = c(sprintf("MLE = %.3f", M_mle),
                                    sprintf("95%% CI: %.3f - %.3f", M_lo95, M_hi95),
                                    "Chi-sq (1.92)"),
             lty = c(1, 3, 2), col = c("black", "grey40", "red"),
             bty = "n", cex = 0.8)
    } else {
      plot.new()
      text(0.5, 0.5, sprintf("M = %.4f\n(profile skipped)", M_mle),
           cex = 1.1, col = "#185FA5")
    }

    ## ---- Panel 2: Estimated F and Z by year, with delta-method CI on F ----
    ylim_b <- c(0, max(Z_yr, F_hi95, na.rm = TRUE) * 1.15)
    plot(years, F_yr, type = "n",
         xlab = "Year", ylab = "Annual mortality rate",
         main = "Estimated F and Z by Year",
         ylim = ylim_b)
    ## F 95% CI ribbon
    if (any(is.finite(F_lo95)) && any(is.finite(F_hi95))) {
      polygon(c(years, rev(years)),
              c(F_hi95, rev(F_lo95)),
              col = adjustcolor("#E05C1A", alpha.f = 0.18), border = NA)
    }
    lines(years, Z_yr, lwd = 2, col = "black")
    lines(years, F_yr, lwd = 2, col = "#E05C1A")
    abline(h = M_mle,  lty = 3, col = "#185FA5", lwd = 1.5)
    abline(h = meanF,  lty = 2, col = "#E05C1A", lwd = 1.2)
    ## Mean F CI as short horizontal segment at right margin
    usr <- par("usr")
    x_right <- usr[2] - diff(usr[1:2]) * 0.01
    if (is.finite(meanF_lo95) && is.finite(meanF_hi95)) {
      segments(x_right, meanF_lo95, x_right, meanF_hi95,
               col = "#E05C1A", lwd = 2.5, lend = "butt")
      points(x_right, meanF, pch = 21, bg = "#E05C1A", col = "white", cex = 1.1)
    }
    mean_F_label <- if (is.finite(meanF_lo95))
      sprintf("mean F = %.3f (%.3f-%.3f)", meanF, meanF_lo95, meanF_hi95)
    else
      sprintf("mean F = %.3f", meanF)
    legend("topright",
           legend = c("Z = M+F", "F  (95% CI shaded)", sprintf("M = %.3f", M_mle),
                      mean_F_label),
           lty  = c(1, 1, 3, 2),
           lwd  = c(2, 2, 1.5, 1.2),
           col  = c("black", "#E05C1A", "#185FA5", "#E05C1A"),
           bty  = "n", cex = 0.75)

    ## ---- Panel 3: Obs vs Pred recaptures by year ----
    obs_yr  <- colSums(obs_mat[, as.character(years), drop = FALSE])
    pred_yr <- colSums(pred_n[,  as.character(years), drop = FALSE])
    ylim_c  <- c(0, max(obs_yr, pred_yr, na.rm = TRUE) * 1.15)
    plot(years, obs_yr, pch = 16, col = "black",
         xlab = "Year", ylab = "Number of recaptures",
         main = "Obs vs Pred: Recaptures by Year", ylim = ylim_c)
    lines(years, pred_yr, col = "#185FA5", lwd = 2)
    legend("topright", legend = c("Observed", "Predicted"),
           pch = c(16, NA), lty = c(NA, 1),
           col = c("black", "#185FA5"), bty = "n", cex = 0.8)

    ## ---- Panel 4: Obs vs Pred recovery by years at liberty ----
    max_lib <- 10
    obs_lib <- pred_lib <- n_lib_vec <- numeric(max_lib)
    for (ri in seq_len(nrel)) {
      rel_t <- yr_idx(rel_years[ri])
      n_rel <- n_rel_vec[ri]
      for (lib in 1:max_lib) {
        t <- rel_t + lib - 1
        if (t > nyears) break
        obs_lib[lib]   <- obs_lib[lib]   + obs_mat[ri, t]
        pred_lib[lib]  <- pred_lib[lib]  + pred_n[ri, t]
        n_lib_vec[lib] <- n_lib_vec[lib] + n_rel
      }
    }
    valid  <- n_lib_vec > 0
    obs_r  <- obs_lib[valid]  / n_lib_vec[valid]
    pred_r <- pred_lib[valid] / n_lib_vec[valid]
    libs   <- (1:max_lib)[valid]
    ylim_d <- c(0, max(obs_r, pred_r, na.rm = TRUE) * 1.2)
    plot(libs, obs_r, pch = 16,
         xlab = "Years at liberty", ylab = "Recovery rate",
         main = "Obs vs Pred: Recovery by Liberty", ylim = ylim_d)
    lines(libs, pred_r, col = "#185FA5", lwd = 2)
    legend("topright", legend = c("Observed", "Predicted"),
           pch = c(16, NA), lty = c(NA, 1),
           col = c("black", "#185FA5"), bty = "n", cex = 0.8)

    ## ---- Panel 5 (wide): Pearson residuals, all cohorts, colour by release year ----
    ## Build residual data frame
    resid_rows <- vector("list", nrel)
    for (ri in seq_len(nrel)) {
      ry    <- rel_years[ri]
      rel_t <- yr_idx(ry)
      n_r   <- n_rel_vec[ri]
      rows  <- lapply(rel_t:nyears, function(t) {
        o <- obs_mat[ri, t]
        p <- pred_n[ri, t]
        if (p < 0.5) return(NULL)   # skip near-empty cells
        data.frame(relyr   = ry,
                   recyr   = years[t],
                   lib     = t - rel_t + 1L,
                   obs     = as.numeric(o),
                   pred    = p,
                   pearson = (o - p) / sqrt(p),
                   n_rel   = n_r)
      })
      resid_rows[[ri]] <- do.call(rbind, Filter(Negate(is.null), rows))
    }
    resid_df <- do.call(rbind, Filter(Negate(is.null), resid_rows))

    ry_labs  <- sort(unique(resid_df$relyr))
    n_cohort <- length(ry_labs)
    pal      <- colorRampPalette(c("#4575B4", "#91BFDB", "#74C476",
                                   "#FD8D3C", "#D73027"))(n_cohort)
    col_map  <- setNames(pal, as.character(ry_labs))

    pr_lim <- max(abs(resid_df$pearson), na.rm = TRUE) * 1.15
    pr_lim <- max(pr_lim, 2.5)   # always show ±2 lines with headroom

    par(mar = c(4, 4.5, 3, 1))
    plot(resid_df$recyr, resid_df$pearson,
         type = "n",
         xlim = range(resid_df$recyr),
         ylim = c(-pr_lim, pr_lim),
         xlab = "Recapture year",
         ylab = expression("Pearson residual  " * (obs - pred) / sqrt(pred)),
         main = "Pearson Residuals by Release Year")
    abline(h =  0, col = "grey40", lwd = 1.2)
    abline(h =  2, col = "grey60", lwd = 0.8, lty = 2)
    abline(h = -2, col = "grey60", lwd = 0.8, lty = 2)

    jit <- 0.15
    for (ry in ry_labs) {
      sub   <- resid_df[resid_df$relyr == ry, ]
      sub   <- sub[order(sub$recyr), ]
      col_i <- col_map[as.character(ry)]
      points(jitter(sub$recyr, amount = jit), sub$pearson,
             pch = 21,
             bg  = adjustcolor(col_i, alpha.f = 0.75),
             col = adjustcolor(col_i, alpha.f = 0.95),
             cex = 0.8 + 0.014 * sqrt(sub$n_rel))
    }

    legend("topright",
           legend = as.character(ry_labs),
           col    = pal,
           pt.bg  = adjustcolor(pal, alpha.f = 0.75),
           pch    = 21,
           pt.cex = 0.85,
           bty    = "n",
           cex    = 0.62,
           ncol   = 2,
           title  = "Release yr")

    ## ---- Overall title ----
    mtext(sprintf("MLL=%g mm  tag_loss=%.2f  report_rate=%.2f  sigma_F=%.2f  |  n=%d  n_recap=%d",
                  MLL, tag_loss[1], report_rate[1], sigma_F, n_legal, n_recap),
          side = 3, outer = TRUE, line = 0.3, cex = 0.72)

    ## ---- Bootstrap histogram (new page if requested) ----
    if (!is.null(M_boot) && length(M_boot) > 10) {
      par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
      hist(M_boot, breaks = 30, col = "#378ADD80", border = "white",
           xlab = "M (bootstrap)", main = "Bootstrap Distribution of M",
           freq = FALSE)
      abline(v = M_mle, col = "black", lwd = 2)
      abline(v = quantile(M_boot, c(0.025, 0.975)), col = "red", lty = 2)
    }
  }

  invisible(list(M = M_mle, M_lo95 = M_lo95, M_hi95 = M_hi95,
                 M_boot = M_boot, F_yr = F_yr, F_lo95 = F_lo95, F_hi95 = F_hi95,
                 F_cv = F_cv, meanF = meanF, meanF_lo95 = meanF_lo95,
                 meanF_hi95 = meanF_hi95, meanF_cv = meanF_cv,
                 Z_yr = Z_yr, nll = nll_mle,
                 n_legal = n_legal, n_recap = n_recap, MLL = MLL,
                 tag_loss = tag_loss, report_rate = report_rate,
                 sigma_F = sigma_F, obs_mat = obs_mat, pred_mat = pred_n,
                 profile = profile_df))
}
