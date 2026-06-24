#' Estimate Natural Mortality Using Brownie Dead-Recovery Model
#'
#' Estimates natural mortality (M) from tag-recapture data using a Brownie
#' dead-recovery model. Follows a matrix-based formulation where rows are
#' release cohorts (by year) and columns are recapture years, with a final
#' 'nsa' (not seen again) column. Year-specific fishing mortality is estimated
#' directly from the recovery pattern — no catch or effort data required.
#' Multiple release cohorts sharing the same F in a given year provide the
#' cross-cohort signal that identifies M.
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
#' @param n_boot Integer. Number of bootstrap resamples for uncertainty.
#'   Set to 0 to skip bootstrapping and use likelihood profiling only.
#'   Default 0.
#'
#' @param plot Logical. Produce diagnostic plots. Default TRUE.
#'
#' @return Invisibly returns a list with:
#'   \itemize{
#'     \item \code{M}         - MLE of natural mortality (annual rate)
#'     \item \code{M_lo95}    - lower 95\% profile CI
#'     \item \code{M_hi95}    - upper 95\% profile CI
#'     \item \code{M_boot}    - bootstrap distribution of M (if n_boot > 0)
#'     \item \code{F_yr}      - named vector of annual F estimates
#'     \item \code{Z_yr}      - named vector of annual Z = M + F
#'     \item \code{nll}       - NLL at MLE
#'     \item \code{n_legal}   - animals included
#'     \item \code{n_recap}   - recapture events
#'     \item \code{obs_mat}   - observed recovery matrix (cohort x year + nsa)
#'     \item \code{pred_mat}  - predicted recovery matrix
#'     \item \code{profile}   - likelihood profile data frame
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
#' @examples
#' \dontrun{
#' result <- EstimateMBrownie(obs, MLL=110, tag_loss=0.05)
#' cat("M:", result$M, "(", result$M_lo95, "-", result$M_hi95, ")\n")
#'
#' # Sensitivity to sigma_F
#' lapply(c(0.2, 0.5, 1.0, Inf), function(sf) {
#'   r <- EstimateMBrownie(obs, MLL=110, tag_loss=0.05, sigma_F=sf, plot=FALSE)
#'   data.frame(sigma_F=sf, M=round(r$M,3), lo=round(r$M_lo95,3),
#'              hi=round(r$M_hi95,3))
#' }) |> do.call(what=rbind)
#'
#' # Sensitivity to tag loss
#' lapply(c(0, 0.05, 0.10, 0.20), function(tl) {
#'   r <- EstimateMBrownie(obs, MLL=110, tag_loss=tl, plot=FALSE)
#'   data.frame(tag_loss=tl, M=round(r$M,3), lo=round(r$M_lo95,3),
#'              hi=round(r$M_hi95,3))
#' }) |> do.call(what=rbind)
#'
#' # Sensitivity to reporting rate
#' lapply(c(0.7, 0.8, 0.9, 1.0), function(rr) {
#'   r <- EstimateMBrownie(obs, MLL=110, tag_loss=0.05, report_rate=rr,
#'                         plot=FALSE)
#'   data.frame(report_rate=rr, M=round(r$M,3), lo=round(r$M_lo95,3),
#'              hi=round(r$M_hi95,3))
#' }) |> do.call(what=rbind)
#' }
#'
#' @export
EstimateMBrownie <- function(obs,
                              MLL         = 110,
                              tag_loss    = 0.05,
                              report_rate = 1,
                              sigma_F     = 0.5,
                              M_init      = 0.1,
                              n_boot      = 0,
                              plot        = TRUE) {

  ## ── 0. Input checks ──────────────────────────────────────────────────────
  required <- c("tag", "LCl", "relyr", "recyr", "isrecap")
  missing  <- setdiff(required, names(obs))
  if (length(missing) > 0)
    stop("obs missing columns: ", paste(missing, collapse=", "))
  if (any(tag_loss < 0) || any(tag_loss >= 1))
    stop("tag_loss must be in [0, 1)")
  if (any(report_rate <= 0) || any(report_rate > 1))
    stop("report_rate must be in (0, 1]")

  ## ── 1. Filter to legal-sized animals ─────────────────────────────────────
  obs     <- obs[!is.na(obs$LCl) & obs$LCl >= MLL, ]
  n_legal <- nrow(obs)
  n_recap <- sum(obs$isrecap, na.rm=TRUE)
  if (n_legal == 0) stop("No animals with LCl >= MLL = ", MLL)
  if (n_recap == 0) stop("No recaptured animals with LCl >= MLL = ", MLL)
  message(sprintf("Using %d animals (%d recaptured) with LCl >= %g mm",
                  n_legal, n_recap, MLL))

  ## ── 2. Build year structure ───────────────────────────────────────────────
  rel_years <- sort(unique(obs$relyr))
  rec_years <- sort(unique(obs$recyr[!is.na(obs$recyr) & obs$isrecap == 1]))
  years     <- sort(unique(c(rel_years, rec_years)))
  yr_min    <- min(years); yr_max <- max(years)
  years     <- yr_min:yr_max
  nyears    <- length(years)
  nrel      <- length(rel_years)
  yr_idx    <- function(y) match(y, years)
  ry_idx    <- function(y) match(y, rel_years)

  ## ── 3. Build observed recovery matrix ─────────────────────────────────────
  # rows = release cohorts, cols = recapture years + 'nsa'
  obs_mat <- matrix(0L,
                    nrow     = nrel,
                    ncol     = nyears + 1L,
                    dimnames = list(rel_years, c(years, "nsa")))

  n_rel_vec <- setNames(integer(nrel), rel_years)

  for (ry in rel_years) {
    sub   <- obs[obs$relyr == ry, ]
    ri    <- ry_idx(ry)
    n_rel_vec[as.character(ry)] <- nrow(sub)

    # Recaptured animals
    recs  <- sub[sub$isrecap == 1 & !is.na(sub$recyr), ]
    for (i in seq_len(nrow(recs))) {
      cy <- as.character(recs$recyr[i])
      if (cy %in% colnames(obs_mat))
        obs_mat[ri, cy] <- obs_mat[ri, cy] + 1L
    }
    # Not seen again
    obs_mat[ri, "nsa"] <- nrow(sub) - sum(obs_mat[ri, as.character(years)])
  }

  ## ── 4. Build tag-loss and reporting-rate matrices ─────────────────────────
  # phi_mat[r, t] = reporting rate for cohort r recaptured in year t
  # lambda_mat[r, t] = cumulative tag retention for cohort r at year t
  #                  = (1 - annual_tag_loss)^(years_at_liberty)

  phi_mat    <- matrix(1,    nrow=nrel, ncol=nyears,
                       dimnames=list(rel_years, years))
  lambda_mat <- matrix(1,    nrow=nrel, ncol=nyears,
                       dimnames=list(rel_years, years))

  # Expand report_rate to vector over years
  phi_vec <- if (length(report_rate) == 1) rep(report_rate, nyears) else {
    if (length(report_rate) != nyears)
      stop("report_rate must be scalar or length nyears = ", nyears)
    report_rate
  }
  for (ri in seq_len(nrel)) phi_mat[ri, ] <- phi_vec

  # Expand tag_loss: scalar -> annual rate; matrix -> use directly
  if (is.matrix(tag_loss)) {
    if (!identical(dim(tag_loss), c(nrel, nyears)))
      stop("tag_loss matrix must be nrel x nyears = ", nrel, " x ", nyears)
    lambda_mat <- tag_loss
  } else {
    # Cumulative tag retention: (1 - rate)^years_at_liberty
    tl_annual <- tag_loss[1]
    for (ri in seq_len(nrel)) {
      rel_t <- yr_idx(rel_years[ri])
      for (t in seq_len(nyears)) {
        lib <- t - rel_t        # years at liberty (0 = year of release)
        lambda_mat[ri, t] <- if (lib < 0) 0 else (1 - tl_annual)^lib
      }
    }
  }

  ## ── 5. u(F,M) — probability of dying from fishing ────────────────────────
  u_fm <- function(F, M) {
    val <- (F / (F + M)) * (1 - exp(-(F + M)))
    pmax(val, 1e-10)
  }

  ## ── 6. Core NLL function (matrix formulation) ────────────────────────────
  nll_fn <- function(log_M, log_F) {

    M   <- exp(log_M)
    F_t <- exp(log_F)           # length nyears

    ## Predicted probability matrix
    pred <- matrix(0, nrow=nrel, ncol=nyears,
                   dimnames=list(rel_years, years))

    for (ri in seq_len(nrel)) {
      rel_t   <- yr_idx(rel_years[ri])
      cum_F   <- 0              # cumulative F from release to start of year t
      cum_M   <- 0              # cumulative M (= M * years_at_liberty)

      for (t in rel_t:nyears) {
        lib <- t - rel_t        # years at liberty at start of year t
        if (lib > 0) {
          cum_F <- sum(F_t[rel_t:(t-1)])
          cum_M <- M * lib
        } else {
          cum_F <- 0; cum_M <- 0
        }
        pred[ri, t] <- phi_mat[ri, t] *
                       lambda_mat[ri, t] *
                       u_fm(F_t[t], M) *
                       exp(-cum_F - cum_M)
      }
    }

    ## p_nsa = 1 - sum of all recovery probabilities for that cohort
    p_nsa <- pmax(1 - rowSums(pred), 1e-10)

    ## Log-likelihood: recaptures + not seen again
    ll_recap <- sum(obs_mat[, as.character(years)] *
                    log(pmax(pred, 1e-10)), na.rm=TRUE)
    ll_nsa   <- sum(obs_mat[, "nsa"] * log(p_nsa), na.rm=TRUE)

    nll <- -(ll_recap + ll_nsa)

    ## Random walk penalty on log(F)
    if (is.finite(sigma_F) && sigma_F > 0 && nyears > 1)
      nll <- nll + sum(diff(log_F)^2) / (2 * sigma_F^2)

    return(nll)
  }

  ## ── 7. Optimise ───────────────────────────────────────────────────────────
  par_init <- c(log(M_init), rep(log(0.1), nyears))
  message(sprintf("Optimising M + %d year-specific F values...", nyears))

  opt <- nlminb(
    par_init,
    function(p) nll_fn(p[1], p[-1]),
    lower   = rep(log(1e-4), 1 + nyears),
    upper   = rep(log(5.0),  1 + nyears),
    control = list(iter.max=2000, eval.max=5000, rel.tol=1e-9)
  )

  M_mle   <- exp(opt$par[1])
  F_yr    <- setNames(exp(opt$par[-1]), years)
  Z_yr    <- setNames(M_mle + F_yr,    years)
  nll_mle <- opt$objective

  message(sprintf("MLE: M = %.4f, mean(F) = %.4f, NLL = %.2f",
                  M_mle, mean(F_yr), nll_mle))

  ## Compute predicted matrix at MLE for output and plotting
  pred_mat <- matrix(0, nrow=nrel, ncol=nyears+1,
                     dimnames=list(rel_years, c(years, "nsa")))
  for (ri in seq_len(nrel)) {
    rel_t <- yr_idx(rel_years[ri])
    for (t in rel_t:nyears) {
      lib   <- t - rel_t
      cum_F <- if (lib > 0) sum(F_yr[rel_t:(t-1)]) else 0
      cum_M <- M_mle * lib
      pred_mat[ri, t] <- phi_mat[ri, t] * lambda_mat[ri, t] *
                         u_fm(F_yr[t], M_mle) * exp(-cum_F - cum_M)
    }
    pred_mat[ri, "nsa"] <- max(1 - sum(pred_mat[ri, as.character(years)]),
                               1e-10)
  }
  # Convert to expected numbers
  pred_n <- pred_mat * n_rel_vec   # element-wise, recycled by row

  ## ── 8. Likelihood profile for M ──────────────────────────────────────────
  message("Computing likelihood profile for M...")
  # Use 40 points (sufficient for CI) on log scale, log-spaced around MLE
  M_grid   <- exp(seq(log(max(0.005, M_mle * 0.1)),
                      log(min(3.0,   M_mle * 10)),
                      length.out = 40))
  prof_nll <- numeric(length(M_grid))

  # Warm-start: sweep left from MLE then right, carrying forward F estimates
  mle_k    <- which.min(abs(M_grid - M_mle))   # closest grid point to MLE
  lF_warm  <- opt$par[-1]                       # start from MLE F

  # Sweep upward from MLE
  for (k in mle_k:length(M_grid)) {
    opt_k <- tryCatch(
      nlminb(lF_warm,
             function(lF) nll_fn(log(M_grid[k]), lF),
             lower   = rep(log(1e-4), nyears),
             upper   = rep(log(5.0),  nyears),
             control = list(iter.max=200, eval.max=400, rel.tol=1e-6)),
      error = function(e) list(objective=Inf, par=lF_warm))
    prof_nll[k] <- opt_k$objective
    if (is.finite(opt_k$objective)) lF_warm <- opt_k$par  # warm-start next step
  }

  # Sweep downward from MLE
  lF_warm <- opt$par[-1]
  for (k in mle_k:1) {
    opt_k <- tryCatch(
      nlminb(lF_warm,
             function(lF) nll_fn(log(M_grid[k]), lF),
             lower   = rep(log(1e-4), nyears),
             upper   = rep(log(5.0),  nyears),
             control = list(iter.max=200, eval.max=400, rel.tol=1e-6)),
      error = function(e) list(objective=Inf, par=lF_warm))
    prof_nll[k] <- opt_k$objective
    if (is.finite(opt_k$objective)) lF_warm <- opt_k$par
  }

  chi_thresh <- 1.92
  in_ci      <- (prof_nll - nll_mle) <= chi_thresh
  M_lo95     <- if (any(in_ci)) min(M_grid[in_ci]) else NA_real_
  M_hi95     <- if (any(in_ci)) max(M_grid[in_ci]) else NA_real_
  profile_df <- data.frame(M=M_grid, nll=prof_nll,
                            delta=prof_nll - nll_mle, in_ci95=in_ci)

  message(sprintf("M = %.4f  (95%% CI: %.4f - %.4f)", M_mle, M_lo95, M_hi95))

  ## ── 9. Bootstrap (optional) ───────────────────────────────────────────────
  M_boot <- NULL
  if (n_boot > 0) {
    message(sprintf("Running %d bootstrap resamples...", n_boot))
    # Resample individual animals within each cohort
    M_boot <- numeric(n_boot)
    obs_orig <- obs
    for (b in seq_len(n_boot)) {
      set.seed(b)
      # Resample rows within each release cohort
      obs_b <- do.call(rbind, lapply(rel_years, function(ry) {
        sub <- obs_orig[obs_orig$relyr == ry, ]
        sub[sample(nrow(sub), replace=TRUE), ]
      }))
      # Quick re-fit (no profiling)
      tryCatch({
        # Rebuild obs_mat for bootstrap sample
        obs_mat_b <- matrix(0L, nrow=nrel, ncol=nyears+1,
                            dimnames=list(rel_years, c(years,"nsa")))
        for (ry in rel_years) {
          ri  <- ry_idx(ry)
          sub <- obs_b[obs_b$relyr==ry, ]
          obs_mat_b[ri, "nsa"] <- nrow(sub)
          recs <- sub[sub$isrecap==1 & !is.na(sub$recyr), ]
          for (i in seq_len(nrow(recs))) {
            cy <- as.character(recs$recyr[i])
            if (cy %in% colnames(obs_mat_b)) {
              obs_mat_b[ri, cy] <- obs_mat_b[ri, cy] + 1L
              obs_mat_b[ri, "nsa"] <- obs_mat_b[ri, "nsa"] - 1L
            }
          }
        }
        # Temporarily swap obs_mat
        obs_mat_save <- obs_mat
        obs_mat      <<- obs_mat_b
        opt_b <- nlminb(opt$par,
                        function(p) nll_fn(p[1], p[-1]),
                        lower=rep(log(1e-4), 1+nyears),
                        upper=rep(log(5.0),  1+nyears),
                        control=list(iter.max=500))
        M_boot[b] <- exp(opt_b$par[1])
        obs_mat   <<- obs_mat_save
      }, error=function(e) { M_boot[b] <<- NA_real_ })
    }
    M_boot <- M_boot[!is.na(M_boot)]
    message(sprintf("Bootstrap M: median=%.4f, 95%% CI: %.4f-%.4f",
                    median(M_boot),
                    quantile(M_boot, 0.025),
                    quantile(M_boot, 0.975)))
  }

  ## ── 10. Diagnostic plots ─────────────────────────────────────────────────
  if (plot) {
    old_par <- par(no.readonly=TRUE)
    on.exit(par(old_par), add=TRUE)
    par(mfrow=c(2,3), mar=c(4,4,3,1))

    ## (a) Likelihood profile
    ylim_a <- c(0, min(15, max(profile_df$delta[is.finite(profile_df$delta)])))
    plot(profile_df$M, profile_df$delta, type="l", lwd=2, col="#185FA5",
         xlab="M (natural mortality)", ylab=expression(Delta~NLL),
         main="Likelihood Profile for M", ylim=ylim_a)
    abline(h=chi_thresh, lty=2, col="red")
    abline(v=M_mle,  col="black",  lwd=1.5)
    if (!is.na(M_lo95)) abline(v=M_lo95, lty=3, col="grey40")
    if (!is.na(M_hi95)) abline(v=M_hi95, lty=3, col="grey40")
    legend("topright",
           legend=c(sprintf("MLE = %.3f", M_mle),
                    sprintf("95%% CI: %.3f - %.3f", M_lo95, M_hi95),
                    "Chi-sq (1.92)"),
           lty=c(1,3,2), col=c("black","grey40","red"), bty="n", cex=0.8)

    ## (b) Annual F and Z
    plot(years, F_yr, type="l", lwd=2, col="#E05C1A",
         xlab="Year", ylab="Annual mortality rate",
         main="Estimated F and Z by Year",
         ylim=c(0, max(Z_yr)*1.15))
    lines(years, Z_yr, lwd=2, col="black")
    abline(h=M_mle, lty=3, col="#185FA5", lwd=1.5)
    legend("topright",
           legend=c("Z = M+F", "F", sprintf("M = %.3f", M_mle)),
           lty=c(1,1,3), lwd=c(2,2,1.5),
           col=c("black","#E05C1A","#185FA5"), bty="n", cex=0.8)

    ## (c) Obs vs predicted: total recaptures by calendar year
    obs_yr  <- colSums(obs_mat[, as.character(years), drop=FALSE])
    pred_yr <- colSums(pred_n[, as.character(years), drop=FALSE])
    ylim_c  <- c(0, max(obs_yr, pred_yr, na.rm=TRUE) * 1.15)
    plot(years, obs_yr, pch=16, col="black",
         xlab="Year", ylab="Number of recaptures",
         main="Obs vs Pred: Recaptures by Year", ylim=ylim_c)
    lines(years, pred_yr, col="#185FA5", lwd=2)
    legend("topright", legend=c("Observed","Predicted"),
           pch=c(16,NA), lty=c(NA,1), col=c("black","#185FA5"),
           bty="n", cex=0.8)

    ## (d) Obs vs predicted: pooled by years at liberty
    max_lib  <- 10
    obs_lib  <- pred_lib <- n_lib_vec <- numeric(max_lib)
    for (ri in seq_len(nrel)) {
      rel_t <- yr_idx(rel_years[ri])
      n_rel <- n_rel_vec[ri]
      for (lib in 1:max_lib) {
        t <- rel_t + lib - 1
        if (t > nyears) break
        obs_lib[lib]     <- obs_lib[lib]     + obs_mat[ri, t]
        pred_lib[lib]    <- pred_lib[lib]    + pred_n[ri, t]
        n_lib_vec[lib]   <- n_lib_vec[lib]   + n_rel
      }
    }
    valid   <- n_lib_vec > 0
    obs_r   <- obs_lib[valid]  / n_lib_vec[valid]
    pred_r  <- pred_lib[valid] / n_lib_vec[valid]
    libs    <- (1:max_lib)[valid]
    ylim_d  <- c(0, max(obs_r, pred_r, na.rm=TRUE) * 1.2)
    plot(libs, obs_r, pch=16,
         xlab="Years at liberty", ylab="Recovery rate",
         main="Obs vs Pred: Recovery by Liberty", ylim=ylim_d)
    lines(libs, pred_r, col="#185FA5", lwd=2)
    legend("topright", legend=c("Observed","Predicted"),
           pch=c(16,NA), lty=c(NA,1), col=c("black","#185FA5"),
           bty="n", cex=0.8)

    ## (e) Observed recovery matrix bubble plot (like your original)
    # sqrt-scaled circles, positive = filled dark, shown by release vs recap yr
    obs_vec  <- as.vector(obs_mat[, as.character(years)])
    pred_vec <- as.vector(pred_n[,  as.character(years)])
    rel_grid <- rep(rel_years, times=nyears)
    rec_grid <- rep(years,     each =nrel)
    # Only show cells where release <= recap and obs or pred > 0
    keep <- rel_grid <= rec_grid & (obs_vec > 0 | pred_vec > 0.5)

    par(mar=c(4,4,3,1))
    plot(rec_grid[keep], rel_grid[keep], type="n",
         xlab="Recapture year", ylab="Release year",
         main="Observed (dark) vs Predicted (red) Recaptures",
         xlim=range(years), ylim=range(rel_years))
    # Observed
    idx_obs <- keep & obs_vec > 0
    symbols(rec_grid[idx_obs], rel_grid[idx_obs],
            circles=sqrt(obs_vec[idx_obs]),
            fg=grey(0.2, 0.7), bg=grey(0.2, 0.5),
            inches=0.15, add=TRUE)
    # Predicted
    idx_pred <- keep & pred_vec > 0.5
    symbols(rec_grid[idx_pred] + 0.15, rel_grid[idx_pred],
            circles=sqrt(pred_vec[idx_pred]),
            fg=rgb(1,0,0,0.7), bg=rgb(1,0,0,0.3),
            inches=0.15, add=TRUE)
    legend("topright", legend=c("Observed","Predicted"),
           pch=21, pt.bg=c(grey(0.5,0.5), rgb(1,0,0,0.3)),
           col=c(grey(0.2,0.7), rgb(1,0,0,0.7)),
           bty="n", cex=0.8)

    ## (f) Recovery by cohort — obs (solid) vs pred (dashed)
    good <- rel_years[n_rel_vec >= 40]
    cols <- rainbow(length(good), alpha=0.85)
    plot(NA, xlim=c(0.5, 8.5), ylim=c(0, 0.15),
         xlab="Years at liberty", ylab="Recovery rate",
         main="Recovery by Cohort (n >= 40 released)")
    for (ci in seq_along(good)) {
      ry    <- good[ci]
      ri    <- ry_idx(ry)
      rel_t <- yr_idx(ry)
      n_rel <- n_rel_vec[ri]
      obs_r <- pred_r <- numeric(8)
      for (lib in 1:8) {
        t <- rel_t + lib - 1
        if (t > nyears) break
        obs_r[lib]  <- obs_mat[ri, t] / n_rel
        pred_r[lib] <- pred_n[ri,  t] / n_rel
      }
      lines(1:8, obs_r,  col=cols[ci], lwd=1.5, lty=1)
      lines(1:8, pred_r, col=cols[ci], lwd=1.5, lty=2)
    }
    legend("topright",
           legend=c(as.character(good), "Obs (—)", "Pred (--)"),
           col=c(cols, "black","black"),
           lty=c(rep(1,length(good)), 1, 2),
           lwd=c(rep(1.5,length(good)), 1.5, 1.5),
           bty="n", cex=0.6, ncol=2)

    # Bootstrap histogram if available
    if (!is.null(M_boot) && length(M_boot) > 10) {
      par(mfrow=c(1,1), mar=c(4,4,3,1))
      hist(M_boot, breaks=30, col="#378ADD80", border="white",
           xlab="M (bootstrap)", main="Bootstrap Distribution of M",
           freq=FALSE)
      abline(v=M_mle,             col="black", lwd=2)
      abline(v=quantile(M_boot, c(0.025,0.975)), col="red", lty=2)
    }

    mtext(sprintf(
      "MLL=%g mm  tag_loss=%.2f  report_rate=%.2f  sigma_F=%.2f  |  n=%d  n_recap=%d",
      MLL, tag_loss[1], report_rate[1], sigma_F, n_legal, n_recap),
      side=3, outer=TRUE, line=-1.2, cex=0.72)
  }

  ## ── 11. Return ────────────────────────────────────────────────────────────
  invisible(list(
    M           = M_mle,
    M_lo95      = M_lo95,
    M_hi95      = M_hi95,
    M_boot      = M_boot,
    F_yr        = F_yr,
    Z_yr        = Z_yr,
    nll         = nll_mle,
    n_legal     = n_legal,
    n_recap     = n_recap,
    MLL         = MLL,
    tag_loss    = tag_loss,
    report_rate = report_rate,
    sigma_F     = sigma_F,
    obs_mat     = obs_mat,
    pred_mat    = pred_n,
    profile     = profile_df
  ))
}
