#' Derive Annual Growth Parameters from a Sub-annual Tag2STM Fit
#'
#' Compounds the fitted sub-annual STMs (the \code{goodts} seasons) into a
#' single 12-month STM, \eqn{A}, on the target model's length bins, then finds
#' the annual parameter set whose one-step STM best reproduces \eqn{A}. The
#' objective is a weighted KL divergence summed over from-length columns.
#' Intended for supplying IMuLT with annual growth parameters, and priors for
#' them, that approximate the compounded sub-annual growth.
#'
#' @param stm_fun function(pars, bins) returning an nbin x nbin
#'   \strong{column-stochastic} STM (columns = from-length, as in ClipSTM)
#'   for ONE annual step on \code{bins} (a list with lbinL, lbinU, lbin). This
#'   must be the SAME builder IMuLT / buildfiles uses, including its
#'   floor and top-bin treatment, so that the returned parameters reproduce
#'   \eqn{A} when IMuLT builds its own STM from them.
#' @param par_template Named list of annual parameters in \code{stm_fun}'s
#'   format, e.g. \code{list(Growth_par = c(...), LCV = log(0.5))}. Values are
#'   used as starting values (estimated) or fixed values (not estimated).
#' @param est Character. Names in \code{par_template} to estimate.
#' @param LowLB,UpLB,Gap Lower edge of the first and last target bins, and
#'   target bin width (mm), i.e. IMuLT's length bins.
#' @param start Element of \code{goodts} at which the 12-month product starts.
#'   Set this to match where growth falls in IMuLT's model year. Default
#'   \code{goodts[1]}, as in ClipSTM.
#' @param years NULL (all supported years) or integer year indices to include,
#'   e.g. one era. Ignored when the fit has no year dimension.
#' @param by_year Logical. If TRUE, fit each year's annual STM separately.
#'   If FALSE, fit the average of the per-year annual STMs.
#' @param wts Weights over target from-bins (length nbin). NULL gives uniform
#'   weights; a representative length distribution (survey/catch LF or IMuLT
#'   numbers-at-length) is strongly recommended.
#' @param plus_group Logical. TRUE accumulates mass above UpLB into the top
#'   bin; FALSE drops and renormalises (ClipSTM behaviour). Match IMuLT.
#' @param ndraws Integer. Number of draws from the Tag2STM fit's fixed-effect
#'   covariance to propagate into the annual parameters (0 = none). Each draw
#'   recompounds A and re-projects, warm-started from the MLE.
#' @param nage Years for the mean-length-at-age diagnostic.
#' @param plot Logical. Print diagnostic plots.
#' @param mod,bins,datain,goodts Fitted RTMB object and Tag2STM inputs;
#'   default to the objects of those names in the calling environment.
#'
#' @return A list (one element per fitted STM: "avg", "all", or "yr<i>"),
#'   each containing: pars (relisted parameters), theta, objective,
#'   convergence, A, Ahat, moments (per-bin increment/sd/P(stay) for A and
#'   Ahat), traj (mean length-at-age for both), and, if ndraws > 0, prior
#'   (MLE, draw mean, draw sd per estimated element) and theta_draws.
#'
#' @details
#' Compounding is done on the native Tag2STM bins and only then coarsened to
#' the target bins. Coarsening sums to-bins and averages from-bins (uniform
#' within a target bin), so it is correct whether or not Gap equals the
#' native bin width.
#'
#' With year effects, each supported year's seasonal STMs are multiplied
#' first and the annual products then averaged. This is the
#' population-average 12-month transition (a mixture over years), so its
#' spread includes between-year variation. ClipSTM instead averages each
#' season across years before multiplying.
#'
#' The parameter draws assume no random effects (re = FALSE), which holds
#' for every settled Tag2STM configuration.
#'
#' @export
MakeAnnualPars <- function(stm_fun, par_template, est = names(par_template),
                           LowLB = 41, UpLB = 151, Gap = 2,
                           start = NULL, years = NULL, by_year = FALSE,
                           wts = NULL, plus_group = FALSE,
                           ndraws = 0, nage = 30, plot = TRUE,
                           mod    = get("mod",    envir = parent.frame()),
                           bins   = get("bins",   envir = parent.frame()),
                           datain = get("datain", envir = parent.frame()),
                           goodts = get("goodts", envir = parent.frame())) {

  ## ---- target bins and native -> target coarsening ----
  tL    <- seq(LowLB, UpLB, Gap)
  nb    <- length(tL)
  tbins <- list(lbinL = tL, lbinU = tL + Gap, lbin = tL + Gap / 2)

  brks    <- c(tL, UpLB + Gap)
  to_map  <- findInterval(bins$lbin, brks)          # 0 below, nb+1 above
  if (plus_group) to_map[to_map == nb + 1] <- nb
  to_map[to_map > nb] <- 0
  from_map <- findInterval(bins$lbin, brks)
  from_map[from_map > nb] <- 0

  Agg <- outer(seq_len(nb), to_map, "==") * 1       # nb x nnat: sums to-bins
  Fm  <- outer(from_map, seq_len(nb), "==") * 1     # nnat x nb: averages from-bins
  if (any(colSums(Fm) == 0))
    stop("Some target bins contain no native Tag2STM bin; check LowLB/UpLB/Gap against bins$lbin")
  Fm  <- sweep(Fm, 2, colSums(Fm), "/")

  coarsen <- function(G) {
    C <- Agg %*% G %*% Fm
    sweep(C, 2, colSums(C), "/")
  }

  ## ---- season order for the 12-month product ----
  if (is.null(start)) start <- goodts[1]
  k <- match(start, goodts)
  if (is.na(k)) stop("'start' must be one of goodts")
  ord <- goodts[c(k:length(goodts), seq_len(k - 1))]

  compound <- function(stm_full, yr = NULL) {
    A <- diag(dim(stm_full)[1])
    for (s in ord) {
      G <- if (is.null(yr)) stm_full[, , s] else stm_full[, , s, yr]
      A <- G %*% A                                  # column convention, as ClipSTM
    }
    A
  }

  annual_from_stm <- function(stm_full) {
    if (length(dim(stm_full)) == 3) return(list(all = coarsen(compound(stm_full))))
    yrs <- if (is.null(years)) which(datain$yr_supported) else years
    if (any(!datain$yr_supported[yrs]))
      warning("Unsupported years included: they are fixed at S = 0, not estimated")
    Ay <- lapply(yrs, function(y) coarsen(compound(stm_full, y)))
    names(Ay) <- paste0("yr", yrs)
    if (by_year) Ay else list(avg = Reduce(`+`, Ay) / length(Ay))
  }

  ## ---- weights ----
  if (is.null(wts)) wts <- rep(1, nb)
  if (length(wts) != nb) stop("wts must have one value per target bin (", nb, ")")
  wts <- wts / sum(wts)

  ## ---- parameter flattening ----
  is_est <- names(par_template) %in% est
  if (!any(is_est)) stop("No names in 'est' match par_template")
  theta0 <- unlist(par_template[is_est])
  build  <- function(theta) {
    p <- par_template
    p[is_est] <- utils::relist(theta, par_template[is_est])
    p
  }

  ## ---- projection fit ----
  eps <- 1e-12
  kl  <- function(A, Ahat) sum(wts * colSums(A * (log(A + eps) - log(Ahat + eps))))

  fit_one <- function(A, th_start) {
    obj <- function(th) {
      Ahat <- stm_fun(build(th), tbins)
      if (any(!is.finite(Ahat))) return(1e10)
      kl(A, Ahat)
    }
    ctrl <- list(eval.max = 4000, iter.max = 2000)
    o1 <- nlminb(th_start, obj, control = ctrl)
    o2 <- nlminb(o1$par, obj, control = ctrl)       # cheap restart check
    o  <- if (o2$objective < o1$objective) o2 else o1
    list(theta = o$par, objective = o$objective, convergence = o$convergence,
         message = o$message, pars = build(o$par),
         Ahat = stm_fun(build(o$par), tbins))
  }

  ## ---- diagnostics ----
  moments <- function(M) {
    m  <- colSums(M * tbins$lbin)
    sd <- sqrt(colSums(M * outer(tbins$lbin, m, "-")^2))
    data.frame(L = tbins$lbin, inc = m - tbins$lbin, sd = sd, p_stay = diag(M))
  }
  traj <- function(M) {
    v <- c(1, rep(0, nb - 1)); out <- numeric(nage)
    for (a in seq_len(nage)) { out[a] <- sum(v * tbins$lbin); v <- M %*% v }
    out
  }

  ## ---- fit at the MLE ----
  A_list <- annual_from_stm(mod$report()$stm)
  res <- lapply(names(A_list), function(nm) {
    A <- A_list[[nm]]
    f <- fit_one(A, theta0)
    mA <- moments(A); mH <- moments(f$Ahat)
    f$A <- A
    f$moments <- rbind(cbind(mA, source = "compounded A"),
                       cbind(mH, source = "annual fit"))
    f$traj <- data.frame(age = rep(seq_len(nage), 2),
                         L = c(traj(A), traj(f$Ahat)),
                         source = rep(c("compounded A", "annual fit"), each = nage))
    if (f$convergence != 0) warning(nm, ": nlminb did not converge (", f$message, ")")
    f
  })
  names(res) <- names(A_list)

  ## ---- propagate Tag2STM uncertainty ----
  if (ndraws > 0) {
    if (length(mod$env$random) > 0) stop("Parameter draws assume no random effects")
    lp <- mod$env$last.par.best
    on.exit(invisible(mod$report(lp)), add = TRUE)  # leave mod at its MLE for ClipSTM etc.
    sdr <- RTMB::sdreport(mod)
    th_d <- MASS::mvrnorm(ndraws, sdr$par.fixed, sdr$cov.fixed)
    draw_theta <- lapply(seq_len(ndraws), function(d) {
      Ad <- annual_from_stm(mod$report(th_d[d, ])$stm)
      lapply(names(Ad), function(nm) fit_one(Ad[[nm]], res[[nm]]$theta)$theta)
    })
    for (j in seq_along(res)) {
      TH <- do.call(rbind, lapply(draw_theta, `[[`, j))
      colnames(TH) <- names(theta0)
      res[[j]]$theta_draws <- TH
      res[[j]]$prior <- data.frame(par  = names(theta0),
                                   mle  = res[[j]]$theta,
                                   mean = colMeans(TH),
                                   sd   = apply(TH, 2, sd),
                                   row.names = NULL)
    }
  }

  ## ---- plots ----
  if (plot) {
    for (nm in names(res)) {
      mm <- tidyr::pivot_longer(res[[nm]]$moments, c(inc, sd, p_stay),
                                names_to = "stat", values_to = "value")
      print(
        ggplot(mm, aes(L, value, colour = source)) +
          geom_line(linewidth = 0.8) +
          facet_wrap(~stat, scales = "free_y", ncol = 1,
                     labeller = as_labeller(c(inc = "Mean annual increment (mm)",
                                              sd = "SD of length after one year (mm)",
                                              p_stay = "P(remain in bin)"))) +
          scale_colour_manual(values = c("compounded A" = "#185FA5", "annual fit" = "#D85A30")) +
          labs(x = "Length (mm)", y = NULL, colour = NULL,
               title = paste0("Annual STM projection: ", nm)) +
          theme_minimal(base_size = 12) +
          theme(legend.position = "top", panel.grid.minor = element_blank())
      )
      print(
        ggplot(res[[nm]]$traj, aes(age, L, colour = source)) +
          geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
          scale_colour_manual(values = c("compounded A" = "#185FA5", "annual fit" = "#D85A30")) +
          labs(x = "Annual time step", y = "Mean length (mm)", colour = NULL,
               title = paste0("Mean length-at-age: ", nm)) +
          theme_minimal(base_size = 12) +
          theme(legend.position = "top", panel.grid.minor = element_blank())
      )
    }
  }

  res
}
