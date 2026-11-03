#' Write Selected Tag2STM Growth Fits Out as IMuLT Growth Inputs
#'
#' Turns the selected sub-annual fits into IMuLT's growth parameter and
#' pointer tables, in one of two forms:
#' \describe{
#'   \item{\code{mode = "timestep"}}{One IMuLT matrix per moult window, using
#'     each window's fitted parameters directly. For an IMuLT model run on the
#'     same sub-annual timesteps.}
#'   \item{\code{mode = "annual"}}{The windows are compounded into one
#'     12-month STM, which is projected onto a single annual parameter set
#'     (\code{\link{AnnualFromSeasonPars}}). One IMuLT matrix per sex x area
#'     x period, for an annual-timestep IMuLT model. Optionally plots each
#'     projection against its compounded target
#'     (\code{\link{PlotAnnualProjection}}).}
#' }
#'
#' @param paroutlist Named list of fits, one element per fit name, each with
#'   \code{goodts}, \code{mpy_floor}, \code{Growth_Pmin} and \code{Pars1..n}
#'   (one 7-vector per window, growmodPar order).
#' @param selected Data frame of the selected fits, one row per sex x area x
#'   period, with columns \code{Name}, \code{sex}, \code{area}, \code{period}.
#'   An \code{area} of \code{"all"} marks a fit pooled across areas: it is
#'   pointed at from every area in \code{imult_areas}.
#' @param bins Tag2STM length bins (\code{\link{MakeLbin}}).
#' @param mode \code{"annual"} or \code{"timestep"}.
#' @param imult_areas Character vector of IMuLT area labels. Needed only when a
#'   selected fit has \code{area = "all"}.
#' @param partypes Parameter labels, in growmodPar order.
#' @param lower,upper Bounds per parameter type (IMuLT's, and the annual
#'   projection's).
#' @param ann_LowLB,ann_UpLB Lower edges of IMuLT's first and last length bins
#'   (annual mode). Default: the regular-width part of \code{bins}.
#' @param ann_start Timestep at which IMuLT's growth year begins.
#' @param ann_est Parameter types the annual projection estimates; the rest
#'   stay at the sub-annual fit's values. A character vector, or
#'   \code{function(sex, area, period)} returning one. Default: all.
#' @param ann_cv_cap Logical, default TRUE. Caps the annual CV (LsigGrow) in
#'   the projection at the larger of the per-moult CV and the CV of a year's
#'   growth among animals that moulted (see Details). Has no effect unless
#'   LsigGrow is in \code{ann_est}.
#' @param ann_floor Moult-probability floor for the annual STM.
#' @param ann_wts NULL, a weight vector over IMuLT bins, or
#'   \code{function(sex, area, period)} returning one.
#' @param ann_plus TRUE if IMuLT's top bin is a plus group.
#' @param n_pmoult1,P5 The fits' \code{datain$n_pmoult1} and
#'   \code{datain$Growth_P5_fixed}.
#' @param imult_est Named list by sex of parameter types IMuLT should
#'   re-estimate, e.g. \code{list(F = c("AveGrowth", "Inflection", "Slope1",
#'   "Slope2"))}. These get \code{Phase = est_phase}, \code{UsePrior = 1} and
#'   prior sd from \code{prior_sd_est}. All others are fixed (\code{Phase =
#'   fixed_phase}).
#' @param est_phase,fixed_phase IMuLT phases for freed and fixed parameters.
#' @param prior_sd_est Named vector of absolute prior sds for freed parameter
#'   types, on each parameter's own scale.
#' @param prior_frac,prior_min Prior sd for all other rows: \code{prior_frac}
#'   x |estimate|, floored at \code{prior_min} (recorded, unused while
#'   \code{UsePrior = 0}).
#' @param startseason,endseason Season range written to the pointer table.
#' @param plot Annual mode: draw \code{\link{PlotAnnualProjection}} for each
#'   matrix.
#' @param plot_file Annual mode: if given, also write those plots to this PDF,
#'   one page per matrix.
#' @param write Write \code{Growth_Parameters_IMuLT_<tag>.csv} and
#'   \code{Growth_Pointers_IMuLT_<tag>.csv} to \code{outdir}.
#' @param outdir Output directory.
#' @param tag File tag. Default \code{"annual_<binw>mm"} or
#'   \code{"<n>ts_<binw>mm"}.
#'
#' @details
#' \strong{The annual CV.} A year's growth is the sum of however many moults
#' fell in it, so its spread differs from one moult's: smaller when animals
#' moult several times a year (about 1/sqrt(moults) of it), similar when they
#' moult once. The annual projection therefore estimates the CV. But when
#' annual growth is smaller than a bin (slow growers), the annual STM cannot
#' identify it, and the optimiser runs to a degenerate solution -- Amax near
#' 0 and CV very large -- that reproduces the stay/next-bin split without
#' meaning anything. That failure always pushes the CV up, so the cap blocks
#' it without constraining well-identified cases: the cap is the larger of
#' the per-moult CV and the CV of a year's growth among animals that moulted
#' at least once, computed from the fitted windows (independent moults, one
#' increment per moult), averaged over the IMuLT length range weighted by the
#' probability of moulting. Where growth is informative the estimate sits
#' below the cap; where it is not, it stops at the cap, which is then the
#' best available value. Crabs moulting at most once a year: cap = per-moult
#' CV, which is what a free fit returns anyway.
#'
#' @return Invisibly, a list: \code{Params}, \code{pertab} (as written, without
#'   bookkeeping columns), \code{annlist} (annual mode: one
#'   \code{AnnualFromSeasonPars} result per matrix), \code{tag}.
#'
#' @seealso \code{\link{AnnualFromSeasonPars}}, \code{\link{PlotAnnualProjection}}
#' @export
GrowthToIMuLT <- function(paroutlist, selected, bins, mode = c("annual", "timestep"),
                          imult_areas = NULL,
                          partypes = c('AveGrowth', 'Inflection', 'Slope1', 'Slope2',
                                       'LsigGrow', 'Pmoultinf', 'Pmoultslp'),
                          lower = c(-5,   0, -5, -5, -5, -20, -15),
                          upper = c( 4, 200,  6, 15,  4,  20,   3),
                          ann_LowLB = NULL, ann_UpLB = NULL, ann_start = 1,
                          ann_est = partypes, ann_cv_cap = TRUE,
                          ann_floor = 0, ann_wts = NULL, ann_plus = FALSE,
                          n_pmoult1 = 1, P5 = 0.1,
                          imult_est = list(), est_phase = 3, fixed_phase = -4,
                          prior_sd_est = c(AveGrowth = 0.3, Inflection = 10,
                                           Slope1 = 0.5, Slope2 = 0.5),
                          prior_frac = 0.1, prior_min = 0.01,
                          startseason = 1900, endseason = 2050,
                          plot = TRUE, plot_file = NULL,
                          write = TRUE, outdir = ".", tag = NULL) {
  mode <- match.arg(mode)
  npar <- length(partypes)
  if (length(lower) != npar || length(upper) != npar)
    stop("lower and upper need one value per parameter type (", npar, ").")
  need <- c("Name", "sex", "area", "period")
  if (!all(need %in% names(selected)))
    stop("selected needs columns: ", paste(need, collapse = ", "))
  selected$area <- as.character(selected$area)
  if (any(selected$area == "all") && is.null(imult_areas))
    stop("A selected fit is pooled (area = \"all\"): supply imult_areas.")
  miss <- setdiff(selected$Name, names(paroutlist))
  if (length(miss)) stop("Selected fits not in paroutlist: ", paste(miss, collapse = ", "))

  ## --- Bins: modal width and the regular-width part of the grid -----------
  binw    <- as.numeric(names(which.max(table(round(diff(bins$lbinL), 8)))))
  reg_lbL <- bins$lbinL[abs((bins$lbinU - bins$lbinL) - binw) < 1e-6]
  if (mode == "annual") {
    if (is.null(ann_LowLB)) ann_LowLB <- min(reg_lbL)
    if (is.null(ann_UpLB))  ann_UpLB  <- max(reg_lbL)
    want <- seq(ann_LowLB, ann_UpLB, binw)
    bad  <- want[is.na(match(round(want, 6), round(reg_lbL, 6)))]
    if (length(bad))
      stop("IMuLT range ", ann_LowLB, "-", ann_UpLB, " includes edges not on the ", binw,
           " mm grid: ", paste(bad, collapse = ", "),
           "\nValid lower edges: ", paste(reg_lbL, collapse = ", "))
  }
  if (is.null(tag)) {
    nts_max <- max(vapply(selected$Name, function(nm) length(paroutlist[[nm]]$goodts), integer(1)))
    tag <- paste0(if (mode == "annual") "annual" else paste0(nts_max, "ts"), "_", binw, "mm")
  }

  ## --- Pointer table --------------------------------------------------------
  ## One row per sex x IMuLT area x period x window ("timestep"), or per sex x
  ## IMuLT area x period ("annual"). Timesteps that are not moult windows get
  ## no row. A pooled fit (area "all") gets a row for every IMuLT area, all
  ## pointing at the same matrix.
  pertab <- do.call(rbind, lapply(seq_len(nrow(selected)), function(i) {
    f    <- selected[i, ]
    tst  <- if (mode == "annual") 1L else paroutlist[[f$Name]]$goodts
    base <- paste0(f$sex, "_a", f$area, "_p", f$period)
    ar   <- if (f$area == "all") as.character(imult_areas) else f$area
    pt   <- expand.grid(startseason = startseason, endseason = endseason, period = f$period,
                        sex = f$sex, area = ar, tstep = tst, grow = 1, stringsAsFactors = FALSE)
    pt$matrix  <- if (mode == "annual") paste0(base, "_ann") else paste0(base, "_ts", pt$tstep)
    pt$fit     <- f$Name
    pt$fitarea <- f$area
    pt
  }))
  pertab <- pertab[order(pertab$sex, pertab$area, pertab$period, pertab$tstep), ]
  rownames(pertab) <- NULL

  ## --- Parameter table: one set per matrix ----------------------------------
  mats <- unique(pertab[, c("matrix", "fit", "sex", "fitarea", "period", "tstep")])
  nmat <- nrow(mats)
  Params <- data.frame(lower       = rep(lower, nmat),
                       upper       = rep(upper, nmat),
                       est         = NA_real_,
                       Phase       = fixed_phase,
                       Link        = 0,
                       UsePrior    = 0,
                       Prior       = NA_real_,
                       Priorsd     = NA_real_,
                       description = rep(mats$matrix, each = npar),
                       Partype     = rep(partypes, nmat))

  estflag <- rep(FALSE, nrow(Params))   # rows the annual projection estimated
  annlist <- list()

  for (i in seq_len(nmat)) {
    Id      <- mats$matrix[i]
    fitpars <- paroutlist[[mats$fit[i]]]
    rows    <- (i - 1) * npar + seq_len(npar)

    if (mode == "timestep") {
      pars <- as.numeric(fitpars[[paste0("Pars", match(mats$tstep[i], fitpars$goodts))]])
    } else {
      tst   <- fitpars$goodts
      spars <- lapply(seq_along(tst), function(k) as.numeric(fitpars[[paste0("Pars", k)]]))
      w     <- if (is.function(ann_wts)) ann_wts(mats$sex[i], mats$fitarea[i], mats$period[i]) else ann_wts
      flr   <- if (is.null(fitpars$mpy_floor)) 0 else fitpars$mpy_floor[tst]
      Pmin  <- if (is.null(fitpars$Growth_Pmin)) 0 else fitpars$Growth_Pmin
      est_i <- if (is.function(ann_est)) ann_est(mats$sex[i], mats$fitarea[i], mats$period[i]) else ann_est
      estflag[rows] <- partypes %in% est_i

      ## Cap the annual CV: see @details
      up_i <- upper
      iL   <- which(partypes == "LsigGrow")
      cap  <- NA_real_
      if (ann_cv_cap && length(iL) == 1 && "LsigGrow" %in% est_i) {
        cap <- annual_cv_cap(spars, flr, bins, Pmin = Pmin, P5 = P5, n_pmoult1 = n_pmoult1,
                             LowLB = ann_LowLB, UpLB = ann_UpLB)
        up_i[iL] <- min(upper[iL], log(cap))
      }

      ann <- AnnualFromSeasonPars(spars, tst, floors = flr, stm_fun = STMfromPars,
                                  partypes = partypes, lower = lower, upper = up_i,
                                  est = est_i, bins = bins, LowLB = ann_LowLB, UpLB = ann_UpLB,
                                  start = ann_start, ann_floor = ann_floor, wts = w,
                                  plus_group = ann_plus, n_pmoult1 = n_pmoult1, P5 = P5,
                                  Pmin = Pmin)
      ## nlminb reports "false convergence" when the minimum is very flat and
      ## the objective tiny (slow growth relative to bin width). That is not a
      ## failed fit -- check the plot. Anything else is a real warning.
      if (ann$convergence != 0) {
        if (grepl("false convergence", ann$message) && ann$objective < 0.01) {
          message(Id, ": annual projection ended with nlminb 'false convergence' at objective ",
                  signif(ann$objective, 3), " (flat minimum; check the projection plot)")
        } else {
          warning(Id, ": annual projection did not converge (", ann$message,
                  ", objective ", signif(ann$objective, 3), ")")
        }
      }
      if (!is.na(cap) && abs(ann$pars[iL] - log(cap)) < 1e-3)
        message(Id, ": annual CV at its cap (", round(cap, 2), "). Expected with one moult ",
                "window; otherwise annual growth is too small to identify it, and the cap is used")
      ann$cv_cap <- cap
      annlist[[Id]] <- ann
      pars <- as.numeric(ann$pars)
    }
    if (length(pars) != npar) stop(Id, ": expected ", npar, " parameters, found ", length(pars))
    Params$est[rows] <- round(pars, 5)
  }

  ## --- Priors and freed parameters -----------------------------------------
  Params$Prior   <- Params$est
  Params$Priorsd <- round(pmax(abs(Params$est) * prior_frac, prior_min), 5)
  row_sex <- rep(mats$sex, each = npar)
  free <- vapply(seq_len(nrow(Params)), function(j)
    Params$Partype[j] %in% imult_est[[as.character(row_sex[j])]], logical(1))
  if (any(free)) {
    nosd <- setdiff(unique(Params$Partype[free]), names(prior_sd_est))
    if (length(nosd)) stop("prior_sd_est has no value for: ", paste(nosd, collapse = ", "))
    Params$Phase[free]    <- est_phase
    Params$UsePrior[free] <- 1
    Params$Priorsd[free]  <- unname(prior_sd_est[Params$Partype[free]])
    message("Estimated in IMuLT with priors:\n",
            paste(capture.output(print(Params[free, c("description", "Partype", "est", "Priorsd", "Phase")],
                                       row.names = FALSE)), collapse = "\n"))
  }

  ## --- Checks ----------------------------------------------------------------
  miss <- is.na(Params$est)
  if (any(miss))
    warning("Missing estimates:\n",
            paste(capture.output(print(Params[miss, c("description", "Partype")])), collapse = "\n"))
  oob <- !miss & (Params$est < Params$lower | Params$est > Params$upper)
  if (any(oob))
    warning("Estimates outside bounds:\n",
            paste(capture.output(print(Params[oob, c("description", "Partype", "lower", "est", "upper")])),
                  collapse = "\n"))
  if (mode == "annual") {
    rng <- Params$upper - Params$lower
    atb <- !miss & estflag & (abs(Params$est - Params$lower) < 1e-3 * rng |
                                abs(Params$upper - Params$est) < 1e-3 * rng)
    if (any(atb))
      warning("Annual estimates at a bound (widen bounds or drop from ann_est):\n",
              paste(capture.output(print(Params[atb, c("description", "Partype", "lower", "est", "upper")])),
                    collapse = "\n"))
  }

  ## --- Plots (annual mode) -----------------------------------------------------
  if (mode == "annual" && (plot || !is.null(plot_file))) {
    plots <- lapply(names(annlist), function(nm) PlotAnnualProjection(annlist[[nm]], title = nm))
    if (plot) for (p in plots) print(p)
    if (!is.null(plot_file)) {
      grDevices::pdf(plot_file, width = 9, height = 7.5)
      for (p in plots) print(p)
      grDevices::dev.off()
    }
  }

  ## --- Write -------------------------------------------------------------------
  pertab_out <- pertab[, setdiff(names(pertab), c("period", "fit", "fitarea"))]
  if (write) {
    utils::write.csv(Params, file.path(outdir, paste0("Growth_Parameters_IMuLT_", tag, ".csv")),
                     row.names = FALSE)
    utils::write.csv(pertab_out, file.path(outdir, paste0("Growth_Pointers_IMuLT_", tag, ".csv")),
                     row.names = FALSE)
    if (mode == "annual")
      saveRDS(annlist, file.path(outdir, paste0("Growth_Annual_Projection_", tag, ".rds")))
  }

  invisible(list(Params = Params, pertab = pertab_out, annlist = annlist, tag = tag))
}


## Cap for the annual CV (see GrowthToIMuLT @details): the larger of the
## per-moult CV and the CV of a year's growth among animals that moulted at
## least once. Windows moult independently with the fitted P(moult); each
## moult adds Normal(mu_k, (cv * mu_k)^2) at the starting length. Averaged
## over the IMuLT range, weighted by P(at least one moult). Same growth and
## P(moult) formulas as STMfromPars.
annual_cv_cap <- function(season_pars, floors, bins, Pmin = 2, P5 = 0.1, n_pmoult1 = 1,
                          LowLB = min(bins$lbinL), UpLB = max(bins$lbinL)) {
  L    <- bins$lbin
  keep <- bins$lbinL >= LowLB - 1e-8 & bins$lbinL <= UpLB + 1e-8
  floors <- rep_len(floors, length(season_pars))
  gfun <- function(p) {
    x <- L - p[2]; sw <- 1 / (1 + exp(x / P5))
    exp(p[1]) * (sw / (1 + exp(x / (Pmin + exp(p[3])))) + (1 - sw) / (1 + exp(x / (Pmin + exp(p[4])))))
  }
  pfun <- function(p, fl) {
    pm <- fl + (1 - fl) * plogis(p[6] - exp(p[7]) * L)
    pm[seq_len(min(n_pmoult1, length(L)))] <- 1
    pm
  }
  K   <- length(season_pars)
  mu  <- vapply(season_pars, gfun, numeric(length(L)))
  pm  <- vapply(seq_len(K), function(k) pfun(season_pars[[k]], floors[k]), numeric(length(L)))
  cvs <- vapply(season_pars, function(p) exp(p[5]), numeric(1))
  mu  <- matrix(mu, ncol = K); pm <- matrix(pm, ncol = K)
  sub <- as.matrix(expand.grid(rep(list(0:1), K)))[-1, , drop = FALSE]   # >= 1 moult
  res <- vapply(which(keep), function(i) {
    pr <- apply(sub, 1, function(b) prod(ifelse(b == 1, pm[i, ], 1 - pm[i, ])))
    if (sum(pr) <= 0 || all(mu[i, ] <= 0)) return(c(NA, 0))
    m  <- as.numeric(sub %*% mu[i, ])
    v  <- as.numeric(sub %*% (cvs * mu[i, ])^2)
    pr <- pr / sum(pr)
    M  <- sum(pr * m); V <- sum(pr * (v + m^2)) - M^2
    c(sqrt(max(V, 0)) / M, 1 - prod(1 - pm[i, ]))
  }, numeric(2))
  ok <- is.finite(res[1, ]) & res[2, ] > 0
  cv_year <- if (any(ok)) sum(res[1, ok] * res[2, ok]) / sum(res[2, ok]) else NA
  max(cvs, cv_year, na.rm = TRUE)
}
