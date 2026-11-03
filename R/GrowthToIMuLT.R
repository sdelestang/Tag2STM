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
#'   stay at the sub-annual fit's values. Default fixes the CV, which the
#'   annual STM cannot identify when growth is small relative to bin width.
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
                          ann_est = setdiff(partypes, "LsigGrow"),
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
      estflag[rows] <- partypes %in% ann_est

      ann <- AnnualFromSeasonPars(spars, tst, floors = flr, stm_fun = STMfromPars,
                                  partypes = partypes, lower = lower, upper = upper,
                                  est = ann_est, bins = bins, LowLB = ann_LowLB, UpLB = ann_UpLB,
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
