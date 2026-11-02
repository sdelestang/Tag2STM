#' Compare the 7-parameter (moult hurdle) and 5-parameter (no hurdle) fits
#'
#' Builds and fits \code{\link{growmodPar}} twice on the same \code{datain}:
#' once with \code{MakepinPar(moult_hurdle = TRUE)} (4 growth + CV + 2
#' moult-probability parameters per season) and once with
#' \code{moult_hurdle = FALSE} (4 growth + CV; every animal moults at every
#' \code{goodts} opportunity). Everything else -- data, smoother, sigError,
#' ident_wt, suppression, temporal settings, map options -- is held
#' identical, so the two fits differ only in the hurdle.
#'
#' @param datain Data list from \code{\link{Makedata}}. Assigned into
#'   \code{.GlobalEnv}, as \code{make_growmod_obj} does.
#' @param pin_args Named list of extra arguments to \code{\link{MakepinPar}}
#'   (anything except \code{moult_hurdle}), applied to both fits.
#' @param map_args Named list of extra arguments to \code{\link{MakemapPar}}
#'   (anything except \code{pin}), applied to both fits.
#' @param Like Likelihood form, passed to \code{make_growmod_obj}.
#' @param nstart Number of starts per model. Start 1 is the MakepinPar
#'   values; the rest jitter every free parameter by
#'   \code{rnorm(sd = jitter)}. The best objective is kept. Multiple local
#'   optima are known in the suppression model, so use more than one.
#' @param jitter SD of the start jitter (on the estimation scale).
#' @param seed RNG seed for the jitter.
#' @param control \code{nlminb} control list.
#' @param plot Logical. Draw the comparison plots.
#'
#' @details
#' \strong{What to compare.} The objectives are KL pseudo-likelihoods
#' weighted by \code{nlob}, not true log-likelihoods, so a difference is NOT
#' a likelihood-ratio test and AIC is not strictly valid. They are
#' comparable as a measure of fit because both fits share the same data and
#' settings. Compare \code{NLL_data} (the KL terms plus the identification
#' mixture), not the objective: the 7-parameter objective also carries
#' \code{PenPmoult}, which the 5-parameter model does not have.
#'
#' The more informative comparisons are structural:
#' \itemize{
#'   \item \code{byK}: observed vs predicted mean increment and mean
#'     per-record fit by number of moult opportunities. A hurdle that
#'     matters shows up as the 5-parameter model overpredicting at K = 1
#'     (it forces a moult) and/or mis-scaling at K >= 2.
#'   \item \code{bysize}: mean annual increment, annual P(no growth) and
#'     expected moults per year by size, from the compounded annual STM.
#'     If Pmoult in the 7-parameter fit sits at ~1 everywhere the two
#'     annual curves should overlay and the extra two parameters are
#'     buying nothing.
#' }
#'
#' With more than one \code{goodts} window the 5-parameter model forces a
#' moult in every window at every size; a large animal that moults once a
#' year can then only be represented by one window's increment shrinking
#' toward zero. Keep \code{Growth_shared = FALSE} for that case.
#'
#' @return A list: \code{summary} (one row per model), \code{byK},
#'   \code{bysize}, and \code{fits} (per model: \code{obj}, \code{opt},
#'   \code{rep}, \code{all_starts}).
#'
#' @seealso \code{\link{MakepinPar}}, \code{\link{MakemapPar}},
#'   \code{\link{growmodPar}}
#' @export
CompareMoultHurdle <- function(datain, pin_args = list(), map_args = list(),
                               Like = 1, nstart = 3, jitter = 0.25, seed = 1,
                               control = list(eval.max = 4000, iter.max = 2000),
                               plot = TRUE) {
  assign("datain", datain, envir = .GlobalEnv)
  pin_args$moult_hurdle <- NULL
  map_args$pin <- NULL

  ## --- Fit one mode, multi-start -------------------------------------------
  fit_one <- function(hurdle) {
    pin <- do.call(MakepinPar, c(pin_args, list(moult_hurdle = hurdle)))
    map <- do.call(MakemapPar, c(list(pin = pin), map_args))
    obj <- make_growmod_obj(pin, datain, map = map, Like = Like,
                            model = "growmodPar")
    obj$env$tracemgc <- FALSE

    set.seed(seed)
    p0 <- obj$par
    starts <- vector("list", nstart)
    for (i in seq_len(nstart)) {
      st <- if (i == 1) p0 else p0 + rnorm(length(p0), 0, jitter)
      starts[[i]] <- tryCatch(
        nlminb(st, obj$fn, obj$gr, control = control),
        error = function(e) list(objective = Inf, message = conditionMessage(e))
      )
    }
    objs <- vapply(starts, function(o) as.numeric(o$objective), numeric(1))
    if (all(!is.finite(objs))) stop("All starts failed for moult_hurdle = ", hurdle)
    best <- starts[[which.min(objs)]]

    rep  <- obj$report(best$par)
    pl   <- obj$env$parList(best$par)
    list(obj = obj, opt = best, rep = rep, pl = pl,
         grad = max(abs(obj$gr(best$par))),
         all_starts = objs)
  }

  cat("\n=== Fitting 7-parameter model (moult hurdle) ===\n")
  f7 <- fit_one(TRUE)
  cat("\n=== Fitting 5-parameter model (no hurdle) ===\n")
  f5 <- fit_one(FALSE)
  fits <- list(hurdle = f7, no_hurdle = f5)

  ## --- Summary table --------------------------------------------------------
  summ <- do.call(rbind, lapply(names(fits), function(nm) {
    f <- fits[[nm]]
    data.frame(
      model       = nm,
      n_free      = length(f$opt$par),
      objective   = f$opt$objective,
      NLL_data    = as.numeric(f$rep$NLL_data),
      penalties   = f$opt$objective - as.numeric(f$rep$NLL_data),
      max_grad    = f$grad,
      iterations  = f$opt$iterations,
      convergence = f$opt$convergence,
      starts_at_best = sum(abs(f$all_starts - f$opt$objective) < 0.01),
      stringsAsFactors = FALSE
    )
  }))
  summ$dNLL_data <- summ$NLL_data - min(summ$NLL_data)

  ## --- Per-record diagnostics by moult opportunities (K) --------------------
  ## Same timestep sequence growmodPar builds (the temporal branch gives the
  ## same seasonal cycle, only the year index differs).
  nobs <- length(datain$Rccl)
  K <- vapply(seq_len(nobs), function(r) {
    nts <- datain$tsteps[r]
    if (nts <= 0) return(0L)
    tsv <- c(datain$relts[r]:datain$ntsteps, rep(1:datain$ntsteps, 30))[1:nts]
    sum(tsv %in% datain$goodts)
  }, integer(1))
  Kc  <- ifelse(K >= 5, "5+", as.character(K))
  obs <- datain$Rccl - datain$Rlcl
  w   <- datain$nlob

  pred_inc <- function(f) {
    sigE <- exp(f$pl$LsigError)
    rel  <- datain$Rlcl + f$pl$MerrorRel
    relmean <- vapply(seq_len(nobs), function(r) {
      d <- pnorm(datain$lbinU, rel[r], sigE) - pnorm(datain$lbinL, rel[r], sigE)
      sum(d * datain$lbin) / sum(d)
    }, numeric(1))
    as.numeric(f$rep$EstRecLen %*% datain$lbin) - relmean
  }
  p7 <- pred_inc(f7); p5 <- pred_inc(f5)

  byK <- do.call(rbind, lapply(sort(unique(Kc)), function(k) {
    i <- Kc == k
    wm <- function(x) sum(x[i] * w[i]) / sum(w[i])
    data.frame(
      K          = k,
      n          = sum(i),
      obs_inc    = wm(obs),
      pred_7     = wm(p7),
      pred_5     = wm(p5),
      resid_7    = wm(obs) - wm(p7),
      resid_5    = wm(obs) - wm(p5),
      meanLL_7   = mean(f7$rep$LL[i] / w[i]),
      meanLL_5   = mean(f5$rep$LL[i] / w[i]),
      sumNLL_7   = -sum(f7$rep$LL[i]) - datain$ident_wt * sum(f7$rep$IdentLL[i]),
      sumNLL_5   = -sum(f5$rep$LL[i]) - datain$ident_wt * sum(f5$rep$IdentLL[i]),
      stringsAsFactors = FALSE
    )
  }))
  byK$dNLL_5_minus_7 <- byK$sumNLL_5 - byK$sumNLL_7

  ## --- Annual growth by size from the compounded STM ------------------------
  ## Rebuilt from parameters with STMfromPars (base growth, no year effects),
  ## which works for both modes because growmodPar reports a saturated
  ## Pmoult_par when there is no hurdle.
  bins <- list(lbin = datain$lbin, lbinL = datain$lbinL, lbinU = datain$lbinU)
  P5   <- if (is.null(datain$Growth_P5_fixed)) 0.1 else datain$Growth_P5_fixed
  n1   <- if (is.null(datain$n_pmoult1)) 1 else datain$n_pmoult1
  annual <- function(f) {
    A <- diag(length(bins$lbin))
    for (ns in sort(datain$goodts)) {
      pars <- c(f$pl$Growth_par[ns, ], f$pl$LsigGrow, f$rep$Pmoult_par[ns, ])
      A <- STMfromPars(pars, bins, floor = f$rep$mpy_floor[ns],
                       n_pmoult1 = n1, P5 = P5,
                       Pmin = f$rep$Growth_Pmin) %*% A
    }
    data.frame(
      inc    = colSums(A * bins$lbin) - bins$lbin,
      p_stay = diag(A),
      moults = colSums(f$rep$Pmoult_vec[datain$goodts, , drop = FALSE])
    )
  }
  a7 <- annual(f7); a5 <- annual(f5)
  bysize <- data.frame(L = bins$lbin,
                       inc_7 = a7$inc, inc_5 = a5$inc,
                       pstay_7 = a7$p_stay, pstay_5 = a5$p_stay,
                       moults_yr_7 = a7$moults, moults_yr_5 = a5$moults)

  ## --- Print -----------------------------------------------------------------
  cat("\n=== Summary (compare NLL_data, not objective) ===\n")
  print(summ, digits = 5, row.names = FALSE)
  cat("\n=== By moult opportunities (K) ===\n")
  print(byK[, c("K", "n", "obs_inc", "pred_7", "pred_5", "resid_7", "resid_5",
                "dNLL_5_minus_7")], digits = 4, row.names = FALSE)
  cat("  dNLL_5_minus_7 > 0: the hurdle fits that K group better.\n")

  ## --- Plots -----------------------------------------------------------------
  if (plot) {
    op <- par(mfrow = c(1, 3), mar = c(4, 4, 2.5, 1))
    on.exit(par(op), add = TRUE)
    col7 <- "#1b6ca8"; col5 <- "#c0392b"

    rng <- range(c(bysize$inc_7, bysize$inc_5), na.rm = TRUE)
    plot(bysize$L, bysize$inc_7, type = "l", lwd = 2, col = col7, ylim = rng,
         xlab = "Length (mm)", ylab = "Mean annual increment (mm)",
         main = "Annual growth")
    lines(bysize$L, bysize$inc_5, lwd = 2, lty = 2, col = col5)
    legend("topright", c("7-par (hurdle)", "5-par (no hurdle)"),
           col = c(col7, col5), lty = 1:2, lwd = 2, bty = "n")

    plot(bysize$L, bysize$moults_yr_7, type = "l", lwd = 2, col = col7,
         ylim = c(0, max(1, bysize$moults_yr_5)),
         xlab = "Length (mm)", ylab = "Expected moults per year",
         main = "Moult frequency")
    lines(bysize$L, bysize$moults_yr_5, lwd = 2, lty = 2, col = col5)

    rr <- range(c(byK$resid_7, byK$resid_5, 0))
    x  <- seq_len(nrow(byK))
    plot(x - 0.1, byK$resid_7, pch = 19, col = col7, ylim = rr, xaxt = "n",
         xlab = "Moult opportunities (K)", ylab = "Obs - pred mean increment (mm)",
         main = "Residual by K", xlim = c(0.5, max(x) + 0.5))
    points(x + 0.1, byK$resid_5, pch = 17, col = col5)
    axis(1, at = x, labels = byK$K)
    abline(h = 0, lty = 3)
  }

  invisible(list(summary = summ, byK = byK, bysize = bysize, fits = fits))
}
