#' Project Per-timestep Growth Parameters onto One Annual Parameter Set
#'
#' Rebuilds each timestep's STM from its fitted growmodPar parameters on the
#' Tag2STM length bins, compounds them into a 12-month STM, clips to the
#' IMuLT length range, then fits a single annual parameter set whose
#' one-step STM best reproduces the compounded matrix. The objective is a
#' weighted KL divergence summed over from-length columns.
#'
#' Parameter-only counterpart of MakeAnnualPars(): needs no fitted model
#' object, so it runs inside the timestep-selection script from paroutlist.
#' Works from mean growth (no year effects), like-for-like with the
#' per-timestep parameters written to IMuLT.
#'
#' @param season_pars List of length-7 parameter vectors, one per timestep,
#'   in growmodPar order: log(Amax), P2, log(P1), log(P3), LsigGrow,
#'   Pmoult intercept, log(-Pmoult slope).
#' @param tsteps Integer timestep index of each element of season_pars.
#' @param floors Per-timestep moult-probability floors (growmodPar's
#'   mpy_floor for those timesteps). Default 0.
#' @param stm_fun STM builder, function(pars, bins, floor, ...). Default
#'   STMfromPars.
#' @param partypes Parameter names, in order (labels only).
#' @param lower,upper Bounds for each parameter.
#' @param est Names of parameters to estimate; the rest stay at start values.
#' @param bins Tag2STM length bins (from MakeLbin), on IMuLT's bin width.
#'   Their wider span is used for compounding.
#' @param LowLB,UpLB Lower edges of IMuLT's first and last bins.
#' @param start Timestep at which IMuLT's growth year begins.
#' @param ann_floor Moult-probability floor for the annual STM. Default 0.
#' @param wts Weights over IMuLT from-bins; NULL = uniform.
#' @param plus_group TRUE accumulates mass above UpLB into the top IMuLT bin;
#'   FALSE drops and renormalises (ClipSTM behaviour).
#' @param ... Passed to stm_fun (e.g. n_pmoult1, P5).
#'
#' @return list(pars, objective, convergence, message, A, Ahat, moments)
#' @export
AnnualFromSeasonPars <- function(season_pars, tsteps, floors = 0,
                                 stm_fun = STMfromPars, partypes,
                                 lower, upper, est = partypes, bins,
                                 LowLB, UpLB, start = min(tsteps),
                                 ann_floor = 0, wts = NULL, plus_group = FALSE, ...) {
  if (length(season_pars) != length(tsteps)) stop("One parameter vector per timestep required")
  floors <- rep_len(floors, length(tsteps))

  ## ---- IMuLT bins: a subset of the Tag2STM bins ----
  Gap    <- as.numeric(names(which.max(table(round(diff(bins$lbinL), 8)))))
  tL     <- seq(LowLB, UpLB, Gap)
  tokeep <- match(tL, bins$lbinL)
  if (anyNA(tokeep)) stop("LowLB/UpLB do not fall on the bins grid (width ", Gap, " mm)")
  if (any(abs(diff(bins$lbinL[c(tokeep, max(tokeep) + 1)]) - Gap) > 1e-8, na.rm = TRUE))
    stop("bins are not a constant ", Gap, " mm within the IMuLT range")
  nb <- length(tL)
  Lm <- bins$lbin[tokeep]

  clip <- function(G) {
    C <- G[tokeep, tokeep]
    if (plus_group) {
      above <- bins$lbinL > UpLB
      if (any(above)) C[nb, ] <- C[nb, ] + colSums(G[above, tokeep, drop = FALSE])
    }
    sweep(C, 2, colSums(C), "/")
  }

  ## ---- compound the timesteps in calendar order from 'start' ----
  o <- order(tsteps)
  tsteps <- tsteps[o]; season_pars <- season_pars[o]; floors <- floors[o]
  k <- which(tsteps >= start)[1]
  if (is.na(k)) k <- 1
  ord <- c(k:length(tsteps), seq_len(k - 1))
  A <- diag(length(bins$lbinL))
  for (s in ord) A <- stm_fun(season_pars[[s]], bins, floor = floors[s], ...) %*% A
  A <- clip(A)

  ## ---- starting values ----
  ## Shape from the season with the largest Amax; annual Amax = sum of the
  ## seasonal Amax (AveGrowth is log(Amax), so sum on the natural scale).
  lower <- setNames(lower, partypes); upper <- setNames(upper, partypes)
  a   <- sapply(season_pars, `[`, 1)
  p0  <- setNames(as.numeric(season_pars[[which.max(a)]]), partypes)
  p0[1] <- log(sum(exp(a)))
  pad <- 1e-4 * (upper - lower)
  p0  <- pmin(pmax(p0, lower + pad), upper - pad)

  ## ---- projection fit ----
  if (is.null(wts)) wts <- rep(1, nb)
  if (length(wts) != nb) stop("wts must have one value per IMuLT bin (", nb, ")")
  wts <- wts / sum(wts)

  ie   <- partypes %in% est
  eps  <- 1e-12
  full <- function(th) { p <- p0; p[ie] <- th; p }
  ahat <- function(p) clip(stm_fun(p, bins, floor = ann_floor, ...))
  obj  <- function(th) {
    Ahat <- ahat(full(th))
    if (any(!is.finite(Ahat))) return(1e10)
    sum(wts * colSums(A * (log(A + eps) - log(Ahat + eps))))
  }
  ctrl <- list(eval.max = 4000, iter.max = 2000)
  o1 <- nlminb(p0[ie], obj, lower = lower[ie], upper = upper[ie], control = ctrl)
  o2 <- nlminb(o1$par, obj, lower = lower[ie], upper = upper[ie], control = ctrl)
  op <- if (o2$objective < o1$objective) o2 else o1

  pars <- full(op$par)
  Ahat <- ahat(pars)

  moments <- function(M, src) {
    m  <- colSums(M * Lm)
    sd <- sqrt(colSums(M * outer(Lm, m, "-")^2))
    data.frame(L = Lm, inc = m - Lm, sd = sd, p_stay = diag(M), source = src)
  }

  list(pars = pars, objective = op$objective, convergence = op$convergence,
       message = op$message, A = A, Ahat = Ahat,
       moments = rbind(moments(A, "compounded"), moments(Ahat, "annual fit")))
}
