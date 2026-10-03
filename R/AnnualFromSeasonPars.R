#' Project Per-timestep Growth Parameters onto One Annual Parameter Set
#'
#' Rebuilds each timestep's STM from its fitted parameters with the IMuLT
#' builder on the Tag2STM length bins, compounds them into a 12-month STM,
#' clips to the IMuLT length range, then fits a single annual parameter set
#' whose one-step STM best reproduces the compounded matrix. The objective is
#' a weighted KL divergence summed over from-length columns.
#'
#' This is the parameter-only counterpart of MakeAnnualPars(). It needs no
#' fitted model object, so it can run inside the timestep-selection script
#' from paroutlist alone. It works from mean growth (no year effects), which
#' is like-for-like with the per-timestep parameters written to IMuLT.
#'
#' @param season_pars List of numeric parameter vectors, one per timestep,
#'   each in \code{partypes} order.
#' @param tsteps Integer timestep index of each element of season_pars.
#' @param stm_fun function(pars, bins) -> column-stochastic STM (columns =
#'   from-length) for one step, where pars is a named vector in partypes
#'   order and bins a list(lbinL, lbinU, lbin). Must be IMuLT's own builder.
#' @param partypes Parameter names, in order.
#' @param lower,upper Bounds for each parameter (length = length(partypes)).
#' @param est Names of parameters to estimate; the rest stay at start values.
#' @param bins Tag2STM length bins (from MakeLbin). Must be on IMuLT's bin
#'   width; their wider span is used for compounding so growth near the
#'   edges of the IMuLT range is not truncated mid-year.
#' @param LowLB,UpLB Lower edges of IMuLT's first and last bins. Must lie on
#'   the bins grid.
#' @param start Timestep at which IMuLT's growth year begins. The product
#'   starts at the first fitted timestep >= start, wrapping around.
#' @param wts Weights over IMuLT from-bins; NULL = uniform.
#' @param plus_group TRUE accumulates mass above UpLB into the top IMuLT bin;
#'   FALSE drops and renormalises (ClipSTM behaviour). Match how stm_fun
#'   treats its top bin.
#'
#' @return list(pars, objective, convergence, message, A, Ahat, moments)
#' @export
AnnualFromSeasonPars <- function(season_pars, tsteps, stm_fun, partypes,
                                 lower, upper, est = partypes, bins,
                                 LowLB, UpLB, start = min(tsteps),
                                 wts = NULL, plus_group = FALSE) {
  if (length(season_pars) != length(tsteps)) stop("One parameter vector per timestep required")

  ## ---- IMuLT bins: a subset of the Tag2STM bins ----
  ## Modal width: tolerates irregular end bins outside the IMuLT range, but
  ## every bin INSIDE the range must be on this width
  Gap    <- as.numeric(names(which.max(table(round(diff(bins$lbinL), 8)))))
  tL     <- seq(LowLB, UpLB, Gap)
  tokeep <- match(tL, bins$lbinL)
  if (anyNA(tokeep)) stop("LowLB/UpLB do not fall on the bins grid (width ", Gap, " mm)")
  if (any(abs(diff(bins$lbinL[c(tokeep, max(tokeep) + 1)]) - Gap) > 1e-8, na.rm = TRUE))
    stop("bins are not a constant ", Gap, " mm within the IMuLT range")
  nb    <- length(tL)
  tbins <- list(lbinL = tL, lbinU = tL + Gap, lbin = bins$lbin[tokeep])

  clip <- function(G) {
    C <- G[tokeep, tokeep]
    if (plus_group) {
      above <- bins$lbinL > UpLB
      if (any(above)) C[nb, ] <- C[nb, ] + colSums(G[above, tokeep, drop = FALSE])
    }
    sweep(C, 2, colSums(C), "/")
  }

  ## ---- compound the timesteps in calendar order from 'start' ----
  o <- order(tsteps); tsteps <- tsteps[o]; season_pars <- season_pars[o]
  k <- which(tsteps >= start)[1]
  if (is.na(k)) k <- 1
  ord <- c(k:length(tsteps), seq_len(k - 1))
  A <- diag(length(bins$lbinL))
  for (s in ord) A <- stm_fun(setNames(season_pars[[s]], partypes), bins) %*% A
  A <- clip(A)

  ## ---- starting values: largest-growth season's shape, summed AveGrowth ----
  lower <- setNames(lower, partypes); upper <- setNames(upper, partypes)
  big <- which.max(sapply(season_pars, `[`, 1))
  p0  <- setNames(season_pars[[big]], partypes)
  if ("AveGrowth" %in% partypes)
    p0["AveGrowth"] <- sum(sapply(season_pars, function(p) p[match("AveGrowth", partypes)]))
  pad <- 1e-4 * (upper - lower)
  p0  <- pmin(pmax(p0, lower + pad), upper - pad)

  ## ---- projection fit: annual STM built on the full bins, then clipped the
  ## ---- same way, so both sides of the comparison see identical edge handling
  if (is.null(wts)) wts <- rep(1, nb)
  if (length(wts) != nb) stop("wts must have one value per IMuLT bin (", nb, ")")
  wts <- wts / sum(wts)

  ie   <- partypes %in% est
  eps  <- 1e-12
  full <- function(th) { p <- p0; p[ie] <- th; p }
  ahat <- function(p) clip(stm_fun(p, bins))
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
    m  <- colSums(M * tbins$lbin)
    sd <- sqrt(colSums(M * outer(tbins$lbin, m, "-")^2))
    data.frame(L = tbins$lbin, inc = m - tbins$lbin, sd = sd, p_stay = diag(M), source = src)
  }

  list(pars = pars, objective = op$objective, convergence = op$convergence,
       message = op$message, A = A, Ahat = Ahat,
       moments = rbind(moments(A, "compounded"), moments(Ahat, "annual fit")))
}
