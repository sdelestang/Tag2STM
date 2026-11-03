#' Build One Size Transition Matrix from growmodPar Parameters
#'
#' Plain-R (non-AD) copy of the STM construction in \code{growmodPar}
#' (\code{Pmoult_fn}, the double-logistic growth-at-length block and
#' \code{make_stm}) for a single timestep, with no year effect. Used to
#' rebuild per-timestep STMs from stored parameters and to fit annual
#' parameters in \code{AnnualFromSeasonPars}.
#'
#' \strong{Keep in sync with growmodPar.} This is a second copy of the
#' construction. Check it against a fitted model whenever growmodPar
#' changes:
#' \preformatted{
#' rp <- mod$report()
#' ns <- goodts[1]
#' p  <- c(rp$Growth_par[ns, ], rp$LsigGrow, rp$Pmoult_par[ns, ])
#' range(STMfromPars(p, bins, floor = rp$mpy_floor[ns]) - rp$stm[, , ns])  # ~0
#' }
#' (For a TemporalGrowth fit compare against a year with S = 0.)
#'
#' @param pars Numeric length 7, in the order growmodPar stores them:
#'   log(Amax), P2, log(P1), log(P3), LsigGrow, Pmoult intercept,
#'   log(-Pmoult slope). Names are ignored; order is what matters.
#' @param bins List with lbinL, lbinU, lbin (as from MakeLbin).
#' @param floor Minimum moult probability for this timestep (growmodPar's
#'   \code{mpy_floor[ns]}). 0 = no floor.
#' @param Pmin Floor on the logistic scales P1 and P3, as growmodPar's
#'   \code{datain$Growth_Pmin}: P = Pmin + exp(par). Default NULL = one
#'   width of \code{bins}, matching growmodPar's default when the bins are
#'   the fitted bins. Pass \code{rep$Growth_Pmin} when they are not, and
#'   \code{0} for parameters from fits made before the floor existed.
#' @param max_moults,moult_tau As \code{datain$max_moults} and
#'   \code{datain$moult_tau} in the fit. With \code{max_moults > 1},
#'   parameters 6-7 describe the expected number of moults per window, not a
#'   probability. Default 1 (the single-moult hurdle).
#' @param n_pmoult1 Number of smallest bins with Pmoult fixed at exactly 1
#'   (\code{datain$n_pmoult1}, default 1).
#' @param P5 Fixed swap width (\code{datain$Growth_P5_fixed}, default 0.1).
#'
#' @return nlbin x nlbin column-stochastic matrix (columns = from-length).
#' @export
STMfromPars <- function(pars, bins, floor = 0, n_pmoult1 = 1, P5 = 0.1,
                        Pmin = NULL, max_moults = 1, moult_tau = 0.01) {
  pars  <- as.numeric(pars)
  lbin  <- bins$lbin; lbinL <- bins$lbinL; lbinU <- bins$lbinU
  n     <- length(lbin)
  ## Floor on P1/P3, as growmodPar (datain$Growth_Pmin). The default, one
  ## width of THESE bins, matches growmodPar's default only when bins are
  ## the fit's bins -- otherwise pass Pmin = rep$Growth_Pmin.
  if (is.null(Pmin)) Pmin <- median(diff(lbin))

  ## growth-at-length: double logistic, as growmodPar
  Amax  <- exp(pars[1])
  P2    <- pars[2]
  P1    <- Pmin + exp(pars[3])
  P3    <- Pmin + exp(pars[4])
  xdev  <- lbin - P2
  swap1 <- 1 / (1 + exp(xdev / P5))
  g     <- Amax * (1 / (1 + exp(xdev / P1)) * swap1 + 1 / (1 + exp(xdev / P3)) * (1 - swap1))
  sdg   <- exp(pars[5]) * g

  ## moult probability (max_moults = 1) or expected moults per window
  ## (max_moults > 1): floor + logistic. Smallest bins fixed at 1 only in
  ## the single-moult case, as growmodPar.
  pm <- floor + (max_moults - floor) * plogis(pars[6] - exp(pars[7]) * lbin)
  if (max_moults == 1) pm[seq_len(min(n_pmoult1, n))] <- 1

  ## one-moult matrix: columns truncated at the from-bin, plus group at the
  ## top bin
  G <- matrix(0, n, n)
  for (fm in seq_len(n)) {
    k  <- fm:n
    mu <- lbin[fm] + g[fm]
    lo <- pnorm(lbinL[k], mu, sdg[fm])
    p  <- pnorm(lbinU[k], mu, sdg[fm]) - lo
    p[length(k)] <- 1 - lo[length(k)]
    G[k, fm] <- p / sum(p)
  }
  moult_mix(G, pm, max_moults, moult_tau)
}

## Numeric copy of growmodPar's moult_wts(): P(N = 0..max_moults moults in a
## window) given the expected number lam. Keep the two in sync.
moult_count_wts <- function(lam, max_moults = 1, moult_tau = 0.01) {
  if (max_moults == 1) return(c(1 - lam, lam))
  ramp <- function(x) moult_tau * (log1p(exp(x / moult_tau)) -
                                   log1p(exp((x - 1) / moult_tau)))
  ge <- c(1, ramp(lam - seq_len(max_moults) + 1), 0)   # P(N >= 0..max_moults+1)
  ge[-length(ge)] - ge[-1]
}

## Window STM from the one-moult matrix G and per-column lambda: column fm is
## sum_k P(N = k) * (G^k)[, fm]. Numeric mirror of growmodPar's win_stm().
moult_mix <- function(G, lam, max_moults = 1, moult_tau = 0.01) {
  n  <- ncol(G)
  Gk <- vector("list", max_moults)
  Gk[[1]] <- G
  if (max_moults > 1) for (k in 2:max_moults) Gk[[k]] <- G %*% Gk[[k - 1]]
  A <- matrix(0, n, n)
  for (fm in seq_len(n)) {
    w <- moult_count_wts(lam[fm], max_moults, moult_tau)
    col <- numeric(n); col[fm] <- w[1]
    for (k in seq_len(max_moults)) col <- col + w[k + 1] * Gk[[k]][, fm]
    A[, fm] <- col
  }
  A
}
