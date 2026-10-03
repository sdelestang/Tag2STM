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
#' @param n_pmoult1 Number of smallest bins with Pmoult fixed at exactly 1
#'   (\code{datain$n_pmoult1}, default 1).
#' @param P5 Fixed swap width (\code{datain$Growth_P5_fixed}, default 0.1).
#'
#' @return nlbin x nlbin column-stochastic matrix (columns = from-length).
#' @export
STMfromPars <- function(pars, bins, floor = 0, n_pmoult1 = 1, P5 = 0.1) {
  pars  <- as.numeric(pars)
  lbin  <- bins$lbin; lbinL <- bins$lbinL; lbinU <- bins$lbinU
  n     <- length(lbin)

  ## growth-at-length: double logistic, as growmodPar
  Amax  <- exp(pars[1])
  P2    <- pars[2]
  P1    <- exp(pars[3])
  P3    <- exp(pars[4])
  xdev  <- lbin - P2
  swap1 <- 1 / (1 + exp(xdev / P5))
  g     <- Amax * (1 / (1 + exp(xdev / P1)) * swap1 + 1 / (1 + exp(xdev / P3)) * (1 - swap1))
  sdg   <- exp(pars[5]) * g

  ## moult probability: floor + logistic, smallest bins fixed at 1
  pm <- floor + (1 - floor) * plogis(pars[6] - exp(pars[7]) * lbin)
  pm[seq_len(min(n_pmoult1, n))] <- 1

  ## columns: truncated at the from-bin, plus group at the top bin
  A <- matrix(0, n, n)
  for (fm in seq_len(n)) {
    k  <- fm:n
    mu <- lbin[fm] + g[fm]
    lo <- pnorm(lbinL[k], mu, sdg[fm])
    p  <- pnorm(lbinU[k], mu, sdg[fm]) - lo
    p[length(k)] <- 1 - lo[length(k)]
    p  <- p / sum(p)
    A[k, fm]  <- pm[fm] * p
    A[fm, fm] <- A[fm, fm] + (1 - pm[fm])
  }
  A
}
