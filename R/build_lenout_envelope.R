#' Build Worst/Average/Best-Year Growth Trajectory Envelope (Post-Fit)
#'
#' Companion function to \code{\link{growmod}}. Reconstructs the
#' \code{scenario_S <- c(0, min(S), max(S))} envelope that cannot be built
#' inside RTMB's AD-traced objective (\code{min()}/\code{max()} on \code{S},
#' derived from the estimated \code{Sraw}, are comparisons on an AD type).
#' This runs entirely on the plain numeric output of \code{mod$rep()} taken
#' after fitting, where those operations are unproblematic.
#'
#' @param rep_out List. Result of \code{mod$rep()} (or
#'   \code{mod$report(mout$par)}) for a fitted
#'   \code{growmod(..., TemporalGrowth = TRUE)} model. Must contain
#'   \code{S}, \code{growthmat}, \code{LsigGrow}, \code{Pmoult_par}, and
#'   \code{mpy_floor}.
#' @param datain List. The same \code{datain} used to fit -- needs
#'   \code{nlbin}, \code{ntsteps}, \code{goodts}, \code{lbin}, \code{lbinL},
#'   \code{lbinU}, and (optionally) \code{n_pmoult1}.
#'
#' @return A 3D array, \code{(30*ntsteps) x nlbin x 3}. Third dimension is
#'   (1) S=0 average year, (2) worst fitted year (min S), (3) best fitted
#'   year (max S). Returns \code{NULL} if \code{rep_out$S} is absent (plain
#'   \code{growmod} fit), so callers such as \code{plotfit} degrade to a
#'   single-trajectory plot.
#'
#' @section Keeping Pmoult in sync with growmod:
#' This function reconstructs STMs itself rather than reusing
#' \code{growmod}'s, so its Pmoult construction MUST mirror
#' \code{growmod}'s \code{Pmoult_fn} exactly. The previous version did not:
#' it used a bare \code{plogis(Pmoult_par[1] + Pmoult_par[2] * lbin[fm])},
#' predating both the \code{mpy} floor and the \code{n_pmoult1} hard-fix,
#' and indexed \code{Pmoult_par} as a flat vector rather than the
#' \code{ntsteps x 2} matrix it now is. With a fitted logistic whose
#' midpoint falls below the smallest length bin -- which is exactly what
#' happens when the \code{mpy} floor is binding, since the logistic
#' collapses and the floor carries all the probability -- that stale
#' formula returned Pmoult ~ 0 at every bin, so nothing ever moulted and
#' the whole envelope came out flat at the first bin for all 30 years.
#' The flat-vector indexing also happened to work only at
#' \code{ntsteps = 1}; with more than one season it silently read the
#' wrong elements.
#'
#' If \code{Pmoult_fn} in \code{growmod} changes again, this block must
#' change with it.
#'
#' @seealso \code{\link{growmod}}
#'
#' @export
build_lenout_envelope <- function(rep_out, datain) {
  if (is.null(rep_out$S)) return(NULL)

  nlbin   <- datain$nlbin
  ntsteps <- datain$ntsteps
  goodts  <- datain$goodts
  lbin    <- datain$lbin
  lbinL   <- datain$lbinL
  lbinU   <- datain$lbinU

  # n_pmoult1: number of smallest bins with Pmoult hard-fixed at exactly 1.
  # Mirrors growmod's own default when absent from datain.
  n_pmoult1 <- if (is.null(datain$n_pmoult1)) 1 else datain$n_pmoult1

  growthmat  <- rep_out$growthmat
  LsigGrow   <- rep_out$LsigGrow          # scalar log-CV
  Pmoult_par <- rep_out$Pmoult_par        # ntsteps x 2 matrix
  mpy_floor  <- rep_out$mpy_floor         # length ntsteps, 0 outside goodts
  S          <- rep_out$S                 # plain numeric here -> min/max fine

  if (is.null(Pmoult_par)) {
    # Fit predating the moult hurdle entirely: Pmoult = 1 everywhere, which
    # is what such a fit implicitly assumed.
    Pmoult_par <- matrix(rep(c(1e6, 0), each = ntsteps), nrow = ntsteps, ncol = 2)
  } else if (!is.matrix(Pmoult_par)) {
    # Fit from the old length-2-vector parameterization: same curve for
    # every season.
    Pmoult_par <- matrix(rep(Pmoult_par, each = ntsteps), nrow = ntsteps, ncol = 2)
  }
  if (is.null(mpy_floor)) mpy_floor <- rep(0, ntsteps)   # pre-mpy fit

  ## Exact mirror of growmod's Pmoult_fn -- see @section above.
  Pmoult_fn <- function(ns, fm) {
    if (fm <= n_pmoult1) return(1)
    fl <- mpy_floor[ns]
    fl + (1 - fl) * plogis(Pmoult_par[ns, 1] + Pmoult_par[ns, 2] * lbin[fm])
  }

  scenario_S <- c(0, min(S), max(S))
  lenout_scenarios <- array(0, c(30 * ntsteps, nlbin, 3))

  for (sc in 1:3) {
    Ssc <- scenario_S[sc]
    lenout_scenarios[1, 1, sc] <- 1
    cnt <- 0
    for (y in 1:(dim(lenout_scenarios)[1] / ntsteps)) {
      for (ts in 1:ntsteps) {
        cnt <- cnt + 1
        if (cnt > 1) lenout_scenarios[cnt, , sc] <- lenout_scenarios[cnt - 1, , sc]
        if (ts %in% goodts) {
          scenario_stm <- matrix(0, nlbin, nlbin)
          for (fm in 1:nlbin) {
            mn_growth <- growthmat[ts, fm] * exp(Ssc)   # matches growmod's exp(S)
            sd_growth <- exp(LsigGrow) * mn_growth      # proportional CV spread
            sd_growth <- max(sd_growth, 1e-4)           # plain max() -- not AD here
            Pmoult    <- Pmoult_fn(ts, fm)
            probs <- rep(0, nlbin)
            for (k in fm:(nlbin - 1)) {
              probs[k] <- pnorm(lbinU[k], lbin[fm] + mn_growth, sd_growth) -
                pnorm(lbinL[k], lbin[fm] + mn_growth, sd_growth)
            }
            probs[nlbin] <- 1 - pnorm(lbinL[nlbin], lbin[fm] + mn_growth, sd_growth)
            probs_norm <- probs[fm:nlbin] / sum(probs[fm:nlbin])
            scenario_stm[fm:nlbin, fm] <- Pmoult * probs_norm
            scenario_stm[fm, fm] <- scenario_stm[fm, fm] + (1 - Pmoult)
          }
          lenout_scenarios[cnt, , sc] <- scenario_stm %*% lenout_scenarios[cnt, , sc]
        }
      }
    }
  }

  lenout_scenarios
}
