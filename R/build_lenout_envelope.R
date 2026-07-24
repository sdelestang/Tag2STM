#' Build Worst/Average/Best-Year Growth Trajectory Envelope (Post-Fit)
#'
#' Companion function to \code{\link{growmod}}. Reconstructs the same
#' \code{scenario_S <- c(0, min(S), max(S))} envelope that an earlier version
#' of \code{growmod(..., TemporalGrowth = TRUE)} tried to build internally — but that approach fails
#' inside RTMB's AD-traced objective function, because \code{S} is derived
#' from the estimated parameter \code{Sraw} and is therefore an AD type while
#' the model is being fit; \code{min()}/\code{max()} require comparisons,
#' and RTMB raises "Comparison is generally unsafe for AD types" for any
#' comparison on an AD value. This function instead runs entirely on the
#' plain numeric output of \code{mod$rep()}, taken \emph{after} fitting is
#' complete — at that point \code{S}, \code{growthmat}, etc. are ordinary
#' doubles, and \code{min()}/\code{max()} are unproblematic.
#'
#' @param rep_out List. The result of \code{mod$rep()} for a fitted
#'   \code{growmod(..., TemporalGrowth = TRUE)} model. Must contain \code{S}, \code{growthmat}, and
#'   \code{LsigGrow} (all reported by the current version of
#'   \code{growmod(..., TemporalGrowth = TRUE)}).
#' @param datain List. The same \code{datain} object used to fit the model —
#'   needed for \code{nlbin}, \code{ntsteps}, \code{goodts}, \code{lbin},
#'   \code{lbinL}, \code{lbinU}.
#'
#' @return A 3D array, \code{(30*ntsteps) x nlbin x 3}, matching the shape
#'   the earlier internal \code{lenout_scenarios} object had. The third
#'   dimension is (1) S=0 average year, (2) worst fitted year (min S),
#'   (3) best fitted year (max S).
#'
#' @details
#' If \code{rep_out$S} is \code{NULL} (i.e. the fit was a plain
#' \code{growmod}, not \code{growmod(..., TemporalGrowth = TRUE)}), this function returns \code{NULL}
#' rather than erroring, so callers (e.g. \code{plotfit}) can check for that
#' and fall back to a single-trajectory plot.
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

  growthmat     <- rep_out$growthmat
  LsigGrow <- rep_out$LsigGrow   # single scalar log-CV, not a 2-vector
  Pmoult_par    <- rep_out$Pmoult_par
  if (is.null(Pmoult_par)) {
    # older fit predating the moult-probability hurdle -- fall back to
    # Pmoult = 1 everywhere so this still runs against a growmod fit predating the hurdle,
    # just without the hurdle mixture it wouldn't have applied either
    Pmoult_par <- c(1e6, 0)   # plogis(1e6 + 0*x) == 1 for any finite x
  }
  S             <- rep_out$S            # plain numeric here -> min/max are fine

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
            mn_growth <- growthmat[ts, fm] * exp(Ssc)   # matches growmod's exp(S) parameterization
            sd_growth <- exp(LsigGrow) * mn_growth   # proportional (CV-type) spread, CONDITIONAL on moult
            sd_growth <- max(sd_growth, 1e-4)   # plain max() -- fine, not AD here
            Pmoult <- plogis(Pmoult_par[1] + Pmoult_par[2] * lbin[fm])
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
