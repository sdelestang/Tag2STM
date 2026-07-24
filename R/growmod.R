#' Tag-Recapture Growth Model with Size Transition Matrices
#'
#' Estimates flexible, data-driven growth patterns from tag-recapture data by
#' constructing size transition matrices (STMs) that describe how animals grow
#' from one length bin to another over time. Includes a moult-probability
#' hurdle (see Details) and, optionally, year-specific proportional growth
#' scaling. The resulting STMs can be directly incorporated into length-based
#' stock assessment models.
#'
#' @param pin List. Parameter list created by \code{\link{Makepin}}, containing
#'   initial values for all model parameters including growth parameters,
#'   the moult-probability hurdle, measurement error terms, and (when
#'   \code{TemporalGrowth = TRUE}) the year-effect parameters.
#' @param Like Integer. Likelihood form to use for comparing projected and
#'   observed length distributions. \code{1} (default): asymmetric KL
#'   divergence; \code{2}: symmetric KL divergence. This is a structural
#'   switch passed via closure, not differentiated by RTMB.
#' @param TemporalGrowth Logical. If \code{TRUE}, fits year-specific
#'   proportional growth scaling (the model previously implemented as a
#'   separate function; if \code{FALSE} (default), fits
#'   the base model with a single average STM and no interannual variation.
#'   This is also a structural switch, decided once when \code{MakeADFun}
#'   traces the function -- like \code{Like}, it is plain data at trace time,
#'   not derived from any AD parameter, so branching on it is safe (unlike,
#'   e.g., \code{min()}/\code{max()} on an AD-derived quantity -- see the
#'   year-effect section below for why that distinction matters here).
#'   \code{pin} must be built with the matching
#'   \code{Makepin(TemporalGrowth = ...)} value (\code{Sraw} is only present
#'   in \code{pin} when \code{TRUE}); \code{make_growmod_obj} checks this
#'   automatically via \code{attr(pin, "TemporalGrowth")}.
#'
#' @details
#' **Model Overview:**
#'
#' This function implements a tag-recapture growth model that:
#' \enumerate{
#'   \item Constructs size transition matrices (STMs) from flexible growth parameters
#'   \item Projects each tagged animal forward through time using appropriate STMs
#'   \item Compares projected size distributions with observed recapture lengths
#'   \item Maximizes a multinomial-like likelihood weighted by release cohort sizes
#' }
#'
#' **Required Data (from \code{datain}):**
#'
#' \itemize{
#'   \item \code{Rlcl}, \code{Rccl}: release/recapture lengths (carapace length)
#'   \item \code{tsteps}: time at liberty (number of time steps) for each animal
#'   \item \code{relts}: release time step (within-season position) for each animal
#'   \item \code{nlob}: number of lobsters/crabs in each release cohort (for weighting)
#'   \item \code{lbin}, \code{lbinL}, \code{lbinU}, \code{nlbin}: length bin structure
#'   \item \code{ntsteps}, \code{goodts}: time step structure and which have adequate data
#'   \item \code{smoother}: smoothness penalty weight for growth parameters
#'   \item When \code{TemporalGrowth = TRUE}, additionally: \code{nyears},
#'     \code{relyr} (release year, 1-based index into 1..nyears -- NOT a raw
#'     calendar year, see \code{\link{add_year_support}}), and
#'     \code{yr_supported} (logical, from \code{\link{add_year_support}}).
#' }
#'
#' **Growth Model Structure:**
#'
#' Growth is modeled using a random walk approach: \code{growth_vecpar} are
#' reshaped into a matrix with \code{nlbin} columns and \code{ntsteps} rows;
#' for each \code{goodts} time step, parameters are converted to cumulative
#' growth starting from the largest length bin (minimal growth) and adding
#' increments for each smaller bin. A smoothness penalty encourages adjacent
#' length bins to have similar growth.
#'
#' **Growth Variance -- Proportional (CV-type), Not Saturating:**
#'
#' \code{sd_growth = exp(LsigGrow) * mn_growth} -- growth spread is directly
#' proportional to mean growth increment at each size/season, tapering to
#' near-zero as growth itself tapers toward zero (e.g. near maximum size).
#' \code{LsigGrow} is a single scalar (log coefficient of variation), not a
#' 2-parameter saturating-ceiling curve as in an earlier version of this
#' model -- stratifying raw single-moult recapture increments by release size
#' showed both mean AND sd shrinking together toward the largest size
#' classes, which a saturating-ceiling formula could not reproduce (it
#' predicts close to the same near-maximal sd regardless of mn_growth once
#' past a small threshold).
#'
#' **Moult-Probability Hurdle:**
#'
#' Raw single-timestep recapture increments showed ~65% essentially zero
#' growth -- not a tail of one distribution but a genuinely separate "did not
#' moult" outcome (consistent with episodic moulting rather than continuous
#' growth). Each STM column is therefore a mixture: with probability
#' \code{1 - Pmoult(fm)}, the animal stays exactly in bin \code{fm} (no
#' moult); with probability \code{Pmoult(fm)}, growth follows the
#' proportional-CV distribution above, CONDITIONAL on a moult having
#' occurred. \code{Pmoult(fm) = plogis(Pmoult_par[1] + Pmoult_par[2] *
#' lbin[fm])} -- a simple 2-parameter logistic in size, not a flexible
#' per-bin parameter like \code{growth_vecpar}, since the observed size trend
#' (moult fraction ~0.41 at small sizes to ~0.25 at large sizes) was smooth
#' and didn't obviously need bin-level flexibility. This hurdle applies
#' identically whether or not \code{TemporalGrowth} is used. \code{Pmoult} is
#' currently constant across years and seasons -- worth revisiting if
#' diagnostics suggest moult probability itself varies by year, not just
#' growth-given-moult.
#'
#' **Year-Specific Growth Scaling (\code{TemporalGrowth = TRUE} only):**
#'
#' Each year t gets a multiplicative scalar \code{exp(S_t)} applied to the
#' mean growth increment (and hence, through the proportional formula above,
#' to sd_growth too) at every length bin and season. \code{S_t} is on the log
#' scale (\code{S_t = 0} => average growth). This is deliberately NOT an
#' additive \code{(1 + S_t)} scaling: that form has a hard wall at
#' \code{S_t = -1} where \code{mn_growth} hits exactly 0 and \code{sd_growth}
#' goes negative just past it (\code{pnorm()} returns \code{NaN}); a profile
#' likelihood run against that parameterization showed the objective still
#' improving monotonically all the way to the wall with no turnover -- i.e.
#' the optimizer was running into an undefined region, not converging on a
#' genuine minimum. \code{exp(S_t)} is strictly positive for any finite
#' \code{S_t}, so "very slow growth" is expressed smoothly as
#' \code{S_t -> -infinity}, with no wall to hit.
#'
#' \code{S} is constrained to sum to zero on the log scale, among
#' \strong{supported} years only (see \code{\link{add_year_support}}):
#' unsupported years are fixed at exactly \code{S = 0} rather than
#' participating in the sum-to-zero constraint, where an unsupported year
#' could otherwise be forced to an arbitrary value by construction (this bit
#' the very first version of this model, in which every unsupported year
#' collapsed onto whatever the automatically-derived last position required,
#' regardless of whether that position happened to be the unsupported one).
#'
#' Because \code{S} is derived from the parameter \code{Sraw}, it is an AD
#' type while the model is being traced -- \code{min()}/\code{max()} (or any
#' comparison) on \code{S} inside this function would raise "Comparison is
#' generally unsafe for AD types". The worst-year/best-year scenario envelope
#' is therefore NOT built here; it is a separate post-fit step, in plain R,
#' from \code{mod$rep()}'s output where \code{S} is already numeric -- see
#' \code{\link{build_lenout_envelope}}.
#'
#' \strong{Design choices worth revisiting:}
#' \itemize{
#'   \item \code{PenS} (only when \code{TemporalGrowth = TRUE}) is a ridge
#'     penalty (weight 1.0) shrinking supported years' \code{S} toward zero
#'     (it also sums over hard-fixed unsupported years, but those contribute
#'     exactly 0 regardless). The weight is arbitrary and should be tuned or
#'     profiled, not left at 1.0 by default.
#'   \item Animals whose liberty period runs past \code{nyears} have their
#'     internal year index clamped to \code{nyears} rather than erroring --
#'     check this is the behaviour you want.
#' }
#'
#' @return Scalar total negative log-likelihood. \code{REPORT()} objects:
#' \describe{
#'   \item{LL}{Vector of individual log-likelihoods for each tagged animal}
#'   \item{lenout}{Matrix, \code{(30*ntsteps) x nlbin} -- average growth
#'     trajectory (S=0 when \code{TemporalGrowth = TRUE}), for visualizing
#'     growth trajectories}
#'   \item{EstCapLen}{Observed recapture length distributions (with
#'     measurement error) for each animal}
#'   \item{EstRecLen}{Model-projected length distributions for each animal}
#'   \item{stm}{Size transition matrices: \code{nlbin x nlbin x ntsteps} when
#'     \code{TemporalGrowth = FALSE}, or \code{nlbin x nlbin x ntsteps x
#'     nyears} when \code{TRUE}}
#'   \item{sigGrowvec}{Growth SDs as a function of growth increment, evaluated
#'     at \code{seq(0, 5, 0.1)} (used by \code{plotfit}'s ribbon)}
#'   \item{growthmat}{Matrix of mean growth by length bin and time step}
#'   \item{MerrorRel}, \item{MerrorRec}{Estimated individual measurement errors}
#'   \item{LsigGrow}{The fitted scalar log-CV parameter}
#'   \item{Pmoult_vec}{P(moult) at each size bin, for diagnostic plotting}
#'   \item{Pmoult_par}{The raw fitted logistic coefficients (intercept, slope),
#'     reported so \code{build_lenout_envelope} can use the actual fitted
#'     hurdle rather than falling back to an "always moult" default}
#'   \item{S}{Only present when \code{TemporalGrowth = TRUE}: length
#'     \code{nyears}, the fitted year effects (sums to 0 among supported years)}
#' }
#'
#' @seealso \code{\link{Makepin}}, \code{\link{Makemap}},
#'   \code{\link{make_growmod_obj}}, \code{\link{add_year_support}},
#'   \code{\link{build_lenout_envelope}}
#'
#' @export
growmod <- function(pin, Like = 1, TemporalGrowth = FALSE) {
  getAll(datain, pin, warn = FALSE)
  npar <- length(names(pin))
  nobs <- length(Rccl)

  if (TemporalGrowth) {
    # relyr is plain data (not AD), so this comparison is safe here
    if (any(relyr < 1) || any(relyr > nyears)) {
      stop("datain$relyr has values outside 1..", nyears, ". Convert to a ",
           "1-based year index before fitting -- see add_year_support().")
    }

    # --- Year-effect vector ---
    # Sraw only covers years with adequate recapture support. Unsupported
    # years are fixed at exactly S = 0 rather than folded into the
    # sum-to-zero constraint -- see @details.
    supported_idx <- which(yr_supported)   # plain data, not AD
    n_sup <- length(supported_idx)
    S <- rep(0, nyears)
    if (n_sup > 1) {
      S[supported_idx] <- c(Sraw, -sum(Sraw))
    }
    ADREPORT(S)
  }

  # --- Initialize output structures ---
  if (TemporalGrowth) {
    stm <- array(0, c(nlbin, nlbin, ntsteps, nyears))
  } else {
    stm <- array(0, c(nlbin, nlbin, ntsteps))
  }
  EstRecLen <- matrix(0, ncol = length(lbin), nrow = nobs)
  EstCapLen <- matrix(0, ncol = length(lbin), nrow = nobs)
  growthmat <- matrix(0, ncol = length(lbin), nrow = ntsteps)
  LL <- rep(0, nobs)

  sigError <- exp(LsigError)

  growth_vecmat <- matrix(growth_vecpar, ncol = nlbin, nrow = ntsteps, byrow = TRUE)

  ## Mean growth curve (shared regardless of TemporalGrowth)
  for (ns in goodts) {
    growth_vec <- rep(0, nlbin)
    growth_vec[nlbin] <- log(1 + exp(growth_vecmat[ns, nlbin]))
    for (i in (nlbin - 1):1) {
      growth_vec[i] <- growth_vec[i + 1] + log(1 + exp(growth_vecmat[ns, i]))
    }
    growthmat[ns, ] <- growth_vec
  }

  ## Build STM(s). Both branches share the proportional-CV sd_growth and
  ## moult-probability hurdle; they differ only in whether growth is scaled
  ## by a year effect and whether stm has a 4th (year) dimension.
  if (TemporalGrowth) {
    for (yr in 1:nyears) {
      for (ns in goodts) {
        for (fm in 1:nlbin) {
          mn_growth <- growthmat[ns, fm] * exp(S[yr])
          sd_growth <- exp(LsigGrow) * mn_growth
          Pmoult <- plogis(Pmoult_par[1] + Pmoult_par[2] * lbin[fm])

          probs <- rep(0, nlbin)
          for (k in fm:(nlbin - 1)) {
            probs[k] <- pnorm(lbinU[k], lbin[fm] + mn_growth, sd_growth) -
              pnorm(lbinL[k], lbin[fm] + mn_growth, sd_growth)
          }
          probs[nlbin] <- 1 - pnorm(lbinL[nlbin], lbin[fm] + mn_growth, sd_growth)
          probs_norm <- probs[fm:nlbin] / sum(probs[fm:nlbin])

          stm[fm:nlbin, fm, ns, yr] <- Pmoult * probs_norm
          stm[fm, fm, ns, yr] <- stm[fm, fm, ns, yr] + (1 - Pmoult)
        }
      }
    }
  } else {
    for (ns in goodts) {
      for (fm in 1:nlbin) {
        mn_growth <- growthmat[ns, fm]
        sd_growth <- exp(LsigGrow) * mn_growth
        Pmoult <- plogis(Pmoult_par[1] + Pmoult_par[2] * lbin[fm])

        probs <- rep(0, nlbin)
        for (k in fm:(nlbin - 1)) {
          probs[k] <- pnorm(lbinU[k], lbin[fm] + mn_growth, sd_growth) -
            pnorm(lbinL[k], lbin[fm] + mn_growth, sd_growth)
        }
        probs[nlbin] <- 1 - pnorm(lbinL[nlbin], lbin[fm] + mn_growth, sd_growth)
        probs_norm <- probs[fm:nlbin] / sum(probs[fm:nlbin])

        stm[fm:nlbin, fm, ns] <- Pmoult * probs_norm
        stm[fm, fm, ns] <- stm[fm, fm, ns] + (1 - Pmoult)
      }
    }
  }

  ## Calculate likelihood for each tagged animal
  for (r in 1:nobs) {
    Relength <- MerrorRel[r] + Rlcl[r]
    lens <- pnorm(lbinU, Relength, sigError) - pnorm(lbinL, Relength, sigError)

    if (TemporalGrowth) {
      # Absolute time indexing -> parallel season and year vectors.
      # Guard tsteps[r] > 0 explicitly: R's 0:(-1) evaluates to c(0, -1), NOT
      # an empty sequence, so without this guard a zero-liberty animal
      # (tsteps[r] == 0, a large fraction of this dataset) would silently get
      # a 2-element tstepsvec/yearvec computed from a bogus negative offset,
      # and the STM would be applied twice when it should be applied zero
      # times. This was a real bug, not just a theoretical edge case, given
      # how many recaptures here have tsteps == 0.
      if (tsteps[r] > 0) {
        abs0 <- (relyr[r] - 1) * ntsteps + relts[r]
        abst <- abs0 + (0:(tsteps[r] - 1))
        tstepsvec <- ((abst - 1) %% ntsteps) + 1
        yearvec   <- ((abst - 1) %/% ntsteps) + 1
        yearvec[yearvec > nyears] <- nyears   # clamp animals outliving the modelled years

        for (ts in 1:length(tstepsvec)) {
          if (tstepsvec[ts] %in% goodts) {
            tmpstm <- stm[, , tstepsvec[ts], yearvec[ts]]
            lens <- tmpstm %*% lens
          }
        }
      }
    } else {
      # Original cycling logic -- no relyr/nyears needed at all
      tstepsvec <- c(relts[r]:ntsteps, rep(1:ntsteps, 30))[1:tsteps[r]]
      if (length(tstepsvec) > 0) {
        for (ts in 1:length(tstepsvec)) {
          if (tstepsvec[ts] %in% goodts) {
            tmpstm <- stm[, , tstepsvec[ts]]
            lens <- tmpstm %*% lens
          }
        }
      }
    }

    EstRecLen[r, ] <- lens

    Reclength <- MerrorRec[r] + Rccl[r]
    EstCapLen[r, ] <- pnorm(lbinU, Reclength, sigError) - pnorm(lbinL, Reclength, sigError)

    eps <- 1e-8
    if (Like == 1) {
      LL[r] <- sum((EstCapLen[r, ] + eps) * log((EstRecLen[r, ] + eps) / (EstCapLen[r, ] + eps))) * nlob[r]
    }
    if (Like == 2) {
      P <- EstRecLen[r, ] + eps
      Q <- EstCapLen[r, ] + eps
      LL[r] <- -sum((P - Q) * log(P / Q)) * nlob[r]
    }
  }

  ## Average growth trajectory. TemporalGrowth = FALSE: reuse stm[,,ts]
  ## directly (already the only STM there is). TemporalGrowth = TRUE:
  ## literal Ssc=0 built fresh (a constant, not a comparison on an AD type,
  ## so safe here -- see @details), since stm[,,,yr] is year-specific and
  ## none of those slices IS the average.
  lenout <- matrix(0, 30 * ntsteps, length(lbin))
  lenout[1, 1] <- 1
  cnt <- 0
  for (y in 1:(nrow(lenout) / ntsteps)) {
    for (ts in 1:ntsteps) {
      cnt <- cnt + 1
      if (cnt > 1) lenout[cnt, ] <- lenout[cnt - 1, ]
      if (ts %in% goodts) {
        if (TemporalGrowth) {
          avg_stm <- matrix(0, nlbin, nlbin)
          for (fm in 1:nlbin) {
            mn_growth <- growthmat[ts, fm]              # Ssc = 0, i.e. * exp(0) = 1x
            sd_growth <- exp(LsigGrow) * mn_growth
            Pmoult <- plogis(Pmoult_par[1] + Pmoult_par[2] * lbin[fm])
            probs <- rep(0, nlbin)
            for (k in fm:(nlbin - 1)) {
              probs[k] <- pnorm(lbinU[k], lbin[fm] + mn_growth, sd_growth) -
                pnorm(lbinL[k], lbin[fm] + mn_growth, sd_growth)
            }
            probs[nlbin] <- 1 - pnorm(lbinL[nlbin], lbin[fm] + mn_growth, sd_growth)
            probs_norm <- probs[fm:nlbin] / sum(probs[fm:nlbin])
            avg_stm[fm:nlbin, fm] <- Pmoult * probs_norm
            avg_stm[fm, fm] <- avg_stm[fm, fm] + (1 - Pmoult)
          }
          lenout[cnt, ] <- avg_stm %*% lenout[cnt, ]
        } else {
          tmpstm <- stm[, , ts]
          lenout[cnt, ] <- tmpstm %*% lenout[cnt, ]
        }
      }
    }
  }

  # gseq/sigGrowvec is a lookup table plotfit() uses (via approx()) to shade
  # the growth-spread ribbon in panel d.
  gseq <- seq(0, 5, 0.1)
  sigGrowvec <- exp(LsigGrow) * gseq

  # Calculate penalties
  PenSigError <- -dnorm(LsigError, log(2.0), 0.5, log = TRUE)
  PenMerrorRel <- -sum(dnorm(0, MerrorRel, exp(LMerrorRelsigma), log = TRUE))
  PenMerrorRec <- -sum(dnorm(0, MerrorRec, exp(LMerrorRecsigma), log = TRUE))
  smooth_penalty <- smoother * sum((growthmat[goodts, 2:nlbin] - growthmat[goodts, 1:(nlbin - 1)])^2)
  drift_penalty <- 0.1 * sum(log(1 + exp(-(growth_vecpar + 10))))

  TLL <- -sum(LL) + PenSigError + PenMerrorRel + PenMerrorRec +
    smooth_penalty + drift_penalty

  if (TemporalGrowth) {
    # Ridge penalty on year effects -- see @details for the weight caveat.
    PenS <- 1.0 * sum(S^2)
    TLL <- TLL + PenS
  }

  # Report objects for post-fit examination
  REPORT(LL)
  REPORT(lenout)
  REPORT(EstCapLen)
  REPORT(EstRecLen)
  REPORT(stm)
  REPORT(sigGrowvec)
  REPORT(growthmat)
  REPORT(MerrorRel)
  REPORT(MerrorRec)
  REPORT(LsigGrow)
  Pmoult_vec <- plogis(Pmoult_par[1] + Pmoult_par[2] * lbin)
  REPORT(Pmoult_vec)
  REPORT(Pmoult_par)
  if (TemporalGrowth) REPORT(S)

  TLL
}
