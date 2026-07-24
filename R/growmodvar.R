#' Tag-Recapture Growth Model with Year-Specific Proportional Growth Scaling
#'
#' Extension of \code{growmod} that adds interannual variation in growth rate.
#' Each year t gets a multiplicative scalar exp(S_t) applied to the mean growth
#' increment (and, since sd_growth is a function of mean growth, to the growth
#' variance as well) at every length bin and season. S_t is on the log scale
#' (S_t = 0 => average growth, i.e. exp(0) = 1x; S_t = log(1.1) ~= 0.0953 =>
#' 10% faster growth that year; S_t = log(0.9) ~= -0.105 => 10% slower). This
#' is deliberately NOT an additive (1 + S_t) scaling -- see @details for why.
#'
#' S is constrained to sum to zero (on the log scale, among supported years)
#' so the "average STM" (i.e. growthmat evaluated with S = 0, exp(0) = 1x)
#' remains a meaningful, estimable reference trajectory rather than an
#' arbitrary reference year. Concretely this constrains the *geometric* mean
#' of the supported years' growth multipliers exp(S_t) to 1, which is the
#' standard convention for multiplicative year effects.
#'
#' @param pin List. As for \code{\link{growmod}}, plus:
#'   \itemize{
#'     \item \code{Sraw}: vector of length \code{nyears - 1}. Free year-effect
#'       parameters; the final year's effect is derived as
#'       \code{-sum(Sraw)} so that all \code{nyears} effects sum to zero.
#'   }
#' @param Like Integer. As for \code{\link{growmod}} (1 = asymmetric KL,
#'   2 = symmetric KL).
#'
#' @details
#' Requires \code{datain} to additionally supply:
#' \itemize{
#'   \item \code{nyears}: number of years modelled (integer)
#'   \item \code{relyr}: release year (1..nyears) for each tagged animal,
#'     parallel to \code{relts}
#'   \item \code{yr_supported}: logical vector, length \code{nyears}, TRUE for
#'     years with adequate recapture support — produced by
#'     \code{\link{add_year_support}}, which must be run (and \code{Sraw} sized
#'     accordingly via \code{Makepin(TempGrowth = TRUE)}) before fitting.
#' }
#'
#' Animals at liberty across a year boundary are handled automatically:
#' the projection loop advances one time step at a time, indexing into
#' \code{stm[, , season, year]}, so an animal recaptured mid-year only
#' passes through that year's scaled STM as many times as it actually
#' moulted before recapture — no separate proration step is needed.
#'
#' \strong{Important RTMB constraint — why the min/max scenario envelope is
#' NOT built here:} \code{S} is derived from the estimated parameter
#' \code{Sraw}, so within this function it is an AD (automatic-differentiation)
#' type, not a plain number. RTMB refuses \code{min()}/\code{max()} (and any
#' other comparison, e.g. \code{>}) on AD types, because differentiating
#' through a branch/comparison is not well-defined — it raises
#' "Comparison is generally unsafe for AD types". This means the worst-year/
#' best-year scenario envelope (previously built here via
#' \code{c(0, min(S), max(S))}) cannot live inside this function at all. It is
#' instead built as a separate post-fit step, in plain R, from
#' \code{mod$rep()}'s output — see \code{\link{build_lenout_envelope}} — where
#' \code{S} is already a plain numeric vector and \code{min()}/\code{max()}
#' are unproblematic. \code{lenout} here is therefore only ever the literal
#' \code{S = 0} average-year trajectory (a fixed constant, not a comparison,
#' so it's safe to compute inside the AD-traced function).
#'
#' \strong{Design choices worth revisiting once this is fitted:}
#' \itemize{
#'   \item \code{PenS} is a ridge penalty (weight 1.0) shrinking the
#'     \emph{supported}-year values of \code{S} toward zero (it also sums over
#'     the hard-fixed unsupported years, but those contribute exactly 0
#'     regardless, so this only actually regularizes supported years). The
#'     weight is arbitrary and should be tuned or profiled, not left at 1.0
#'     by default. Unsupported years no longer risk drifting to extreme
#'     values — they are fixed at \code{S = 0} directly, not merely penalized
#'     — since a year with no data participating in a sum-to-zero constraint
#'     could previously be forced to an arbitrary value by construction,
#'     regardless of how strong the ridge penalty was.
#'   \item Animals whose liberty period runs past \code{nyears} (i.e. beyond
#'     the last modelled year) have their \code{yearvec} clamped to
#'     \code{nyears} — they keep experiencing the last year's S_t rather than
#'     erroring out. Check this is the behaviour you want; the alternative is
#'     to exclude such animals from \code{datain} entirely.
#'   \item \code{mn_growth} and \code{sd_growth} both scale with \code{S_t}
#'     (coupled, per your call) — a fast-growth year is also a
#'     higher-growth-variance year. If post-fit diagnostics show \code{S_t}
#'     doing double duty (explaining both an unusual mean AND an unusual
#'     spread in a given year's recapture lengths), that's this coupling at
#'     work, not necessarily a data problem.
#' }
#'
#' @return Scalar total negative log-likelihood, with the same \code{REPORT()}
#'   objects as \code{growmod} plus:
#' \describe{
#'   \item{S}{Vector of length \code{nyears}, the fitted year effects (sums to 0)}
#'   \item{stm}{Now a 4D array: \code{nlbin x nlbin x ntsteps x nyears}}
#'   \item{lenout}{2D, \code{(30*ntsteps) x nlbin} — the average-year (S=0)
#'     trajectory only, same shape as \code{growmod}'s original \code{lenout},
#'     so existing diagnostic code (e.g. \code{plotfit}) works unmodified.}
#'   \item{LsigGrow}{The fitted \code{LsigGrow} parameter (length 2), reported
#'     so that \code{\link{build_lenout_envelope}} can recompute sd_growth
#'     from mn_growth outside the AD-traced function.}
#' }
#'
#' @seealso \code{\link{build_lenout_envelope}} for the post-fit min/max
#'   scenario envelope this function deliberately does not compute.
#'
#' @export
growmodvar <- function(pin, Like = 1) {
  getAll(datain, pin, warn = FALSE)
  npar <- length(names(pin))
  nobs <- length(Rccl)

  # relyr must be a 1-based index into 1..nyears, not a raw calendar year --
  # see add_year_support() for the full explanation. relyr is plain data
  # (not AD), so this comparison is safe here regardless of the AD-comparison
  # restriction that applies to parameter-derived quantities elsewhere in
  # this function.
  if (any(relyr < 1) || any(relyr > nyears)) {
    stop("datain$relyr has values outside 1..", nyears, ". Convert to a ",
         "1-based year index before fitting -- see add_year_support().")
  }

  # --- Year-effect vector ---
  # Sraw only covers years with adequate recapture support (see
  # add_year_support()/Makepin's TempGrowth sizing). Unsupported years are
  # fixed at exactly S = 0 (average growth) rather than being folded into the
  # sum-to-zero constraint, where an unsupported year would otherwise collapse
  # onto whatever value satisfies the constraint given the supported years --
  # potentially large, and previously landing on whichever position TMB
  # derives automatically (always the last chronological year), regardless of
  # whether that specific year happened to be the unsupported one.
  supported_idx <- which(yr_supported)   # plain data (logical from datain), not AD
  n_sup <- length(supported_idx)
  S <- rep(0, nyears)
  if (n_sup > 1) {
    S[supported_idx] <- c(Sraw, -sum(Sraw))  # sums to 0 among supported years only
  }
  ADREPORT(S)

  # --- Initialize output structures ---
  stm <- array(0, c(nlbin, nlbin, ntsteps, nyears))
  EstRecLen <- matrix(0, ncol = length(lbin), nrow = nobs)
  EstCapLen <- matrix(0, ncol = length(lbin), nrow = nobs)
  growthmat <- matrix(0, ncol = length(lbin), nrow = ntsteps)
  LL <- rep(0, nobs)

  # Transform sigma parameters
  sigError <- exp(LsigError)
  LsigGrow_base <- LsigGrow[1]
  Lsig_rate <- LsigGrow[2]

  # Reshape growth parameters into matrix
  growth_vecmat <- matrix(growth_vecpar, ncol = nlbin, nrow = ntsteps, byrow = TRUE)

  ## Average growth curve (identical to growmod) — this is the S=0 reference
  for (ns in goodts) {
    growth_vec <- rep(0, nlbin)
    growth_vec[nlbin] <- log(1 + exp(growth_vecmat[ns, nlbin]))
    for (i in (nlbin - 1):1) {
      growth_vec[i] <- growth_vec[i + 1] + log(1 + exp(growth_vecmat[ns, i]))
    }
    growthmat[ns, ] <- growth_vec
  }

  ## Build year- AND season-specific STMs: mean (and hence sd) scaled by exp(S[yr])
  ## (multiplicative on the log scale, NOT (1 + S[yr]) -- the additive form has
  ## a hard wall at S = -1 where mn_growth hits exactly 0, sd_growth goes
  ## negative just past it, and pnorm() returns NaN. A profile likelihood run
  ## against that parameterization showed the objective still improving
  ## monotonically all the way to S = -1 with no turnover -- i.e. the
  ## optimizer was running into that undefined region, not converging on a
  ## genuine minimum. exp(S) is strictly positive for any finite S, so
  ## "very slow growth" is expressed as S -> -infinity smoothly, with no wall
  ## to hit. S = 0 still means "average year" exactly as before; the
  ## sum-to-zero constraint on Sraw now constrains supported years' geometric
  ## mean multiplier to 1, rather than their arithmetic mean deviation to 0.
  for (yr in 1:nyears) {
    for (ns in goodts) {
      for (fm in 1:nlbin) {
        mn_growth <- growthmat[ns, fm] * exp(S[yr])
        #sd_growth <- exp(LsigGrow_base) * (1 - exp(-exp(Lsig_rate) * mn_growth))
        sd_growth <- exp(LsigGrow_base) * mn_growth

        probs <- rep(0, nlbin)
        for (k in fm:(nlbin - 1)) {
          probs[k] <- pnorm(lbinU[k], lbin[fm] + mn_growth, sd_growth) -
            pnorm(lbinL[k], lbin[fm] + mn_growth, sd_growth)
        }
        probs[nlbin] <- 1 - pnorm(lbinL[nlbin], lbin[fm] + mn_growth, sd_growth)
        stm[fm:nlbin, fm, ns, yr] <- probs[fm:nlbin] / sum(probs[fm:nlbin])
      }
    }
  }

  ## Calculate likelihood for each tagged animal
  for (r in 1:nobs) {
    # Release length distribution with measurement error
    Relength <- MerrorRel[r] + Rlcl[r]
    lens <- pnorm(lbinU, Relength, sigError) - pnorm(lbinL, Relength, sigError)

    # Absolute time indexing -> parallel season (tstepsvec) and year (yearvec) vectors
    abs0 <- (relyr[r] - 1) * ntsteps + relts[r]
    abst <- abs0 + (0:(tsteps[r] - 1))
    tstepsvec <- ((abst - 1) %% ntsteps) + 1
    yearvec   <- ((abst - 1) %/% ntsteps) + 1
    yearvec[yearvec > nyears] <- nyears   # TODO: decide if clamping is right for your data,
    # vs. excluding animals that outlive the modelled years

    if (length(tstepsvec) > 0) {
      for (ts in 1:length(tstepsvec)) {
        if (tstepsvec[ts] %in% goodts) {
          tmpstm <- stm[, , tstepsvec[ts], yearvec[ts]]
          lens <- tmpstm %*% lens
        }
      }
    }

    EstRecLen[r, ] <- lens

    # Recapture length distribution with measurement error
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

  ## Average-year (S = 0) growth trajectory — literal constant, no comparison
  ## on an AD type involved, so safe to compute here (unlike the min/max
  ## envelope — see @details and build_lenout_envelope()).
  lenout <- matrix(0, 30 * ntsteps, length(lbin))
  lenout[1, 1] <- 1
  cnt <- 0
  for (y in 1:(nrow(lenout) / ntsteps)) {
    for (ts in 1:ntsteps) {
      cnt <- cnt + 1
      if (cnt > 1) lenout[cnt, ] <- lenout[cnt - 1, ]
      if (ts %in% goodts) {
        avg_stm <- matrix(0, nlbin, nlbin)
        for (fm in 1:nlbin) {
          mn_growth <- growthmat[ts, fm]              # Ssc = 0, i.e. * (1 + 0)
          sd_growth <- exp(LsigGrow_base) * (1 - exp(-exp(Lsig_rate) * mn_growth))
          probs <- rep(0, nlbin)
          for (k in fm:(nlbin - 1)) {
            probs[k] <- pnorm(lbinU[k], lbin[fm] + mn_growth, sd_growth) -
              pnorm(lbinL[k], lbin[fm] + mn_growth, sd_growth)
          }
          probs[nlbin] <- 1 - pnorm(lbinL[nlbin], lbin[fm] + mn_growth, sd_growth)
          avg_stm[fm:nlbin, fm] <- probs[fm:nlbin] / sum(probs[fm:nlbin])
        }
        lenout[cnt, ] <- avg_stm %*% lenout[cnt, ]
      }
    }
  }

  gseq <- seq(0, 5, 0.1)
  sigGrowvec <- exp(LsigGrow_base) * (1 - exp(-exp(Lsig_rate) * gseq))

  # Calculate penalties
  sig_at_5 <- exp(LsigGrow[1]) * (1 - exp(-exp(LsigGrow[2]) * 5))
  PensigGrowsd <- 1.0 * log(1 + exp(sig_at_5 - 5))
  PenSigError <- -dnorm(LsigError, log(2.0), 0.5, log = TRUE)
  PenMerrorRel <- -sum(dnorm(0, MerrorRel, exp(LMerrorRelsigma), log = TRUE))
  PenMerrorRec <- -sum(dnorm(0, MerrorRec, exp(LMerrorRecsigma), log = TRUE))
  smooth_penalty <- smoother * sum((growthmat[goodts, 2:nlbin] - growthmat[goodts, 1:(nlbin - 1)])^2)
  drift_penalty <- 0.1 * sum(log(1 + exp(-(growth_vecpar + 10))))

  # Ridge penalty on year effects — keeps Sraw identifiable/stable when some
  # years have thin recapture data. Weight (1.0) is a starting point, not tuned.
  PenS <- 1.0 * sum(S^2)

  # Total negative log-likelihood
  TLL <- -sum(LL) + PenSigError + PenMerrorRel + PenMerrorRec +
    smooth_penalty + drift_penalty + PensigGrowsd + PenS

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
  REPORT(S)
  REPORT(LsigGrow)

  TLL
}
