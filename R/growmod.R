#' Tag-Recapture Growth Model with Size Transition Matrices
#'
#' Estimates flexible, data-driven growth patterns from tag-recapture data by
#' constructing size transition matrices (STMs) that describe how animals grow
#' from one length bin to another over time. Includes a moult-probability
#' hurdle with an optional biological minimum-moults-per-year floor, an
#' internal identification mixture (replacing the old external auxiliary
#' data machinery), and, optionally, year-specific proportional growth
#' scaling.
#'
#' @param pin List. Parameter list from \code{Makepin}.
#' @param Like Integer, 1 (default, asymmetric KL) or 2 (symmetric KL).
#'   Structural switch, plain data at trace time.
#' @param TemporalGrowth Logical. Structural switch, as before.
#'
#' @details
#' **CHANGES FROM THE PREVIOUS VERSION**
#'
#' \strong{1. Minimum-moults-per-year floor (\code{datain$mpy}).} Some
#' species have a known biological lower bound on annual moult frequency --
#' e.g. large western rock lobster must moult at least once a year
#' (\code{mpy = 1}) even though there are two \code{goodts} moulting windows
#' each capable of falling as low as 0.5; deep-sea crab may moult as rarely
#' as once every three years (\code{mpy = 1/3}) with a single annual
#' \code{goodts}; species with a terminal moult have no such floor
#' (\code{mpy = 0}, the default, reproduces the original unclamped
#' logistic exactly).
#'
#' The floor is enforced by construction, not by a penalty: each
#' \code{goodts} row's logistic is rescaled to asymptote at
#' \code{mpy_floor} (as size grows large) rather than at 0, so
#' \code{Pmoult(ns, fm) = mpy_floor + (1 - mpy_floor) * plogis(...)}, always
#' \code{>= mpy_floor}, with no comparison/branch on an AD-derived quantity
#' anywhere.
#'
#' How \code{mpy} is allocated across the year's \code{goodts} seasons is
#' an ESTIMATED parameter (\code{mpy_split_par}, present in \code{pin} only
#' when \code{length(goodts) > 1} -- see \code{Makepin}), not a preset
#' \code{datain} vector: a species may have one moulting window collapsing
#' toward 0 while another stays near 1 (rather than an even split), and
#' which window gets which end is generally not known in advance. The
#' split is built via a softmax over \code{mpy_split_par} with the last
#' \code{goodts} season's log-odds anchored at 0 for identifiability
#' (arbitrary anchor choice, unbiased since the parameter starts at 0 --
#' an exactly even split -- and is free to move either direction). With a
#' single annual \code{goodts} (e.g. deep-sea crab) there is nothing to
#' split and \code{mpy} goes entirely to that one season, same as before.
#'
#' A guard at the top of the function checks \code{mpy <= 1} (plain data,
#' safe to compare) and stops with an informative error otherwise, since a
#' single season's floor can now carry the entire \code{mpy} value (e.g.
#' fully allocated by the fitted split) and so cannot itself exceed 1.
#'
#' \strong{Smallest bins fixed at exactly 1 (\code{datain$n_pmoult1}).}
#' The smallest \code{n_pmoult1} length bins (default 1) have Pmoult
#' hard-fixed at exactly 1 -- juveniles moult with certainty by
#' definition, which a fitted logistic can only ever approach
#' asymptotically, never reach exactly. Set to 0 to disable (reverts to
#' the plain floor+logistic construction for every bin).
#'
#' \strong{2. Internal identification mixture (replaces \code{make_aux_data}/
#' \code{make_aux2_data} entirely).} The old two-step design (external
#' preprocessing building \code{aux_*}/\code{aux2_*} data blocks, restricted
#' to single- or Kmax-capped multi-opportunity records) is gone. Instead,
#' for every animal, the same \code{tstepsvec} already built for STM
#' chaining is reused directly to run a moment-matched forward recursion
#' over its \code{K} \code{goodts} opportunities (arbitrary \code{K}, no
#' cap), giving an exact Poisson-binomial mixture over moult COUNT with a
#' moment-matched Gaussian per count bucket (exact for animals with
#' \code{K <= 1}, an approximation for \code{K > 1} when growth means vary
#' across seasons within the record -- see \code{ident_mixture_ll} below).
#' Controlled by \code{datain$ident_wt} (default 0, reproduces a pure
#' marginal fit exactly -- fit at 0 first to confirm plumbing, then > 0).
#'
#' \code{datain$use_individual_error} (default \code{FALSE}) switches
#' between the population-level \code{sqrt(2)*sigError} measurement-error
#' term (as before) and a version using each animal's own fitted
#' \code{MerrorRel}/\code{MerrorRec}. The latter couples this
#' identification term to the same latent random effects the main
#' likelihood already fits -- compare both before trusting \code{TRUE}.
#'
#' All \code{aux_*} REPORT objects (\code{AuxLL}, \code{aux_pmoult},
#' \code{aux_post}) are gone; downstream code (\code{plotfit()} etc.) that
#' referenced them will need updating to use \code{IdentLL} instead. There
#' is no longer a "double counting" question or a \code{drop_from_main} /
#' seasonal-selection caveat to manage, since every animal contributes to
#' \code{IdentLL} through the same one mechanism regardless of liberty
#' length.
#'
#' @return Scalar total negative log-likelihood. Key \code{REPORT()}
#'   objects: \code{LL} (per-animal marginal log-likelihoods), \code{IdentLL}
#'   (per-animal identification log-densities, 0 where an animal had no
#'   \code{goodts} opportunity), \code{lenout}, \code{EstCapLen},
#'   \code{EstRecLen}, \code{stm}, \code{sigGrowvec}, \code{growthmat},
#'   \code{MerrorRel}, \code{MerrorRec}, \code{LsigGrow}, \code{Pmoult_vec},
#'   \code{Pmoult_par}, \code{mpy_floor} (the fitted floor per \code{goodts}
#'   row, for diagnostics), and \code{S} when \code{TemporalGrowth = TRUE}.
#'
#' @export
growmod <- function(pin, Like = 1, TemporalGrowth = FALSE) {

  ## --- Backwards compatibility -------------------------------------------
  ## Old aux_* fields are no longer consumed at all. mpy/ident_wt/
  ## use_individual_error/n_pmoult1 are new; all plain data, safe to
  ## default here.
  if (is.null(datain$mpy))                 datain$mpy <- 0
  if (is.null(datain$ident_wt))            datain$ident_wt <- 0
  if (is.null(datain$use_individual_error)) datain$use_individual_error <- FALSE
  if (is.null(datain$n_pmoult1))            datain$n_pmoult1 <- 1
  ## Tagging-induced moult suppression, as a recovery function of the
  ## animal's TOTAL time at liberty. Plain data / structural switch,
  ## decided at trace time. Default reproduces the previous behaviour
  ## exactly (no suppression, no extra parameters).
  if (is.null(datain$suppress))            datain$suppress <- FALSE
  ## Period mode: S is common within period rather than unique per year.
  ## Plain data / structural switch, decided at trace time.
  if (is.null(datain$period_mode))         datain$period_mode <- FALSE
  ## Weak prior on the Pmoult logistic, so the curve rests on a sensible
  ## value rather than drifting to an unidentified boundary when the data
  ## do not constrain it. Centred on a flat curve near 1 (slope 0). Set
  ## Pmoult_prior_sd to c(Inf, Inf) to switch it off entirely.
  if (is.null(datain$Pmoult_prior_mean))   datain$Pmoult_prior_mean <- c(qlogis(0.95), 0)
  if (is.null(datain$Pmoult_prior_sd))     datain$Pmoult_prior_sd   <- c(3, 0.1)
  ## PenSigError's centre and width. Previously hardcoded to
  ## log(2.0) / 0.5, which silently pulled toward 2 mm for any species
  ## whose measurement error was not 2 mm. Now supplied by add_sigError()
  ## via Makedata; the old values remain the fallback for a datain built
  ## before that existed.
  if (is.null(datain$LsigError_prior_mean)) datain$LsigError_prior_mean <- log(2.0)
  if (is.null(datain$LsigError_prior_sd))   datain$LsigError_prior_sd   <- 0.5

  getAll(datain, pin, warn = FALSE)
  npar <- length(names(pin))
  nobs <- length(Rccl)

  ## mpy is plain data; the feasibility bound is just mpy <= 1 now (a
  ## single season could in principle carry the whole floor -- see below),
  ## not mpy/length(goodts) <= 1 as in the old fixed-even-split version.
  n_goodts <- length(goodts)
  if (mpy > 1) {
    stop("datain$mpy = ", mpy, " exceeds 1 -- not a valid probability floor.")
  }

  ## How mpy is allocated across the year's goodts seasons is now an
  ## ESTIMATED parameter, not a preset datain vector -- because which
  ## season carries more of the floor (e.g. one lobster moulting window
  ## descending to 0 while the other stays near 1, vs. an even 0.5/0.5
  ## split) is exactly the kind of thing you don't know in advance and
  ## want the data to decide.
  ##
  ## mpy_split_par has length (n_goodts - 1) when n_goodts > 1 (present in
  ## pin only in that case -- see Makepin), and is transformed via a
  ## softmax with the LAST goodts season's log-odds anchored at 0 for
  ## identifiability (standard multinomial-logit-to-simplex construction,
  ## arbitrary which category is anchored -- doesn't bias the fit since
  ## the parameter starts at exactly 0, i.e. a perfectly even split, and is
  ## free to move either direction from there).
  ##
  ## When n_goodts <= 1 there's nothing to split -- mpy goes entirely to
  ## the single season, same as before.
  if (n_goodts > 1) {
    split_raw <- c(mpy_split_par, 0)     # length n_goodts, last entry anchored
    split_exp <- exp(split_raw)
    split     <- split_exp / sum(split_exp)   # sums to 1 by construction
    mpy_floor_vec <- mpy * split
  } else {
    mpy_floor_vec <- rep(mpy, n_goodts)  # length 0 or 1
  }
  ## Index by raw timestep number (1..ntsteps), not by position in goodts,
  ## so Pmoult_fn(ns, fm) can look it up directly with ns as given elsewhere.
  mpy_floor <- rep(0, ntsteps)
  mpy_floor[goodts] <- mpy_floor_vec

  if (TemporalGrowth) {
    if (any(relyr < 1) || any(relyr > nyears)) {
      stop("datain$relyr has values outside 1..", nyears, ". Convert to a ",
           "1-based year index before fitting -- see add_year_support().")
    }
    if (period_mode) {
      ## One growth scalar per PERIOD, shared by every year in it. Sraw has
      ## length nperiods - 1 and the last period is the sum-to-zero anchor,
      ## exactly as the annual version does across supported years.
      ##
      ## year_period covers ALL years including those with no records, so an
      ## animal whose liberty spans a gap still gets a defined S -- which is
      ## the main practical advantage over annual effects on a series with a
      ## hole in it, where gap years are pinned at 0 by construction.
      ##
      ## No yr_supported logic here: support is a property of the period,
      ## which Makedata checks, not of the individual year.
      Sper <- rep(0, nperiods)
      if (nperiods > 1) Sper <- c(Sraw, -sum(Sraw))
      S <- Sper[year_period]
      Spen <- Sper   # penalise once per period, not once per year
    } else {
      supported_idx <- which(yr_supported)
      n_sup <- length(supported_idx)
      S <- rep(0, nyears)
      if (n_sup > 1) {
        S[supported_idx] <- c(Sraw, -sum(Sraw))
      }
      Spen <- S
    }
    ADREPORT(S)
  }

  ## --- Initialize output structures --------------------------------------
  if (TemporalGrowth) {
    stm <- array(0, c(nlbin, nlbin, ntsteps, nyears))
  } else {
    stm <- array(0, c(nlbin, nlbin, ntsteps))
  }
  EstRecLen <- matrix(0, ncol = length(lbin), nrow = nobs)
  EstCapLen <- matrix(0, ncol = length(lbin), nrow = nobs)
  growthmat <- matrix(0, ncol = length(lbin), nrow = ntsteps)
  LL       <- rep(0, nobs)
  IdentLL  <- rep(0, nobs)

  sigError <- exp(LsigError)

  ## Pmoult_fn: the smallest n_pmoult1 length bins are hard-fixed at exactly
  ## 1 (juveniles moult with certainty, by definition -- not merely "very
  ## likely", which is all a fitted logistic asymptote can ever give you,
  ## since plogis() approaches but never reaches exactly 1 for a finite
  ## argument). fm and n_pmoult1 are both plain data (fm is a loop index,
  ## n_pmoult1 is a datain scalar), so this comparison is safe -- same class
  ## as the goodts membership checks elsewhere in this function. Returning
  ## a plain numeric 1 from the true branch is fine even though the false
  ## branch returns an AD value: RTMB/TMB arithmetic between AD types and
  ## plain doubles is completely ordinary (already relied on throughout,
  ## e.g. "1 - Pmoult").
  ##
  ## Beyond n_pmoult1, mpy_floor[ns] is plain data (computed above from data
  ## only), so the floor rescaling of plogis(...) is a fixed transform of an
  ## AD quantity, not a branch -- safe. mpy_floor[ns] = 0 (default, e.g.
  ## non-goodts rows or mpy = 0) reproduces the original unclamped logistic.
  Pmoult_fn <- function(ns, fm) {
    if (fm <= n_pmoult1) return(1)
    fl <- mpy_floor[ns]
    fl + (1 - fl) * plogis(Pmoult_par[ns, 1] + Pmoult_par[ns, 2] * lbin[fm])
  }

  ## --- Tagging-induced moult suppression ---------------------------------
  ## Handling and tag insertion reduce an animal's chance of moulting soon
  ## after release; this is documented in other crustaceans. In the
  ## deep-sea crab data it showed as ~2.4 mm growth among animals with one
  ## moult opportunity against ~5.5 mm per opportunity at two and three,
  ## with no detectable shortfall by two opportunities.
  ##
  ## Modelled as a recovery function of the animal's TOTAL time at
  ## liberty, one value per animal, multiplying Pmoult in every timestep
  ## that animal occupies:
  ##
  ##   r_i = r0 + (1 - r0) * plogis((liberty_i - lib50) / slope)
  ##
  ## r0 in (0, 1) is the short-term impact of release, lib50 the liberty
  ## (in DAYS) at which half the recovery has occurred, slope the width of
  ## the transition in days. r rises from r0 at release toward 1.
  ##
  ## \strong{What this is and is not.} r * Pmoult is pure loss -- a
  ## suppressed moult is gone, with no mechanism for it to occur later. So
  ## this does NOT model catch-up; it asserts that animals recaptured after
  ## a long liberty were never suppressed. Since r depends on liberty,
  ## which is an outcome (when the animal happened to be recaptured) rather
  ## than a covariate of the animal, this is a correction for "observed too
  ## soon" rather than a biological rate: an animal at liberty long enough
  ## may have fitted two moults in and the shortfall is no longer visible.
  ## Report r0 as "animals recaptured shortly after release show a fraction
  ## r0 of the expected moult probability", not as "handling suppresses
  ## moulting by (1 - r0)".
  ##
  ## This replaces an earlier opportunity-indexed pair (a multiplier on the
  ## first goodts opportunity plus a compensating boost in the second).
  ## That version could not depress a single-opportunity animal without
  ## also depressing the first opportunity of every longer-liberty animal,
  ## because at any ntsteps the first opportunity sits at the same elapsed
  ## time regardless of how long the animal ends up at liberty.
  ##
  ## Note r scales Pmoult uniformly across size bins, so it does NOT need
  ## its own STM array -- see the blend identity at the chaining step.
  ## r deliberately overrides both the n_pmoult1 hard-fix and the mpy
  ## floor: a freshly tagged animal is neither certain to moult nor bound
  ## by a floor that describes untagged animals.
  if (suppress) {
    r0_p    <- plogis(r0_par)
    lib50_p <- exp(lib50_par)
    slope_p <- exp(slope_par)
    ## liberty is plain data (days); r_vec is AD via the parameters.
    r_vec <- r0_p + (1 - r0_p) * plogis((liberty - lib50_p) / slope_p)
  } else {
    r_vec <- rep(1, nobs)
  }

  growth_vecmat <- matrix(growth_vecpar, ncol = nlbin, nrow = ntsteps, byrow = TRUE)

  for (ns in goodts) {
    growth_vec <- rep(0, nlbin)
    growth_vec[nlbin] <- log(1 + exp(growth_vecmat[ns, nlbin]))
    for (i in (nlbin - 1):1) {
      growth_vec[i] <- growth_vec[i + 1] + log(1 + exp(growth_vecmat[ns, i]))
    }
    growthmat[ns, ] <- growth_vec
  }

  ## --- Build STM(s) -------------------------------------------------------
  ## Factored into a helper taking the Pmoult function to use, so the normal
  ## / suppressed / compensated arrays are built by identical code and
  ## cannot drift apart. (Divergence between two hand-maintained copies of
  ## this construction is exactly what produced a flat growth trajectory in
  ## build_lenout_envelope() previously.)
  make_stm <- function(pm_fn) {
    if (TemporalGrowth) {
      A <- array(0, c(nlbin, nlbin, ntsteps, nyears))
      for (yr in 1:nyears) {
        for (ns in goodts) {
          for (fm in 1:nlbin) {
            mn_growth <- growthmat[ns, fm] * exp(S[yr])
            sd_growth <- exp(LsigGrow) * mn_growth
            Pmoult    <- pm_fn(ns, fm)

            probs <- rep(0, nlbin)
            for (k in fm:(nlbin - 1)) {
              probs[k] <- pnorm(lbinU[k], lbin[fm] + mn_growth, sd_growth) -
                pnorm(lbinL[k], lbin[fm] + mn_growth, sd_growth)
            }
            probs[nlbin] <- 1 - pnorm(lbinL[nlbin], lbin[fm] + mn_growth, sd_growth)
            probs_norm <- probs[fm:nlbin] / sum(probs[fm:nlbin])

            A[fm:nlbin, fm, ns, yr] <- Pmoult * probs_norm
            A[fm, fm, ns, yr] <- A[fm, fm, ns, yr] + (1 - Pmoult)
          }
        }
      }
    } else {
      A <- array(0, c(nlbin, nlbin, ntsteps))
      for (ns in goodts) {
        for (fm in 1:nlbin) {
          mn_growth <- growthmat[ns, fm]
          sd_growth <- exp(LsigGrow) * mn_growth
          Pmoult    <- pm_fn(ns, fm)

          probs <- rep(0, nlbin)
          for (k in fm:(nlbin - 1)) {
            probs[k] <- pnorm(lbinU[k], lbin[fm] + mn_growth, sd_growth) -
              pnorm(lbinL[k], lbin[fm] + mn_growth, sd_growth)
          }
          probs[nlbin] <- 1 - pnorm(lbinL[nlbin], lbin[fm] + mn_growth, sd_growth)
          probs_norm <- probs[fm:nlbin] / sum(probs[fm:nlbin])

          A[fm:nlbin, fm, ns] <- Pmoult * probs_norm
          A[fm, fm, ns] <- A[fm, fm, ns] + (1 - Pmoult)
        }
      }
    }
    A
  }

  ## stm is the ordinary (untagged) transition matrix and is what gets
  ## REPORTed -- downstream consumers (ClipSTM, Vb2STM, the stock
  ## assessment) want the biological STM, not the tagging-affected one.
  ##
  ## No suppressed variant is built. Because r scales Pmoult uniformly
  ## across every size bin, and each STM column is
  ## stm(P) = P * A + (1 - P) * I (with A holding the normalised growth
  ## columns), scaling P -> r*P gives exactly
  ##
  ##   stm(r*P) = r * stm(P) + (1 - r) * I
  ##
  ## so a suppressed step is an exact linear blend of the base STM and the
  ## identity, applied at the chaining step below. That holds column-wise
  ## for any P, including the n_pmoult1 hard-fixed columns. It costs no
  ## extra arrays and no extra memory, which is what makes a fine-timestep
  ## configuration (monthly, say) feasible.
  stm <- make_stm(Pmoult_fn)

  ## --- Internal identification mixture (local function) -------------------
  ## Moment-matched forward recursion over K goodts opportunities. No
  ## comparisons on AD-derived quantities: the only branches (on m, mi, K)
  ## are over plain-integer loop indices, never over a weight/probability
  ## value itself; division-by-zero is avoided with an additive eps instead
  ## of a guarded branch.
  ident_mixture_ll <- function(inc, fm, ns_seq, yr_seq, r_i) {
    K <- length(ns_seq)
    if (K == 0) return(0)   # plain-data branch: fine

    p_k  <- rep(0, K)
    mu_k <- rep(0, K)
    vr_k <- rep(0, K)
    for (k in 1:K) {
      ns_k <- ns_seq[k]
      ## r_i is this animal's recovery factor, constant across its
      ## timesteps (it is a function of total liberty), so it scales every
      ## opportunity's moult probability identically -- matching what the
      ## STM chain does.
      p_k[k] <- r_i * Pmoult_fn(ns_k, fm)
      mg <- growthmat[ns_k, fm]
      if (TemporalGrowth) mg <- mg * exp(S[yr_seq[k]])
      mu_k[k] <- mg
      vr_k[k] <- (exp(LsigGrow) * mg)^2
    }

    eps <- 1e-10
    W  <- rep(0, K + 1); W[1] <- 1
    Mn <- rep(0, K + 1)
    Vr <- rep(0, K + 1)

    for (k in 1:K) {
      newW  <- rep(0, K + 1)
      newMn <- rep(0, K + 1)
      newVr <- rep(0, K + 1)

      for (mi in 1:(K + 1)) {
        m <- mi - 1   # plain integer, safe to compare below

        wA  <- W[mi] * (1 - p_k[k])
        muA <- Mn[mi]
        vrA <- Vr[mi]

        if (m >= 1) {
          wB  <- W[mi - 1] * p_k[k]
          muB <- Mn[mi - 1] + mu_k[k]
          vrB <- Vr[mi - 1] + vr_k[k]
        } else {
          wB <- 0; muB <- 0; vrB <- 0
        }

        wtot <- wA + wB
        newW[mi]  <- wtot
        newMn[mi] <- (wA * muA + wB * muB) / (wtot + eps)
        ex2       <- (wA * (vrA + muA^2) + wB * (vrB + muB^2)) / (wtot + eps)
        newVr[mi] <- ex2 - newMn[mi]^2
      }
      W <- newW; Mn <- newMn; Vr <- newVr
    }

    dens <- 0
    for (mi in 1:(K + 1)) {
      sd_tot <- sqrt(Vr[mi] + 2 * sigError^2)
      dens <- dens + W[mi] * dnorm(inc, Mn[mi], sd_tot)
    }
    log(dens + 1e-300)
  }

  ## --- Main (marginal) likelihood + identification term, one animal at a time
  for (r in 1:nobs) {
    Relength <- MerrorRel[r] + Rlcl[r]
    lens <- pnorm(lbinU, Relength, sigError) - pnorm(lbinL, Relength, sigError)

    if (TemporalGrowth) {
      if (tsteps[r] > 0) {
        abs0 <- (relyr[r] - 1) * ntsteps + relts[r]
        abst <- abs0 + (0:(tsteps[r] - 1))
        tstepsvec <- ((abst - 1) %%  ntsteps) + 1
        yearvec   <- ((abst - 1) %/% ntsteps) + 1
        yearvec[yearvec > nyears] <- nyears
      } else {
        tstepsvec <- integer(0)
        yearvec   <- integer(0)
      }
    } else {
      tstepsvec <- if (tsteps[r] > 0) {
        c(relts[r]:ntsteps, rep(1:ntsteps, 30))[1:tsteps[r]]
      } else integer(0)
      yearvec <- rep(1L, length(tstepsvec))   # unused when TemporalGrowth = FALSE
    }

    ## --- Identification term for this animal ---
    good_idx <- which(tstepsvec %in% goodts)
    if (length(good_idx) > 0) {
      fm_r <- which.min(abs(lbin - Rlcl[r]))   # release bin, fixed-bin approx
      inc_r <- if (use_individual_error) {
        (Rccl[r] + MerrorRec[r]) - (Rlcl[r] + MerrorRel[r])
      } else {
        Rccl[r] - Rlcl[r]
      }
      IdentLL[r] <- ident_mixture_ll(
        inc = inc_r, fm = fm_r,
        ns_seq = tstepsvec[good_idx],
        yr_seq = if (TemporalGrowth) yearvec[good_idx] else NULL,
        r_i = r_vec[r]
      )
    }

    ## --- STM chaining -------------------------------------------------
    ## Suppression enters as the exact blend stm(r*P) = r*stm(P) + (1-r)*I
    ## derived above, so no suppressed STM array is needed. r_vec[r] = 1
    ## when suppress is FALSE, which reduces this to the plain product.
    if (length(tstepsvec) > 0) {
      for (ts in 1:length(tstepsvec)) {
        if (tstepsvec[ts] %in% goodts) {
          tmpstm <- if (TemporalGrowth) stm[, , tstepsvec[ts], yearvec[ts]]
                    else stm[, , tstepsvec[ts]]
          lens <- (1 - r_vec[r]) * lens + r_vec[r] * (tmpstm %*% lens)
        }
      }
    }

    EstRecLen[r, ] <- lens

    Reclength <- MerrorRec[r] + Rccl[r]
    EstCapLen[r, ] <- pnorm(lbinU, Reclength, sigError) - pnorm(lbinL, Reclength, sigError)

    eps <- 1e-8
    if (Like == 1) {
      LL[r] <- sum((EstCapLen[r, ] + eps) *
                     log((EstRecLen[r, ] + eps) / (EstCapLen[r, ] + eps))) * nlob[r]
    }
    if (Like == 2) {
      P <- EstRecLen[r, ] + eps
      Q <- EstCapLen[r, ] + eps
      LL[r] <- -sum((P - Q) * log(P / Q)) * nlob[r]
    }
  }

  TIdentLL <- sum(IdentLL) * ident_wt

  ## --- Average growth trajectory (unchanged) -------------------------------
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
            mn_growth <- growthmat[ts, fm]
            sd_growth <- exp(LsigGrow) * mn_growth
            Pmoult    <- Pmoult_fn(ts, fm)
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

  gseq <- seq(0, 5, 0.1)
  sigGrowvec <- exp(LsigGrow) * gseq

  ## --- Penalties ------------------------------------------------------------
  PenSigError    <- -dnorm(LsigError, LsigError_prior_mean,
                           LsigError_prior_sd, log = TRUE)
  PenMerrorRel   <- -sum(dnorm(0, MerrorRel, exp(LMerrorRelsigma), log = TRUE))
  PenMerrorRec   <- -sum(dnorm(0, MerrorRec, exp(LMerrorRecsigma), log = TRUE))
  smooth_penalty <- smoother * sum((growthmat[goodts, 2:nlbin] -
                                      growthmat[goodts, 1:(nlbin - 1)])^2)
  drift_penalty  <- 0.1 * sum(log(1 + exp(-(growth_vecpar + 10))))

  ## --- Weak prior on the Pmoult logistic ----------------------------------
  ## Pmoult = 1 sits at the EDGE of the parameter space in logit terms
  ## (plogis(a) -> 1 only as a -> Inf), so whenever the data want a flat
  ## curve at 1 the likelihood is flat in (a, b) and the optimiser wanders
  ## to arbitrary values -- lobster fitted a = 0.77, b = 1.04, which puts
  ## the logistic argument above 50 at the smallest bin. Every large
  ## positive combination fits identically, so the parameters carry no
  ## information and their standard errors are meaningless.
  ##
  ## This prior gives the curve somewhere to rest instead: intercept
  ## centred on qlogis(0.95) and slope on 0, i.e. flat and near 1. It is
  ## deliberately weak (default sd 3 on the intercept, 0.1 per mm on the
  ## slope, roughly 10 logit units across a 100 mm size range), so data
  ## that genuinely want a declining curve override it easily -- the
  ## deep-sea crab slope of -0.166 is 1.7 sd, costing ~1.4 against the
  ## several hundred points that fit prefers. It is an ANCHOR for the
  ## unidentified case, not a constraint on the identified one.
  ##
  ## Summed over goodts rows only: other rows are never evaluated by
  ## Pmoult_fn and are fixed by Makemap, so penalising them would add a
  ## constant and nothing else.
  ##
  ## Note the prior is on the logistic BEFORE the mpy floor is applied, so
  ## with mpy > 0 the resting Pmoult is above 0.95 accordingly.
  PenPmoult <- 0
  for (ns in goodts) {
    PenPmoult <- PenPmoult -
      dnorm(Pmoult_par[ns, 1], Pmoult_prior_mean[1], Pmoult_prior_sd[1], log = TRUE) -
      dnorm(Pmoult_par[ns, 2], Pmoult_prior_mean[2], Pmoult_prior_sd[2], log = TRUE)
  }

  TLL <- -sum(LL) - TIdentLL + PenSigError + PenMerrorRel + PenMerrorRec +
    smooth_penalty + drift_penalty + PenPmoult

  if (TemporalGrowth) {
    ## Spen is the period vector in period mode and the year vector
    ## otherwise. Using S directly in period mode would penalise a long
    ## period more than a short one purely for containing more years.
    PenS <- 1.0 * sum(Spen^2)
    TLL <- TLL + PenS
  }

  ## --- Report objects ---------------------------------------------------
  REPORT(LL)
  REPORT(IdentLL)
  REPORT(lenout)
  REPORT(EstCapLen)
  REPORT(EstRecLen)
  REPORT(stm)
  REPORT(sigGrowvec)
  REPORT(growthmat)
  REPORT(MerrorRel)
  REPORT(MerrorRec)
  REPORT(LsigGrow)
  Pmoult_vec <- matrix(0, ntsteps, nlbin)
  for (ns in 1:ntsteps) {
    for (fm in 1:nlbin) Pmoult_vec[ns, fm] <- Pmoult_fn(ns, fm)
  }
  REPORT(Pmoult_vec)
  REPORT(Pmoult_par)
  REPORT(PenPmoult)
  REPORT(mpy_floor)
  if (suppress) {
    ## r0_p    : recovery multiplier at zero liberty (1 = no effect)
    ## lib50_p : liberty in DAYS at half recovery
    ## slope_p : transition width in days
    ## r_vec   : per-animal multiplier actually applied
    ## r_curve : r evaluated on a day grid, for plotting
    REPORT(r0_p)
    REPORT(lib50_p)
    REPORT(slope_p)
    REPORT(r_vec)
    rgrid   <- seq(0, 1500, 10)
    r_curve <- r0_p + (1 - r0_p) * plogis((rgrid - lib50_p) / slope_p)
    REPORT(rgrid)
    REPORT(r_curve)
  }
  if (TemporalGrowth) {
    REPORT(S)
    if (period_mode) {
      ## Sper is the estimated quantity; S is it broadcast over years.
      REPORT(Sper)
    }
  }

  TLL
}
