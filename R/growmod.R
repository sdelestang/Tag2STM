#' Tag-Recapture Growth Model with Size Transition Matrices
#'
#' Estimates flexible, data-driven growth patterns from tag-recapture data by
#' constructing size transition matrices (STMs) that describe how animals grow
#' from one length bin to another over time. The model accounts for measurement
#' error at both release and recapture, and can include individual-level variation
#' in measurement precision. The resulting STMs can be directly incorporated into
#' length-based stock assessment models.
#'
#' @param pin List. Parameter list created by \code{\link{Makepin}}, containing
#'   initial values for all model parameters including growth parameters,
#'   measurement error terms, and random effects.
#'
#' @return Scalar numeric value representing the total negative log-likelihood
#'   (TLL) to be minimized. Lower values indicate better model fit. The function
#'   also populates numerous objects via \code{REPORT()} for post-fit diagnostics
#'   (see Details).
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
#' The function uses \code{getAll(datain, pin)} to load data and parameters, requiring:
#' \itemize{
#'   \item \code{Rlcl}: Release lengths (carapace length) for each tagged animal
#'   \item \code{Rccl}: Recapture lengths for each tagged animal
#'   \item \code{tsteps}: Time at liberty (number of time steps) for each animal
#'   \item \code{relts}: Release time step for each animal
#'   \item \code{nlob}: Number of lobsters in each release cohort (for weighting)
#'   \item \code{lbin}: Vector of length bin midpoints
#'   \item \code{lbinL}: Lower bounds of length bins
#'   \item \code{lbinU}: Upper bounds of length bins
#'   \item \code{nlbin}: Number of length bins
#'   \item \code{ntsteps}: Total number of time steps (e.g., 12 for monthly)
#'   \item \code{goodts}: Vector of time steps with sufficient data for estimation
#'   \item \code{smoother}: Smoothness penalty weight for growth parameters
#' }
#'
#' **Growth Model Structure:**
#'
#' Growth is modeled using a random walk approach:
#' \enumerate{
#'   \item Growth parameters (\code{growth_vecpar}) are reshaped into a matrix
#'     with \code{nlbin} columns and \code{ntsteps} rows
#'   \item For each time step in \code{goodts}, parameters are converted to
#'     cumulative growth starting from the largest length bin (minimal growth)
#'     and adding increments for each smaller bin
#'   \item Growth can vary flexibly with both size and season/time
#'   \item A smoothness penalty encourages adjacent length bins to have similar growth
#' }
#'
#' **Size Transition Matrix (STM) Construction:**
#'
#' For each time step in \code{goodts}:
#' \itemize{
#'   \item For each "from" length bin \code{fm}, animals grow according to
#'     \code{growthmat[timestep, fm]}
#'   \item Growth variability is modeled with constant variance \code{sigGrowsd}
#'   \item Transition probabilities to "to" length bins are calculated using
#'     truncated normal distributions: \code{pnorm(lbinU, mean, sd) - pnorm(lbinL, mean, sd)}
#'   \item This creates an \code{nlbin × nlbin} matrix for each time step
#' }
#'
#' **Measurement Error Model:**
#'
#' The model accounts for two sources of measurement error:
#' \enumerate{
#'   \item **Baseline error** (\code{sigError}): Applied to all measurements,
#'     representing instrument precision and rounding
#'   \item **Individual random effects** (\code{MerrorRel}, \code{MerrorRec}):
#'     Individual-specific deviations at release and recapture, useful for
#'     capturing differences between taggers or measurement conditions
#' }
#'
#' At both release and recapture, the observed length is treated as uncertain,
#' and this uncertainty is propagated through the model by creating probability
#' distributions across length bins.
#'
#' **Likelihood Calculation:**
#'
#' For each tagged animal:
#' \enumerate{
#'   \item Create release length distribution including measurement error
#'   \item Project forward through time using STMs for each time step at liberty
#'   \item Create recapture length distribution including measurement error
#'   \item Calculate Kullback-Leibler divergence between projected and observed
#'     distributions
#'   \item Weight by \code{nlob} (number of animals in release cohort)
#' }
#'
#' The likelihood is:
#' \deqn{LL_r = \sum_{l} (p_{obs,l} + \epsilon) \times \log\frac{p_{proj,l} + \epsilon}{p_{obs,l} + \epsilon} \times n_{lob,r}}
#' where \code{epsilon = 1e-8} prevents numerical issues.
#'
#' **Time Step Cycling:**
#'
#' The model handles cyclical time steps (e.g., annual cycles):
#' \itemize{
#'   \item If an animal is at liberty longer than \code{ntsteps}, time steps
#'     cycle back to 1
#'   \item This allows modeling seasonal growth patterns with multi-year recaptures
#'   \item Growth in time steps not in \code{goodts} is held constant (no STM applied)
#' }
#'
#' **Penalties and Priors:**
#'
#' The model includes several penalty terms to regularize parameter estimates:
#' \itemize{
#'   \item \strong{Measurement error prior}: Normal prior on \code{LsigError}
#'     centered at log(2) with SD = 0.5, providing weak regularization
#'   \item \strong{Random effects}: Normal priors on individual measurement errors
#'     \code{MerrorRel} and \code{MerrorRec} centered at 0 with SDs controlled
#'     by \code{LMerrorRelsigma} and \code{LMerrorRecsigma}
#'   \item \strong{Smoothness penalty}: Penalizes differences between adjacent
#'     growth increments within the first time step, weighted by \code{smoother}.
#'     This prevents unrealistic zigzag growth patterns across length bins.
#' }
#'
#' **Reported Objects:**
#'
#' The following objects are made available via \code{REPORT()} for post-fit
#' examination using \code{sdreport()}:
#' \describe{
#'   \item{LL}{Vector of individual log-likelihoods for each tagged animal}
#'   \item{lenout}{Matrix showing growth progression over 30 years starting
#'     from the smallest length bin, useful for visualizing growth trajectories}
#'   \item{EstCapLen}{Matrix of observed recapture length distributions
#'     (including measurement error) for each animal}
#'   \item{EstRecLen}{Matrix of model-projected length distributions for each
#'     animal after growth}
#'   \item{stm}{3D array (\code{nlbin × nlbin × ntsteps}) of size transition
#'     matrices - the primary output for use in stock assessment models}
#'   \item{sigGrowvec}{Vector of growth SDs (for legacy compatibility)}
#'   \item{growthmat}{Matrix of mean growth by length bin and time step}
#'   \item{MerrorRel}{Estimated individual measurement errors at release}
#'   \item{MerrorRec}{Estimated individual measurement errors at recapture}
#' }
#'
#' @section Typical Workflow:
#' \preformatted{
#' # 1. Prepare data
#' datain <- list(
#'   Rlcl = release_lengths,
#'   Rccl = recapture_lengths,
#'   tsteps = time_at_liberty,
#'   relts = release_timesteps,
#'   nlob = cohort_sizes,
#'   lbin = length_bins,
#'   lbinL = bin_lower_bounds,
#'   lbinU = bin_upper_bounds,
#'   nlbin = length(length_bins),
#'   ntsteps = 12,  # monthly
#'   goodts = c(1, 2, 3, 10, 11, 12),  # moulting seasons
#'   smoother = 1.0
#' )
#'
#' # 2. Create initial parameters
#' pin <- Makepin(LsigError = log(2), LsigGrow = log(0.15))
#'
#' # 3. Create parameter mapping
#' map <- Makemap(pin, re = FALSE)
#' map$LsigError <- NULL  # Estimate measurement error
#'
#' # 4. Build and fit model
#' obj <- MakeADFun(growmod, pin, map = map)
#' opt <- nlminb(obj$par, obj$fn, obj$gr)
#'
#' # 5. Extract results
#' rep <- obj$report()
#' stm_array <- rep$stm  # Use in stock assessment
#' growth_curves <- rep$growthmat
#' }
#'
#' @section Model Assumptions:
#' \itemize{
#'   \item Growth variability is constant across length bins (not proportional)
#'   \item Measurement error is normally distributed with constant variance
#'   \item Growth is deterministic within a time step (individual variation
#'     captured by STM spread, not stochastic growth)
#'   \item Animals in the same release cohort are exchangeable
#'   \item Time steps represent biologically meaningful periods (e.g., months, seasons)
#' }
#'
#' @section Computational Notes:
#' \itemize{
#'   \item Uses automatic differentiation via RTMB/TMB for efficient optimization
#'   \item Individual random effects can be integrated out using Laplace approximation
#'     when \code{random = c("MerrorRel", "MerrorRec")} is specified in \code{MakeADFun()}
#'   \item Smoothness penalty improves numerical stability and prevents overfitting
#'   \item The \code{eps = 1e-8} term prevents log(0) in likelihood calculation
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic model fit
#' datain <- PrepareTagData(tagdata)  # Hypothetical data prep function
#' pin <- Makepin()
#' map <- Makemap(pin, re = FALSE)
#' obj <- MakeADFun(growmod, pin, map = map)
#' opt <- nlminb(obj$par, obj$fn, obj$gr)
#'
#' # Extract STMs for stock assessment
#' rep <- obj$report()
#' stm_monthly <- rep$stm
#'
#' # Example 2: Model with estimated measurement error and random effects
#' pin <- Makepin(LMerrorRelsigma = log(0.5), LMerrorRecsigma = log(0.5))
#' map <- Makemap(pin, re = TRUE)
#' map$LsigError <- NULL  # Estimate
#' map$LMerrorRelsigma <- NULL
#' map$LMerrorRecsigma <- NULL
#'
#' obj <- MakeADFun(
#'   growmod,
#'   pin,
#'   map = map,
#'   random = c("MerrorRel", "MerrorRec")
#' )
#' opt <- nlminb(obj$par, obj$fn, obj$gr)
#' sdr <- sdreport(obj)
#'
#' # Example 3: Visualize growth patterns
#' rep <- obj$report()
#'
#' # Plot growth by length and season
#' library(ggplot2)
#' growth_df <- as.data.frame(rep$growthmat)
#' names(growth_df) <- datain$lbin
#' growth_df$timestep <- 1:nrow(growth_df)
#' growth_long <- reshape2::melt(growth_df, id.vars = "timestep")
#'
#' ggplot(growth_long, aes(x = as.numeric(as.character(variable)),
#'                         y = value, color = factor(timestep))) +
#'   geom_line() +
#'   labs(x = "Carapace Length (mm)",
#'        y = "Growth Increment (mm)",
#'        color = "Time Step") +
#'   theme_minimal()
#'
#' # Example 4: Examine fit quality
#' plot_df <- data.frame(
#'   observed = datain$Rccl,
#'   projected = apply(rep$EstRecLen, 1, function(x) sum(x * datain$lbin)),
#'   residual = datain$Rccl - apply(rep$EstRecLen, 1, function(x) sum(x * datain$lbin))
#' )
#'
#' ggplot(plot_df, aes(x = observed, y = projected)) +
#'   geom_point(alpha = 0.3) +
#'   geom_abline(slope = 1, intercept = 0, color = "red") +
#'   labs(x = "Observed Recapture Length",
#'        y = "Model-Projected Length") +
#'   theme_minimal()
#' }
#'
#' @seealso
#' \code{\link{Makepin}} for creating initial parameter values
#' \code{\link{Makemap}} for controlling which parameters are estimated
#'
#'
#' @export
growmod <- function(pin) {
  getAll(datain, pin, warn = FALSE)
  npar <- length(names(pin))
  nobs <- length(Rccl)

  # Initialize output structures
  stm <- array(0, c(nlbin, nlbin, ntsteps))  # Size transition matrices
  EstRecLen <- matrix(0, ncol = length(lbin), nrow = nobs)
  EstCapLen <- matrix(0, ncol = length(lbin), nrow = nobs)
  moultinc <- matrix(0, ncol = length(lbin), nrow = ntsteps)
  growthmat <- matrix(0, ncol = length(lbin), nrow = ntsteps)
  estmnlen <- rep(0, nobs)
  LL <- rep(0, nobs)  # Log-likelihood for each observation

  # Transform sigma parameters
  sigError <- exp(LsigError)   # Measurement error SD
  sigGrowsd <- exp(LsigGrow)   # Growth variability SD

  # Reshape growth parameters into matrix
  growth_vecmat <- matrix(growth_vecpar, ncol = nlbin, nrow = ntsteps, byrow = TRUE)

  ## Build Size Transition Matrices for each time step
  for (ns in goodts) {
    # Convert log-increments to cumulative growth by length
    growth_vec <- rep(0, nlbin)
    growth_vec[nlbin] <- exp(growth_vecmat[ns, nlbin])  # Largest bin (near zero)
    for (i in (nlbin - 1):1) {
      growth_vec[i] <- growth_vec[i + 1] + exp(growth_vecmat[ns, i])
    }

    growthmat[ns, ] <- growth_vec

    # Build STM using normal distribution of growth
    for (fm in 1:nlbin) {
      growth <- growthmat[ns, fm]
      sd_growth <- sigGrowsd

      # Middle bins: normal transitions
      if (fm + 1 <= nlbin - 1) {
        stm[(fm+1):(nlbin-1), fm, ns] <- pnorm(lbinU[(fm+1):(nlbin-1)], lbin[fm] + growth, sd_growth) -
          pnorm(lbinL[(fm+1):(nlbin-1)], lbin[fm] + growth, sd_growth)
      }

      # Floor: current bin absorbs all probability of staying same size or shrinking
      stm[fm, fm, ns] <- pnorm(lbinU[fm], lbin[fm] + growth, sd_growth)

      # Ceiling: last bin absorbs everything above its lower bound
      stm[nlbin, fm, ns] <- 1 - pnorm(lbinL[nlbin], lbin[fm] + growth, sd_growth)
    }
  }

  ## Calculate likelihood for each tagged animal
  for (r in 1:nobs) {
    # Release length distribution with measurement error
    Relength <- MerrorRel[r] + Rlcl[r]
    lens <- pnorm(lbinU, Relength, sigError) - pnorm(lbinL, Relength, sigError)

    # Determine sequence of time steps at liberty (with cycling)
    tstepsvec <- c(relts[r]:ntsteps, rep(1:ntsteps, 30))[1:tsteps[r]]

    # Project forward through time using STMs
    if (length(tstepsvec) > 0) {
      for (ts in 1:length(tstepsvec)) {
        if (tstepsvec[ts] %in% goodts) {
          tmpstm <- stm[, , tstepsvec[ts]]
          lens <- tmpstm %*% lens  # Matrix multiplication to grow
        }
      }
    }

    EstRecLen[r, ] <- lens

    # Recapture length distribution with measurement error
    Reclength <- MerrorRec[r] + Rccl[r]
    EstCapLen[r, ] <- pnorm(lbinU, Reclength, sigError) - pnorm(lbinL, Reclength, sigError)

    # Calculate KL-divergence weighted by cohort size
    eps <- 1e-8
    LL[r] <- sum((EstCapLen[r, ] + eps) * log((EstRecLen[r, ] + eps) / (EstCapLen[r, ] + eps))) * nlob[r]
  }

  ## Create growth trajectory output (30 years starting from smallest bin)
  lenout <- matrix(0, ncol = length(lbin), nrow = 30 * ntsteps)
  lenout[1, 1] <- 1
  tsteps_seq <- 1:ntsteps
  cnt <- 0
  for (y in 1:(nrow(lenout) / ntsteps)) {
    for (ts in 1:ntsteps) {
      cnt <- cnt + 1
      if (cnt > 1) lenout[cnt, ] <- lenout[cnt - 1, ]
      if (ts %in% goodts) {
        tmpstm <- stm[, , tsteps_seq[ts]]
        lenout[cnt, ] <- tmpstm %*% lenout[cnt, ]
      }
    }
  }

  sigGrowvec <- rep(sigGrowsd, length(seq(0, 5, 0.1)))

  # Calculate penalties
  PensigGrowsd <- -dnorm(sigGrowsd, log(2.0), 0.5, log = TRUE)
  PenSigError <- -dnorm(LsigError, log(2.0), 0.5, log = TRUE)
  PenMerrorRel <- -sum(dnorm(0, MerrorRel, exp(LMerrorRelsigma), log = TRUE))
  PenMerrorRec <- -sum(dnorm(0, MerrorRec, exp(LMerrorRecsigma), log = TRUE))
  smooth_penalty <- smoother * sum((growth_vecpar[2:nlbin] - growth_vecpar[1:(nlbin - 1)])^2)

  # Total negative log-likelihood
  TLL <- -sum(LL) + PensigGrowsd + PenSigError + PenMerrorRel + PenMerrorRec + smooth_penalty

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

  TLL
}
