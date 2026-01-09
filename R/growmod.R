#' Growth Model Objective Function for Tag-Recapture Data
#'
#' Core RTMB objective function that constructs size transition matrices (STMs)
#' from asymmetric inverse logistic growth curves and calculates negative
#' log-likelihood for tag-recapture observations with measurement error.
#'
#' @param pin Named list of parameters
#'
#' @details
#' This function is designed to be passed to `RTMB::MakeADFun()` and expects
#' the following objects in the calling environment:
#'
#' **From `datain` (via `getAll()`):**
#' \describe{
#'   \item{nlbin}{Number of length bins}
#'   \item{ntsteps}{Number of time steps (typically months) per year}
#'   \item{goodts}{Integer vector of time steps for which to estimate growth}
#'   \item{lbin, lbinL, lbinU}{Length bin midpoints, lower, and upper bounds}
#'   \item{Tlbin, Clbin}{Release and recapture length bin indices}
#'   \item{Rlcl, Rccl}{Release and recapture carapace lengths}
#'   \item{relts}{Release time step for each observation}
#'   \item{tsteps}{Time steps at liberty for each observation}
#'   \item{nlob}{Number of lobsters per observation (for weighting)}
#' }
#'
#' **From `pin` (parameters to estimate):**
#' \describe{
#'   \item{mxpin, mnpin}{Log-scale maximum and minimum growth increments by time step}
#'   \item{ipin}{Inflection point of growth curve by time step}
#'   \item{spin_left, spin_right}{Log-scale steepness parameters (left/right of inflection)}
#'   \item{LsigError}{Log-scale measurement error standard deviation}
#'   \item{sigGrow}{Log-scale growth variability}
#'   \item{MerrorRel, MerrorRec}{Measurement error random effects (release/recapture)}
#'   \item{LMerrorRelsigma, LMerrorRecsigma}{Log-scale random effect standard deviations}
#' }
#'
#' **Growth Model:**
#' Uses an asymmetric inverse logistic function with smoothly transitioning steepness:
#' \deqn{growth = \frac{max - min}{1 + \exp((L - inflection)/steepness)} + min}
#' where steepness transitions from `spin_left` to `spin_right` across the inflection point.
#'
#' **Likelihood:**
#' Kullback-Leibler divergence between predicted and observed recapture length
#' distributions, with penalty terms for measurement error priors and parameter
#' constraints.
#'
#' @return Scalar value of total negative log-likelihood (to be minimized).
#'   Model outputs are made available via `REPORT()` for later extraction with
#'   `mod$rep()`.
#'
#' @examples
#' \dontrun{
#' # Prepare data and parameters
#' datain <- list(...)  # Tag-recapture data
#' pin <- list(...)     # Initial parameter values
#'
#' # Fit model
#' mod <- MakeADFun(data = datain, parameters = pin,
#'                  map = Mapfunc(pin), DLL = "growmod")
#' opt <- nlminb(mod$par, mod$fn, mod$gr)
#' }
#'
#' @export
growmod <- function(pin){
  getAll(datain, pin, warn=FALSE)
  npar <- length(names(pin))
  nobs <- length(Tlbin)
  IsMin <- max(grepl('mnpin',names(pin)))
  stm <- array(0, c(nlbin, nlbin,ntsteps)) # for growth with sigma scaler
  EstRecLen <- matrix(0, ncol=length(lbin), nrow=nobs)
  EstCapLen <- matrix(0, ncol=length(lbin), nrow=nobs)
  moultinc <- matrix(0, ncol=length(lbin), nrow=ntsteps)
  growthmat <- matrix(0, ncol=length(lbin), nrow=ntsteps)
  estmnlen <- rep(0, nobs)
  LL <- rep(0, nobs)  ## final LL of each obs
  ##Sigmas
  sigError <- exp(LsigError)  ## Single sigma for measurement error - spreading the lbins
  sigGrowsd <- exp(sigGrow)

  ## Make STMs
  for(ns in goodts){
    # Transition centered at inflection point
    transition_width <- 5  # Width of transition region (in mm)
    weight <- 1 / (1 + exp(-(lbin - ipin[ns]) / transition_width))

    # Weighted average of left and right steepness
    s <- exp(spin_left[ns]) * (1 - weight) + exp(spin_right[ns]) * weight

    # Calculate growth with smoothly-varying steepness
    growthmat[ns,] <- (exp(mxpin[ns])-exp(mnpin[ns])) / (1+exp((lbin-ipin[ns])/s)) + exp(mnpin[ns])

    # Build STM (same as before but with vector on second loop)
    for(fm in 1:nlbin){
      growth <- growthmat[ns, fm]
      # Probability of growing from fm to each larger bin
      stm[fm:nlbin, fm, ns] <- pnorm(lbinU[fm:nlbin], lbin[fm] + growth, sigGrowsd) -
        pnorm(lbinL[fm:nlbin], lbin[fm] + growth, sigGrowsd)
    }
  }

  #Run through each observation and grow to calculate LL including measurement error
  for(r in 1:nobs){
    # make structures and release lobster with error
    Relength <- MerrorRel[r]+Rlcl[r]
    lens <- pnorm(lbinU, Relength, sigError) - pnorm(lbinL, Relength, sigError)
    ## Make recapture time structure
    tstepsvec <- c(relts[r]:ntsteps,rep(1:ntsteps,30))[1:tsteps[r]]
    # Run through time at liberty and grow when appropriate tstep occurs
    if(length(tstepsvec)>0) {
      for(ts in 1:length(tstepsvec)){
        cnt <- 0
        if(tstepsvec[ts]%in%goodts){
          tmpstm <- stm[,,tstepsvec[ts]]
          cnt <- cnt + 1
          lens <- tmpstm %*% lens}
      }}
    EstRecLen[r,] <- lens  ## Record estimate length distribution
    Reclength <- MerrorRec[r]+Rccl[r]
    EstCapLen[r,] <- pnorm(lbinU, Reclength, sigError) - pnorm(lbinL,Reclength, sigError)


    ## Determine likelihood
    LL[r] = sum((EstCapLen[r,]+1e-5)*log((EstRecLen[r,]+1e-5)/(EstCapLen[r,]+1e-5)))*nlob[r]
  }

  ## Make output to examine the growth without Measurement error
  lenout <- matrix(0, ncol=length(lbin), nrow=30*ntsteps); lenout[1,1] <- 1
  tsteps <- 1:ntsteps; cnt <- 0
  for(y in 1:(nrow(lenout)/ntsteps)){
    for(ts in 1:ntsteps){
      cnt <- cnt + 1
      if(cnt>1) lenout[cnt,] <- lenout[cnt-1,]
      if(ts%in%goodts){
        tmpstm <- stm[,,tsteps[ts]]
        lenout[cnt,] <- tmpstm %*% lenout[cnt,]
      }
    }
  }

  sigGrowvec <- rep(exp(sigGrow), length(seq(0,5,0.1)))
  #sigGrowvec <- sigGrowvec + min(sigGrowvec)#+exp(sigGrow[2])

  # Measurement error priors
  PenSigError <- -dnorm(LsigError, log(2), 0.5, log=TRUE)

  ## Stop min going too small
  minPen <- -sum(dnorm(-5, mnpin, 1, T))
  spin_leftPen <- -sum(dnorm(-5, spin_left, 5, T))

  ## Keep RE on release and recapture  close to zero
  PenMerrorRel <- -sum(dnorm(0,MerrorRel,exp(LMerrorRelsigma),log = TRUE))
  PenMerrorRec <- -sum(dnorm(0,MerrorRec,exp(LMerrorRecsigma),log = TRUE))

  TLL <- -sum(LL) + PenSigError + PenMerrorRel + PenMerrorRec + spin_leftPen + minPen
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
