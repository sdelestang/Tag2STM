#' Convert Von Bertalanffy Growth Parameters to Size Transition Matrix
#'
#' Creates a size transition matrix (STM) from Von Bertalanffy growth parameters
#' with no-shrinkage constraint. The function accounts for truncation bias when
#' preventing backward transitions using truncated normal distribution theory.
#'
#' @param lbinmin Numeric. Minimum length bin (mm or cm)
#' @param lbinmax Numeric. Maximum length bin (mm or cm)
#' @param gap Numeric. Width of length bins (mm or cm)
#' @param Linf Numeric. Von Bertalanffy asymptotic length parameter
#' @param k Numeric. Von Bertalanffy growth rate parameter (per year)
#' @param to Numeric. Von Bertalanffy theoretical age at length 0
#' @param growsd Numeric. Standard deviation of growth variability (measurement error
#'   and individual variation). Larger values increase uncertainty in transitions.
#' @param maxage Numeric. Maximum age for diagnostic plot (default = 30 years)
#'
#' @return A size transition matrix where each column sums to 1, representing
#'   transition probabilities from one length bin (column) to another (row) over
#'   one year. Matrix is truncated to the range [lbinmin, lbinmax].
#'
#' @details
#' The function implements a no-shrinkage constraint by:
#' \itemize{
#'   \item Truncating the normal growth distribution at the current bin
#'   \item Adjusting the mean growth using inverse Mills ratio to correct for
#'     truncation bias
#'   \item Adding the truncated probability mass to the diagonal (staying in same bin)
#' }
#'
#' This approach maintains consistency with the Von Bertalanffy growth curve even
#' with large growth variability (growsd), unlike naive truncation which causes
#' upward bias in mean length-at-age.
#'
#' A diagnostic plot is produced showing STM predictions (black points with 95% CI)
#' against the Von Bertalanffy curve (red line).
#'
#' @examples
#' # Create STM for Western Rock Lobster
#' stm <- vb2STM(lbinmin = 20, lbinmax = 140, gap = 2,
#'               Linf = 120, k = 0.1, to = 0,
#'               growsd = 3, maxage = 30)
#'
#' # Check that columns sum to 1
#' colSums(stm)
#'
#' @seealso \code{\link{pnorm}}, \code{\link{dnorm}}
#'
#' @export
vb2STM <- function(lbinmin,lbinmax,gap, Linf, k, to, growsd, maxage=30)  {
  ## Make vb data
  dat <- data.frame(age=seq(0, 200, 0.01))
  dat$len <- Linf*(1-exp(-k*(dat$age-to)))

  # Make lbin structure to handle conversion - it must start below age 1 of vb
  lbinL <- seq(lbinmin-(gap*30),lbinmax,gap)
  begin <- max(match(floor(dat$len[dat$age<1]),lbinL), na.rm=T)
  lbinL<- lbinL[begin:length(lbinL)]

  #Make other bin constraints
  lbinU <- lbinL+gap
  lbin <- (lbinU+lbinL)/((gap-1)*2)
  lbinU[length(lbin)] <- lbinU[length(lbin)]+gap*20
  lbinL[1] <- abs(lbinL[1]-gap*10)

  nlbin <- length(lbin)
  #Function to match the lengthbin
  finflbin <- function(x) {
    which.min(abs(dat$len - x)) }
  pos <- apply(as.matrix(lbin), 1, finflbin)
  age1 <- (dat$age[pos])
  age2 <- age1+1

  ## Makes growth vector, increase over one year for each age
  fingwth <- function(x) {
    x <- round(x,2)
    growth <- dat$len[round(dat$age,2)==x]-dat$len[round(dat$age,2)==round((x-1),1)]
    if(length(growth)==0) growth <- 0
    return(growth)
  }

  growthvec <- apply(as.matrix(age2), 1, fingwth)
  # Blank stm to fill
  stm <- matrix(0, ncol=nlbin, nrow=nlbin)

  # Build STM with adjusted growth to account for truncation
  for(fm in 1:nlbin){
    growth <- growthvec[fm]

    # Calculate truncation point (no shrinkage below current bin)
    alpha <- (lbinL[fm] - (lbin[fm] + growth)) / growsd

    # Probability of truncation (would-be shrinkage)
    shrink_prob <- pnorm(alpha)

    # Adjusted mean to account for truncation bias
    # This is the mean of truncated normal: μ + σ * λ(α)
    # where λ(α) = φ(α)/(1-Φ(α)) is the inverse Mills ratio
    lambda <- dnorm(alpha) / (1 - pnorm(alpha))
    growth_adjusted <- growth - growsd * lambda

    # Calculate forward movement probabilities with adjusted growth
    stm[fm:nlbin, fm] <- pnorm(lbinU[fm:nlbin], lbin[fm] + growth_adjusted, growsd) -
      pnorm(lbinL[fm:nlbin], lbin[fm] + growth_adjusted, growsd)

    # Add shrinkage probability to diagonal
    stm[fm, fm] <- stm[fm, fm] + shrink_prob

    # Normalize column to sum to 1
    stm[fm:nlbin, fm] <- stm[fm:nlbin, fm]/sum(stm[fm:nlbin, fm])
  }

  ## Make output to examine the growth without Measurement error
  lenout <- matrix(0, ncol=length(lbin), nrow=30);

  ## Redo vb with same ages we will plot
  tdat <- data.frame(age=seq(1, maxage, 1))
  tdat$len <- Linf*(1-exp(-k*(tdat$age-to)))
  lenat1 <- tdat$len[tdat$age==1]
  pos <- which(abs(lbin-lenat1)==min(abs(lbin-lenat1)))
  lenout[pos ,1] <- 1
  for(y in 2:(nrow(lenout))){
    lenout[y,] <- stm %*% lenout[y-1,]
  }

  qfunc <- function(x,y=lbin) {
    x <- rep(y,x*1000)
    return(quantile(x, probs=c(0.025,0.5,0.975)))}
  mn <- apply(lenout,1,qfunc)
  par(las=1)
  plot(1:maxage, mn[2,], ylim=c(0,Linf*1.1), xlab='Age', ylab='Length', pch=16)
  polygon(c(tdat$age,rev(tdat$age)), c(mn[1,],rev(mn[3,])), border=F, col=rgb(0,1,0,0.3))
  lines(tdat$age, tdat$len, col=2, lwd=2)
  legend('top', col=c(1,2), pch=c(16,-1), lty=c(0,1), lwd=2, legend=c('STM', 'VonBert'))

  tokeep <- (lbinL>=lbinmin&lbinL<=lbinmax)
  stm2 <- stm[tokeep,tokeep]
  funcsum <- function(x) {  x/sum(x,na.rm=T)  }
  stm2 <- apply(stm2, 2, funcsum)

  return(stm2)
}
