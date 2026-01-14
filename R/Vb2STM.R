#' Convert Von Bertalanffy Growth Parameters to Size Transition Matrix with No Shrinkage
#'
#' Creates a size transition matrix (STM) from Von Bertalanffy growth parameters
#' with a strict no-shrinkage constraint. The function optimizes growth increments
#' to ensure the resulting no-shrinkage STM matches the target Von Bertalanffy
#' growth trajectory when projected forward in time.
#'
#' @param lbinmin Numeric. Minimum length bin (mm or cm)
#' @param lbinmax Numeric. Maximum length bin (mm or cm)
#' @param gap Numeric. Width of length bins (mm or cm)
#' @param Linf Numeric. Von Bertalanffy asymptotic length parameter
#' @param k Numeric. Von Bertalanffy growth rate parameter (per year)
#' @param to Numeric. Von Bertalanffy theoretical age at length 0
#' @param growsd Numeric. Standard deviation of growth variability representing
#'   measurement error and individual variation. Larger values increase uncertainty
#'   in size transitions.
#' @param maxage Numeric. Maximum age for diagnostic plot (default = 30 years)
#'
#' @return A lower-triangular size transition matrix where each column sums to 1,
#'   representing transition probabilities from one length bin (column) to another
#'   (row) over one year. Matrix has zeros above the diagonal (no backward growth).
#'   Matrix is trimmed to the range [lbinmin, lbinmax].
#'
#' @details
#' This function addresses the challenge of creating a no-shrinkage size transition
#' matrix that accurately reproduces Von Bertalanffy growth trajectories. The
#' approach works in three steps:
#'
#' \strong{1. Initial Growth Vector:}
#' Calculates initial annual growth increments from the Von Bertalanffy curve for
#' each length bin.
#'
#' \strong{2. Optimization:}
#' Uses numerical optimization to find growth multipliers that, when applied to the
#' initial growth vector and processed through the no-shrinkage constraint, produce
#' an STM that matches the target Von Bertalanffy trajectory when projected forward.
#'
#' \strong{3. No-Shrinkage Constraint:}
#' Builds the STM allowing full forward and backward movement, then post-processes
#' to remove backward transitions (sets upper triangle to zero) and adds the lost
#' probability mass to the diagonal (staying in same bin).
#'
#' The optimization is essential because naively preventing shrinkage creates
#' systematic upward bias in projected length-at-age, particularly with high growth
#' rates (large k) or high variability (large growsd). By optimizing growth
#' parameters specifically for the no-shrinkage constraint, the function ensures
#' the final STM accurately represents the intended growth dynamics.
#'
#' A diagnostic plot is produced showing STM-projected median length-at-age (black
#' points) with 95% confidence intervals (green shading) against the target Von
#' Bertalanffy curve (red line). Close agreement indicates successful optimization.
#'
#' @examples
#' # Western Rock Lobster with slow growth (k=0.1)
#' stm1 <- vb2STM(lbinmin = 20, lbinmax = 140, gap = 2,
#'                Linf = 120, k = 0.1, to = 0,
#'                growsd = 2, maxage = 30)
#'
#' # Fast-growing species (k=0.3)
#' stm2 <- vb2STM(lbinmin = 20, lbinmax = 130, gap = 2,
#'                Linf = 120, k = 0.3, to = -1,
#'                growsd = 2, maxage = 30)
#'
#' # Check that columns sum to 1
#' colSums(stm1)
#'
#' # Verify no backward transitions
#' any(stm1[upper.tri(stm1)] > 0)  # Should be FALSE
#'
#' @note
#' The optimization prints progress messages showing the final sum of squared
#' errors (SSE) and range of growth multipliers applied. Lower SSE indicates
#' better fit to the target trajectory.
#'
#' @seealso \code{\link{pnorm}}, \code{\link{optim}}
#'
#' @export
vb2STM <- function(lbinmin, lbinmax, gap, Linf, k, to, growsd, maxage = 30) {

  dat <- data.frame(age = seq(0, 200, 0.01))
  dat$len <- Linf * (1 - exp(-k * (dat$age - to)))

  # Make lbin structure
  lbinL <- seq(lbinmin - (gap * 30), lbinmax, gap)
  begin <- max(match(floor(dat$len[dat$age < 1]), lbinL), na.rm = T)
  lbinL <- lbinL[begin:length(lbinL)]

  lbinU <- lbinL + gap
  lbin <- lbinL + (gap/2)
  lbinU[length(lbin)] <- lbinU[length(lbin)] + gap * 20
  lbinL[1] <- abs(lbinL[1] - gap * 10)

  nlbin <- length(lbin)

  # Get initial growth vector from VonBert
  finflbin <- function(x) { which.min(abs(dat$len - x)) }
  pos <- apply(as.matrix(lbin), 1, finflbin)
  age1 <- dat$age[pos]
  age2 <- age1 + 1

  fingwth <- function(x) {
    x <- round(x, 2)
    growth <- dat$len[round(dat$age, 2) == x] -
      dat$len[round(dat$age, 2) == round((x - 1), 1)]
    if (length(growth) == 0) growth <- 0
    return(growth)
  }

  growthvec_init <- apply(as.matrix(age2), 1, fingwth)
  growthvec_init[growthvec_init <= 0] <- 0.01  # Ensure positive

  # Target VonBert trajectory
  tdat <- data.frame(age = seq(1, maxage, 1))
  tdat$len <- Linf * (1 - exp(-k * (tdat$age - to)))

  # Function to build NO-SHRINKAGE STM and test it
  build_noshrink_stm <- function(growthvec) {
    stm <- matrix(0, ncol = nlbin, nrow = nlbin)

    # Build with full forward/backward
    for (fm in 1:nlbin) {
      growth <- growthvec[fm]
      stm[, fm] <- pnorm(lbinU, lbin[fm] + growth, growsd) -
        pnorm(lbinL, lbin[fm] + growth, growsd)
      stm[, fm] <- stm[, fm] / sum(stm[, fm])
    }

    # POST-PROCESS: Remove backward transitions
    orig_sums <- colSums(stm)
    stm[upper.tri(stm)] <- 0
    new_sums <- colSums(stm)
    diag(stm) <- diag(stm) + (orig_sums - new_sums)

    # Project THIS no-shrinkage STM forward
    lenout <- matrix(0, ncol = nlbin, nrow = maxage)
    lenat1 <- tdat$len[tdat$age == 1]
    pos <- which.min(abs(lbin - lenat1))
    lenout[1, pos] <- 1

    for (y in 2:maxage) {
      lenout[y, ] <- stm %*% lenout[y - 1, ]
    }

    # Calculate median length at age from THIS no-shrinkage STM
    median_len <- apply(lenout, 1, function(x) sum(x * lbin))

    # Return SSE and STM
    sse <- sum((median_len - tdat$len)^2)
    return(list(sse = sse, stm = stm, median_len = median_len))
  }

  # Optimize using multipliers (more stable than log)
  cat("Optimizing growth parameters for no-shrinkage STM...\n")
  opt_result <- optim(
    par = rep(1, nlbin),  # Start with multipliers of 1
    fn = function(pars) {
      # Bounds check
      if (any(pars < 0.1) || any(pars > 3)) return(1e10)

      growthvec <- growthvec_init * pars
      result <- tryCatch({
        build_noshrink_stm(growthvec)$sse
      }, error = function(e) {
        return(1e10)
      })
      return(result)
    },
    method = "L-BFGS-B",
    lower = rep(0.1, nlbin),
    upper = rep(3, nlbin),
    control = list(maxit = 100)
  )

  # Build final optimized no-shrinkage STM
  growthvec_opt <- growthvec_init * opt_result$par
  final <- build_noshrink_stm(growthvec_opt)
  stm <- final$stm

  cat("Optimization complete. SSE:", round(opt_result$value, 2), "\n")
  cat("Growth multiplier range:", round(range(opt_result$par), 3), "\n")

  # Project final STM for plotting
  lenout <- matrix(0, ncol = nlbin, nrow = maxage)
  lenat1 <- tdat$len[tdat$age == 1]
  pos <- which.min(abs(lbin - lenat1))
  lenout[1, pos] <- 1

  for (y in 2:maxage) {
    lenout[y, ] <- stm %*% lenout[y - 1, ]
  }

  x <- lenout[2,]
  # Plot
  qfunc <- function(x, y = lbin) {
    mult <- ifelse((x * 100000)<0, 0,(x * 100000))
    x <- rep(y, mult)
    return(quantile(x, probs = c(0.025, 0.5, 0.975)))
  }
  mn <- apply(lenout, 1, qfunc)

  par(las = 1)
  plot(1:maxage, mn[2, ], ylim = c(0, Linf * 1.4),
       xlab = "Age", ylab = "Length", pch = 16)
  polygon(c(tdat$age, rev(tdat$age)), c(mn[1, ], rev(mn[3, ])),
          border = F, col = rgb(0, 1, 0, 0.3))
  lines(tdat$age, tdat$len, col = 2, lwd = 2)
  legend("top", col = c(1, 2), pch = c(16, -1), lty = c(0, 1), lwd = 2,
         legend = c("STM", "VonBert"))

  # Return trimmed no-shrinkage matrix
  tokeep <- (lbinL >= lbinmin & lbinL <= lbinmax)
  stm2 <- stm[tokeep, tokeep]
  funcsum <- function(x) { x / sum(x, na.rm = T) }
  stm2 <- apply(stm2, 2, funcsum)

  return(stm2)
}

