#' Estimate Annual Tag Loss Rate from Double-Tagged Recaptures
#'
#' Estimates the annual tag loss rate from animals that were released with two
#' tags and subsequently recaptured. Animals recaptured with only one tag
#' provide direct evidence of tag loss independent of natural or fishing
#' mortality, since recapture confirms the animal survived. A likelihood
#' profile is used to derive confidence intervals.
#'
#' @param double_tags Data frame of double-tagged animals that were recaptured,
#'   containing columns:
#'   \itemize{
#'     \item \code{years_at_liberty} - Time between release and recapture
#'       (numeric, in years or fractional years). Must be > 0.
#'     \item \code{tags_at_recapture} - Number of tags present at recapture
#'       (integer, 1 or 2 only). Animals recaptured with 0 tags should be
#'       excluded as tag identity cannot be confirmed.
#'   }
#'   Only recaptured animals should be included — animals never recaptured
#'   cannot contribute to this likelihood since we cannot distinguish tag
#'   loss from mortality.
#'
#' @param plot Logical. If TRUE (default), produces a likelihood profile plot
#'   for the tag loss rate with MLE and 95\% CI marked.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{tag_loss}    - MLE annual tag loss rate (per tag per year)
#'     \item \code{tl_lo95}    - lower 95\% profile CI
#'     \item \code{tl_hi95}    - upper 95\% profile CI
#'     \item \code{n}          - number of recaptured double-tagged animals
#'     \item \code{n_lost_one} - number retaining only 1 tag at recapture
#'     \item \code{n_both}     - number retaining both tags at recapture
#'     \item \code{profile}    - data frame of tag loss values and profile NLL
#'   }
#'
#' @details
#' **Model:**
#'
#' For an animal released with 2 tags and at liberty for \eqn{t} years,
#' the per-tag annual loss rate is \eqn{\delta}. Tag losses are assumed
#' independent between the two tags. The probabilities of retaining both
#' tags or exactly one tag at recapture are:
#'
#' \deqn{p_{both} = (1 - \delta)^{2t}}
#' \deqn{p_{one}  = 2 \cdot \left[1-(1-\delta)^t\right] \cdot (1-\delta)^t}
#'
#' where \eqn{(1-\delta)^t} is the probability of retaining a single tag
#' over \eqn{t} years, and the factor of 2 accounts for either tag being
#' the one lost.
#'
#' Conditional on recapture (i.e. excluding animals that lost both tags and
#' could not be identified), the normalised probabilities are:
#'
#' \deqn{p_{both|recap}  = \frac{p_{both}}{p_{both} + p_{one}}}
#' \deqn{p_{one|recap}   = \frac{p_{one}}{p_{both} + p_{one}}}
#'
#' The log-likelihood is:
#'
#' \deqn{\ell(\delta) = \sum_i \left[ \mathbf{1}(n_i=2) \log(p_{both,i|recap}) +
#'                                     \mathbf{1}(n_i=1) \log(p_{one,i|recap}) \right]}
#'
#' where the sum is over all recaptured double-tagged animals \eqn{i}.
#'
#' **Liberty time:** The model correctly handles varying liberty times —
#' an animal at liberty longer has had more opportunity to lose a tag, so
#' its contribution to the likelihood accounts for the expected cumulative
#' tag retention over its specific liberty period.
#'
#' **Independence assumption:** Tag losses are assumed independent between
#' the two tags on the same animal. If tags interact (e.g. both attached
#' near the same wound site), this assumption may be violated and tag loss
#' will be underestimated.
#'
#' **Use in Brownie model:** The estimated \code{tag_loss} should be passed
#' directly as the \code{tag_loss} argument to \code{EstimateMBrownie()}.
#' Since tag loss and M are inversely confounded in the Brownie model
#' (both explain missing animals), having an independent tag loss estimate
#' from double-tags is essential for unbiased M estimation.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' result <- EstimateTagLoss(double_tags)
#' cat("Annual tag loss:", result$tag_loss, "\n")
#' cat("95% CI:", result$tl_lo95, "-", result$tl_hi95, "\n")
#'
#' # Pass result directly to EstimateMBrownie
#' m_result <- EstimateMBrownie(obs, MLL=115, tag_loss=result$tag_loss)
#'
#' # Sensitivity: how does M change across the tag loss CI?
#' lapply(c(result$tl_lo95, result$tag_loss, result$tl_hi95), function(tl) {
#'   r <- EstimateMBrownie(obs, MLL=115, tag_loss=tl, plot=FALSE)
#'   data.frame(tag_loss=round(tl,3), M=round(r$M,3),
#'              lo=round(r$M_lo95,3), hi=round(r$M_hi95,3))
#' }) |> do.call(what=rbind)
#' }
#'
#' @seealso \code{\link{EstimateMBrownie}} for the Brownie dead-recovery model
#'   that uses the tag loss estimate.
#'
#' @export
EstimateTagLoss <- function(double_tags, plot = TRUE) {

  ## ── 0. Input checks ──────────────────────────────────────────────────────
  required <- c("years_at_liberty", "tags_at_recapture")
  missing  <- setdiff(required, names(double_tags))
  if (length(missing) > 0)
    stop("double_tags missing columns: ", paste(missing, collapse=", "))
  if (any(double_tags$years_at_liberty <= 0, na.rm=TRUE))
    stop("years_at_liberty must be > 0 for all animals.")
  if (!all(double_tags$tags_at_recapture %in% c(1L, 2L)))
    stop("tags_at_recapture must be 1 or 2 only. Remove animals with 0 tags.")
  if (any(is.na(double_tags$years_at_liberty) |
          is.na(double_tags$tags_at_recapture)))
    stop("No NAs allowed in years_at_liberty or tags_at_recapture.")

  n          <- nrow(double_tags)
  n_lost_one <- sum(double_tags$tags_at_recapture == 1)
  n_both     <- sum(double_tags$tags_at_recapture == 2)

  message(sprintf("Double-tag recaptures: %d total  (%d both tags, %d one tag)",
                  n, n_both, n_lost_one))

  ## ── 1. NLL function ───────────────────────────────────────────────────────
  # Parameterise on log scale to enforce positivity
  nll <- function(log_tl) {
    tl <- exp(log_tl)
    t  <- double_tags$years_at_liberty

    # Per-tag retention over t years
    retain_t <- (1 - tl)^t

    # Joint probabilities for a double-tagged animal
    p_both <- retain_t^2                    # both tags retained
    p_one  <- 2 * (1 - retain_t) * retain_t  # exactly one tag retained

    # Conditional on recapture (i.e. at least one tag retained)
    p_total <- p_both + p_one
    p_total <- pmax(p_total, 1e-15)         # numerical floor

    ll <- sum(
      ifelse(double_tags$tags_at_recapture == 2,
             log(pmax(p_both / p_total, 1e-15)),
             log(pmax(p_one  / p_total, 1e-15)))
    )
    -ll
  }

  ## ── 2. Optimise ───────────────────────────────────────────────────────────
  opt     <- optimise(nll, interval=c(log(1e-4), log(0.99)), tol=1e-8)
  tl_mle  <- exp(opt$minimum)
  nll_mle <- opt$objective

  message(sprintf("MLE tag loss: %.4f (%.1f%% per tag per year)", tl_mle, tl_mle*100))

  ## ── 3. Likelihood profile for CI ─────────────────────────────────────────
  tl_grid  <- exp(seq(log(max(1e-4, tl_mle * 0.05)),
                      log(min(0.99,  tl_mle * 20)),
                      length.out = 200))
  prof_nll <- sapply(log(tl_grid), nll)

  chi_thresh <- 1.92
  in_ci      <- (prof_nll - nll_mle) <= chi_thresh
  tl_lo95    <- if (any(in_ci)) min(tl_grid[in_ci]) else NA_real_
  tl_hi95    <- if (any(in_ci)) max(tl_grid[in_ci]) else NA_real_

  profile_df <- data.frame(tag_loss  = tl_grid,
                            nll       = prof_nll,
                            delta_nll = prof_nll - nll_mle,
                            in_ci95   = in_ci)

  message(sprintf("95%% CI: %.4f - %.4f", tl_lo95, tl_hi95))

  ## ── 4. Diagnostic plot ────────────────────────────────────────────────────
  if (plot) {
    old_par <- par(no.readonly=TRUE)
    on.exit(par(old_par), add=TRUE)
    par(mfrow=c(1,2), mar=c(4,4,3,1))

    ## (a) Likelihood profile
    ylim_a <- c(0, min(10, max(profile_df$delta_nll[is.finite(profile_df$delta_nll)])))
    plot(profile_df$tag_loss, profile_df$delta_nll,
         type="l", lwd=2, col="#185FA5",
         xlab="Annual tag loss rate (per tag)",
         ylab=expression(Delta~NLL),
         main="Likelihood Profile: Tag Loss Rate",
         ylim=ylim_a)
    abline(h  = chi_thresh, lty=2, col="red")
    abline(v  = tl_mle,  lty=1, col="black",  lwd=1.5)
    if (!is.na(tl_lo95)) abline(v=tl_lo95, lty=3, col="grey40")
    if (!is.na(tl_hi95)) abline(v=tl_hi95, lty=3, col="grey40")
    legend("topright",
           legend=c(sprintf("MLE = %.3f (%.1f%%/yr)", tl_mle, tl_mle*100),
                    sprintf("95%% CI: %.3f - %.3f", tl_lo95, tl_hi95),
                    "Chi-sq threshold (1.92)"),
           lty=c(1,3,2), col=c("black","grey40","red"),
           bty="n", cex=0.85)

    ## (b) Observed vs expected proportion retaining both tags by liberty time
    # Bin liberty times for display
    breaks   <- quantile(double_tags$years_at_liberty,
                         probs=seq(0,1,length.out=6), na.rm=TRUE)
    breaks   <- unique(breaks)
    if (length(breaks) < 3) breaks <- seq(min(double_tags$years_at_liberty),
                                           max(double_tags$years_at_liberty),
                                           length.out=5)
    lib_class <- cut(double_tags$years_at_liberty, breaks=breaks,
                     include.lowest=TRUE)
    midpoints <- tapply(double_tags$years_at_liberty, lib_class, mean, na.rm=TRUE)
    obs_prop  <- tapply(double_tags$tags_at_recapture == 2, lib_class, mean, na.rm=TRUE)
    n_class   <- tapply(double_tags$tags_at_recapture,      lib_class, length)

    # Predicted proportion retaining both tags at MLE
    t_seq     <- seq(min(double_tags$years_at_liberty),
                     max(double_tags$years_at_liberty),
                     length.out=100)
    retain_t  <- (1 - tl_mle)^t_seq
    p_both_t  <- retain_t^2
    p_one_t   <- 2 * (1 - retain_t) * retain_t
    pred_prop <- p_both_t / (p_both_t + p_one_t)

    ylim_b <- c(0, 1)
    plot(as.numeric(midpoints), as.numeric(obs_prop),
         pch=16, cex=sqrt(as.numeric(n_class)/max(n_class))*2.5,
         col="#185FA580",
         xlab="Years at liberty",
         ylab="Proportion retaining both tags",
         main="Observed vs Predicted Tag Retention",
         ylim=ylim_b,
         xlim=range(double_tags$years_at_liberty))
    lines(t_seq, pred_prop, col="#E05C1A", lwd=2)
    legend("topright",
           legend=c("Observed (size ~ n)", "Predicted"),
           pch=c(16, NA), lty=c(NA,1), lwd=c(NA,2),
           col=c("#185FA580","#E05C1A"),
           bty="n", cex=0.85)

    mtext(sprintf("n=%d recaptured  (%d both tags, %d one tag)",
                  n, n_both, n_lost_one),
          side=3, outer=TRUE, line=-1.2, cex=0.8)
  }

  ## ── 5. Return ─────────────────────────────────────────────────────────────
  invisible(list(
    tag_loss    = tl_mle,
    tl_lo95     = tl_lo95,
    tl_hi95     = tl_hi95,
    n           = n,
    n_lost_one  = n_lost_one,
    n_both      = n_both,
    profile     = profile_df
  ))
}
