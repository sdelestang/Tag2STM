#' Plot an Annual Projection Against the Compounded Sub-annual Growth
#'
#' Compares the 12-month STM compounded from the fitted sub-annual windows
#' (\code{ann$A}, the target) with the one-step annual STM built from the
#' projected annual parameters (\code{ann$Ahat}, what IMuLT will use), from
#' \code{\link{AnnualFromSeasonPars}}. Four panels:
#' \enumerate{
#'   \item Mean annual increment by length, with a +/- 1 SD band.
#'   \item Probability of staying in the same length bin over the year.
#'   \item Length distribution after one year from a few starting lengths:
#'     bars = compounded, points = annual fit.
#'   \item Length-at-age from the smallest bin, median and 95\% interval,
#'     which shows whether small per-bin differences accumulate over time.
#' }
#'
#' @param ann One element of the list returned by
#'   \code{\link{AnnualFromSeasonPars}} (needs \code{A}, \code{Ahat} and
#'   \code{moments}).
#' @param title Plot title, typically the IMuLT matrix name.
#' @param from Starting lengths (mm) for panel 3. Default: four lengths
#'   spread across the range.
#' @param nyears Years for the length-at-age panel. Default 25.
#'
#' @return A patchwork object (print it to draw).
#' @export
PlotAnnualProjection <- function(ann, title = NULL, from = NULL, nyears = 25) {
  requireNamespace("ggplot2"); requireNamespace("patchwork")
  A    <- ann$A
  Ahat <- ann$Ahat
  mom  <- ann$moments
  L    <- mom$L[mom$source == "compounded"]
  nb   <- length(L)
  cols <- c(compounded = "#185FA5", `annual fit` = "#D85A30")
  base <- ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "none", panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 9))

  ## 1. Mean increment +/- 1 SD
  p1 <- ggplot2::ggplot(mom, ggplot2::aes(L, inc, colour = source, fill = source)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(inc - sd, 0), ymax = inc + sd),
                         alpha = 0.15, colour = NA) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_colour_manual(values = cols) + ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(x = "Length (mm)", y = "Annual increment (mm)",
                  title = "Mean increment +/- 1 SD") + base +
    ggplot2::guides(colour = "none", fill = "none")

  ## 2. P(stay in bin)
  p2 <- ggplot2::ggplot(mom, ggplot2::aes(L, p_stay, colour = source)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "Length (mm)", y = "P(same bin after 1 yr)",
                  title = "Probability of no growth") + base +
    ggplot2::guides(colour = "none")

  ## 3. One-year length distributions from a few starting lengths
  if (is.null(from)) from <- L[round(seq(0.1, 0.75, length.out = 4) * nb)]
  jj <- unique(vapply(from, function(x) which.min(abs(L - x)), integer(1)))
  dd <- do.call(rbind, lapply(jj, function(j) {
    keep <- which(A[, j] > 1e-3 | Ahat[, j] > 1e-3)
    rbind(data.frame(start = paste0("from ", L[j], " mm"), L = L[keep], p = A[keep, j], source = "compounded"),
          data.frame(start = paste0("from ", L[j], " mm"), L = L[keep], p = Ahat[keep, j], source = "annual fit"))
  }))
  dd$start <- factor(dd$start, levels = unique(dd$start))
  p3 <- ggplot2::ggplot(dd, ggplot2::aes(L, p)) +
    ggplot2::geom_col(data = subset(dd, source == "compounded"), fill = cols["compounded"],
                      alpha = 0.6, width = diff(L)[1] * 0.9) +
    ggplot2::geom_point(data = subset(dd, source == "annual fit"), colour = cols["annual fit"], size = 1.2) +
    ggplot2::facet_wrap(~start, nrow = 1, scales = "free_x") +
    ggplot2::labs(x = "Length after 1 year (mm)", y = "Probability",
                  title = "One-year length distribution (bars = compounded, points = annual fit)") +
    base + ggplot2::theme(strip.text = ggplot2::element_text(size = 7),
                          axis.text.x = ggplot2::element_text(size = 6))

  ## 4. Length-at-age from the smallest bin
  traj <- function(M, src) {
    x <- c(1, rep(0, nb - 1)); out <- NULL
    for (y in 0:nyears) {
      cdf <- cumsum(x)
      q <- function(pr) L[min(which(cdf >= pr))]
      out <- rbind(out, data.frame(age = y, med = q(0.5), lo = q(0.025), hi = q(0.975), source = src))
      x <- as.numeric(M %*% x)
    }
    out
  }
  tj <- rbind(traj(A, "compounded"), traj(Ahat, "annual fit"))
  p4 <- ggplot2::ggplot(tj, ggplot2::aes(age, med, colour = source, fill = source)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    ggplot2::geom_step(linewidth = 0.7) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::scale_fill_manual(values = cols, name = NULL) +
    ggplot2::labs(x = "Years from smallest bin", y = "Length (mm)",
                  title = "Length-at-age, median and 95% interval") +
    base

  ((p1 | p2) / p3 / p4) +
    patchwork::plot_layout(heights = c(1, 0.8, 0.9), guides = "collect") +
    patchwork::plot_annotation(
      title = title,
      subtitle = sprintf("Projection objective (weighted KL) %.4g; convergence %d",
                         ann$objective, ann$convergence),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold"),
                             plot.subtitle = ggplot2::element_text(size = 8))) &
    ggplot2::theme(legend.position = "bottom")
}
