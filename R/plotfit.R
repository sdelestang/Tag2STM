#' Plot Tag-Recapture Model Diagnostics
#'
#' Creates a multi-panel diagnostic plot for fitted tag-recapture growth models,
#' showing residuals, growth curves, length distributions over time, and projected
#' growth trajectories with confidence intervals.
#'
#' @details
#' This function expects the following objects to exist in the calling environment:
#' \describe{
#'   \item{mod}{A fitted RTMB model object from `MakeADFun()`}
#'   \item{tdat}{Data frame with tag-recapture observations (columns: rccl, Clbin, Tlbin)}
#'   \item{lbin}{Numeric vector of length bin midpoints}
#'   \item{datain}{List containing model input data (must include tsteps)}
#'   \item{goodts}{Integer vector of valid time step indices}
#'   \item{ntsteps}{Number of time steps per year}
#'   \item{lbinL}{Numeric vector of length bin limits}
#'   \item{M}{Natural mortality rate}
#' }
#'
#' The function creates six diagnostic panels:
#' \enumerate{
#'   \item Residuals vs. release length bin
#'   \item Residuals vs. liberty period (log months)
#'   \item Growth increment curves by length
#'   \item Growth variability (sigma) across increment sizes
#'   \item Length distributions at 1, 5, 10, 15, and 20 years
#'   \item Projected growth trajectory with 95% confidence intervals
#' }
#'
#' @return Invisibly returns the model report object from `mod$rep()`.
#'
#' @examples
#' \dontrun{
#' # After preparing data and fitting model
#' mod <- MakeADFun(...)
#' out <- plotfit()
#' }
#'
#' @importFrom dplyr mutate
#' @importFrom magrittr %>% %<>%
#' @export
plotfit <- function(datIn = tdat) {
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(dplyr)

  out <- mod$rep()
  lbin  <- bins$lbin
  lbinL <- bins$lbinL
  lenout <- out$lenout

  estlen <- apply(out$EstRecLen, 1, function(x) weighted.mean(lbin, w = x))

  obs <- data.frame(
    tag    = datIn$tag,
    relyr  = datIn$relyr,
    recyr  = datIn$recyr,
    estlen = estlen,
    obslen = datIn$rccl,
    Rlcl   = datIn$rlcl,
    ntstep = datIn$ntstep,
    lat    = datIn$dLat,
    lon    = datIn$dLon
  ) %>% mutate(
    res      = obslen - estlen,
    midyr    = (relyr + recyr) / 2,
    midyr_rd = round(midyr),
    Clbin    = as.numeric(cut(datIn$rccl,
                              breaks = c(bins$lbinL, max(bins$lbinL) + 30),
                              include.lowest = TRUE, right = FALSE))
  )

  pa <- ggplot(obs, aes(x = Rlcl, y = res, colour = midyr)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_colour_viridis_c(option = "viridis", name = "Mid-year") +
    labs(x = "Release size", y = "Residual", tag = "a") +
    theme_bw() +
    theme(legend.position = "none", plot.tag = element_text(face = "bold"))

  pb <- ggplot(obs, aes(x = log(ntstep + 1), y = res, colour = midyr)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_colour_viridis_c(option = "viridis", name = "Mid-year") +
    labs(x = "Liberty (log number timesteps)", y = "Residual", tag = "b") +
    theme_bw() +
    theme(legend.position = "none", plot.tag = element_text(face = "bold"))

  midyr_vals <- sort(unique(obs$midyr_rd))
  nmid <- length(midyr_vals)

  pc <- ggplot(obs, aes(x = factor(midyr_rd), y = res)) +
    geom_boxplot(fill = "grey90", alpha = 0.7, outlier.size = 1,
                 colour = "grey40", linewidth = 0.4) +
    stat_summary(aes(colour = factor(midyr_rd)), fun = median,
                 geom = "crossbar", width = 0.75, linewidth = 1,
                 show.legend = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
    scale_colour_viridis_d(option = "viridis") +
    guides(colour = "none") +
    labs(x = "Mid-year at liberty", y = "Residual", tag = "c") +
    theme_bw() +
    theme(plot.tag = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  growthmat <- out$growthmat[goodts, ]
  if (is.matrix(growthmat)) {
    gdf <- do.call(rbind, lapply(1:nrow(growthmat), function(i)
      data.frame(lbin = lbin, increment = growthmat[i, ], ts = factor(i))))
  } else {
    gdf <- data.frame(lbin = lbin, increment = growthmat, ts = factor(1))
  }

  pd <- ggplot(gdf, aes(x = lbin, y = increment, colour = ts)) +
    geom_line(linewidth = 0.8) +
    scale_colour_viridis_d(option = "plasma") +
    guides(colour = "none") +
    labs(x = "Size", y = "Increment", title = "Growth/Moult", tag = "d") +
    theme_bw() +
    theme(plot.tag = element_text(face = "bold"), plot.title = element_text(hjust = 0.5))

  sigdf <- data.frame(growth = seq(0, 5, 0.1), sigma = out$sigGrowvec)

  pe <- ggplot(sigdf, aes(x = growth, y = sigma)) +
    geom_point(size = 1.5) +
    geom_line() +
    labs(x = "Growth", y = "Growth spread", tag = "e") +
    theme_bw() +
    theme(plot.tag = element_text(face = "bold"))

  weighted.probs <- function(x) quantile(rep(lbin, 1e+05 * x), probs = c(0.025, 0.5, 0.975))
  mx_ci  <- apply(lenout, 1, weighted.probs)
  xs     <- 0:(ntsteps * 30)
  Mxage  <- trunc((3.5 / M) / 5) * 5 + 5
  Mxtstep <- Mxage * ntsteps
  keep   <- xs < Mxtstep

  traj_df <- data.frame(
    ts  = xs[keep],
    med = c(lbin[1], mx_ci[2, ])[keep],
    lo  = c(lbin[1], mx_ci[1, ])[keep],
    hi  = c(lbin[1], mx_ci[3, ])[keep]
  )

  age_breaks <- seq(2, (ntsteps * Mxage) + 2, ntsteps * 2)
  age_labels  <- seq(2, Mxage + 2, 2)

  pf <- ggplot(traj_df, aes(x = ts)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70", alpha = 0.8) +
    geom_line(aes(y = med), linewidth = 0.8) +
    scale_x_continuous(breaks = age_breaks, labels = age_labels) +
    scale_y_continuous(limits = c(0, max(lbinL))) +
    labs(x = "Relative Age (y)", y = "Size", tag = "f") +
    theme_bw() +
    theme(plot.tag = element_text(face = "bold"))

  has_spatial <- any(!is.na(obs$lat) & !is.na(obs$lon))

  if (has_spatial) {
    obs_sp <- obs %>%
      filter(!is.na(lat) & !is.na(lon)) %>%
      mutate(
        direction = ifelse(res >= 0, "Positive", "Negative"),
        dot_size  = sqrt(abs(res))
      )

    pg <- ggplot(obs_sp, aes(x = lon, y = lat, colour = direction, size = dot_size)) +
      geom_point(alpha = 0.4) +
      scale_colour_manual(values = c(Positive = "#3481c8", Negative = "#e84b33"),
                          name = "Residual") +
      scale_size_continuous(range = c(0.5, 3)) +
      guides(size = "none") +                          # suppress size legend
      labs(x = "Longitude", y = "Latitude", tag = "g") +
      theme_bw() +
      theme(plot.tag = element_text(face = "bold"))

    layout <- "
AB
CD
EF
GH
"
    final <- pa + pb + pc + pg + pd + pe + pf  + guide_area() +
      plot_layout(design = layout, guides = "collect") +
      plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2))) &
      theme(legend.position    = "right",
            legend.box         = "vertical",        # stack legends on top of each other
            legend.direction   = "horizontal",      # each legend runs horizontally
            legend.margin      = margin(0, 0, 0, 0),
            legend.box.margin  = margin(0, 0, 0, 0),
            plot.margin        = margin(2, 4, 2, 4))

  } else {
    final <- (pa + pb) / (pc + pd) / (pe + pf) +
      plot_layout(guides = "collect") +
      plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2))) &
      theme(legend.position = "bottom",
            legend.margin    = margin(0, 0, 0, 0),
            legend.box.margin = margin(-5, 0, 0, 0),
            plot.margin      = margin(2, 4, 2, 4))
  }

  suppressWarnings(print(final))
  return(list(out = out, obs = obs))
}
