#' Plot Tag-Recapture Model Diagnostics
#'
#' Creates a multi-panel diagnostic plot for fitted tag-recapture growth models,
#' showing residuals, growth curves, length distributions over time, and projected
#' growth trajectories with confidence intervals.
#'
#' @param datIn Data frame of tag-recapture observations. Default \code{tdat}.
#' @param mout Optional. The result of \code{nlminb(mod$par, mod$fn, mod$gr, ...)}.
#'   If supplied, the model is reported explicitly at \code{mout$par} via
#'   \code{mod$report(mout$par)}, guaranteeing the plot reflects the actual
#'   fitted optimum. If omitted (\code{NULL}, the default), falls back to
#'   \code{mod$rep()}, which reports at \code{mod}'s current internal
#'   parameter state -- usually the optimum right after fitting, but this can
#'   silently go stale if anything else (e.g. \code{sdreport()},
#'   \code{tmbprofile()}) has called \code{mod$fn()}/\code{mod$gr()} since.
#'   Passing \code{mout} explicitly is recommended.
#'
#' @details
#' This function expects the following objects to exist in the calling environment:
#' \describe{
#'   \item{mod}{A fitted RTMB model object from `MakeADFun()` — `growmod()`,
#'     either `TemporalGrowth = FALSE` or `TRUE`. The function auto-detects
#'     which was used by checking whether `mod$rep()$S` is present, and
#'     builds the scenario envelope itself (via `build_lenout_envelope()`)
#'     when it is.}
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
#'   \item Residuals vs. mid-year at liberty (boxplots)
#'   \item Growth increment curves by length, with a growth-spread ribbon
#'     (mean +/- sd_growth, sd_growth interpolated from \code{sigGrowvec} at
#'     each bin's mean increment)
#'   \item If the model was fit with \code{TemporalGrowth = TRUE} (i.e.
#'     \code{mod$rep()$S} exists): year effects (\code{S}) shown as percent
#'     deviation in mean growth from the average year, one bar per modelled
#'     year. Otherwise, falls back to the original growth-spread-vs-increment-size panel.
#'   \item Projected growth trajectory with 95% CI ribbon for the average
#'     year. If \code{mod$rep()$lenout_scenarios} exists, the slowest- and
#'     fastest-fitted-year median trajectories are overlaid as dashed lines.
#'   \item (if spatial columns present) Residuals by release location
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
#' @importFrom dplyr mutate bind_rows
#' @importFrom magrittr %>% %<>%
#' @export
plotfit <- function(datIn = tdat, mout = NULL) {
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(dplyr)

  # mod$rep() alone reports at whatever mod's internal parameter state
  # currently is -- that state gets updated by ANY call to mod$fn()/mod$gr()
  # (e.g. sdreport(), tmbprofile(), or anything else run since fitting), so it
  # is not reliably "the fitted optimum" just because nlminb() ran earlier in
  # the session. Pass mout explicitly for a guaranteed-correct report.
  if (!is.null(mout)) {
    out <- mod$report(mout$par)
  } else {
    out <- mod$rep()
  }
  lbin  <- bins$lbin
  lbinL <- bins$lbinL

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
    scale_colour_viridis_c(option = "viridis", name = "Mid-year",
                            guide = guide_colourbar(barwidth = unit(5, "cm"), barheight = unit(0.4, "cm"))) +
    labs(x = "Release size", y = "Residual") +
    theme_bw() +
    theme(legend.position = "none", plot.tag = element_text(face = "bold"))

  pb <- ggplot(obs, aes(x = log(ntstep + 1), y = res, colour = midyr)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_colour_viridis_c(option = "viridis", name = "Mid-year",
                            guide = guide_colourbar(barwidth = unit(5, "cm"), barheight = unit(0.4, "cm"))) +
    labs(x = "Liberty (log number timesteps)", y = "Residual") +
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
    labs(x = "Mid-year at liberty", y = "Residual") +
    theme_bw() +
    theme(plot.tag = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  ## --- Panel d: growth increment curve(s), with a growth-spread ribbon AND
  ## P(moult)-by-size on a calibrated secondary axis (different units/scale
  ## to increment, so NOT plotted on the shared raw axis) ---
  growthmat <- out$growthmat[goodts, ]
  if (is.matrix(growthmat)) {
    gdf <- do.call(rbind, lapply(1:nrow(growthmat), function(i)
      data.frame(lbin = lbin, increment = growthmat[i, ], ts = factor(i))))
  } else {
    gdf <- data.frame(lbin = lbin, increment = growthmat, ts = factor(1))
  }
  # sd_growth is now proportional (sd = exp(LsigGrow) * increment), and
  # sigGrowvec = exp(LsigGrow) * gseq is exactly linear -- so this
  # interpolation is now an exact readout, not an approximation.
  gdf$sd <- approx(x = seq(0, 5, 0.1), y = out$sigGrowvec, xout = gdf$increment, rule = 2)$y

  # Pmoult_vec is a function of size only (not season), so the same nlbin
  # values repeat once per ts block in gdf's row order.
  gdf$Pmoult <- rep(out$Pmoult_vec, times = length(unique(gdf$ts)))

  # Rescale P(moult) [0,1] onto the ribbon's own [min, max] range, rather than
  # an arbitrary [0, some_constant] -- this makes the dashed P(moult) line
  # occupy exactly the same vertical space as the ribbon, and the secondary
  # axis (0 to 1) lines up directly with where the ribbon actually sits,
  # instead of spanning most of an axis the ribbon itself only occupies a
  # narrow band of.
  ribbon_lo <- min(pmax(gdf$increment - gdf$sd, 0), na.rm = TRUE)
  ribbon_hi <- max(gdf$increment + gdf$sd, na.rm = TRUE)
  ribbon_range <- ribbon_hi - ribbon_lo

  pd <- ggplot(gdf, aes(x = lbin)) +
    geom_ribbon(aes(ymin = pmax(increment - sd, 0), ymax = increment + sd, fill = ts),
                alpha = 0.15, colour = NA) +
    geom_line(aes(y = increment, colour = ts), linewidth = 0.8) +
    geom_line(aes(y = ribbon_lo + Pmoult * ribbon_range), colour = "black", linetype = "dashed", linewidth = 0.7) +
    scale_colour_viridis_d(option = "plasma") +
    scale_fill_viridis_d(option = "plasma") +
    guides(colour = "none", fill = "none") +
    scale_y_continuous(
      name = "Increment",
      sec.axis = sec_axis(~ (. - ribbon_lo) / ribbon_range, name = "P(moult)")
    ) +
    labs(x = "Size", title = "Growth/Moult (dashed = P(moult))") +
    theme_bw() +
    theme(plot.tag = element_text(face = "bold"), plot.title = element_text(hjust = 0.5))

  ## --- Panel e: year effects (TemporalGrowth = TRUE), falling back to growth-spread plot (TemporalGrowth = FALSE) ---
  has_year_effects <- !is.null(out$S)

  if (has_year_effects) {
    # Year 1 of S corresponds to the earliest calendar year spanned by the
    # tagging data — same convention used to derive datain$nyears
    # (diff(range(relyr, recyr)) + 1), so anchor the labels the same way.
    yr0  <- min(c(datIn$relyr, datIn$recyr), na.rm = TRUE)
    yrs  <- yr0 + seq_along(out$S) - 1
    sdf  <- data.frame(year = yrs, pct = (exp(out$S) - 1) * 100)

    pe <- ggplot(sdf, aes(x = year, y = pct, fill = pct > 0)) +
      geom_col(width = 0.7, show.legend = FALSE) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      scale_fill_manual(values = c(`TRUE` = "#3481c8", `FALSE` = "#e84b33")) +
      labs(x = "Year", y = "Relative growth (%)") +
      theme_bw() +
      theme(plot.tag = element_text(face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  } else {
    sigdf <- data.frame(growth = seq(0, 5, 0.1), sigma = out$sigGrowvec)
    pe <- ggplot(sigdf, aes(x = growth, y = sigma)) +
      geom_point(size = 1.5) +
      geom_line() +
      labs(x = "Growth", y = "Growth spread") +
      theme_bw() +
      theme(plot.tag = element_text(face = "bold"))
  }

  ## --- Panel f: projected growth trajectory, with year-scenario overlay if available ---
  # x can contain tiny negative values (e.g. -1e-14) from floating-point rounding
  # in the pnorm(upper) - pnorm(lower) subtraction used to build each STM's tail
  # bins -- clamp before using as rep() counts, or a single such bin among nlbin
  # crashes the whole call ("invalid 'times' argument").
  weighted.probs <- function(x) {
    x <- pmax(x, 0)
    quantile(rep(lbin, round(1e+05 * x)), probs = c(0.025, 0.5, 0.975))
  }
  xs      <- 0:(ntsteps * 30)
  Mxage   <- trunc((3.5 / M) / 5) * 5 + 5
  Mxtstep <- Mxage * ntsteps
  keep    <- xs < Mxtstep

  build_traj <- function(mat) {
    ci <- apply(mat, 1, weighted.probs)
    data.frame(
      ts  = xs[keep],
      med = c(lbin[1], ci[2, ])[keep],
      lo  = c(lbin[1], ci[1, ])[keep],
      hi  = c(lbin[1], ci[3, ])[keep]
    )
  }

  # lenout_scenarios is no longer REPORT()ed directly by growmod (it can't
  # be -- min(S)/max(S) inside the AD-traced function is unsafe). Build it
  # here instead, from the already-numeric mod$rep() output. Returns NULL for
  # a plain growmod fit (no S reported), so this degrades gracefully.
  lenout_scenarios <- build_lenout_envelope(out, datain)
  has_scenarios <- !is.null(lenout_scenarios)

  if (has_scenarios) {
    traj_df <- bind_rows(
      build_traj(lenout_scenarios[, , 1]) %>% mutate(scenario = "Average"),
      build_traj(lenout_scenarios[, , 2]) %>% mutate(scenario = "Slowest year"),
      build_traj(lenout_scenarios[, , 3]) %>% mutate(scenario = "Fastest year")
    )
  } else {
    traj_df <- build_traj(out$lenout) %>% mutate(scenario = "Average")
  }

  age_breaks <- seq(2, (ntsteps * Mxage) + 2, ntsteps * 2)
  age_labels <- seq(2, Mxage + 2, 2)

  pf <- ggplot(traj_df, aes(x = ts)) +
    geom_ribbon(data = subset(traj_df, scenario == "Average"),
                aes(ymin = lo, ymax = hi), fill = "grey70", alpha = 0.8) +
    geom_line(aes(y = med, colour = scenario, linetype = scenario), linewidth = 0.8) +
    scale_colour_manual(values = c(Average = "black", `Slowest year` = "#e84b33",
                                    `Fastest year` = "#3481c8")) +
    scale_linetype_manual(values = c(Average = "solid", `Slowest year` = "dashed",
                                       `Fastest year` = "dashed")) +
    scale_x_continuous(breaks = age_breaks, labels = age_labels) +
    scale_y_continuous(limits = c(0, max(lbinL))) +
    labs(x = "Relative Age (y)", y = "Size",
         colour = "Year scenario", linetype = "Year scenario") +
    theme_bw() +
    theme(plot.tag = element_text(face = "bold"),
          legend.position = if (has_scenarios) "right" else "none")

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
      labs(x = "Longitude", y = "Latitude") +
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
      plot_annotation(tag_levels = "a", theme = theme(plot.margin = margin(2, 2, 2, 2))) &
      theme(legend.position    = "right",
            legend.box         = "vertical",        # stack legends on top of each other
            legend.direction   = "horizontal",      # each legend runs horizontally
            legend.margin      = margin(0, 0, 0, 0),
            legend.box.margin  = margin(0, 0, 0, 0),
            plot.margin        = margin(2, 4, 2, 4))

  } else {
    final <- (pa + pb) / (pc + pd) / (pe + pf) +
      plot_layout(guides = "collect") +
      plot_annotation(tag_levels = "a", theme = theme(plot.margin = margin(2, 2, 2, 2))) &
      theme(legend.position = "bottom",
            legend.margin    = margin(0, 0, 0, 0),
            legend.box.margin = margin(-5, 0, 0, 0),
            plot.margin      = margin(2, 4, 2, 4))
  }

  suppressWarnings(print(final))
  return(list(out = out, obs = obs))
}
