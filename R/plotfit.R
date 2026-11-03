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
#'     year. Otherwise, a calibration (PIT) histogram: for each animal the
#'     fitted probability of a recapture length at or below the observed one.
#'     Flat (within the grey 95\% band) means the predicted distributions have
#'     the right location and spread; a U shape means the predicted spread is
#'     too narrow, a hump too wide, and a tilt a bias.
#'   \item Projected growth trajectory with 95% CI ribbon for the average
#'     year. If \code{mod$rep()$lenout_scenarios} exists, the slowest- and
#'     fastest-fitted-year median trajectories are overlaid as dashed lines.
#'   \item (if spatial columns present) Mean residual by grid cell: each
#'     \code{sp_cell} x \code{sp_cell} cell is filled by its mean residual
#'     (blue = grew more than predicted, red = less), labelled with its
#'     record count, and outlined in black where the mean differs from zero
#'     by more than 2 SE. Cells with fewer than \code{sp_min_n} records are
#'     not drawn.
#' }
#'
#' @param area_col Optional name of a column in \code{datIn} holding the
#'   growth area. Adds each area's overall mean residual (+/- SE) to the
#'   spatial panel's subtitle.
#' @param sp_cell Grid cell size for the spatial panels, in the units of
#'   \code{dLon}/\code{dLat}. Default 1.
#' @param sp_min_n Minimum records for a cell or band to be shown. Default 3.
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
plotfit <- function(datIn = tdat, mout = NULL, run_label = NULL,
                    area_col = NULL, sp_cell = 1, sp_min_n = 3) {
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
    lon    = datIn$dLon,
    area   = if (!is.null(area_col) && area_col %in% names(datIn))
               as.character(datIn[[area_col]]) else "all"
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

  ## Pmoult_vec now genuinely varies by season as well as size (Makemap's
  ## default frees each goodts row of Pmoult_par independently -- see
  ## growmod/Makemap changes), so it must be subset and aligned row-by-row
  ## with growthmat exactly the same way, not broadcast with rep().
  pmoultmat <- out$Pmoult_vec[goodts, , drop = FALSE]

  if (is.matrix(growthmat)) {
    gdf <- do.call(rbind, lapply(1:nrow(growthmat), function(i)
      data.frame(lbin = lbin, increment = growthmat[i, ],
                 Pmoult = pmoultmat[i, ], ts = factor(i))))
  } else {
    gdf <- data.frame(lbin = lbin, increment = growthmat,
                      Pmoult = as.numeric(pmoultmat), ts = factor(1))
  }
  # sd_growth is now proportional (sd = exp(LsigGrow) * increment), and
  # sigGrowvec = exp(LsigGrow) * gseq is exactly linear -- so this
  # interpolation is now an exact readout, not an approximation.
  gdf$sd <- approx(x = seq(0, 5, 0.1), y = out$sigGrowvec, xout = gdf$increment, rule = 2)$y

  # (old `gdf$Pmoult <- rep(out$Pmoult_vec, times = length(unique(gdf$ts)))`
  # line removed -- Pmoult is now built correctly above, aligned per season.)

  # Rescale P(moult) [0,1] onto the ribbon's own [min, max] range, rather than
  # an arbitrary [0, some_constant] -- this makes the dashed P(moult) line(s)
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
    # group = ts added: with season-varying Pmoult, this now draws one dashed
    # line PER SEASON rather than connecting points across seasons in raw
    # x-order (which produced a jagged, meaningless line once Pmoult stopped
    # being constant across ts). colour still fixed at black/dashed so the
    # moult-probability line stays visually distinct from the growth curves
    # above, which already use the viridis palette for colour = ts.
    geom_line(aes(y = ribbon_lo + Pmoult * ribbon_range, group = ts),
              colour = "black", linetype = "dashed", linewidth = 0.7) +
    scale_colour_viridis_d(option = "plasma") +
    scale_fill_viridis_d(option = "plasma") +
    guides(colour = "none", fill = "none") +
    scale_y_continuous(
      name = "Increment",
      sec.axis = sec_axis(~ (. - ribbon_lo) / ribbon_range, name = "P(moult)")
    ) +
    labs(x = "Size", title = "Growth/Moult (dashed = P(moult), one line per season)") +
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
    ## Calibration (PIT) histogram. For each animal, the fitted probability
    ## that its recapture length is at or below the one observed: the
    ## predicted true-length distribution (EstRecLen) convolved with
    ## recapture measurement error. If the predicted distributions are right
    ## -- location AND spread, including the moult/no-moult mixture -- these
    ## values are uniform on (0, 1). This tests what the old growth-spread
    ## panel could not: that panel only redrew sd = CV x increment, a straight
    ## line through the origin whatever the fit.
    pl_now   <- if (!is.null(mout)) mod$env$parList(mout$par) else mod$env$parList()
    sigE     <- exp(pl_now$LsigError)
    erl      <- out$EstRecLen
    pit <- vapply(seq_len(nrow(erl)), function(i)
      sum(erl[i, ] * pnorm(datIn$rccl[i], lbin, sigE)) / sum(erl[i, ]), numeric(1))
    nb   <- 10
    n    <- length(pit)
    cnt  <- as.numeric(table(cut(pit, seq(0, 1, length.out = nb + 1), include.lowest = TRUE)))
    band <- qbinom(c(0.025, 0.975), n, 1 / nb) / (n / nb)
    pitdf <- data.frame(mid = (seq_len(nb) - 0.5) / nb, rel = cnt / (n / nb))

    pe <- ggplot(pitdf, aes(x = mid, y = rel)) +
      annotate("rect", xmin = 0, xmax = 1, ymin = band[1], ymax = band[2],
               fill = "grey85") +
      geom_col(width = 1 / nb, fill = "#3481c8", colour = "white", alpha = 0.8) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      scale_x_continuous(breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
      labs(x = "Predicted P(length <= observed)",
           y = "Frequency / expected",
           subtitle = "Flat = calibrated; U = spread too narrow;\nhump = too wide; tilt = biased") +
      theme_bw() +
      theme(plot.tag = element_text(face = "bold"),
            plot.subtitle = element_text(size = 8))
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
    obs_sp <- obs %>% filter(!is.na(lat) & !is.na(lon))

    ## --- Panel: mean residual by grid cell ----------------------------------
    ## Individual points hide imbalance: positive and negative residuals of
    ## similar size overplot each other. Averaging within cells shows where
    ## the model is systematically high or low, and the 2-SE outline shows
    ## which of those departures are more than noise.
    cells <- obs_sp %>%
      mutate(cx = (floor(lon / sp_cell) + 0.5) * sp_cell,
             cy = (floor(lat / sp_cell) + 0.5) * sp_cell) %>%
      group_by(cx, cy) %>%
      summarise(n = n(), mres = mean(res), se = sd(res) / sqrt(n()), .groups = "drop") %>%
      filter(n >= sp_min_n) %>%
      mutate(sig = abs(mres) > 2 * se)
    lim <- max(abs(cells$mres), na.rm = TRUE)

    ## With growth areas supplied, add each area's overall mean residual
    ## (+/- SE) to the subtitle: one number per area that should sit near 0.
    sub_txt <- "Outlined: |mean| > 2 SE; number = records"
    if (length(unique(obs_sp$area)) > 1) {
      am <- obs_sp %>% group_by(area) %>%
        summarise(m = mean(res), se = sd(res) / sqrt(n()), .groups = "drop")
      sub_txt <- paste0(sub_txt, "\n",
                        paste(sprintf("%s: %+.2f (SE %.2f)", am$area, am$m, am$se),
                              collapse = "\n"))
    }

    pg <- ggplot(cells, aes(cx, cy)) +
      geom_tile(aes(fill = mres), width = sp_cell, height = sp_cell, colour = "grey85") +
      geom_tile(data = subset(cells, sig), width = sp_cell, height = sp_cell,
                fill = NA, colour = "black", linewidth = 0.6) +
      geom_text(aes(label = n), size = 2.3, colour = "grey20") +
      scale_fill_gradient2(low = "#e84b33", mid = "white", high = "#3481c8",
                           midpoint = 0, limits = c(-lim, lim),
                           name = "Mean residual (mm)") +
      ## No coord_equal(): a fixed aspect ratio makes patchwork shrink the
      ## whole right-hand column to fit this panel.
      labs(x = "Longitude", y = "Latitude",
           subtitle = sub_txt) +
      theme_bw() +
      theme(plot.tag = element_text(face = "bold"),
            plot.subtitle = element_text(size = 8))

    layout <- "
AB
CD
EF
GG
"
    final <- pa + pb + pc + pg + pd + pe + pf +
      plot_layout(design = layout, guides = "collect", widths = c(1, 1)) +
      plot_annotation(tag_levels = "a", caption = run_label,
                      theme = theme(plot.margin = margin(2, 2, 2, 2),
                                    plot.caption = element_text(hjust = 0.5, size = 9, face = "italic"))) &
      theme(legend.position    = "bottom",
            legend.box         = "horizontal",
            legend.direction   = "horizontal",
            legend.margin      = margin(0, 0, 0, 0),
            legend.box.margin  = margin(0, 0, 0, 0),
            plot.margin        = margin(2, 4, 2, 4))

  } else {
    final <- (pa + pb) / (pc + pd) / (pe + pf) +
      plot_layout(guides = "collect") +
      plot_annotation(tag_levels = "a", caption = run_label,
                      theme = theme(plot.margin = margin(2, 2, 2, 2),
                                    plot.caption = element_text(hjust = 0.5, size = 9, face = "italic"))) &
      theme(legend.position = "bottom",
            legend.margin    = margin(0, 0, 0, 0),
            legend.box.margin = margin(-5, 0, 0, 0),
            plot.margin      = margin(2, 4, 2, 4))
  }

  suppressWarnings(print(final))
  return(list(out = out, obs = obs))
}
