#' Extract and Save Subset Size Transition Matrices
#'
#' Extracts a size-restricted subset of the fitted model's size transition
#' matrix (or matrices), renormalizes columns to sum to 1, and writes CSV
#' files suitable for downstream stock assessment use.
#'
#' @param LowLB,UpLB,Gap Numeric. Lower bound, upper bound, and step size (mm)
#'   for the subset length bins to keep.
#' @param Annual Logical. If \code{TRUE}, combines all \code{goodts} seasons
#'   into a single within-year-cycle STM via sequential matrix multiplication
#'   (as before), plus a mean-length-at-age diagnostic plot. If \code{FALSE}
#'   (default), writes one STM per \code{goodts} season, uncombined.
#' @param AvgYears Logical. Only relevant when the fitted model used
#'   \code{TemporalGrowth = TRUE} (i.e. \code{mod$report()$stm} is 4D:
#'   \code{nlbin x nlbin x ntsteps x nyears}). If \code{TRUE} (default),
#'   averages the STM across \strong{supported} years
#'   (\code{datain$yr_supported}) for each \code{goodts}, giving the same
#'   single-STM-per-season output shape as a model with no year effects. If
#'   \code{FALSE}, produces one STM per \code{(goodts, year)} pair instead of
#'   averaging, with the year appended to each output filename. Ignored
#'   entirely (no effect) when \code{stm} is still the older 3D shape (i.e.
#'   the fit used \code{TemporalGrowth = FALSE}).
#' @param return Logical. If \code{TRUE}, returns the (first) subset STM
#'   instead of writing it to a file. Only meaningful for a single result --
#'   with \code{AvgYears = FALSE} producing several STMs, only the first is
#'   returned; the rest are still written to disk regardless.
#'
#' @details
#' Unsupported years (\code{!datain$yr_supported}) are excluded from both the
#' averaging and the per-year output entirely -- they are hard-fixed at
#' average growth (\code{S = 0}) inside \code{growmod}, not independently
#' estimated, so including them would either dilute a genuine average with
#' duplicate copies of it, or misleadingly present a non-informative year as
#' if it carried real interannual signal.
#'
#' @export
ClipSTM <- function(LowLB = 41, UpLB = 151, Gap = 2, Annual = FALSE,
                     AvgYears = TRUE, return = FALSE) {
  lbinL  <- bins$lbinL
  mlbinL <- seq(LowLB, UpLB, Gap)
  tokeep <- match(mlbinL, lbinL)
  funcsum <- function(x) x / sum(x, na.rm = TRUE)

  stm_full  <- mod$report()$stm
  has_years <- length(dim(stm_full)) == 4   # TemporalGrowth = TRUE -> 4D array
  yr_idx    <- if (has_years) which(datain$yr_supported) else NULL

  finalize_stm <- function(stm) {
    stm2 <- stm[tokeep, tokeep]
    stm2 <- apply(stm2, 2, funcsum)
    stm2[stm2 < 1e-7] <- 0
    stm2
  }

  ## ---- Annual = TRUE: combine goodts seasons into one within-year STM ----
  if (Annual) {
    combine_seasons <- function(stm_seasons) {
      # stm_seasons: nlbin x nlbin x length(goodts), for ONE year (or the
      # year-averaged case). Sequentially multiply seasons together, exactly
      # as the original fixed dim3==2/3/4 cases did, but for any length.
      dim3 <- dim(stm_seasons)[3]
      out <- stm_seasons[, , 1]
      if (dim3 >= 2) for (i in 2:dim3) out <- stm_seasons[, , i] %*% out
      out
    }

    if (!has_years) {
      annual_list <- list(avg = combine_seasons(stm_full[, , goodts]))
    } else if (AvgYears) {
      # Average each season across supported years FIRST, then combine
      # seasons -- gives one averaged annual STM.
      avg_seasons <- array(0, c(dim(stm_full)[1], dim(stm_full)[2], length(goodts)))
      for (i in seq_along(goodts)) {
        avg_seasons[, , i] <- apply(stm_full[, , goodts[i], yr_idx, drop = FALSE], c(1, 2), mean)
      }
      annual_list <- list(avg = combine_seasons(avg_seasons))
    } else {
      # One combined annual STM PER supported year, no averaging.
      annual_list <- setNames(
        lapply(yr_idx, function(yr) combine_seasons(stm_full[, , goodts, yr])),
        paste0("yr", yr_idx)
      )
    }

    for (nm in names(annual_list)) {
      stm <- annual_list[[nm]]

      # Diagnostic mean-length-at-age plot -- averaged case only (one plot
      # per supported year would be excessive by default; extend this loop
      # if per-year diagnostic plots are ever wanted).
      if (nm == "avg") {
        lenout <- matrix(0, ncol = ncol(stm), nrow = 30)
        lenout[1, 1] <- 1
        for (y in 2:30) lenout[y, ] <- stm %*% lenout[y - 1, ]
        templen <- apply(lenout, 1, function(x) {
          wm  <- weighted.mean(bins$lbin, x)
          wv  <- sum(x * (bins$lbin - wm)^2) / sum(x)
          wsd <- sqrt(wv)
          n   <- sum(x > 0)
          se  <- wsd / sqrt(n)
          c(mean = wm, lo95 = wm - 1.96 * se, hi95 = wm + 1.96 * se)
        }) %>% t() %>% as.data.frame()
        templen$step <- seq_len(nrow(templen))

        print(
          ggplot(templen, aes(x = step)) +
            geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "#378ADD", alpha = 0.15) +
            geom_line(aes(y = mean), colour = "#185FA5", linewidth = 0.8) +
            geom_point(aes(y = mean), colour = "#185FA5", size = 1.8) +
            labs(x = "Annual Time step", y = "Mean length (mm)",
                 title = "Mean length-at-age with 95% CI") +
            theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  plot.title = element_text(size = 13, face = "plain"))
        )
      }

      stm2 <- finalize_stm(stm)

      Sex <- ifelse(unique(tdat$Lsex) == 'F', 'Fem', 'Male')
      if (exists('a') & !exists('l')) { l <- a }

      yr_suffix <- if (nm == "avg") "" else paste0("_", nm)

      if (exists('p') & exists('l') & exists('Sex')) {
        Fname <- paste0('STM_', Sex, '_L', l, '_p', p, '_Annual', yr_suffix, '.csv') }
      if (exists('p') & !exists('l') & !exists('Sex')) {
        Fname <- paste0('STM_p', p, '_Annual', yr_suffix, '.csv') }
      if (!exists('p') & exists('l') & !exists('Sex')) {
        Fname <- paste0('STM_L', l, '_Annual', yr_suffix, '.csv') }
      if (!exists('p') & !exists('l') & exists('Sex')) {
        Fname <- paste0('STM_', Sex, '_Annual', yr_suffix, '.csv') }
      if (exists('p') & !exists('l') & exists('Sex')) {
        Fname <- paste0('STM_', Sex, '_p', p, '_Annual', yr_suffix, '.csv') }
      if (!exists('p') & exists('l') & exists('Sex')) {
        Fname <- paste0('STM_', Sex, '_L', l, '_Annual', yr_suffix, '.csv') }

      if (return == FALSE) {
        print(paste("Saving:", Fname))
        write.csv(stm2, Fname, row.names = FALSE)
      }
      if (return == TRUE) return(stm2)
    }
  }

  ## ---- Annual = FALSE: one STM per goodts season, uncombined ----
  if (!Annual) {
    for (tt in goodts) {
      if (!has_years) {
        stm_list <- list(avg = stm_full[, , tt])
      } else if (AvgYears) {
        stm_list <- list(avg = apply(stm_full[, , tt, yr_idx, drop = FALSE], c(1, 2), mean))
      } else {
        stm_list <- setNames(lapply(yr_idx, function(yr) stm_full[, , tt, yr]),
                              paste0("yr", yr_idx))
      }

      for (nm in names(stm_list)) {
        stm2 <- finalize_stm(stm_list[[nm]])

        Sex <- ifelse(unique(tdat$Lsex) == 'F', 1, 2)
        if (exists('a') & !exists('l')) { l <- a }

        yr_suffix <- if (nm == "avg") "" else paste0("_", nm)

        if (exists('p') & exists('l') & exists('Sex')) {
          Fname <- paste0('STM_s', Sex, '_L', l, '_ts', tt, '_p', p, yr_suffix, '.csv') }
        if (exists('p') & !exists('l') & !exists('Sex')) {
          Fname <- paste0('STM_p', p, '_ts', tt, yr_suffix, '.csv') }
        if (!exists('p') & exists('l') & !exists('Sex')) {
          Fname <- paste0('STM_L', l, '_ts', tt, yr_suffix, '.csv') }
        if (!exists('p') & !exists('l') & exists('Sex')) {
          Fname <- paste0('STM_s', Sex, '_ts', tt, yr_suffix, '.csv') }
        if (exists('p') & !exists('l') & exists('Sex')) {
          Fname <- paste0('STM_s', Sex, '_ts', tt, '_p', p, yr_suffix, '.csv') }
        if (!exists('p') & exists('l') & exists('Sex')) {
          Fname <- paste0('STM_s', Sex, '_L', l, '_ts', tt, yr_suffix, '.csv') }

        print(paste("Saving:", Fname))
        write.csv(stm2, Fname, row.names = FALSE)
      }
    }
  }
  invisible(NULL)
}
