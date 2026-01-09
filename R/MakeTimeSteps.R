#' Create Time Step Mapping Template
#'
#' Generates a complete time step lookup table for tag-recapture data,
#' mapping month and half-month periods to sequential time steps. Fills
#' missing time steps circularly and calculates fractional allocation
#' when multiple periods map to the same time step.
#'
#' @param x Data frame with columns defining time step structure:
#' \describe{
#'   \item{Column 1}{Month (1-12)}
#'   \item{Column 2}{Half-month indicator (0 = first half, 1 = second half)}
#'   \item{Column 3 (optional)}{Time step number (created if missing)}
#' }
#'
#' @return A data frame with 24 rows (12 months × 2 half-months) and columns:
#' \describe{
#'   \item{month}{Month number (1-12)}
#'   \item{hm}{Half-month indicator (0 or 1)}
#'   \item{tstep}{Sequential time step number}
#'   \item{tfrac}{Fractional allocation (1/n when n periods share a time step)}
#'   \item{id}{Character identifier "month hm"}
#' }
#'
#' @details
#' This function creates a complete 24-row template covering all month/half-month
#' combinations, even if some are missing from the input. Missing time steps are
#' filled forward, with circular wrapping (time steps before the first defined
#' one are filled from the last defined time step).
#'
#' The `tfrac` column accounts for cases where multiple half-month periods map
#' to the same time step, useful for allocating catch or effort data.
#'
#' @examples
#' \dontrun{
#' # Define seasonal time steps (e.g., 4 quarters)
#' tstep_def <- data.frame(
#'   month = c(11, 2, 5, 8),
#'   hm = c(0, 0, 0, 0)
#' )
#' time_map <- Maketimes(tstep_def)
#' }
#'
#' @importFrom dplyr mutate arrange group_by ungroup if_else row_number last n
#' @importFrom tidyr fill
#' @importFrom magrittr %>% %<>%
#' @export
MakeTsteps <- function(x) {
  x$tstep <- 1:nrow(x)

  # template
  tmp <- expand.grid(month = 1:12, hm = 0:1) %>%
    arrange(month, hm) %>%
    mutate(tstep = NA_integer_)

  # map tstep from x onto template
  tmp$tstep <- x[match(
    paste(tmp$month, tmp$hm),
    paste(x[[1]], x[[2]])
  ), 3]

  # ---- fill missing tstep, circular ----
  tmp <- tmp %>%
    mutate(
      tstep = if_else(
        is.na(tstep) & row_number() < which(!is.na(tstep))[1],
        last(na.omit(tstep)),
        tstep
      )
    ) %>%
    fill(tstep)

  # ---- fractional allocation per timestep ----
  tmp <- tmp %>%
    group_by(tstep) %>%
    mutate(tfrac = 1 / n()) %>%
    ungroup()

  tmp %<>% mutate(id = paste(month, hm))
  return(tmp)
}
