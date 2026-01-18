#' Extract Subset of Size Transition Matrix for Stock Assessment
#'
#' Extracts and saves a smaller size transition matrix (STM) from a fitted growth
#' model for use in stock assessment models that may operate on a different or
#' smaller length bin structure. The function subsets the full STM to specified
#' length bins, renormalizes transition probabilities, and exports results as
#' CSV files suitable for stock assessment software.
#'
#' @param LowLB Numeric. Lower bound of the smallest length bin to include in
#'   the subset STM (in mm). Should match a length bin lower bound from the
#'   original growth model. Default is 41 mm.
#' @param UpLB Numeric. Upper bound of the largest length bin to include (in mm).
#'   Should be less than or equal to the maximum length in the growth model.
#'   Default is 151 mm.
#' @param Gap Numeric. Width of length bins (in mm) for the subset STM. This
#'   defines the resolution of the output matrix and should typically match the
#'   bin width used in the stock assessment model. Default is 2 mm.
#'
#' @return NULL (invisibly). The function is called for its side effect of writing
#'   CSV files to the current working directory.
#'
#' @details
#' **Purpose:**
#'
#' Growth models are often fitted across a wide size range to capture all available
#' tag-recapture data, but stock assessment models may focus on a narrower size
#' range relevant to the fishery (e.g., legal-size animals only). This function
#' bridges the gap by extracting appropriately-sized STMs from the growth model
#' for direct use in stock assessment.
#'
#' **Algorithm:**
#'
#' For each time step in \code{goodts}:
#' \enumerate{
#'   \item Define target length bins: \code{seq(LowLB, UpLB, Gap)}
#'   \item Extract the full STM for that time step from the fitted model
#'   \item Find which rows/columns of the full STM correspond to target bins
#'   \item Subset the STM to create a smaller matrix
#'   \item Renormalize each column so probabilities sum to 1 (critical step!)
#'   \item Set very small probabilities (<1e-7) to exactly 0
#'   \item Save as CSV with informative filename
#' }
#'
#' **Renormalization:**
#'
#' When subsetting an STM, some probability mass may be lost (e.g., transitions
#' to length bins outside the new range). The renormalization step ensures each
#' column still represents a valid probability distribution that sums to 1. This
#' is essential for stock assessment models that assume STM columns are proper
#' probability distributions.
#'
#' **File Naming Convention:**
#'
#' Output files are named as:
#' \itemize{
#'   \item If \code{l} exists in environment: \code{STM_s[Sex]_L[l]_ts[timestep].csv}
#'   \item Otherwise: \code{STM_s[Sex]_ts[timestep].csv}
#' }
#'
#' Where:
#' \itemize{
#'   \item \code{[Sex]}: 1 = female, 2 = male (from \code{tdat$Lsex})
#'   \item \code{[l]}: Optional locality/region identifier if defined
#'   \item \code{[timestep]}: Time step number from \code{goodts}
#' }
#'
#' **Required Objects in Environment:**
#'
#' The function expects the following objects to exist:
#' \itemize{
#'   \item \code{mod}: Fitted RTMB/TMB model object from \code{MakeADFun()}
#'     containing the growth model results
#'   \item \code{lbinL}: Vector of length bin lower bounds from the original
#'     growth model (from \code{datain})
#'   \item \code{tdat}: Tag-recapture data frame containing \code{Lsex} column
#'     with sex information ('F' for female, 'M' for male)
#'   \item \code{goodts}: Vector of time steps with estimated growth (from \code{datain})
#'   \item \code{l}: (Optional) Locality or region identifier for file naming
#' }
#'
#' @section Important Notes:
#' \itemize{
#'   \item The specified length bins (\code{LowLB}, \code{UpLB}, \code{Gap})
#'     must align with bins in the original growth model, otherwise matching will fail
#'   \item Very small probabilities (<1e-7) are set to 0 to prevent numerical
#'     issues in stock assessment models
#'   \item Files are saved to the current working directory - use \code{setwd()}
#'     to control output location
#'   \item All animals in \code{tdat} must have the same sex (function uses
#'     \code{unique(tdat$Lsex)})
#'   \item The function contains a typo in the code: \code{LobLB} should be
#'     \code{LowLB} in the first line
#' }
#'
#' @section Stock Assessment Integration:
#' The exported CSV files can be directly imported into length-based stock
#' assessment models such as:
#' \itemize{
#'   \item Integrated Size-Structured Assessment (ISSA)
#'   \item Catch-at-Length Analysis (CALA)
#'   \item Custom length-based population models
#' }
#'
#' Each CSV contains an \code{nlbin × nlbin} matrix where element (i,j) represents
#' the probability of an animal in length bin j transitioning to length bin i
#' during that time step.
#'
#' @examples
#' \dontrun{
#' # Typical workflow after fitting growth model:
#'
#' # 1. Fit the growth model
#' datain <- PrepareTagData(tagdata)
#' pin <- Makepin()
#' map <- Mapfunc(pin, re = FALSE)
#' mod <- MakeADFun(growmod, pin, map = map)
#' opt <- nlminb(mod$par, mod$fn, mod$gr)
#'
#' # 2. Extract STMs for stock assessment
#' # Stock assessment uses 2mm bins from 60-140mm
#' setwd("output/stms")
#' ClipSTM(LowLB = 60, UpLB = 140, Gap = 2)
#'
#' # Output files: STM_s1_ts1.csv, STM_s1_ts2.csv, etc. (if female)
#'
#' # 3. Example with locality identifier
#' l <- "ZoneA"  # Define locality
#' ClipSTM(LowLB = 60, UpLB = 140, Gap = 2)
#' # Output files: STM_s1_LZoneA_ts1.csv, etc.
#'
#' # 4. Verify output
#' stm_test <- read.csv("STM_s1_ts1.csv")
#' colSums(stm_test)  # Should all be 1.0 (or very close)
#'
#' # 5. Use different resolution for sensitivity analysis
#' # Stock assessment with 5mm bins
#' ClipSTM(LowLB = 60, UpLB = 140, Gap = 5)
#'
#' # 6. Extract for legal-size range only
#' # If minimum legal size is 76mm
#' ClipSTM(LowLB = 76, UpLB = 140, Gap = 2)
#' }
#'
#' @section Common Issues:
#' \describe{
#'   \item{Mismatched bins}{If \code{match(mlbinL, lbinL)} returns NAs, the
#'     specified length bins don't exist in the growth model. Adjust \code{LowLB},
#'     \code{UpLB}, or \code{Gap} to match the original bin structure.}
#'   \item{Column sums ≠ 1}{After subsetting, column sums may be slightly less
#'     than 1 due to truncation. The renormalization step fixes this.}
#'   \item{Multiple sexes in data}{The function assumes \code{unique(tdat$Lsex)}
#'     returns a single value. Filter \code{tdat} by sex before calling if needed.}
#' }
#'
#' @seealso
#' \code{\link{growmod}} for fitting the growth model
#' \code{\link{Makepin}} for initial parameters
#' \code{\link{Mapfunc}} for parameter mapping
#'
#' @export
ClipSTM <- function(LowLB = 41, UpLB = 151, Gap = 2) {
  # Create sequence of length bin lower bounds for subset STM
  mlbinL <- seq(LowLB, UpLB, Gap)  # Note: Code has typo "LobLB"

  # Process each time step with estimated growth
  for (tt in goodts) {
    # Extract STM for this time step from fitted model
    stm <- mod$report()$stm[, , tt]

    # Find which rows/columns correspond to desired length bins
    tokeep <- match(mlbinL, lbinL)

    # Subset to create smaller STM
    stm2 <- stm[tokeep, tokeep]

    # Renormalize columns to sum to 1 (critical for stock assessment)
    funcsum <- function(x) { x / sum(x, na.rm = TRUE) }
    stm2 <- apply(stm2, 2, funcsum)

    # Set very small probabilities to exactly zero
    stm2[stm2 < 1e-7] <- 0

    # Determine sex code (1 = female, 2 = male)
    Sex <- ifelse(unique(tdat$Lsex) == 'F', 1, 2)

    # Construct filename
    if (exists('l')) {
      Fname <- paste0('STM_s', Sex, '_L', l, '_ts', tt, '.csv')
    } else {
      Fname <- paste0('STM_s', Sex, '_ts', tt, '.csv')
    }

    print(paste("Saving:", Fname))
    write.csv(stm2, Fname, row.names = FALSE)
  }

  invisible(NULL)
}
