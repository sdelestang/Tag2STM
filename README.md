# Tag2STM - Tag-recapture data to Size Transition Matrices

A package to convert tag-recapture data into size transition matrices for use in crustacean stock assessment models like IMuLT.

## Installation

Tag2STM requires several packages. Install directly from GitHub using devtools:
```r
# Install devtools if you don't have it
install.packages("devtools")

# Install CRAN dependencies
install.packages("RTMB")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("sp")

# Install GitHub dependency
devtools::install_github("kassambara/ggpubr")

# Install Tag2STM
devtools::install_github("sdelestang/Tag2STM")
```

## Required Dependencies

The following packages will be automatically installed if not present:
- RTMB (for model fitting)
- dplyr (for data manipulation)
- ggplot2 (for plotting)
- ggpubr (for publication-ready plots)
- sp (for spatial data handling)

## Usage
```r
library(Tag2STM)
library(RTMB)

# 1. Set up length bins
bins <- MakeLbin(start = 20, stop = 200, gap = 2)
lbin <- bins$lbin
lbinL <- bins$lbinL
lbinU <- bins$lbinU

# 2. Define time step structure (e.g., seasonal growth periods)
tstep_def <- data.frame(
  month = c(11, 2, 5, 8),  # Nov, Feb, May, Aug
  hm = c(0, 0, 0, 0)        # First half of each month
)
times <- Maketimes(tstep_def)

# 3. Prepare your tag-recapture data
# Data should have: release length, recapture length, release/recapture dates
# Filter for species/sex/location as needed before processing

# 4. Build data structure for model
datain <- list(
  lbin = lbin,
  lbinL = lbinL, 
  lbinU = lbinU,
  ntsteps = max(times$tstep),
  goodts = unique(times$tstep),
  # ... add your tag data here
)

# 5. Set initial parameters
pin <- list(
  mxpin = rep(0, length(unique(times$tstep))),
  mnpin = rep(-5, length(unique(times$tstep))),
  ipin = rep(80, length(unique(times$tstep))),
  spin_left = rep(0, length(unique(times$tstep))),
  spin_right = rep(0, length(unique(times$tstep))),
  LsigError = log(2),
  sigGrow = log(5),
  # ... add measurement error terms if needed
)

# 6. Create parameter map
map <- Mapfunc(pin, re = FALSE)

# 7. Fit the model
mod <- MakeADFun(growmod, pin, map = map)
opt <- nlminb(mod$par, mod$fn, mod$gr)

# 8. View results
mout <- opt
plotfit()      # Diagnostic plots
plotpars()     # Growth curves by time step

# 9. Extract size transition matrices for use in IMuLT
out <- mod$rep()
stm <- out$stm  # 3D array: [from_bin, to_bin, timestep]
```

## Key Functions

- `MakeLbin()` - Create length bin structure
- `Maketimes()` - Map months to time steps
- `cntTsteps()` - Calculate time at liberty
- `growmod()` - RTMB objective function for growth model
- `Mapfunc()` - Create parameter mapping for estimation control
- `plotfit()` - Diagnostic plots for model fit
- `plotpars()` - Plot estimated growth parameters

## Model Features

- Asymmetric inverse logistic growth curves with smooth steepness transition
- Time-varying growth parameters by season
- Measurement error on both release and recapture lengths
- Flexible time step structure (monthly, seasonal, etc.)
- Output ready for integration with IMuLT population models

## Author

Simon de Lestang (DPIRD, Western Australia)

## See Also

- [IMuLT](https://github.com/sdelestang/IMuLT) - Integrated Model using Length Transition for crustacean stock assessment
