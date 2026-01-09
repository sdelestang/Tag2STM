# Tag2STM

R package for analyzing tag-recapture data using size transition matrices (STMs) with asymmetric inverse logistic growth curves.

## Installation
```r
# Install from GitHub
devtools::install_github("your-username/Tag2STM")
```

## Usage
```r
library(RTMB)
library(Tag2STM)

# Prepare your tag-recapture data
datain <- list(
  Tlbin = your_data$tag_length_bin,
  Rccl = your_data$recapture_length,
  Rlcl = your_data$release_length,
  relts = your_data$release_timestep,
  tsteps = your_data$time_at_liberty,
  nlob = your_data$number_of_lobsters,
  nlbin = 72,                    # Number of length bins
  ntsteps = 6,                   # Number of time steps
  lbinL = bins$lbinL,            # Length bin lower bounds
  lbinU = bins$lbinU,            # Length bin upper bounds
  lbin = bins$lbin,              # Length bin midpoints
  goodts = 1                     # Time steps to model
)

# Set up parameters
pin <- list(
  mxpin = rep(log(40), ntsteps),           # Max growth by timestep
  mnpin = rep(log(5), ntsteps),            # Min growth by timestep
  ipin = rep(80, ntsteps),                 # Inflection point
  spin_left = rep(log(5), ntsteps),        # Left steepness
  spin_right = rep(log(5), ntsteps),       # Right steepness
  LsigError = log(2),                      # Measurement error
  sigGrow = log(5),                        # Growth variability
  MerrorRel = rep(0, nrow(data)),          # Release error (RE)
  MerrorRec = rep(0, nrow(data)),          # Recapture error (RE)
  LMerrorRelsigma = log(2),                # RE sigma for release
  LMerrorRecsigma = log(2)                 # RE sigma for recapture
)

# Optional: fix some parameters
map <- list(
  mxpin = factor(c(1, NA, NA, NA, NA, NA)),      # Only estimate first
  mnpin = factor(c(1, NA, NA, NA, NA, NA)),
  ipin = factor(c(1, NA, NA, NA, NA, NA)),
  spin_left = factor(c(1, NA, NA, NA, NA, NA)),
  spin_right = factor(c(1, NA, NA, NA, NA, NA))
)

# Create model
mod <- make_growmod_obj(
  pin = pin, 
  map = map, 
  random = c("MerrorRel", "MerrorRec")
)

# Optimize
opt <- nlminb(
  start = mod$par,
  objective = mod$fn,
  gradient = mod$gr,
  control = list(eval.max = 2000, iter.max = 1000)
)

# Check convergence
opt$convergence  # 0 = success

# Get results
sdr <- sdreport(mod)
summary(sdr, "fixed")
```

## Growth Model

The package uses an asymmetric inverse logistic growth model with smoothly transitioning steepness:

$$growth = \frac{\max - \min}{1 + \exp\left(\frac{L - inflection}{steepness}\right)} + \min$$

where steepness transitions from `spin_left` to `spin_right` across the inflection point.

## Key Functions

- `growmod()` - Core RTMB objective function for tag-recapture growth model
- `make_growmod_obj()` - Helper to create RTMB objective with proper environment handling

## Requirements

- R >= 4.0.0
- RTMB
- dplyr
- magrittr

## Citation

If you use this package, please cite...

## License

[Your license]
