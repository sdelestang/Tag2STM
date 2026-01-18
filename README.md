# Tag2STM

An R package for estimating size transition matrices (STMs) from tag-recapture data for use in length-based stock assessment models. Implements flexible, data-driven growth models using RTMB/TMB with automatic differentiation.

## Overview

`Tag2STM` fits tag-recapture growth models that:
- Estimate growth patterns that vary by both size and season/time
- Account for measurement error at release and recapture
- Handle individual-level variation in measurement precision
- Produce size transition matrices directly usable in stock assessment models

The package is specifically designed for crustacean stock assessment (e.g., Western Rock Lobster) but is applicable to any species with tag-recapture data.

## Installation
```r
# Install from GitHub
devtools::install_github("sdelestang/Tag2STM")

# Load package
library(Tag2STM)
```

## Key Features

- **Flexible Growth Modeling**: Random walk structure allows growth to vary smoothly across length bins and time steps
- **Measurement Error**: Accounts for both baseline measurement error and individual-level random effects
- **Model Selection**: Tools for comparing models with different seasonal growth patterns
- **Direct Stock Assessment Integration**: Exports STMs in formats ready for length-based models
- **Efficient Computation**: Uses RTMB/TMB automatic differentiation for fast optimization

## Typical Workflow

### 1. Data Preparation
```r
library(Tag2STM)
library(tidyr)
library(dplyr)

# Load tag-recapture data
obs <- read.csv("Tag.data.csv")

# Set up length bins (31-173 mm in 2mm increments)
bins <- MakeLbin(31, 173, 2)

# Clean data: remove missing values, negative liberties
obs <- obs %>% 
  filter(!is.na(Ldate), !is.na(sex)) %>% 
  mutate(Ldays = as.numeric(as.Date(date) - as.Date(Ldate))) %>%  
  mutate(Lmnth = floor(Ldays/30)) %>% 
  filter(Ldays > 0 & floor(Ccl) %in% 40:200 & floor(LCl) %in% 40:200) %>%
  filter(!is.na(Ccl)) %>%  
  mutate(growth = Ccl - LCl) %>% 
  filter(sex %in% c('F','M')) %>%  
  filter(toupper(Lsex) == toupper(sex)) %>% 
  filter(Lmnth >= 0)

# Remove duplicates
obs <- obs %>% 
  mutate(id = paste(tag, Ldate, dLon, LdLon, Ccl, LCl)) %>% 
  filter(!duplicated(id)) %>% 
  select(-id)
```

### 2. Assign Locations and Time Steps
```r
# Assign spatial locations (example: north/south of 30°)
obs <- obs %>% 
  mutate(
    rcloc = ifelse(dLat < 30, 2, 1),
    rlloc = ifelse(LdLat < 30, 2, 1),
    sex = ifelse(toupper(sex) == 'M', 2, 1)
  )

# Define time steps (e.g., quarterly)
timesteps <- data.frame(month = c(1, 4, 7, 10), halfmonth = c(0, 0, 0, 0))
ntsteps <- nrow(timesteps)
times <- MakeTsteps(timesteps)

# Add time step information to data
obs <- obs %>% 
  mutate(
    relmn = as.month(Ldate),
    relhm = ifelse(as.day(Ldate) > 14, 1, 0),
    recmn = as.month(date),
    rechm = ifelse(as.day(date) > 14, 1, 0),
    relts = times$tstep[match(paste(relmn, relhm), paste(times$month, times$hm))],
    rects = times$tstep[match(paste(recmn, rechm), paste(times$month, times$hm))],
    relyr = as.year(Ldate),
    recyr = as.year(date),
    rccl = Ccl,
    rlcl = LCl
  )

# Calculate time steps at liberty
obs$ntstep <- apply(
  as.matrix(obs[, c('relyr', 'recyr', 'relmn', 'relhm', 'recmn', 'rechm')]), 
  1, 
  cntTsteps
)
```

### 3. Prepare Modeling Data
```r
# Filter to specific sex and location, aggregate by release cohort
s <- 'F'  # Female
l <- 2    # Location

tdat <- obs %>% 
  filter(Lsex == s & rlloc == l) %>% 
  group_by(Lsex, rlloc, rccl, rlcl, relts, ntstep) %>% 
  summarise(num = length(Lsex))

# Set natural mortality (for plotting)
M <- 0.14
```

### 4. Model Selection: Find Best Time Steps

Test all combinations of time steps to find which periods have sufficient data for growth estimation:
```r
# Generate all possible combinations of time steps
vec <- 1:ntsteps
totest <- unlist(lapply(1:length(vec), function(i) {
  combn(vec, i, FUN = function(x) paste(x, collapse = ","), simplify = TRUE)
}))

# Create results dataframe
out <- data.frame(id = totest, nLL = NA, AIC = NA, BIC = NA, converge = NA)

# Optionally limit to combinations with <= 3 timesteps
out <- out[apply(as.matrix(out[, 1]), 1, nchar) <= 4, ]

# Fit each model combination
for (i in 1:nrow(out)) {
  goodts <- as.numeric(strsplit(out$id[i], ',')[[1]])
  print(paste("Testing time steps:", paste(goodts, collapse = ", ")))
  
  # Create initial parameters
  pin <- Makepin()
  
  # Prepare data input
  datain <- list(
    Rccl = tdat$rccl,
    Rlcl = tdat$rlcl,
    relts = tdat$relts,
    tsteps = as.integer(tdat$ntstep),
    nlob = tdat$num,
    nlbin = length(lbin),
    ntsteps = ntsteps,
    lbinL = lbinL,
    lbinU = lbinU,
    lbin = lbin,
    goodts = goodts,
    M = M,
    smoother = 5
  )
  
  # Fit base model (no random effects)
  map <- Makemap(pin, re = FALSE)
  mod <- make_growmod_obj(pin = pin, map = map)
  mout <- nlminb(mod$par, mod$fn, mod$gr, 
                 control = list(eval.max = 2000, iter.max = 2000))
  
  # Calculate fit statistics
  dout <- plotfit()
  logLike <- sum(dout$LL)
  k <- length(mout$par)
  n <- nrow(tdat)
  
  out$nLL[i] <- logLike
  out$AIC[i] <- 2 * k - (2 * logLike)
  out$BIC[i] <- log(n) * k - (2 * logLike)
  out$converge[i] <- mout$convergence
}

# View results
out[order(out$AIC), ]

# Select best model (lowest AIC/BIC)
best_model <- out[which.min(out$AIC), ]
goodts <- as.numeric(strsplit(best_model$id, ',')[[1]])
print(paste("Best time steps:", paste(goodts, collapse = ", ")))
```

### 5. Fit Final Model with Random Effects
```r
# Prepare data with selected time steps
datain <- list(
  Rccl = tdat$rccl,
  Rlcl = tdat$rlcl,
  relts = tdat$relts,
  tsteps = as.integer(tdat$ntstep),
  nlob = tdat$num,
  nlbin = length(lbin),
  ntsteps = ntsteps,
  lbinL = lbinL,
  lbinU = lbinU,
  lbin = lbin,
  goodts = goodts,
  M = M,
  smoother = 5
)

# Create initial parameters
pin <- Makepin(LMerrorRelsigma = log(0.5), LMerrorRecsigma = log(0.5))

# Fit model with individual random effects
map <- Makemap(pin, re = TRUE)
mod <- make_growmod_obj(
  pin = pin, 
  map = map, 
  random = c("MerrorRel", "MerrorRec")
)

mout <- nlminb(mod$par, mod$fn, mod$gr, 
               control = list(eval.max = 1000, iter.max = 1000))

# Check convergence and examine fit
print(mout)
dout <- plotfit()
```

### 6. Extract and Save STMs for Stock Assessment
```r
# Extract subset STM matching stock assessment length bins
# Example: 41-151mm in 2mm bins for IMuLT model
ClipSTM(LowLB = 41, UpLB = 151, Gap = 2)

# Files saved as: STM_s1_ts1.csv, STM_s1_ts2.csv, etc.
# where s1 = female, s2 = male
```

### 7. Verify and Use STMs
```r
# Verify that STM columns sum to 1
stm_check <- read.csv("STM_s1_ts1.csv")
colSums(stm_check)  # Should all be ~1.0

# Import into stock assessment model (e.g., IMuLT, ISSA, custom model)
```

## Key Functions

### Data Preparation
- `MakeLbin()`: Create length bin structure
- `MakeTsteps()`: Define time step structure
- `cntTsteps()`: Calculate time at liberty in time steps

### Model Fitting
- `Makepin()`: Create initial parameter list
- `Makemap()`: Control which parameters are estimated
- `growmod()`: Main growth model function (RTMB/TMB objective)
- `make_growmod_obj()`: Wrapper to create model object

### Output Processing
- `plotfit()`: Visualize model fit and diagnostics
- `ClipSTM()`: Extract subset STM for stock assessment

## Model Selection: AIC vs BIC

### When to use AIC:
- **Prediction focus**: AIC is designed to minimize prediction error
- **Moderate sample sizes**: Works well with typical tag-recapture datasets
- **Exploratory analysis**: Better for identifying potentially important time steps
- **Use AIC when**: You want to capture seasonal growth patterns that might be real, even if complex

### When to use BIC:
- **Many parameters**: BIC penalizes complexity more heavily (log(n) vs 2)
- **Large sample sizes**: BIC penalty increases with n, preventing overfitting
- **Parsimony priority**: Better when you want the simplest adequate model
- **Use BIC when**: You have many time steps and want to avoid overfitting

### Recommendation for Tag2STM:
For tag-recapture growth models with many parameters (nlbin × ntsteps growth parameters), **consider using both**:
```r
# Compare both criteria
out <- out[order(out$AIC), ]
print("Top 5 models by AIC:")
print(head(out, 5))

out <- out[order(out$BIC), ]
print("Top 5 models by BIC:")
print(head(out, 5))

# If AIC and BIC agree, use that model
# If they disagree, BIC is likely more conservative (fewer time steps)
```

**General guideline**: 
- If sample size (n) > 40 × number of parameters (k), use BIC
- Otherwise, use AIC or compare both
- For Western Rock Lobster with thousands of tag returns, **BIC is recommended**

## Parameters

### Growth Parameters
- `growth_vecpar`: Log-scale increment parameters defining growth curves
- `LsigGrow`: Log SD for growth variability in STM

### Measurement Error
- `LsigError`: Log SD for baseline measurement error
- `MerrorRel`: Individual random effects at release
- `MerrorRec`: Individual random effects at recapture
- `LMerrorRelsigma`: Log SD for release random effects
- `LMerrorRecsigma`: Log SD for recapture random effects

## Model Assumptions

1. Growth variability is constant across length bins (not proportional)
2. Measurement error is normally distributed
3. Growth is deterministic within time steps (variation captured by STM spread)
4. Animals in the same release cohort are exchangeable
5. Time steps represent biologically meaningful periods

## Output Files

STM CSV files contain nlbin × nlbin matrices where element [i,j] is the probability of an animal in length bin j transitioning to length bin i during that time step.

## Citation

If you use this package, please cite:
```
de Lestang, S. (2025). Tag2STM: Size Transition Matrices from Tag-Recapture Data. 
R package. https://github.com/sdelestang/Tag2STM
```

## References

- Wang, Y.G., and Thomas, M.R. (1995). Accounting for individual variability in the von Bertalanffy growth model. Canadian Journal of Fisheries and Aquatic Sciences, 52(7), 1368-1375.

## Support

For questions or issues, please open an issue on GitHub or contact simon.delestang@dpird.wa.gov.au

## License

[Add your license here]
