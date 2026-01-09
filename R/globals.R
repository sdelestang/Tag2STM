#' @importFrom stats na.omit quantile setNames weighted.mean
#' @importFrom graphics abline axis layout lines mtext par polygon plot
#' @importFrom grDevices grey
#' @importFrom RTMB getAll REPORT MakeADFun
NULL

# Tell R CMD check these are intentionally not imported
# (RTMB provides its own versions for automatic differentiation)
utils::globalVariables(c("pnorm", "dnorm"))

# Variables used via getAll() or non-standard evaluation
utils::globalVariables(c(
  "mod", "tdat", "lbin", "lbinL", "lbinU", "datain", "goodts",
  "ntsteps", "M", "times", "mout",
  "Tlbin", "nlbin", "LsigError", "sigGrow", "ipin", "spin_left",
  "spin_right", "mxpin", "mnpin", "MerrorRel", "Rlcl", "relts",
  "MerrorRec", "Rccl", "nlob", "LMerrorRelsigma", "LMerrorRecsigma",
  "obslen", "Clbin", "month", "hm", "tstep",
  "sigGrowsd", "nobs", "sigError"  # Add these for growmod2
))
