#' @importFrom stats pnorm dnorm na.omit quantile setNames weighted.mean
#' @importFrom graphics abline axis layout lines mtext par polygon plot
#' @importFrom grDevices grey
#' @importFrom RTMB getAll REPORT MakeADFun
NULL

#' Global variables used in package functions
#' @noRd
utils::globalVariables(c(
  "mod", "tdat", "lbin", "lbinL", "lbinU", "datain", "goodts",
  "ntsteps", "M", "times", "mout",
  "Tlbin", "nlbin", "LsigError", "sigGrow", "ipin", "spin_left",
  "spin_right", "mxpin", "mnpin", "MerrorRel", "Rlcl", "relts",
  "MerrorRec", "Rccl", "nlob", "LMerrorRelsigma", "LMerrorRecsigma",
  "obslen", "Clbin", "month", "hm", "tstep"
))
