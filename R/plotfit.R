#' Plot Tag-Recapture Model Diagnostics
#'
#' Creates a multi-panel diagnostic plot for fitted tag-recapture growth models,
#' showing residuals, growth curves, length distributions over time, and projected
#' growth trajectories with confidence intervals.
#'
#' @details
#' This function expects the following objects to exist in the calling environment:
#' \describe{
#'   \item{mod}{A fitted RTMB model object from `MakeADFun()`}
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
#'   \item Growth increment curves by length
#'   \item Growth variability (sigma) across increment sizes
#'   \item Length distributions at 1, 5, 10, 15, and 20 years
#'   \item Projected growth trajectory with 95% confidence intervals
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
#' @importFrom dplyr mutate
#' @importFrom magrittr %>% %<>%
#' @export
plotfit <- function(){
  out <- mod$rep()
  lbin <- bins$lbin
  lbinL <- bins$lbinL
  lbinU <- bins$lbinU
  lenout <- out$lenout
  par(las=1, mar=c(5,4,1,1))
  lout <- matrix(c(1,2,3,4,5,6), ncol=2, byrow = T)
  layout(lout)
  let <- 1;
  estlen <- apply((out$EstRecLen),1,function(x) weighted.mean(lbin, w = x))
  obs <- data.frame(tag=tdat$tag, estlen=estlen, obslen=tdat$rccl) %>% mutate(res=obslen-estlen)
  estlbin <- apply((out$EstRecLen),1,function(x) which.max(x)[1])
  obs %<>% mutate(Rlcl=datain$Rlcl, Clbin=as.numeric(cut(tdat$rccl, breaks = c(bins$lbinL,(max(bins$lbinL)+30)), include.lowest = T, right = F)),  ntstep=datain$tsteps)
str(obs)
  plot(obs$Rlcl, obs$res, pch=16, col=grey(0.2,0.1), xlab='Release length bin', ylab='Residual')
  abline(h=0,lty=3)
  mtext(letters[let], 3, adj=0); let <- let + 1
  plot(log(obs$ntstep), obs$res, pch=16, col=grey(0.2,0.1), xlab='Liberty (log number timesteps)', ylab='Residual')
  abline(h=0,lty=3)
  mtext(letters[let], 3, adj=0); let <- let + 1

  growthmat <- out$growthmat[goodts,]
  mx <- ifelse(max(growthmat)>1,max(growthmat),1)
  if(is.matrix(growthmat)){
    plot(lbin,growthmat[1,],type='l',ylim=c(0,mx),ylab='Increment',xlab='Length bin',main="Growth/Moult")
    for(i in 2:nrow(growthmat)){ lines(lbin,growthmat[i,],col=i)  }}
  if(!is.matrix(growthmat)){
    plot(lbin,growthmat,type='l',ylim=c(0,mx),ylab='Increment',xlab='Length bin',main="Growth/Moult")}

  plot(  seq(0,5,0.1), out$sigGrowvec, type='o',ylab="Growth spread (sigma)", xlab='growth')

  plot(lbin, lenout[ntsteps,], type='o',pch=16, col=1, ylab='Proportion', xlab='Carapace length', ylim=c(0,max(lenout[ntsteps*1:15,])), xlim=c(min(lbinL),max(lbinL)), main="1,5,10,15,20 years")
  lines(lbin, lenout[ntsteps*5,]/sum(lenout[ntsteps*5,]), type='o',pch=16, col=2)
  lines(lbin, lenout[ntsteps*10,]/sum(lenout[ntsteps*10,]), type='o',pch=16, col=3)
  lines(lbin, lenout[ntsteps*15,]/sum(lenout[ntsteps*15,]), type='o',pch=16, col=4)
  lines(lbin, lenout[ntsteps*20,]/sum(lenout[ntsteps*20,]), type='o',pch=16, col=5)
  mtext(letters[let], 3, adj=0); let <- let + 1

  mxorig <- apply((lenout),1,function(x) weighted.mean(lbin, w = x))
  weighted.probs <- function(x) quantile(rep(lbin, 100000*x), probs=c(0.025, 0.5, 0.975))
  mx <- apply(lenout,1,weighted.probs)
  xs <- (0:(ntsteps*30))
  Mxage <-  trunc((3.5 / M)/5)*5 +5
  Mxtstep <- Mxage * ntsteps

  plot(xs[xs<Mxtstep],c(lbin[1],mx['50%',])[xs<Mxtstep], type='o', cex=0.1, xlim=c(0,ntsteps*30), ylim=c(0, max(lbinL)), axes=F, xlab='Relative Age (y)', ylab='Carapace length')
  y0 <- c(lbin[1],mx[2,])[xs<Mxtstep];y1 <- c(lbin[1],mx[1,])[xs<Mxtstep];y2 <- c(lbin[1],mx[3,])[xs<Mxtstep]
  polygon(c(xs[xs<Mxtstep],rev(xs[xs<Mxtstep])),c(y1[xs<Mxtstep],rev(y2[xs<Mxtstep])),col='grey70',border=NA)
  lines(xs[xs<Mxtstep],y0[xs<Mxtstep])
  axis(1,seq(2,((ntsteps*Mxage)+2),ntsteps), 2:(Mxage+2));
  axis(2,las=1)
  mtext(letters[let], 3, adj=0); let <- let + 1
  return(list(out, obs)) #obs2
}
