circadianDrawing <- function(tod, expr, apar, labels, specInfo=NULL){	
  getPred <- function(parS, xx) {	
    parS$A * sin(2*pi/24 * (xx + parS$phase)) + parS$offset
  }
  
  geneName <- apar$genes
  peak <- round(apar$peak)
  if(peak==18) peak <- -6
  #pvalue <- signif(apar$pvalue,3)
  
  # amain <- paste('PV L3 healthy\n',geneName,':',probeName,'\n','p-value =',apvalue,sep='')
  amain <- paste(specInfo,', ',geneName,': ','; peak = ',peak,sep='')
  
  times <- seq(-6,18,0.1)
  pred <- getPred(apar,times)
  
  labelColor <- as.numeric(factor(labels))
  
  plot(tod,expr,col=labelColor, pch=16,cex=2,
       main=amain,xlim=c(-6,18),
       xaxt="n", yaxt="n",
       xlab='',ylab='')
  ytick = pretty(par("usr")[3:4])
  yl = formatC(ytick, format="f", digits=1) 
  axis(2,cex.axis=2.2,at=ytick,labels = yl)
  mtext("Expression", side=2, line=2.6, cex=2.2)
  xtick = pretty(par("usr")[1:2])
  xl = formatC(xtick, format="f", digits=0) 
  axis(1,cex.axis=2.2,at=xtick,labels = xl)
  mtext("TOD", side=1, line=2.6, cex=2.2)
  
  smoothingSpline = smooth.spline(times, pred, spar=0.35)
  lines(smoothingSpline,col='red',lwd=4)
  box(which = "plot", lty = "solid",lwd=3)	
  #legend('topright',legend=unique(labels),col=unique(labelColor),pch=16,cex=2)
}