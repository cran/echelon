#-------------------------------------------------------------------------------
# Axis for dendrogram
#-------------------------------------------------------------------------------
e.dendrogram.axis <- function(main, ylab, yaxes, ylim, xaxes, xdper, dmai, peaks, x){

  dev.new()

  xleft <- min(xdper)*(sum(peaks) - 1)
  xright <- max(xdper)*(sum(peaks) - 1)

  if(is.null(ylim)) ylim <- c(min(x), max(x))

  if(is.null(ylab)) ylab <- ""

  if(is.null(dmai)) par(mai = c(0.4, 0.8, 0.3, 0.01))
  else par(mai = dmai)

  if(!length(peaks) == 1)
    plot(1, 1, type = "n", main = main, xlab = "", ylab = ylab, 
         xlim = c(xleft - 0.3, xright + 0.2), ylim = ylim, axes = F)

  else  plot(1, 1, type = "n", main = main, xlab = "", ylab = ylab,
             xlim = c(xleft - 0.3, xright), ylim = ylim, axes = F)

  if(yaxes){
    measure <- NULL
    for(i in 1:4){
      measure <- c(measure, quantile(ylim)[i])
      measure <- c(measure, median(c(quantile(ylim)[i], quantile(ylim)[i+1])))
    }
    measure <- c(measure, max(ylim))
    axis(side = 2, at = round(measure, digits = 3))
  }
  if(xaxes){
    measure <- c(min(xdper), median(c(min(xdper), max(xdper))), max(xdper))
    axis(side = 1, at = c(xleft - 0.3, median(c(xleft - 0.3, xright)), xright),
         labels = measure, lty = 3)
  }
}
