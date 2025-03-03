\name{echebin}
\alias{echebin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Echelon spatial scan statistic based on Binomial model
}
\description{
The \code{echebin} function detects spatial clusters using the echelon spatial scan statistic with a Binomial model.
}
\usage{
echebin(echelon.obj, cas, ctl, K = length(cas)/2, Kmin = 1, n.sim = 99,
        cluster.type = "high", cluster.legend.pos = "bottomleft",
        dendrogram = TRUE, cluster.info = FALSE, coo = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{echelon.obj}{
An object of class \code{echelon}. For details, see \code{\link{echelon}}.
}
  \item{cas}{
A numeric (integer) vector of case counts. \code{NA} values are not allowed.
}
  \item{ctl}{
A numeric (integer) vector of control counts. \code{NA} values are not allowed.
}
  \item{K}{
Maximum cluster size. If \code{K} >= 1 (integer), the cluster size is limited to \code{K} regions. If 0 < \code{K} < 1, the cluster size is limited to \code{K} * 100\% of the total population.
}
  \item{Kmin}{
Minimum cluster size. 
}
  \item{n.sim}{
The number of Monte Carlo replications used for significance testing of detected clusters. If set to 0, significance is not assessed.
}
  \item{cluster.type}{
A character string specifying the cluster type. If \code{"high"}, the detected clusters have high rates (hotspot). If \code{"low"}, the detected clusters have low rates (coldspot).
}
  \item{cluster.legend.pos}{
The location of the legend on the dendrogram. (See \code{\link{legend}} for details.)
}
  \item{dendrogram}{
Logical. If TRUE, draws an echelon dendrogram with the detected clusters.
}
  \item{cluster.info}{
Logical. If TRUE, returns detailed results of the detected clusters.
}
  \item{coo}{
An array of (x, y) coordinates for the region centroids to plot a cluster map.
}
  \item{\dots}{
Related to dendrogram drawing. (See the help for \code{\link{echelon}})
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{clusters}{Each detected cluster.}
\item{scanned.regions}{A region list of all scanning processes.}
\item{simulated.LLR}{Monte Carlo samples of the log-likelihood ratio.}
%% ...
}
\references{
[1] Kulldorff M, Nagarwalla N. (1995). Spatial disease clusters: Detection and inference. \emph{Statistics in Medicine}, \strong{14}, 799--810.

[2] Kulldorff M. (1997). A spatial scan statistic. \emph{Communications in Statistics: Theory and Methods}, \strong{26}, 1481--1496.
}
\author{
Fumio Ishioka
}
\note{
The function \code{echebin} requires either \code{cas} or \code{ctl}.

Population is defined as the sum of \code{cas} and \code{ctl}.

Typical values of \code{n.sim} are 99, 999, 9999, ...
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{echelon}} for the echelon analysis.

\code{\link{echepoi}} for cluster detection based on echelons using Poisson model.
}

\examples{
##Hotspot detection for non-white birth of North Carolina using echelon scan

#Non-white birth from 1974 to 1984 (case data)
library(spData)
data("nc.sids")
nwb <- nc.sids$NWBIR74 + nc.sids$NWBIR79

#White birth from 1974 to 1984 (control data)
wb <- (nc.sids$BIR74 - nc.sids$NWBIR74) + (nc.sids$BIR79 - nc.sids$NWBIR79)

#Hotspot detection based on Binomial model
nwb.echelon <- echelon(x = nwb/wb, nb = ncCR85.nb, name = row.names(nc.sids))
echebin(nwb.echelon, cas = nwb, ctl = wb, K = 20,
  main = "Hgih rate clusters", ens = FALSE)
text(nwb.echelon$coord, labels = nwb.echelon$regions.name,
  adj = -0.1, cex = 0.7)

#Detected clusters and neighbors map
#XY coordinates of each polygon centroid point
NC.coo <- cbind(nc.sids$lon, nc.sids$lat)
echebin(nwb.echelon, cas = nwb, ctl = wb, K = 20,
  coo = NC.coo, dendrogram = FALSE)

\donttest{
##Detected clusters map
#Here is an example using the sf class "sf"
nwb.clusters <- echebin(nwb.echelon, cas = nwb,
   ctl = wb, K = 20, dendrogram = FALSE)
MLC <- nwb.clusters$clusters[[1]]
Secondary <- nwb.clusters$clusters[[2]]
cluster.col <- rep(0,times=length(nwb))
cluster.col[MLC$regionsID] <- 2
cluster.col[Secondary$regionsID] <- 3

library(sf)
nc <- st_read(system.file("shape/nc.shp", package = "sf"))
plot(nc$geometry, col = cluster.col,
main = "Detected high rate clusters")
text(st_coordinates(st_centroid(st_geometry(nc))),
  labels = nc$CRESS_ID, cex =0.75)
legend("bottomleft",
c(paste("1- p-value:", MLC$p),
  paste("2- p-value:", Secondary$p)),
  text.col = c(2,3))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ echelon analysis }% use one of  RShowDoc("KEYWORDS")
\keyword{ spatial scan statistic }% __ONLY ONE__ keyword per line
\keyword{ spatial cluster deteciotn }% __ONLY ONE__ keyword per line
