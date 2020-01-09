\name{echepoi}
\alias{echepoi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Echelon spatial scan statistic based on Poisson model
}
\description{
\code{echepoi} detects spatial clusters using echelon spatial scan statistic based on Poisson model.
}
\usage{
echepoi(echelon.obj, cas, pop = NULL, ex = NULL, K = length(cas)/2, n.sim = 99,
        cluster.type = "high", cluster.legend.pos = "bottomleft",
        dendrogram = TRUE, cluster.info = FALSE, coo = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{echelon.obj}{
An object of class \code{echelon}. See \code{\link{echelon}}.
}
  \item{cas}{
A numeric (integer) vector of case counts.
\code{NA}s are not allowed.
}
  \item{pop}{
A numeric (integer) vector for population.
\code{NA}s are not allowed.
}
  \item{ex}{
A numeric vector for expected cases.
\code{NA}s are not allowed.
}
  \item{K}{
Maximum cluster size. if \code{K} >= 1 (integer), the cluster size is limit to less than or equal to number of regions \code{K}. On the other hand, if 0 < \code{K} < 1, the cluster size is limit to less than or equal to \code{K *} 100\% of the total population.
}
  \item{n.sim}{
Number of Monte Carlo replications used for significance testing of detected clusters. If 0, the significance is not assessed.
}
  \item{cluster.type}{
A character string specifying the cluster type. If \code{"high"}, the detected clusters have high rates (hotspot clusters). On the other hand, If \code{"low"}, the detected clusters have low rates (coldspot cluster).
}
  \item{cluster.legend.pos}{
A location of the legend on the dendrogram. (See the help for \code{\link{legend}})
}
  \item{dendrogram}{
Logical. if TRUE, draw an echelon dendrogram with detected clusters.
}
  \item{cluster.info}{
Logical. if TRUE, return the result of detected clusters for detail.
}
  \item{coo}{
An array of (x,y)-coordinates of the region centroid to draw a cluster map.
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
[1] Kulldorff M. (1997). A spatial scan statistic. \emph{Communications in Statistics: Theory and Methods}, \strong{26}, 1481--1496.

[2] Ishioka F, Kawahara J, Mizuta M, Minato S, and Kurihara K. (2019) Evaluation of hotspot cluster detection using spatial scan statistic based on exact counting. \emph{Japanese Journal of Statistics and Data Science}, \strong{2}, 241--262.
}
\author{
Fumio Ishioka
}
\note{
\code{echepoi} requires either \code{pop} or \code{ex}.

Typical values of \code{n.sim} are 99, 999, 9999, ...
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{echelon}} for the echelon analysis.

\code{\link{echebin}} for cluster detection based on echelons using Binomial model.
}

\examples{
##Hotspot detection for SIDS data of North Carolina using echelon scan

#Mortality rate per 1,000 live births from 1974 to 1984
library(spData)
data("nc.sids")
SIDS.cas <- nc.sids$SID74 + nc.sids$SID79
SIDS.pop <- nc.sids$BIR74 + nc.sids$BIR79
SIDS.rate <- SIDS.cas * 1000 / SIDS.pop

#Hotspot detection based on Poisson model
SIDS.echelon <- echelon(x = SIDS.rate, nb = ncCR85.nb, name = row.names(nc.sids))
echepoi(SIDS.echelon, cas = SIDS.cas, pop = SIDS.pop, K = 20,
  main = "Hgih rate clusters", ens = FALSE)
text(SIDS.echelon$coord, labels = SIDS.echelon$regions.name,
  adj = -0.1, cex = 0.7)

#Detected clusters and neighbors map
#XY coordinates of each polygon centroid point
NC.coo <- cbind(nc.sids$lon, nc.sids$lat)
echepoi(SIDS.echelon, cas = SIDS.cas, pop = SIDS.pop, K = 20,
  coo = NC.coo, dendrogram = FALSE)

\donttest{
##Detected clusters map
#Here is an example using the sf class "sf"
SIDS.clusters <- echepoi(SIDS.echelon, cas = SIDS.cas,
  pop = SIDS.pop, K = 20, dendrogram = FALSE)
MLC <- SIDS.clusters$clusters[[1]]
Secondary <- SIDS.clusters$clusters[[2]]
cluster.col <- rep(0,times=length(SIDS.rate))
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
