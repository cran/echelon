\name{echelon}
\alias{echelon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Echelon analysis for spatial data
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The \code{echelon} function divides the study area into structural entities, called \emph{'echelons'}, based on neighbor information and draws a dendrogram.
}
\usage{
echelon(x, nb, dendrogram = TRUE, name = NULL,
      main = NULL, ylab = NULL, yaxes = TRUE, ylim = NULL,
      xaxes = FALSE, xdper = c(0, 1), dmai = NULL,
      col = 1, lwd = 1, symbols = 4, cex.symbols = 1, col.symbols = 4,
      ens = TRUE, adj.ens = 1, cex.ens = 0.8, col.ens = 1,
      profiles = FALSE, nb.check = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A numeric vector containing data values.}
  \item{nb}{
Neighbor information data: an object of class \code{nb} or a weights matrix.
}
  \item{name}{
Region names. if NULL, it is assigned \code{seq_along(x)}.
}
  \item{dendrogram}{
Logical. if TRUE, draws an echelon dendrogram.
}
  \item{main}{
Related to dendrogram drawing. The main title for the dendrogram.
}
  \item{ylab}{
Related to dendrogram drawing. The title for the y-axis.
}
  \item{yaxes}{
Related to dendrogram drawing. Logical. if TRUE, draws the y-axis.
}
  \item{ylim}{
Related to dendrogram drawing. If not specified, the y-axis scale is set to \code{c(min, max)}.
}
  \item{xaxes}{
Related to dendrogram drawing. Logical. if TRUE, draws the x-axis.
}
  \item{xdper}{
Related to dendrogram drawing. The percentage of the x-axis to display, specified in \emph{[0, 1]}.
}
  \item{dmai}{
Related to dendrogram drawing. A numeric vector of the form \code{c(bottom, left, top, right)} specifying margin sizes in inches. Default is \code{c(0.4, 0.8, 0.3, 0.01)}.
}
\item{col}{
Related to dendrogram drawing. The line color of the dendrogram.
}
  \item{lwd}{
Related to dendrogram drawing. The line width of the dendrogram.
}
  \item{symbols}{
Related to dendrogram drawing. An integer specifying a symbol or a single character. If integer, it corresponds to \code{\link{pch}} in \code{\link{par}}.
}
  \item{cex.symbols}{
Related to dendrogram drawing. A magnification factor for the plotting symbols.
}
  \item{col.symbols}{
Related to dendrogram drawing. The color for the plotting symbols.
}
  \item{ens}{
Related to dendrogram drawing. Logical. if TRUE, draw the labels of echelon numbers.
}
  \item{adj.ens}{
Related to dendrogram drawing. Adjusts the position of echelon number labels (see \code{\link{text}} for 'adj').
}
  \item{cex.ens}{
Related to dendrogram drawing. A magnification factor for the echelon number labels.
}
  \item{col.ens}{
Related to dendrogram drawing. The color for the echelon number labels.
}
  \item{profiles}{
Logical. If TRUE, returns the echelon profiles result (see [2] for details).
}
  \item{nb.check}{
Logical. if TRUE, checks for errors in the neighbor information data.
}
}
\value{
The \code{echelon} function returns an object of class \code{echelon}, which contains the following components: \item{Table}{A summary of each echelon.} \item{Echelons}{The regions that make up each echelon.}
%% ...
}
\references{
[1] Myers, W.L., Patil, G.P. and Joly, K. (1997). Echelon approach to areas of concern in synoptic regional monitoring. \emph{Environmental and Ecological Statistics}, \strong{4}, 131--152.

[2] Kurihara, K., Myers, W.L. and Patil, G.P. (2000) Echelon analysis of the relationship between population and land cover patter based on remote sensing data. \emph{Community ecology}, \strong{1}, 103--122.
}
\author{
Fumio Ishioka
}
\note{
Any \code{NA} values in \code{x} are replaced with the minimum value of \code{x}.

The functions \code{Sf::st_read} and \code{spdep::poly2nb} are helpful for creating the object specified in the \code{nb} argument.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{echepoi}} and \code{\link{echebin}} for cluster detection based on echelons.
}

\examples{
##Echelon analysis for one-dimensional data with 25 regions
#A weights matrix
one.nb <- matrix(0,25,25)
one.nb[1,2] <- 1
for(i in 2:24) one.nb[i,c(i-1,i+1)] <- c(1,1)
one.nb[25,24] <- 1

#25 random values
one.dat <- runif(25) * 10

#Echelon analysis
echelon(x = one.dat, nb = one.nb)


##Echelon analysis for SIDS data for North Carolina
#Mortality rate per 1,000 live births from 1974 to 1984
library(spData)
data("nc.sids")
SIDS.cas <- nc.sids$SID74 + nc.sids$SID79
SIDS.pop <- nc.sids$BIR74 + nc.sids$BIR79
SIDS.rate <- SIDS.cas * 1000 / SIDS.pop

#Echelon analysis
SIDS.echelon <- echelon(x = SIDS.rate, nb = ncCR85.nb, name = row.names(nc.sids),
  symbols = 12, cex.symbols = 1.5, ens = FALSE)
text(SIDS.echelon$coord, labels = SIDS.echelon$regions.name,
  adj = -0.1, cex = 0.7)

#Echelon Profiles
echelon(x = SIDS.rate, nb = ncCR85.nb, profiles = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ echelon analysis }% use one of  RShowDoc("KEYWORDS")
\keyword{ spatial scan statistic }% __ONLY ONE__ keyword per line
\keyword{ spatial cluster deteciotn }% __ONLY ONE__ keyword per line
