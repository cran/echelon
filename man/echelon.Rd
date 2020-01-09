\name{echelon}
\alias{echelon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Echelon analysis for spatial data
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{echelon} divides study area into structural entities consisting of peaks or foundations, say \emph{'echelons'}, based on neighbor information, and draw its dendrogram.
}
\usage{
echelon(x, nb, dendrogram = TRUE, name = NULL,
      main = NULL, ylab = NULL, yaxes = TRUE, ylim = NULL,
      xaxes = FALSE, xdper = c(0, 1), dmai = NULL,
      col = 1, lwd = 1, symbols = 4, cex.symbols = 1, col.symbols = 4,
      ens = TRUE, adj.ens = 1, cex.ens = 0.8, col.ens = 1,
      profiles = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A numeric vector of data values.}
  \item{nb}{
Neighbor information data. An object of class \code{nb} or a weights matrix.
}
  \item{name}{
The region names. if NULL, it is assigned \code{1:length(x)}.
}
  \item{dendrogram}{
Logical. if TRUE, draw an echelon dendrogram.
}
  \item{main}{
Related to dendrogram drawing. An overall title for the dendrogram.
}
  \item{ylab}{
Related to dendrogram drawing. A title for the y axis.
}
  \item{yaxes}{
Related to dendrogram drawing. Logical. if TRUE, draw the y axis.
}
  \item{ylim}{
Related to dendrogram drawing. A scale of y axis given by \code{c(min, max)}.
}
  \item{xaxes}{
Related to dendrogram drawing. Logical. if TRUE, draw the x axis.
}
  \item{xdper}{
Related to dendrogram drawing. A display percentage of x axis. The full display percentage is given in \emph{[0, 1]}.
}
  \item{dmai}{
Related to dendrogram drawing. A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches. Default is set to \code{c(0.4, 0.8, 0.3, 0.01)}.
}
\item{col}{
Related to dendrogram drawing. A line color of the dendrogram.
}
  \item{lwd}{
Related to dendrogram drawing. A line width of the dendrogram.
}
  \item{symbols}{
Related to dendrogram drawing. Either an integer specifying a symbol or a single character. If integer, it is synonymous with \code{\link{pch}} in \code{\link{par}}.
}
  \item{cex.symbols}{
Related to dendrogram drawing. A magnification to be used for the plotting symbols.
}
  \item{col.symbols}{
Related to dendrogram drawing. A color to be used for the plotting symbols.
}
  \item{ens}{
Related to dendrogram drawing. Logical. if TRUE, draw the labels of echelon numbers.
}
  \item{adj.ens}{
Related to dendrogram drawing. Adjustment of the labels of the echelon numbers.(see the help for \code{\link{text}}("adj")).
}
  \item{cex.ens}{
Related to dendrogram drawing. A magnification to be used for the labels of echelon numbers.
}
  \item{col.ens}{
Related to dendrogram drawing. A color to be used for the labels of echelon numbers.
}
  \item{profiles}{
Logical. if TRUE, return the result of echelon profiles. (See [2] for the details of echelon profiles)
}
}
\value{
The function \code{echelon} returns an object of class \code{echelon}.
An object of class \code{echelon} contains the following components:
\item{Table}{Summary of each echelon.}
\item{Echelons}{Regions that composes each echelon.}
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
If there are \code{NA} in \code{x}, then that is set the minimum value of \code{x}.

\code{Sf::st_read} and \code{spdep::poly2nb} are useful for creating the object specified in the argument \code{nb}.
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
