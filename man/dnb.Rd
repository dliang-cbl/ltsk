\name{dnb}
\alias{dnb}
\title{
Search Neighbours in Time and Space Within Specified Ranges
}
\description{
A neighbor search implementation based on ANN trees to identify observed data points within a given distance around location and time interval.
}
\usage{
dnb(query, obs, th, xcoord='x',ycoord='y',tcoord='t',
	future=TRUE,by=NULL,nbin=NULL,Large=2000,
	cl=NULL,cluster=c(1,3))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{query}{
	a vector; the x, y coordinates and the time stamp of the query point
}
  \item{obs}{
	a matrix; the x, y coordinates and time stamps of the spatiotemporal locations
}
  \item{th}{
	a vector; the distance threshold and time lag
}
  \item{xcoord}{
a character constant, the field name for x coordinate in both \code{query} and \code{obs}
}
  \item{ycoord}{
a character constant, the field name for y coordinate in both \code{query} and \code{obs}
}
  \item{tcoord}{
a character constant, the field name for time coordinate in both \code{query} and \code{obs}
}
	\item{future}{
	logical, whether include observed spatiotemporal points future in time relative to the query spatiotemporal location.
}
  \item{by}{
a vector of \code{xcoord}, \code{ycoord} and \code{zcoord} resolutions used in neighbor search
}
  \item{nbin}{
if \code{by} is \code{NULL}, a vector of bins to divide the \code{xcoord}, \code{ycoord} and \code{zcoord} data range, Default is 200 for space and unique time units.
}
  \item{Large}{
a numeric constant, upper limit of neighbor points, beyond which subsampling is performance
}
  \item{cl}{
	a \pkg{parallel} cluster object (default NULL means single core)
}
  \item{cluster}{
	if cl is a \pkg{parallel} object, the time and space domain in terms of \code{th} to divide the query and send to nodes (default (2,10))
}
}
\details{

Implementation involves first calculating the time lags between query point and observed data (with locational coordinates and time); for observed locations within time lag of query, the function calculates the Euclidean distances between query location and all potential neighbors and select those within specified distance threshold.

The future argument can be used to exclude data in the future in neighbor search. This is useful in an extrapolation application.
}
\value{
A list of vectors, row numbers in the observed data matrix, that are within the given distance threshold and time lag of the query location.
}
\author{
Dong Liang (dliang@umces.edu)
}
\note{
For large dataset, use ANN (for spatial kriging) and Range Tree for spatiotemporal Kriging.
}

\seealso{
\code{get.knn} in \pkg{FNN}
}

\examples{
data(epa_cl)
coords <- c('x','y','t')
ii <- dnb(query[1,coords],obs[,coords],c(0.1,10))
}
