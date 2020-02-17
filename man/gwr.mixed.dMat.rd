\name{gwr.mixed.dMat}
\title{Mixed GWR}
\description{
This function implements mixed (semiparametric) GWR
}
\usage{
gwr.mixed.dMat(formula, data, fixed.vars,
                     intercept.fixed=FALSE, bw, diagnostic=T, kernel="bisquare", 
                     adaptive=FALSE, dMat)
}
\arguments{
  \item{formula}{Regression model formula of a \link{formula} object }
  \item{data}{a Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{fixed.vars}{independent variables that appeared in the formula that are to be treated as global}
  \item{intercept.fixed}{logical, if TRUE the intercept will be treated as global}
  \item{bw}{bandwidth used in the weighting function, possibly calculated by \link{bw.gwr};fixed (distance) or adaptive bandwidth(number of nearest neighbours)}
  \item{diagnostic}{logical, if TRUE the diagnostics will be calculated}
  \item{kernel}{function chosen as follows:
  
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                
                exponential: wgt = exp(-vdist/bw);
                
                bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
                
                tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise; 
                
                boxcar: wgt=1 if dist < bw, wgt=0 otherwise}
  \item{adaptive}{if TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)}
  \item{dMat}{a pre-specified distance matrix, it can be calculated by the function \code{\link{gw.dist}}}
}
\value{
A list of class \dQuote{mgwr}:
  \item{GW.arguments}{a list class object including the model fitting parameters for generating the report file}
  \item{aic}{AICc value from this calibration}
  \item{df.used}{ effective degree of freedom}
  \item{rss}{residual sum of squares}
  \item{SDF}{a SpatialPointsDataFrame (may be gridded) or 
             SpatialPolygonsDataFrame object (see package \dQuote{sp}) integrated with coefficient estimates in its "data" slot.}
  \item{timings}{starting and ending time.}
  \item{this.call}{the function call used.}
}
\references{
Fotheringham S, Brunsdon, C, and Charlton, M (2002),
Geographically Weighted Regression: The Analysis of Spatially Varying Relationships, Chichester: Wiley.

Brunsdon C, Fotheringham AS, Charlton ME (1999) Some notes on parametric signficance 
tests for geographically weighted regression. Journal of Regional Science 39(3):497-524

Mei L-M, He S-Y, Fang K-T (2004) A note on the mixed geographically weighted regression 
model. Journal of regional science 44(1):143-157

Mei L-M, Wang N, Zhang W-X (2006) Testing the importance of the explanatory variables 
in a mixed geographically weighted regression model. Environment and Planning A 38:587-598

Nakaya T, Fotheringham AS, Brunsdon C, Charlton M (2005) Geographically Weighted Poisson Regression for Disease Association Mapping,
Statistics in Medicine 24: 2695-2717

Nakaya T et al. (2011) GWR4.0, \url{http://gwr.nuim.ie/}.

}
\note{
This is a modification of function gwr.mixed that uses a C++ implementation to allow faster run
times. It requires the distance matrix to be specified. 
}
\author{Binbin Lu \email{binbinlu@whu.edu.cn}, Fiona Evans}
\keyword{mixed GWR}
\concept{mixed GWR}

