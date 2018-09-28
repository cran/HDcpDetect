\name{wild.binary.segmentation}
\alias{wild.binary.segmentation}

\title{
wild.binary.segmentation
}
\description{
Detects change points in time series data using the wild binary segmentation algorithm from Fryzlewicz (2014).
}
\usage{
wild.binary.segmentation(data_M,minsize=15,num_intervals=1250,M_threshold=0.05)
}
\arguments{
  \item{data_M}{
  An nxp matrix representing a times series of length n with p dimensions.
  }
  \item{minsize}{
  The minimum interval length.
  }
  \item{num_intervals}{
  The number of random intervals to be generated and tested for change points.
  }
  \item{M_threshold}{
  Value used as a threshold to estimate temporal dependence by determining how small of a standardized difference is indistinguishable from zero. }
  }

\details{
Increasing the minimum interval length will generally reduce type I error while increasing type II error.
}
\value{
The returned value is a list of the estimated change point locations.
}
\references{
Li, J., Li, L., Xu, M., Zhong, P (2018). Change Point Detection in the Mean of
High-Dimensional Time Series Data under Dependence. Manuscript.
Fryzlewicz, P. (2014). Wild Binary Segmentation for Multiple Change-point Detection. The Annals of Statistics.
}
\author{
Jun Li, Jeffrey Okamoto, and Natasha Stewart}

\examples{
library(HDcpDetect)
HAPT2 <- as.matrix(HAPT[1:35,])
wild.binary.segmentation(HAPT2)
}

\keyword{time-series}
\keyword{change-point}
