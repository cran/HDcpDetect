\name{binary.segmentation}
\alias{binary.segmentation}

\title{
binary.segmentation
}
\description{
Detects change points in time series data using a binary segmentation algorithm.
}
\usage{
binary.segmentation(data_M,alpha=.05,power_enhancement=TRUE,M_threshold=0.05)
}

\arguments{
  \item{data_M}{
     An nxp matrix representing a times series of length n with p dimensions.
}
\item{alpha}{
     The critical value for the hypothesis testing procedure.
}
\item{power_enhancement}{
     Indicates whether to add a power enhancement term to the test statistic.
}
\item{M_threshold}{
     Value used as a threshold to estimate temporal dependence by determining how small of a standardized difference is indistinguishable from zero. 
}
}
\details{
    The power enhancement term reduces type II error but slows the algorithm.
}
\value{
  The returned value is a list with the following components
  \item{Foundlist}{The estimated locations of the change points}
  \item{pvalues}{The p values corresponding to each change point estimate}
}
\references{
Li, J., Li, L., Xu, M., Zhong, P (2018). Change Point Detection in the Mean of
High-Dimensional Time Series Data under Dependence. Manuscript.

}
\author{
Jun Li, Jeffrey Okamoto, and Natasha Stewart
}
\examples{
library(HDcpDetect)
HAPT2 <- as.matrix(HAPT[1:35,])
binary.segmentation(data_M=HAPT2,power_enhancement=FALSE)

}
\keyword{time-series}
\keyword{change-point}
