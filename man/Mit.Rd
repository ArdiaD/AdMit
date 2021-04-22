\name{Mit}
\alias{dMit}
\alias{rMit}
\concept{Mixture of Student-t Distributions}
\concept{Adaptive mixture}
\concept{Importance sampling}
\title{Mixture of Student-t Distributions}
\description{
  Density function or random generation for an adaptive mixture of
  Student-t distributions
}
\usage{%
dMit(theta, mit = list(), log = TRUE)
rMit(N = 1, mit = list())
}
\arguments{%
  \item{theta}{matrix (of size \eqn{N \times d}{Nxd}, where
    \eqn{N,d \geq 1}{N,d>=1}) of real values.}
  \item{mit}{list containing information on the mixture approximation (see *Details*).}
  \item{log}{logical; if \code{log = TRUE}, returns (natural) logarithm
    values of the density. Default: \code{log = TRUE}.}
  \item{N}{number of draws (positive integer number).}
}
\details{%
  \code{dMit} returns the density values while \code{rMit} generates
  draws from a mixture of Student-t distributions.
   
  The argument \code{mit} is a list containing information on the
  adaptive mixture of Student-t distributions. The following components must
  be provided:
  
  \describe{%
    \item{\code{p}}{vector (of length \eqn{H}) of mixture probabilities.}
    \item{\code{mu}}{matrix (of size \eqn{H \times d}{Hxd}) containing
      the vectors of modes (in row) of the mixture components.}
    \item{\code{Sigma}}{matrix (of size \eqn{H \times d^2}{Hxd*d})
      containing the scale matrices (in row) of the mixture components.}
    \item{\code{df}}{degrees of freedom parameter of the Student-t
      components (integer number not smaller than one).}
  }
  
  where \eqn{H (\geq 1)}{H (>=1)} is the number of components and
  \eqn{d (\geq 1)}{d (>=1)} is
  the dimension of the mixture approximation. Typically,
  \code{mit} is estimated by the function \code{\link{AdMit}}. If the
  \code{mit = list()}, a Student-t distribution located
  at \code{rep(0,d)} with scale matrix \code{diag(d)} and one
  degree of freedom parameter is used.
}
\value{%
  Vector (of length \eqn{N} of density values, or matrix (of size
  \code{N}x\eqn{d}) of random draws, where \eqn{d (\geq 1)}{d (>=1)} is the
  dimension of the mixture approximation.
}
\note{%
  Further details and examples of the \R package \code{AdMit}
  can be found in Ardia, Hoogerheide, van Dijk (2009a,b). See also
  the package vignette by typing \code{vignette("AdMit")}.

  Please cite the package in publications. Use \code{citation("AdMit")}.
}
\author{David Ardia}
\references{%
Ardia, D., Hoogerheide, L.F., van Dijk, H.K. (2009a).
  AdMit: Adaptive Mixture of Student-t Distributions.
  \emph{R Journal} \bold{1}(1), pp.25-30.
  \doi{10.32614/RJ-2009-003}

  Ardia, D., Hoogerheide, L.F., van Dijk, H.K. (2009b).
  Adaptive Mixture of Student-t Distributions as a Flexible Candidate
  Distribution for Efficient Simulation: The R Package AdMit.
  \emph{Journal of Statistical Software} \bold{29}(3), pp.1-32.
  \doi{10.18637/jss.v029.i03}
}
\seealso{%
  \code{\link{AdMit}} for fitting an adaptive mixture of
  Student-t distributions to a given function \code{KERNEL},
  \code{\link{AdMitIS}} for importance sampling using an adaptive
  mixture of Student-t distributions as the importance density,
  \code{\link{AdMitMH}} for the independence chain Metropolis-Hastings
  using an adaptive mixture of Student-t distributions as the
  candidate density.
}
\examples{%
  ## NB : Low number of draws for speedup. Consider using more draws!
  ## One dimensional two components mixture of Student-t distributions
  mit <- list(p = c(0.5, 0.5),
              mu = matrix(c(-2.0, 0.5), 2, 1, byrow = TRUE),
              Sigma = matrix(0.1, 2),
              df = 10)
  ## Generate draws from the mixture
  hist(rMit(1e4, mit = mit), nclass = 100, freq = FALSE)
  x <- seq(from = -5.0, to = 5.0, by = 0.01)
  ## Add the density to the histogram
  lines(x, dMit(x, mit = mit, log = FALSE), col = "red", lwd = 2)

  ## Two dimensional (one component mixture) Student-t distribution
  mit <- list(p = 1,
              mu = matrix(0.0, 1.0, 2.0),
              Sigma = matrix(c(1.0, 0.0, 0.0, 1.0), 1, 4),
              df = 10)
  ## Function used to plot the mixture in two dimensions
  dMitPlot <- function(x1, x2, mit = mit)
  {
    dMit(cbind(x1, x2), mit = mit, log = FALSE)
  }
  x1 <- x2 <- seq(from = -10.0, to = 10.0, by = 0.1)
  thexlim <- theylim <- range(x1)
  z <- outer(x1, x2, FUN = dMitPlot, mit = mit)
  ## Contour plot of the mixture
  contour(x1, x2, z, nlevel = 20, las = 1, 
          col = rainbow(20),
          xlim = thexlim, ylim = theylim)
  par(new = TRUE)
  ## Generate draws from the mixture
  plot(rMit(1e4, mit = mit), pch = 20, cex = 0.3, 
            xlim = thexlim, ylim = theylim, col = "red", las = 1)

}
\keyword{distribution}
