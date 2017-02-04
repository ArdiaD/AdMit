\name{AdMitIS}
\alias{AdMitIS}
\concept{Mixture of Student-t Distributions}
\concept{Adaptive mixture}
\concept{Importance sampling}
\title{Importance Sampling using an Adaptive Mixture of Student-t
  Distributions as the Importance Density}
\description{
  Performs importance sampling using an adaptive
  mixture of Student-t distributions as the importance density
}
\usage{%
AdMitIS(N = 1e5, KERNEL, G = function(theta){theta}, mit = list(), ...)
}
\arguments{%
  \item{N}{number of draws used in importance sampling (positive
    integer number). Default: \code{N = 1e5}.}
  \item{KERNEL}{kernel function of the target density on which the
    adaptive mixture of Student-t distributions is fitted. This
    function should be vectorized for speed purposes (i.e., its first
    argument should be a matrix and its output a vector). Moreover, the function must contain
    the logical argument \code{log}. If \code{log = TRUE}, the function
    returns (natural) logarithm values of the kernel function. \code{NA}
    and \code{NaN} values are not allowed. (See the function
    \code{\link{AdMit}} for examples of \code{KERNEL} implementation.)}
  \item{G}{function of interest used in importance
    sampling (see *Details*).}
  \item{mit}{list containing information on the mixture approximation (see *Details*).}
  \item{\ldots}{further arguments to be passed to \code{KERNEL} and/or \code{G}.}
}
\details{%
  The \code{AdMitIS} function estimates
  \eqn{E_p[g(\theta)]}{E_p[g(theta)]}, where \eqn{p} is the target
  density, \eqn{g} is an (integrable w.r.t. \eqn{p}) function and \eqn{E}{E} denotes
  the expectation operator, by importance sampling using an adaptive
  mixture of Student-t distributions as the importance density.
  
  By default, the function \code{G} is given by:
  \preformatted{%
    G <- function(theta)
    {
      theta
    } 
  }

  and therefore, \code{AdMitIS} estimates the mean of
  \code{theta} by importance sampling. For other definitions of
  \code{G}, see *Examples*.

  The argument \code{mit} is a list containing information on the
  mixture approximation. The following components must be provided:
  
  \describe{%
    \item{\code{p}}{vector (of length \eqn{H}) of mixing probabilities.}
    \item{\code{mu}}{matrix (of size \eqn{H \times d}{Hxd}) containing
      the vectors of modes (in row) of the mixture components.}
    \item{\code{Sigma}}{matrix (of size \eqn{H \times d^2}{Hxd*d})
      containing the scale matrices (in row) of the mixture components.}
    \item{\code{df}}{degrees of freedom parameter of the Student-t
      components (real number not smaller than one).}
  }
  
  where \eqn{H (\geq 1)}{H (>=1)} is the number of components of the
  adaptive mixture of Student-t distributions and
  \eqn{d (\geq 1)}{d (>=1)} is the dimension of the first argument in \code{KERNEL}. Typically,
  \code{mit} is estimated by the function \code{\link{AdMit}}.
}
\value{%
  A list with the following components:\cr
  
  \code{ghat}: a vector containing the importance sampling estimates.
  \code{NSE}: a vector containing the numerical standard error of the components of \code{ghat}.
  \code{RNE}: a vector containing the relative numerical efficiency of the
  components of \code{ghat}.
}
\note{%
  Further details and examples of the \R package \code{AdMit}
  can be found in Ardia, Hoogerheide, van Dijk (2009a,b). See also
  the package vignette by typing \code{vignette("AdMit")}.

  Further information on importance sampling can be found
  in Geweke (1989) or Koop (2003).

  Please cite the package in publications. Use \code{citation("AdMit")}.
}
\author{David Ardia}
\references{%
  Ardia, D., Hoogerheide, L.F., van Dijk, H.K. (2009a).
  AdMit: Adaptive Mixture of Student-t Distributions.
  \emph{The R Journal} \bold{1}(1), pp.25-30.
  \url{https://journal.r-project.org/archive/2009-1/}

  Ardia, D., Hoogerheide, L.F., van Dijk, H.K. (2009b).
  Adaptive Mixture of Student-t Distributions as a Flexible Candidate
  Distribution for Efficient Simulation: The R Package AdMit.
  \emph{Journal of Statistical Software} \bold{29}(3), pp.1-32.
  \doi{10.18637/jss.v029.i03}

  Geweke, J.F. (1989).
  Bayesian Inference in Econometric Models Using Monte Carlo Integration.
  \emph{Econometrica} \bold{57}(6), pp.1317-1339. 

  Koop, G. (2003). 
  \emph{Bayesian Econometrics}. 
  Wiley-Interscience (London, UK). ISBN: 0470845678.
}
\seealso{%
  \code{\link{AdMit}} for fitting an adaptive mixture of Student-t
  distributions to a target density through its \code{KERNEL} function,
  \code{\link{AdMitMH}} for the independence chain Metropolis-Hastings
  algorithm using an adaptive mixture of Student-t distributions
  as the candidate density.
}
\examples{%
  ## NB : Low number of draws for speedup. Consider using more draws!
  ## Gelman and Meng (1991) kernel function
  GelmanMeng <- function(x, A = 1, B = 0, C1 = 3, C2 = 3, log = TRUE)
  {
    if (is.vector(x))
      x <- matrix(x, nrow = 1)
    r <- -.5 * (A * x[,1]^2 * x[,2]^2 + x[,1]^2 + x[,2]^2
              - 2 * B * x[,1] * x[,2] - 2 * C1 * x[,1] - 2 * C2 * x[,2])
    if (!log)
      r <- exp(r)
    as.vector(r)
  }

  ## Run the AdMit function to fit the mixture approximation
  set.seed(1234)
  outAdMit <- AdMit(KERNEL = GelmanMeng, 
                    mu0 = c(0.0, 0.1), control = list(Ns = 1e4))

  ## Use importance sampling with the mixture approximation as the
  ## importance density
  outAdMitIS <- AdMitIS(N = 1e4, KERNEL = GelmanMeng, mit = outAdMit$mit)
  print(outAdMitIS)
}
\keyword{htest}
