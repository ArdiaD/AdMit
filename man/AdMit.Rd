\name{AdMit}
\alias{AdMit}
\concept{Mixture of Student-t Distributions}
\concept{Adaptive mixture}
\title{Adaptive Mixture of Student-t Distributions}
\description{
  Function which performs the fitting of an adaptive mixture of
  Student-t distributions to approximate a target density through its
  kernel function}
\usage{
AdMit(KERNEL, mu0, Sigma0 = NULL, control = list(), ...)
}
\arguments{%
  \item{KERNEL}{kernel function of the target density on which the adaptive mixture is fitted. This
    function should be vectorized for speed purposes (i.e., its first
    argument should be a matrix and its output a vector). Moreover, the function must contain
    the logical argument \code{log}. If \code{log = TRUE}, the function
    returns (natural) logarithm values of the kernel function. \code{NA} and
    \code{NaN} values are not allowed. (See *Details* for examples
    of \code{KERNEL} implementation.)}
  \item{mu0}{initial value in the first stage optimization (for the location of
    the first Student-t component) in the adaptive mixture, or
    location of the first Student-t component if \code{Sigma0} is not \code{NULL}.}
  \item{Sigma0}{scale matrix of the first Student-t component (square, symmetric and positive definite). Default:
    \code{Sigma0 = NULL}, i.e., the scale matrix of the first Student-t
    component is estimated by the function \code{AdMit}.}
  \item{control}{control parameters (see *Details*).}
  \item{\ldots}{further arguments to be passed to \code{KERNEL}.}
}
\details{%
  The argument \code{KERNEL} is the kernel function of the target
  density, and it should be vectorized for speed purposes.

  As a first example, consider the kernel function proposed by Gelman-Meng (1991):
  \deqn{
    k(x_1,x_2) = \exp\left( -\frac{1}{2} [A x_1^2 x_2^2 + x_1^2 + x_2^2
    - 2 B x_1 x_2 - 2 C_1 x_1 - 2 C_2 x_2] \right)
  }{
    k(x1,x2) = exp( 0.5*[A*x1^2*x2^2 + x1^2 + x2^2
    - 2*B*x1*x2 - 2*C1*x1 - 2*C2*x2] )
  }
  where commonly used values
  are \eqn{A=1}, \eqn{B=0}, \eqn{C_1=3}{C1=3} and \eqn{C_2=3}{C2=3}.

  A vectorized implementation of this function might be:
  \preformatted{
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
  }
  This way, we may supply a point \eqn{(x_1,x_2)}{(x1,x2)}
  for \code{x} and the function will output a single value (i.e., the kernel
  estimated at this point). But the function is vectorized, in the sense
  that we may supply a \eqn{(N \times 2)}{Nx2} matrix
  of values for \code{x}, where rows of \code{x} are
  points \eqn{(x_1,x_2)}{(x1,x2)} and the output will be a vector of
  length \eqn{N}, containing the kernel values for these points.
  Since the \code{AdMit} procedure evaluates \code{KERNEL} for a
  large number of points, a vectorized implementation is important. Note
  also the additional argument \code{log = TRUE} which is used for
  numerical stability.

  As a second example, consider the following (simple) econometric model:
  \deqn{
    y_t \sim \, i.i.d. \, N(\mu,\sigma^2) \quad t=1,\ldots,T
  }{
    y_t ~ i.i.d. N(mu,sigma^2) t=1,...,T
  }
  where \eqn{\mu}{mu} is the mean value and \eqn{\sigma}{sigma} is the
  standard deviation. Our purpose is to estimate
  \eqn{\theta = (\mu,\sigma)}{theta=(mu,sigma)} within a Bayesian
  framework, based on a vector \eqn{y} of \eqn{T} observations; the
  kernel thus consists of the product of the
  prior and the likelihood function. As previously mentioned, the
  kernel function should be vectorized, i.e., treat a \eqn{(N \times 2)}{(Nx2)} matrix of points
  \eqn{\theta}{theta} for which the kernel will be evaluated.
  Using the common (Jeffreys) prior \eqn{p(\theta) = \frac{1}{\sigma}}{p(theta)=1/sigma}
  for \eqn{\sigma > 0}{sigma>0}, a vectorized implementation of the
  kernel function might be:
   \preformatted{
     KERNEL <- function(theta, y, log = TRUE)
     {
       if (is.vector(theta))
         theta <- matrix(theta, nrow = 1)

       ## sub function which returns the log-kernel for a given
       ## thetai value (i.e., a given row of theta)
       KERNEL_sub <- function(thetai)
       {
         if (thetai[2] > 0) ## check if sigma>0
	 { ## if yes, compute the log-kernel at thetai
           r <- - log(thetai[2])
	         + sum(dnorm(y, thetai[1], thetai[2], TRUE))
	 }
	 else
	 { ## if no, returns -Infinity
	   r <- -Inf
	 }
	 as.numeric(r)
       }

       ## 'apply' on the rows of theta (faster than a for loop)
       r <- apply(theta, 1, KERNEL_sub)
       if (!log)
         r <- exp(r)
       as.numeric(r)
     }
   }

   Since this kernel function also depends on the vector \eqn{y}, it
   must be passed to \code{KERNEL} in the \code{AdMit} function. This is
   achieved via the argument \eqn{\ldots}{...}, i.e., \code{AdMit(KERNEL, mu = c(0, 1), y = y)}.
   
   To gain even more speed, implementation of \code{KERNEL} might rely on C or Fortran
   code using the functions \code{.C} and \code{.Fortran}. An example is
   provided in the file \file{AdMitJSS.R} in the package's folder.

  The argument \code{control} is a list that can supply any of
  the following components:
  
  \describe{%
    \item{\code{Ns}}{number of draws used in the evaluation of the
      importance sampling weights (integer number not smaller than 100). Default: \code{Ns = 1e5}.}
    \item{\code{Np}}{number of draws used in the optimization of the mixing
      probabilities (integer number not smaller than 100 and not larger
      than \code{Ns}). Default: \code{Np = 1e3}.}
    \item{\code{Hmax}}{maximum number of Student-t components in the
      adaptive mixture (integer number not smaller than one). Default: \code{Hmax = 10}.}
    \item{\code{df}}{degrees of freedom parameter of the
      Student-t components (real number not smaller than one). Default: \code{df = 1}.}
    \item{\code{CVtol}}{tolerance for the relative change of the coefficient of
      variation (real number in [0,1]). The
      adaptive algorithm stops if the new
      component leads to a relative change in the coefficient of
      variation that is smaller or equal than
      \code{CVtol}. Default: \code{CVtol = 0.1}, i.e., 10\%.}
    \item{\code{weightNC}}{weight assigned to the new
      Student-t component of the adaptive mixture as
      a starting value in the optimization of the mixing probabilities
      (real number in [0,1]). Default: \code{weightNC = 0.1}, i.e., 10\%.}
    \item{\code{trace}}{tracing information on
      the adaptive fitting procedure (logical). Default:
      \code{trace = FALSE}, i.e., no tracing information.}
    \item{\code{IS}}{should importance sampling be used to estimate the
      mode and the scale matrix of the Student-t components (logical). Default: \code{IS = FALSE},
      i.e., use numerical optimization instead.}
    \item{\code{ISpercent}}{vector of percentage(s) of largest weights used to
      estimate the mode and the scale matrix of the Student-t
      components of the adaptive mixture by importance
      sampling (real number(s) in [0,1]). Default:
      \code{ISpercent = c(0.05, 0.15, 0.3)}, i.e., 5\%, 15\% and 30\%.}
    \item{\code{ISscale}}{vector of scaling factor(s) used to rescale the
      scale matrix of the mixture components (real positive numbers).
      Default: \code{ISscale = c(1, 0.25, 4)}.}
    \item{\code{trace.mu}}{Tracing information on
      the progress in the optimization of the mode of the mixture
      components (non-negative integer number). Higher values
      may produce more tracing information (see the source code
      of the function \code{\link{optim}} for further details).
      Default: \code{trace.mu = 0}, i.e., no tracing information.}
    \item{\code{maxit.mu}}{maximum number of iterations used
      in the optimization of the modes of the mixture components
      (positive integer). Default: \code{maxit.mu = 500}.}
    \item{\code{reltol.mu}}{relative convergence tolerance
      used in the optimization of the modes of the mixture components
      (real number in [0,1]). Default: \code{reltol.mu = 1e-8}.}
    \item{\code{trace.p}, \code{maxit.p}, \code{reltol.p}}{the same as for
      the arguments above, but for the optimization of the mixing
      probabilities of the mixture components.}
  }
}
\value{%
  A list with the following components:\cr

  \code{CV}: vector (of length \eqn{H}) of coefficients of variation of
  the importance sampling weights.\cr

  \code{mit}: list (of length 4) containing information on the fitted mixture of
  Student-t distributions, with the following components:\cr
  
  \code{p}: vector (of length \eqn{H}) of mixing probabilities.
  \code{mu}: matrix (of size \eqn{H \times d}{Hxd}) containing the
  vectors of modes (in row) of the mixture components.
  \code{Sigma}: matrix (of size \eqn{H \times d^2}{Hxd*d}) containing the scale
  matrices (in row) of the mixture components.
  \code{df}: degrees of freedom parameter of the Student-t components.\cr
  
  where \eqn{H (\geq 1)}{H (>=1)} is the number of components in the adaptive
  mixture of Student-t distributions and \eqn{d (\geq 1)}{d (>=1)} is
  the dimension of the first argument in \code{KERNEL}.\cr
  
  \code{summary}: data frame containing information on the optimization
  procedures. It returns for each component of the adaptive mixture of
  Student-t distribution: 1. the method used to estimate the mode
  and the scale matrix of the Student-t component (`USER' if \code{Sigma0} is
  provided by the user; numerical optimization: `BFGS', `Nelder-Mead';
  importance sampling: `IS', with percentage(s) of importance weights
  used and scaling factor(s)); 2. the time required for this optimization;
  3. the method used to estimate the mixing probabilities
  (`NLMINB', `BFGS', `Nelder-Mead', `NONE'); 4. the time required for this
  optimization; 5. the coefficient of variation of the importance
  sampling weights.  
}
\note{%
  Further details and examples of the \R package \code{AdMit}
  can be found in Ardia, Hoogerheide, van Dijk (2009a,b). See also
  the package vignette by typing \code{vignette("AdMit")} and the
  files \file{AdMitJSS.txt} and \file{AdMitRnews.txt} in the \file{/doc} package's folder.

  Further details on the core algorithm are given in
  Hoogerheide (2006), Hoogerheide, Kaashoek, van Dijk (2007) and
  Hoogerheide, van Dijk (2008).

  The adaptive mixture \code{mit} returned by the function \code{AdMit} is used by the
  function \code{\link{AdMitIS}} to perform importance sampling using
  \code{mit} as the importance density or by the function \code{\link{AdMitMH}} to perform
  independence chain Metropolis-Hastings sampling using \code{mit} as the
  candidate density.

  Please cite the package in publications. Use \code{citation("AdMit")}.
}
\author{%
  David Ardia for the \R port,
  Lennart F. Hoogerheide and Herman K. van Dijk for the \code{AdMit} algorithm.
}
\references{%
  Ardia, D., Hoogerheide, L.F., van Dijk, H.K. (2009a).
  AdMit: Adaptive Mixture of Student-t Distributions.
  \emph{The R Journal} \bold{1}(1), pp.25--30.
  \url{http://journal.r-project.org/2009-1/}

  Ardia, D., Hoogerheide, L.F., van Dijk, H.K. (2009b).
  Adaptive Mixture of Student-t Distributions as a Flexible Candidate
  Distribution for Efficient Simulation: The R Package AdMit.
  \emph{Journal of Statistical Software} \bold{29}(3), pp.1--32.
  \url{http://www.jstatsoft.org/v29/i03/}

  Gelman, A., Meng, X.-L. (1991). 
  A Note on Bivariate Distributions That Are Conditionally Normal. 
  \emph{The American Statistician} \bold{45}(2), pp.125--126.

  Hoogerheide, L.F. (2006). 
  \emph{Essays on Neural Network Sampling Methods and Instrumental Variables}. 
  PhD thesis, Tinbergen Institute, Erasmus University Rotterdam (NL). 
  ISBN: 9051708261. 
  (Book nr. 379 of the Tinbergen Institute Research Series.)

  Hoogerheide, L.F., Kaashoek, J.F., van Dijk, H.K. (2007). 
  On the Shape of Posterior Densities and Credible Sets in Instrumental Variable Regression Models with Reduced
  Rank: An Application of Flexible Sampling Methods using Neural Networks.
  \emph{Journal of Econometrics} \bold{139}(1), pp.154--180. 
  doi: 10.1016/j.jeconom.2006.06.009.

  Hoogerheide, L.F., van Dijk, H.K. (2008). 
  Possibly Ill-Behaved Posteriors in Econometric Models: On the Connection between Model
  Structures, Non-elliptical Credible Sets and Neural Network
  Simulation Techniques. 
  \emph{Tinbergen Institute discussion paper} \bold{2008-036/4}. 
  \url{http://www.tinbergen.nl/discussionpapers/08036.pdf}
}
\seealso{%
  \code{\link{AdMitIS}} for importance sampling using an
  adaptive mixture of Student-t distributions as the importance density,
  \code{\link{AdMitMH}} for the independence chain Metropolis-Hastings
  algorithm using an adaptive mixture of Student-t distributions as
  the candidate density.
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

  ## Run AdMit (with default values)
  set.seed(1234)
  outAdMit <- AdMit(GelmanMeng, mu0 = c(0.0, 0.1), control = list(Ns = 1e4))
  print(outAdMit)

  ## Run AdMit (using importance sampling to estimate
  ## the modes and the scale matrices)
  set.seed(1234)
  outAdMit <- AdMit(KERNEL = GelmanMeng, 
                    mu0 = c(0.0, 0.1), 
                    control = list(IS = TRUE, Ns = 1e4))
  print(outAdMit)
}
\keyword{htest}
