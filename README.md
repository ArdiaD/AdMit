# AdMit
Adaptive mixture of Student-t distributions

[![DOI](https://zenodo.org/badge/59887530.svg)](https://zenodo.org/badge/latestdoi/59887530)
[![CRAN](http://www.r-pkg.org/badges/version/AdMit)](https://cran.r-project.org/package=AdMit) 
[![Downloads](http://cranlogs.r-pkg.org/badges/AdMit?color=brightgreen)](http://www.r-pkg.org/pkg/AdMit)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/AdMit?color=brightgreen)](http://www.r-pkg.org/pkg/AdMit)

`AdMit` ([Ardia et al., 2009a](http://dx.doi.org/10.18637/jss.v029.i03)) is an R package which provides 
flexible functions to approximate a certain target distribution and to efficiently generate a sample of 
random draws from it, given only a kernel of the target density function. The core 
algorithm fits an adaptive mixture of Student-t distributions to the density of interest, and then, 
importance sampling or the independence chain Metropolis-Hastings algorithm is used to obtain 
quantities of interest for the target density, using the fitted mixture as the importance or 
candidate density. The estimation procedure is fully automatic and thus avoids the 
time-consuming and difficult task of tuning a sampling algorithm.
Full description of the algorithm and numerous applications are available in [Ardia et al. (2009a)](http://dx.doi.org/10.18637/jss.v029.i03) and [Ardia et al. (2009b)](http://journal.r-project.org/archive/2009-1/).

Please cite `AdMit` in publications:

Ardia, D., Hoogerheide, L., van Dijk, H.K. (2009a)  
Adaptive mixture of Student-t distributions as a flexible candidate 
distribution for efficient simulation: The R package AdMit.  
_Journal of Statistical Software_ **29**(3), pp.1--32.   
http://dx.doi.org/10.18637/jss.v029.i03


Ardia, D., Hoogerheide, L., van Dijk, H.K. (2009b)  
AdMit: Adaptive mixture of Student-t distributions. 
_R Journal_ **1**(1), pp.25--30.   
http://journal.r-project.org/archive/2009-1/