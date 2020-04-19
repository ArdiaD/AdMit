# AdMit
[![CRAN](http://www.r-pkg.org/badges/version/AdMit)](https://cran.r-project.org/package=AdMit) 
[![Downloads](http://cranlogs.r-pkg.org/badges/AdMit?color=brightgreen)](http://www.r-pkg.org/pkg/AdMit)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/AdMit?color=brightgreen)](http://www.r-pkg.org/pkg/AdMit)
[![Pending Pull-Requests](http://githubbadges.herokuapp.com/ArdiaD/AdMit/pulls.svg?style=flat)](https://github.com/ArdiaD/AdMit/pulls)
[![Github Issues](http://githubbadges.herokuapp.com/ArdiaD/AdMit/issues.svg)](https://github.com/ArdiaD/AdMit/issues)

`AdMit` ([Ardia et al., 2009a](https://doi.org/10.18637/jss.v029.i03)) is an R package which provides 
flexible functions to approximate a certain target distribution and to efficiently generate a sample of 
random draws from it, given only a kernel of the target density function. The core 
algorithm fits an adaptive mixture of Student-t distributions to the density of interest, and then, 
importance sampling or the independence chain Metropolis-Hastings algorithm is used to obtain 
quantities of interest for the target density, using the fitted mixture as the importance or 
candidate density. The estimation procedure is fully automatic and thus avoids the 
time-consuming and difficult task of tuning a sampling algorithm.
Full description of the algorithm and numerous applications are available in [Ardia et al. (2009a)](https://doi.org/10.18637/jss.v029.i03) and [Ardia et al. (2009b)](https://doi.org/10.32614/RJ-2009-003).

## Please cite the package in publications!

By using `AdMit` you agree to the following rules: 

1) You must cite Ardia et al. (2009a) in working papers and published papers that use `AdMit`.
2) You must place the following URL in a footnote to help others find `AdMit`: https://CRAN.R-project.org/package=AdMit. 
3) You assume all risk for the use of `AdMit`.

Ardia, D., Hoogerheide, L., van Dijk, H.K. (2009a).    
Adaptive mixture of Student-t distributions as a flexible candidate 
distribution for efficient simulation: The R package AdMit.    
_Journal of Statistical Software_ **29**(3), pp.1-32.     
[https://doi.org/10.18637/jss.v029.i03](https://doi.org/10.18637/jss.v029.i03)  

Ardia, D., Hoogerheide, L., van Dijk, H.K. (2009b).    
AdMit: Adaptive mixture of Student-t distributions.   
_R Journal_ **1**(1), pp.25-30.     
[https://doi.org/10.32614/RJ-2009-003](https://doi.org/10.32614/RJ-2009-003)  
