---
title: 'AdMit: Adaptive mixture of Student-t distributions in R'
bibliography: AdMit.bib
date: "21 December 2016"
tags:
- adaptive mixture
- Student-t distributions
- importance sampling
- independence chain Metropolis-Hasting algorithm
- R software
authors:
- affiliation: Institute of Financial Analysis - University of Neuchâtel, Switzerland
- affiliation: Département de Finance, Assurance et Immobilier - Université Laval, Québec, Canada
  name: David Ardia
  orcid: 0000-0003-2823-782X
- affiliation: Department of Econometrics - Vrije Universiteit Amsterdam, The Netherlands
  name: Lennart F. Hoogerheide
  orcid: null
- affiliation: Econometric Institute - Erasmus University Rotterdam, The Netherlands
  name: Herman K. van Dijk
  orcid: null
---

# Summary

`AdMit` is an R package (@R) which provides fexible functions to approximate a
certain target distribution and to effciently generate a sample of random draws from it, given only
a kernel of the target density function. The core algorithm fits an adaptive mixture of Student-t
distributions to the density of interest, and then, importance sampling or the independence chain
Metropolis-Hastings algorithm is used to obtain quantities of interest for the target density, using
the fitted mixture as the importance or candidate density. The estimation procedure is fully
automatic and thus avoids the time-consuming and diffiult task of tuning a sampling algorithm.
Full description of the algorithm and numerous application are available in @ArdiaEtAl2009a
and @ArdiaEtAl2009b. The latest version of the package is available at 'https://github.com/ArdiaD/AdMit'.

# References
