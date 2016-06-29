# RisingAshes

Data and code to the to reproduce the results in [Muñoz et al.
(2016)](http://dx.doi.org/10.1094/PHYTO-11-15-0284-R) [1] ([bioRXiv
preprint](http://dx.doi.org/10.1101/031393)), which analyzes the genetic
variation in susceptibility to a pathogen using a spatio-temporal hierarchical
Bayesian model implemented with [INLA](http://www.r-inla.org/).

This package contains three vignettes with all the code necessary to reproduce
the results from the paper (and more results that did not make it to the paper).
You can browse the code 
[online](https://github.com/famuvie/RisingAshes/tree/master/inst/doc).

To browse offline or to play with the data yourself, install the package locally
as follows.

```r 
devtools:::install_github('famuvie/RisingAshes') 
library(RisingAshes) 
browseVignettes(package = 'RisingAshes') 
data(Devecey) 
```

Note that the package includes some [helper functions](R/helpers.R) in a
separate file.

Since the computations can take several hours, all the expensive results are 
lazy-loaded with the package, so that you can explore the results without the 
need of recomputing everything. The code for producing the results is included 
in the vignettes, but is not executed by default at compilation time. For 
reproducing the results you might want to setup the computation on a server.

[1] F. Muñoz, B. Marçais, J. Dufour and A. Dowkiw (2016). Rising Out of the 
Ashes: Additive Genetic Variation for Crown and Collar Resistance to 
Hymenoscyphus fraxineus in Fraxinus excelsior. Phytopathology. In Press. DOI: 
10.1094/PHYTO-11-15-0284-R. A preprint is available at
[bioRXiv](http://dx.doi.org/10.1101/031393).
