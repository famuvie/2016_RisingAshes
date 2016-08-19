# RisingAshes

Data and code to the to reproduce the results in [Muñoz et al. 
(2016)](http://dx.doi.org/10.1094/PHYTO-11-15-0284-R) [1], which analyzes the
genetic variation in susceptibility to a pathogen using a spatio-temporal
hierarchical Bayesian model implemented with [INLA](http://www.r-inla.org/).

This package contains three vignettes with all the code necessary to reproduce
the results from the paper (and more results that did not make it to the paper).

You can use this repository in one of several ways:

1. You can peek at the code 
[online](https://github.com/famuvie/RisingAshes/tree/master/inst/doc).

2. To browse offline or to explore the data yourself, install the package
locally
as follows.

    ```r 
    if(!require(devtools)) install.packages('devtools')
    devtools:::install_github('famuvie/2016_RisingAshes') 
    library(RisingAshes) 
    browseVignettes('RisingAshes') 
    data(Devecey) 
    str(Devecey)
    ```

    Note that the package includes some [helper functions](R/helpers.R) in a
    separate file.

3. To reproduce the vignettes or to explore the results, download the cached
results file and then reload the package

    ```r
    url <- released file
    destfile <- system.file('R', package = 'RisingAshes')
    download.file(url, destfile)
    library(RisingAshes)
    ```
  Since the computations can take several hours, all the expensive results are 
  lazy-loaded with the package (184 Mb file), so that you can explore the results
  without the need of recomputing everything. The code for producing the results
  is included in the vignettes, but is not executed by default at compilation
  time. 
  
4. If you want to fit the models and fully reproducing the results for yourself,
install `INLA` with:

    ```r
    install.packages('INLA', repos = 'http://www.math.ntnu.no/inla/R/testing')
    ```
    
    For refitting the models you might want to setup the computation on a
    server.


[1] **Muñoz, F., Marçais, B., Dufour, J., and Dowkiw, A. 2016**. Rising Out of 
the Ashes: Additive Genetic Variation for Crown and Collar Resistance to 
Hymenoscyphus fraxineus in Fraxinus excelsior. *Phytopathology* 106:1-10. DOI: 
[10.1094/PHYTO-11-15-0284-R](http://dx.doi.org/10.1094/PHYTO-11-15-0284-R). A
preprint is available at [bioRXiv](http://dx.doi.org/10.1101/031393).
