#' RisingAshes
#'
#' Data and code to reproduce the results in Muñoz et al. (2016), which analyzes
#' the genetic variation in susceptibility to a pathogen using a spatio-temporal
#' hierarchical Bayesian model.
#'
#' This package contains three vignettes with all the code necessary to
#' reproduce the results from the paper (and more results that did not make it
#' to the paper). You can browse the code with \code{browseVignettes(package =
#' 'RisingAshes')} and work with the data with \code{data(Devecey) }
#'
#' Since the computations can take several hours, all the expensive results are
#' lazy-loaded with the package, so that you can explore the results without the
#' need of recomputing everything. The code for producing the results is
#' included in the vignettes, but is not executed by default at compilation
#' time. For reproducing the results you might want to setup the computation on
#' a server.

#' @references F. Muñoz, B. Marçais, J. Dufour and A. Dowkiw (2016). Rising Out
#'   of the Ashes: Additive Genetic Variation for Crown and Collar Resistance to
#'   Hymenoscyphus fraxineus in Fraxinus excelsior. Phytopathology.
#'   http://dx.doi.org/10.1094/PHYTO-11-15-0284-R
#'
#'   Preprint bioRXiv http://dx.doi.org/10.1101/031393
#' @name RisingAshes
#' @docType package
NULL


#' Devecey
#'
#' Forest field trial consisting on 23 open-pollinated common ash (Fraxinus
#' excelsior) maternal progenies from 3 North-Eastern French provenances.
#'
#' Variables are:
#' \describe{
#'  \item{Seq}{tree id}
#'  \item{X, Y}{coordinates (row/cols, spaced 4 m x 4 m)}
#'  \item{BLC}{(randomized incomplete) block}
#'  \item{SBLC}{sub-block}
#'  \item{FAM}{OP Family (23 levels)}
#'  \item{PROV}{Provenance (3 levels)}
#'  \item{BF98}{Bud-Flush precocity ranking (1 - late; 5 - early)}
#'  \item{CD10..CD14}{Crown Dieback proportion (classes 0, 0.05, 0.30, 0.65, 0.90, 1) in 2010..2014}
#'  \item{CL12..CL13}{Collar Lesion proportion (continuous measure from 0 to 1) in 2012 and 2013}
#'  \item{BC12..BC13}{Basal Circumference (in cm) in 2012 and 2013}
#' }
#'
#' @references F. Muñoz, B. Marçais, J. Dufour and A. Dowkiw (2016). Rising Out
#'   of the Ashes: Additive Genetic Variation for Crown and Collar Resistance to
#'   Hymenoscyphus fraxineus in Fraxinus excelsior. Phytopathology.
#'   http://dx.doi.org/10.1094/PHYTO-11-15-0284-R
#'
#'   Preprint bioRXiv http://dx.doi.org/10.1101/031393
#' @name Devecey
#' @docType data
NULL
