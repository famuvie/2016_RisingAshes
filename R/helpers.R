#' Additive-genetic matrices for a family structure
#'
#' Computes the relevant matrices to set-up an additive-genetic effect given a
#' family structure
#'
#' @param id numeric vector of tree ids
#' @param fam factor with family codes
#'
#' @return a list with a \code{recoding} of the original \code{id} codes, the
#'   Inverse relationship matrix (\code{Ainv}) in sparse format, and the
#'   family-wise sum-to-zero \code{constraint.matrix} for use in conjuntion with
#'   an explicit family effect
#' @export
fam_additive <- function(id, fam) {

  ## Pedigree implicit in the family structure
  ped <- with(Devecey, famped(id, fam))

  ## Inverse relationship matrix
  ## Need the package AnimalINLA
  ## http://www.r-inla.org/related-projects/animalinla
  Ainv <- AnimalINLA::compute.Ainverse(ped)

  ## Ainverse in sparse format
  Cmat <- with(Ainv,
               sparseMatrix(i = Ainverse[,1],
                            j = Ainverse[,2],
                            x = Ainverse[,3]))


  ## Family-wise sum-to-zero constraint matrix
  # If explicit, the FAM effect is confounded with the genetic effect
  # So I need to use the linear constrains family-wise
  # Watch out: this need to be in the recoded space

  ## which individuals in the (recoded) pedigree have mum with code x
  whose_mum <- function(x)
    as.integer(ped$mum[order(match(ped$self, Ainv$map[, 1]))] == x)
  const.mat <- t(sapply(
    sort(unique(ped$mum)), whose_mum))

  return(
    list(
      recoding = match(id, Ainv$map[, 1]),
      Ainv     = Cmat,
      constraint.matrix = const.mat
    )
  )
}


#' Pedigree for a family structure
#'
#' Compute the pedigree from tree id and family codes
#'
#' Note that renumbering might take place.
#'
#' @param id numeric vector of tree ids
#' @param fam factor with family codes
#'
#' @return A three-column data.frame with numeric codes for the individual trees
#'   and their progenitors
#'
#' @export
famped <- function(id, fam) {
  ## Build the pedigree
  ## Members of the same family descend from the same mother
  ## Watch out! renumbering might take place!!
  n.fam <- nlevels(fam)
  ped <- data.frame(self = c(id, max(id) + 1:n.fam),
                    mum  = c(max(id) + as.numeric(fam), rep(0, n.fam)),
                    dad  = 0)
  ped$mum[is.na(ped$mum)] <- 0

  return(ped)
}


#' SPDE structure for the Devecey dataset
#'
#' @return A list with the observation locations (loc) a prediction mesh and an
#'   SPDE object which includes the prior distributions for its hyperparameters
#' @export
Devecey_spde <- function(coord = Devecey[,c('X', 'Y')]) {
  loc <- as.matrix(unique(coord))
  rownames(loc) <- NULL
  mesh <- inla.mesh.2d(loc = loc, offset = c(2, 5), max.edge = 5)

  # plot(mesh)
  size = min(apply(mesh$loc[,1:2], 2, function(x) diff(range(x))))
  sigma0 = 0.5 # prior median sd
  range0 = size/2 # prior median range
  kappa0 = sqrt(8)/range0
  tau0   = 1/(sqrt(4*pi)*kappa0*sigma0)



  ## Code for tuning the prior range
  # ## the range can be between 5 and size = 38, so:
  # sqrt(8)/c(5, size)  # kappa0 \in (.566, .074), and
  # log(sqrt(8)/c(5, size))  # theta2 \in (-.57, -2.60), with median
  # log(kappa0)   # -1.9.
  # ## This gives a semi-support of
  # max(abs(log(sqrt(8)/c(5, size)) - log(kappa0) )) # 1.3 units in the log scale
  # ## so we use a gaussian prior for theta2 with mean log(kappa0)
  # ## and sd = 1.3/3 (as 3 sd is the effective support of a Normal dist.)
  # ## Being generous, we use an sd = 1
  rho0.mar <-
    inla.tmarginal(function(x) sqrt(8)/exp(x),
                   transform(data.frame(x = log(kappa0) + seq(-1.5, 5,
                                                              length=2021)),
                             y = dnorm(x, mean = log(kappa0),
                                       sd = 1)))
  # plot(rho0.mar, type = 'l', xlab = 'distance', ylab = '')

  sigma20.mar <-
    inla.tmarginal(function(x) 1/(4*pi*kappa0^2*exp(2*x)),
                   transform(data.frame(x = log(tau0) + seq(-0.5, 5,
                                                            length=2021)),
                             y = dnorm(x, mean = log(tau0),
                                       sd = 1)))
  # plot(sigma20.mar, type = 'l', xlab = 'variance', ylab = '')

  sigma_gt1 <- function(x, kappa) 1/sqrt(4*pi)/kappa/exp(x)
  # sigma_gt1(log(tau0)+3, kappa0)

  spde = inla.spde2.matern(
    mesh,
    B.tau = cbind(log(tau0), -1, 1),
    B.kappa = cbind(log(kappa0), 0, -1),
    theta.prior.mean = c(0, 0),
    theta.prior.prec = c(1, 1),
    constr = TRUE
  )

  return(list(loc = loc, mesh = mesh, spde = spde,
              priors = list(rho = rho0.mar,
                            sigma = sigma20.mar)))
}


# Plotting fitted vs. obs.
#' @export
plot.fitvsobs <- function(x, y, color = NULL, shape = NULL, data, ...) {
  eval(substitute(
    ggplot(data,
           aes(x, y, color = color, shape = shape)) +
      geom_jitter(...) +
      geom_abline(int = 0, sl = 1, col = 'darkgray') +
      guides(shape = 'none',
             color = guide_legend(keyheight = 0.5))
  ))
}


#' @export
sufix <- function(x, sufix, makeNA = FALSE) {
  z <- structure(x, names = paste(names(x), sufix, sep = '.'))
  if(makeNA) z[] <- NA
  z
}


## perform a shift on a (list of) marginals
#' @export
shift_marginal <- function(mar, mu) {
  if (is.list(mar))
    mar <- lapply(mar, shift_marginal, mu)
  else
    mar[, 1] <- mar[, 1] + mu
  mar
}
