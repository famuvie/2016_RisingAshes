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
  Ainv <- compute_Ainverse(ped)

  ## Ainverse in sparse format
  Cmat <- with(Ainv,
               Matrix::sparseMatrix(i = Ainverse[,1],
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
#' @param coord a two-column matrix with coordinates. By default it takes the
#'   coordinates of the Devecey dataset.
#'
#' @return A list with the observation locations (loc) a prediction mesh and an
#'   SPDE object which includes the prior distributions for its hyperparameters
#' @import INLA
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


#' Include a suffix to variable names
#'
#' Used with \code{\link[INLA]{inla.stack}} for building effect stacks.
#'
#' @param x a data.frame
#' @param suffix character string to be appended to variable names
#' @param makeNA logical. If \code{TRUE}, returns the suffixed data.frame filled
#'   with \code{NA}.
#' @export
suffix <- function(x, suffix, makeNA = FALSE) {
  z <- structure(x, names = paste(names(x), suffix, sep = '.'))
  if(makeNA) z[] <- NA
  z
}


#' Perform a shift on a (list of) marginals
#'
#' Add a constant to a \code{INLA}'s marginal distribution or list of marginals.
#'
#' @param mar either a single marginal or a list of marginals
#' @param mu numeric. Constant to be added.
#'
#' \code{INLA}'s marginals are two-column numeric matrices with abscissas and
#' ordinates. This function adds a constant \code{mu} to the ordinates (i.e.,
#' second column).
#'
#' @export
shift_marginal <- function(mar, mu) {
  if (is.list(mar))
    mar <- lapply(mar, shift_marginal, mu)
  else
    mar[, 1] <- mar[, 1] + mu
  mar
}

#' Inverse logit function
#'
#' Compute the inverse logit
#'
#' @param x value(s) to be transformed
#'
#' \deqn{p = exp(x)/(1 + exp(x))}
#'
#' @export
inv.logit <- function(x) {
  p <- exp(x)/(1 + exp(x))
  p <- ifelse(is.na(p) & !is.na(x), 1, p)
}


#' Pedigree functions
#'
#' These functions are borrowed from the package \code{AnimalINLA}
#'
#' @param pedigree 3-col numerical matrix
#' @references http://www.r-inla.org/related-projects/animalinla
#' @name pedigree
NULL


#' @describeIn pedigree Check, recode and compute the inverse relationship matrix
#' @export
compute_Ainverse <- function (pedigree)
{
  # print("Checking pedigree...")
  # CheckPedigree = CheckPedigree(pedigree)
  # print("Sorting the pedigree chronologically...")
  sorted.pedigree = SortPedigree(pedigree)
  xx = RenamePedigree(sorted.pedigree)
  # print("Computing the A matrix...")
  Amatrix = Ainverse(xx$pedigree)
  ret = list(Ainverse = Amatrix, map = xx$map)
  class(ret) = "ped"
  return(ret)
}

#' @describeIn pedigree Sort a pedigree chronologically
SortPedigree <- function (pedigree)
{
  pedigree = as.matrix(pedigree)
  n <- dim(pedigree)[1]
  g <- rep(0, n)
  s <- rep(1, n)
  t <- matrix(c(NA, NA, NA), 1, 3)
  p <- pedigree
  continue = TRUE
  iteration = 0
  ID = pedigree[, 1]
  while (continue == TRUE) {
    iteration = iteration + 1
    for (i in 1:n) {
      if (!(p[i, 2]) %in% t[, 1]) {
        t <- rbind(t, p[p[, 1] == p[i, 2], ])
        ii = which(ID == p[i, 2])
        g[ii] <- g[ii] + 1
      }
      if (!(p[i, 3]) %in% t[, 1]) {
        t <- rbind(t, p[p[, 1] == p[i, 3], ])
        ii = which(ID == p[i, 3])
        g[ii] <- g[ii] + 1
      }
    }
    if (dim(t)[1] == 1) {
      continue = FALSE
    }
    else {
      nn = dim(t)[1]
      t = t[-1, ]
      p = t
      {
        if (nn == 2) {
          n = 1
          p = matrix(p, nrow = 1, ncol = 3)
        }
        else n <- dim(p)[1]
      }
      t = matrix(c(NA, NA, NA), 1, 3)
    }
  }
  oo = order(g, decreasing = TRUE)
  sorted.pedigree = pedigree[oo, ]
  return(sorted.pedigree)
}


#' @describeIn pedigree Recode a pedigree
RenamePedigree <- function (pedigree)
{
  n = dim(pedigree)[1]
  map = cbind(pedigree[, 1], 1:n)
  ped1 = pedigree
  for (i in 1:n) ped1[pedigree == map[i, 1]] = i
  red = list(pedigree = ped1, map = map)
  return(red)
}

#' @describeIn pedigree Compute the inverse relationship matrix
Ainverse <- function (pedigree)
{
  individ = pedigree[, 1]
  parent1 = pedigree[, 2]
  parent2 = pedigree[, 3]
  n = dim(pedigree)[1]
  u = numeric(n)
  v = numeric(n)
  p = pmin(parent1, parent2)
  q = pmax(parent1, parent2)
  Asave <- rep(NA, 3)
  for (k in 1:n) {
    if (p[k] == 0 & q[k] == 0)
      v[k] = 1
    if (p[k] == 0 & q[k] != 0)
      v[k] = sqrt(1 - 0.25 * u[q[k]])
    if (p[k] != 0 & q[k] != 0)
      v[k] = sqrt(1 - 0.25 * u[q[k]] - 0.25 * u[p[k]])
    if (k < n) {
      for (i in (k + 1):n) {
        if (p[i] >= k)
          v[i] = 0.5 * v[p[i]] + 0.5 * v[q[i]]
        else if (p[i] < k & q[i] >= k)
          v[i] = 0.5 * v[q[i]]
        else if (q[i] < k)
          v[i] = 0
        else stop("this should not happen")
      }
      u[k:n] = u[k:n] + v[k:n]^2
    }

    d = v[k]^(-2)
    if (p[k] > 0) {
      xx = cbind(c(k, p[k], q[k], p[k], p[k], q[k]), c(k,
                                                       k, k, p[k], q[k], q[k]), c(d, -0.5 * d, -0.5 *
                                                                                    d, 0.25 * d, 0.25 * d, 0.25 * d))
    }
    if (p[k] == 0 & q[k] > 0) {
      xx = cbind(c(k, q[k], q[k]), c(k, k, q[k]), c(d,
                                                    -0.5 * d, 0.25 * d))
    }
    if (q[k] == 0 & p[k] == 0) {
      xx = c(k, k, d)
    }
    Asave = rbind(Asave, xx)
  }
  Asave = Asave[-1, ]
  nA <- dim(Asave)[1]
  # print("Summing up.....")
  for (i in 1:(nA - 1)) {
    idx = i - 1 + which((Asave[i, 1] == Asave[i:nA, 1]) &
                          (Asave[i, 2] == Asave[i:nA, 2]))
    Asave[i, 3] = sum(Asave[idx, 3])
    Asave[idx[-1], 3] = NA
  }
  Asave = Asave[!is.na(Asave[, 3]), ]
  row.names(Asave) <- NULL
  return(Asave)
}
