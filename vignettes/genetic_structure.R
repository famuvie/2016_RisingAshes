## This script returns:
##
## - the Inverse relationship matrix for the Devecey dataset (gen_str$Ainv)
## - the recoding map (gen_str$recoding)
## - sets-up the family-wise sum-to-zero constraint matrix (gen_str$constraint.matrix)
## - sets-up the prior distribution for the additive-genetic precision (prec.A)

gen_str <- with(Devecey, fam_additive(Seq, FAM))

## check
const.mat.811 <- gen_str$constraint.matrix[24, ] > 0  # (recoded) desc. from last mum
code.811.idx <- match(which(const.mat.811), gen_str$recoding)  # original codes
stopifnot(all(as.numeric(Devecey$FAM[code.811.idx]) == 23))  # all from last family


## Prior Precision for the additive precision
## this is: shape (\alpha) and rate (\beta) of the Gamma distribution
## such that mean = \alpha/\beta and variance = \alpha/\beta^2
prec.A <- list(param = c(.5, 5e-01), fixed = FALSE)

## Visual check
ggplot(transform(data.frame(x = seq(0, 10, length = 1001)),
                 y = with(prec.A,
                          dgamma(x, shape = param[1], rate = param[2]))),
       aes(x, y)) + geom_line()

## The variance follows an Inverse-gamma with shape \alpha and
## **scale** \beta, and thus has a mode at \beta/(\alpha+1) = 1/3
ggplot(transform(data.frame(x = seq(1e-03, .5, length = 1001)),
                 y = with(prec.A,
                          dgamma(1/x, shape = param[1], rate =
                                   param[2]))/x**2),
       aes(x, y)) + geom_line()
# Note that this precludes a priori additive variances close to zero
# this might be inconvenient.
# The 80% of the density is concentrated between .05 and 15:
diff(pgamma(1/c(15, .05), .5, .5))

