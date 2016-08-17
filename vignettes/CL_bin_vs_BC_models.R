
#### MODEL FIT ####

inla.setOption(inla.call = 'submit')  # Remote computing
# inla.setOption(inla.call = NULL)        # Local computing


## Basic logistic model with categorical covariates only
ml$bin_BC$res$ref <-
  inla(formula = CL0 ~ 0 + BC.cat + Year2013,
       data = transform(ml$bin_BC$dat,
                        Year2013 = as.numeric(Year == 2013)),
       family = "binomial",
       Ntrials = 1,
       control.fixed = list(
         expand.factor.strategy = 'inla'), # account for NAs in factor
       control.compute = list(
         dic = TRUE,
         waic = TRUE,
         cpo = TRUE,
         po = TRUE,
         config = TRUE)
  )


## Reference model: without BC

ml$bin_BC$res$noBC <-
  inla(formula = ml$bin_BC$fml,
       data = inla.stack.data(ml$bin_BC$stk),
       family = "binomial",
       Ntrials = 1,
       ## Correction of the Laplace Approximation
       ## for Binomial data with many zeros
       ## http://arxiv.org/abs/1503.07307
       control.inla = list(
         correct = TRUE,
         correct.factor = 10),
       control.predictor = list(
         A = inla.stack.A(ml$bin_BC$stk),
         compute = TRUE, link = 1),
       control.fixed = list(
         expand.factor.strategy = 'inla'), # account for NAs in factor
       control.compute = list(
         dic = TRUE,
         waic = TRUE,
         cpo = TRUE,
         po = TRUE,
         config = TRUE)
  )


## BC as a categorical covariate
ml$bin_BC$res$cat <-
  inla(formula = update(ml$bin_BC$fml, ~ BC.cat.bin + Year2013.bin + . -Year.bin),
       data = inla.stack.data(ml$bin_BC$stk),
       family = "binomial",
       Ntrials = 1,
       ## Correction of the Laplace Approximation
       ## for Binomial data with many zeros
       ## http://arxiv.org/abs/1503.07307
       control.inla = list(
         correct = TRUE,
         correct.factor = 10),
       control.predictor = list(
         A = inla.stack.A(ml$bin_BC$stk),
         compute = TRUE, link = 1),
       control.fixed = list(
         expand.factor.strategy = 'inla'), # account for NAs in factor
       control.compute = list(
         dic = TRUE,
         waic = TRUE,
         cpo = TRUE,
         po = TRUE,
         config = TRUE)
  )


## BC2: only two categories of BC
ml$bin_BC$res$cat2 <-
  inla(formula = update(ml$bin_BC$fml, ~ . -Year.bin + cut(BC.bin, breaks = c(0, 30, 150)):Year.bin),
       data = inla.stack.data(ml$bin_BC$stk),
       family = "binomial",
       Ntrials = 1,
       ## Correction of the Laplace Approximation
       ## for Binomial data with many zeros
       ## http://arxiv.org/abs/1503.07307
       control.inla = list(
         correct = TRUE,
         correct.factor = 10),
       control.predictor = list(
         A = inla.stack.A(ml$bin_BC$stk),
         compute = TRUE, link = 1),
       control.fixed = list(
         expand.factor.strategy = 'inla'), # account for NAs in factor
       control.compute = list(
         dic = TRUE,
         waic = TRUE,
         cpo = TRUE,
         po = TRUE,
         config = TRUE)
  )


## BC*Year as a categorical covariate
ml$bin_BC$res$catxyear <-
  inla(formula = update(ml$bin_BC$fml, ~ . -Year.bin + Year.bin:BC.cat.bin),
       data = inla.stack.data(ml$bin_BC$stk),
       family = "binomial",
       Ntrials = 1,
       ## Correction of the Laplace Approximation
       ## for Binomial data with many zeros
       ## http://arxiv.org/abs/1503.07307
       control.inla = list(
         correct = TRUE,
         correct.factor = 10),
       control.predictor = list(
         A = inla.stack.A(ml$bin_BC$stk),
         compute = TRUE, link = 1),
       control.fixed = list(
         expand.factor.strategy = 'inla'), # account for NAs in factor
       control.compute = list(
         dic = TRUE,
         waic = TRUE,
         cpo = TRUE,
         po = TRUE,
         config = TRUE)
  )

## Linear regression on the scaled BC

## For some reason, I needed to get rid of NAs in BC for this model to run
## For identifiability, I need to code the Year with respect to the baseline of 2012

ml$bin_BC$res$lin <-
  inla(formula = update(ml$bin_BC$fml, ~ BC.sc.bin + . ),
       data = inla.stack.data(ml$bin_BC$stk),
       family = "binomial",
       Ntrials = 1,
       ## Correction of the Laplace Approximation
       ## for Binomial data with many zeros
       ## http://arxiv.org/abs/1503.07307
       control.inla = list(
         correct = TRUE,
         correct.factor = 10),
       control.predictor = list(
         A = inla.stack.A(ml$bin_BC$stk),
         compute = TRUE, link = 1),
       control.fixed = list(
         expand.factor.strategy = 'inla'), # account for NAs in factor
       control.compute = list(
         dic = TRUE,
         waic = TRUE,
         cpo = TRUE,
         po = TRUE,
         config = TRUE)
  )

## Quadratic regression on the scaled BC

ml$bin_BC$res$quad <-
  inla(formula = update(ml$bin_BC$fml, ~ BC.sc.bin + I(BC.sc.bin^2) + .),
       data = inla.stack.data(ml$bin_BC$stk),
       family = "binomial",
       Ntrials = 1,
       ## Correction of the Laplace Approximation
       ## for Binomial data with many zeros
       ## http://arxiv.org/abs/1503.07307
       control.inla = list(
         correct = TRUE,
         correct.factor = 10),
       control.predictor = list(
         A = inla.stack.A(ml$bin_BC$stk),
         compute = TRUE, link = 1),
       control.fixed = list(
         expand.factor.strategy = 'inla'), # account for NAs in factor
       control.compute = list(
         dic = TRUE,
         waic = TRUE,
         cpo = TRUE,
         po = TRUE,
         config = TRUE)
  )



## Non-parametric regression on the scaled BC
ml$bin_BC$res$np <-
  inla(formula = update(ml$bin_BC$fml,
                        ~ f(cut(BC.sc.bin, 51),
                            model = 'rw2',
                            scale.model = TRUE,
                            hyper = list(theta = prec.A)) +
                          .),
       data = inla.stack.data(ml$bin_BC$stk),
       family = "binomial",
       Ntrials = 1,
       ## Correction of the Laplace Approximation
       ## for Binomial data with many zeros
       ## http://arxiv.org/abs/1503.07307
       control.inla = list(
         correct = TRUE,
         correct.factor = 10),
       control.predictor = list(
         A = inla.stack.A(ml$bin_BC$stk),
         compute = TRUE, link = 1),
       control.fixed = list(
         expand.factor.strategy = 'inla'), # account for NAs in factor
       control.compute = list(
         dic = TRUE,
         waic = TRUE,
         cpo = TRUE,
         po = TRUE,
         config = TRUE)
  )

