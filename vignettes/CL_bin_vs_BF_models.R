
#### MODEL FIT ####

inla.setOption(inla.call = 'submit')  # Remote computing
# inla.setOption(inla.call = NULL)        # Local computing


## Basic logistic model with categorical covariates only
ml$bin_BF$res$ref <- 
  inla(formula = CL0 ~ 0 + BF98 + Year2013,
       data = ml$bin_BF$dat,
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


## Reference model: without BF

ml$bin_BF$res$noBF <- 
  inla(formula = ml$bin_BF$fml,
       data = inla.stack.data(ml$bin_BF$stk),
       family = "binomial",
       Ntrials = 1,
       ## Correction of the Laplace Approximation
       ## for Binomial data with many zeros
       ## http://arxiv.org/abs/1503.07307
       control.inla = list(
         correct = TRUE,
         correct.factor = 10),
       control.predictor = list(
         A = inla.stack.A(ml$bin_BF$stk),
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


## BF as a categorical covariate
ml$bin_BF$res$cat <- 
  inla(formula = update(ml$bin_BF$fml, ~ BF98.bin + Year2013.bin + . -Year.bin),
       data = inla.stack.data(ml$bin_BF$stk),
       family = "binomial",
       Ntrials = 1,
       ## Correction of the Laplace Approximation
       ## for Binomial data with many zeros
       ## http://arxiv.org/abs/1503.07307
       control.inla = list(
         correct = TRUE,
         correct.factor = 10),
       control.predictor = list(
         A = inla.stack.A(ml$bin_BF$stk),
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


## BF*Year as a categorical covariate
ml$bin_BF$res$catxyear <- 
  inla(formula = update(ml$bin_BF$fml, ~ . -Year.bin + Year.bin:BF98.bin),
       data = inla.stack.data(ml$bin_BF$stk),
       family = "binomial",
       Ntrials = 1,
       ## Correction of the Laplace Approximation
       ## for Binomial data with many zeros
       ## http://arxiv.org/abs/1503.07307
       control.inla = list(
         correct = TRUE,
         correct.factor = 10),
       control.predictor = list(
         A = inla.stack.A(ml$bin_BF$stk),
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

