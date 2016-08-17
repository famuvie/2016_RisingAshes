
#### MODEL FIT ####

inla.setOption(inla.call = 'submit')  # Remote computing
# inla.setOption(inla.call = NULL)        # Local computing



## Basic logistic model with categorical covariates only
ml$bin_PROV$res$ref <- 
  inla(formula = CL0 ~ 0 + PROV + Year2013,
       data = ml$bin_PROV$dat,
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


## Reference model: without PROV

ml$bin_PROV$res$noPROV <- 
  inla(formula = ml$bin_PROV$fml,
       data = inla.stack.data(ml$bin_PROV$stk),
       family = "binomial",
       Ntrials = 1,
       ## Correction of the Laplace Approximation
       ## for Binomial data with many zeros
       ## http://arxiv.org/abs/1503.07307
       control.inla = list(
         correct = TRUE,
         correct.factor = 10),
       control.predictor = list(
         A = inla.stack.A(ml$bin_PROV$stk),
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


## PROV as a categorical covariate
ml$bin_PROV$res$cat <- 
  inla(formula = update(ml$bin_PROV$fml, ~ PROV.bin + Year2013.bin + . -Year.bin),
       data = inla.stack.data(ml$bin_PROV$stk),
       family = "binomial",
       Ntrials = 1,
       ## Correction of the Laplace Approximation
       ## for Binomial data with many zeros
       ## http://arxiv.org/abs/1503.07307
       control.inla = list(
         correct = TRUE,
         correct.factor = 10),
       control.predictor = list(
         A = inla.stack.A(ml$bin_PROV$stk),
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


## PROV*Year as a categorical covariate
ml$bin_PROV$res$catxyear <- 
  inla(formula = update(ml$bin_PROV$fml, ~ . -Year.bin + Year.bin:PROV.bin),
       data = inla.stack.data(ml$bin_PROV$stk),
       family = "binomial",
       Ntrials = 1,
       ## Correction of the Laplace Approximation
       ## for Binomial data with many zeros
       ## http://arxiv.org/abs/1503.07307
       control.inla = list(
         correct = TRUE,
         correct.factor = 10),
       control.predictor = list(
         A = inla.stack.A(ml$bin_PROV$stk),
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

