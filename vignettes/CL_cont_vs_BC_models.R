
#### MODEL FIT ####

inla.setOption(inla.call = 'submit')  # Remote computing
# inla.setOption(inla.call = NULL)        # Local computing

## Basic logistic model with categorical covariates only
ml$cont_BC$res$ref <-
  inla(formula = yt ~ 0 + BC.cat + Year2013,
       data = ml$cont$dat,
       family = "gamma",
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

ml$cont_BC$res$noBC <-
  inla(formula = ml$cont$fml,
       data = inla.stack.data(ml$cont$stk),
       family = "gamma",
       control.predictor = list(
         A = inla.stack.A(ml$cont$stk),
         compute = TRUE, link = 1),
       control.family = list(
         hyper = list(theta = prec.od)
         # hyper = list(theta = list(fixed = TRUE, initial = 10))
       ),
       control.fixed = list(
         expand.factor.strategy = 'inla'), # account for NAs in factor
       control.compute = list(
         dic = TRUE,
         waic = TRUE,
         cpo = TRUE,
         po = TRUE,
         config = TRUE)
, control.inla = list(h = 0.005)
)


## BC as a categorical covariate
ml$cont_BC$res$cat <-
  inla(formula = update(ml$cont$fml, ~ BC.cat.cont + Year2013.cont + . -Year.cont),
       data = inla.stack.data(ml$cont$stk),
       family = "gamma",
       control.family = list(
         hyper = list(theta = prec.od)
         # hyper = list(theta = list(fixed = TRUE, initial = 10))
       ),
       control.predictor = list(
         A = inla.stack.A(ml$cont$stk),
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
ml$cont_BC$res$cat2 <-
  inla(formula = update(ml$cont$fml, ~ . -Year.cont + cut(BC.cont, breaks = c(0, 30, 150)):Year.cont),
       data = inla.stack.data(ml$cont$stk),
       family = "gamma",
       control.family = list(
         hyper = list(theta = prec.od)
         # hyper = list(theta = list(fixed = TRUE, initial = 10))
       ),
       control.predictor = list(
         A = inla.stack.A(ml$cont$stk),
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
ml$cont_BC$res$catxyear <-
  inla(formula = update(ml$cont$fml, ~ . -Year.cont + Year.cont:BC.cat.cont),
       data = inla.stack.data(ml$cont$stk),
       family = "gamma",
       control.family = list(
         hyper = list(theta = prec.od)
         # hyper = list(theta = list(fixed = TRUE, initial = 10))
       ),
       control.predictor = list(
         A = inla.stack.A(ml$cont$stk),
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

## For some reason, I needed to get rid of NAs in BC.cat for this model to run
## For identifiability, I need to code the Year with respect to the baseline of 2012

ml$cont_BC$res$lin <-
  inla(formula = update(ml$cont$fml, ~ BC.sc.cont + . ),
       data = inla.stack.data(ml$cont$stk),
       family = "gamma",
       control.family = list(
         hyper = list(theta = prec.od)
         # hyper = list(theta = list(fixed = TRUE, initial = 10))
       ),
       control.predictor = list(
         A = inla.stack.A(ml$cont$stk),
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

ml$cont_BC$res$quad <-
  inla(formula = update(ml$cont$fml, ~ BC.sc.cont + I(BC.sc.cont^2) + .),
       data = inla.stack.data(ml$cont$stk),
       family = "gamma",
       control.family = list(
         hyper = list(theta = prec.od)
         # hyper = list(theta = list(fixed = TRUE, initial = 10))
       ),
       control.predictor = list(
         A = inla.stack.A(ml$cont$stk),
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
ml$cont_BC$res$np <-
  inla(formula = update(ml$cont$fml,
                        ~ f(cut(BC.sc.cont, 51),
                            model = 'rw2',
                            scale.model = TRUE,
                            hyper = list(theta = prec.A)) +
                          .),
       data = inla.stack.data(ml$cont$stk),
       family = "gamma",
       control.family = list(
         hyper = list(theta = prec.od)
         # hyper = list(theta = list(fixed = TRUE, initial = 10))
       ),
       control.predictor = list(
         A = inla.stack.A(ml$cont$stk),
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

