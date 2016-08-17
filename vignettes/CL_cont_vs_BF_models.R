
#### MODEL FIT ####

inla.setOption(inla.call = 'submit')  # Remote computing
# inla.setOption(inla.call = NULL)        # Local computing



## Basic logistic model with categorical covariates only
ml$cont_BF$res$ref <- 
  inla(formula = yt ~ 0 + BF98 + Year2013,
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


## Reference model: without BF

ml$cont_BF$res$noBF <- 
  inla(formula = ml$cont$fml,
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


## BF as a categorical covariate
ml$cont_BF$res$cat <- 
  inla(formula = update(ml$cont$fml, ~ BF98.cont + Year2013.cont + . -Year.cont),
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


## BF*Year as a categorical covariate
ml$cont_BF$res$catxyear <- 
  inla(formula = update(ml$cont$fml, ~ . -Year.cont + Year.cont:BF98.cont),
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

