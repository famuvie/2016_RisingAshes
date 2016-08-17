
#### MODEL FIT ####

inla.setOption(inla.call = 'submit')  # Remote computing
# inla.setOption(inla.call = NULL)        # Local computing

## Full model with selected covariates
ml$bin_msel$res$full <- 
  inla(formula = update(ml$bin_msel$fml,
                        ~ BF_2.bin +
                          BF_3.bin +
                          BF_4.bin +
                          BF_5.bin +
                          Year.bin:BC.bin +
                          . -Year.bin),
       data = inla.stack.data(ml$bin_msel$stk),
       family = "binomial",
       Ntrials = 1,
       ## Correction of the Laplace Approximation
       ## for Binomial data with many zeros
       ## http://arxiv.org/abs/1503.07307
       control.inla = list(
         correct = TRUE,
         correct.factor = 10),
       control.predictor = list(
         A = inla.stack.A(ml$bin_msel$stk),
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
