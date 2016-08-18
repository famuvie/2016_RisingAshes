## ----setup, echo=FALSE, message=FALSE, cache = FALSE---------------------
## Required packages
library(RisingAshes)
library(lme4)     # for reproducing the reference model
library(lattice)  # some lme4 plotting functions
library(boot)     # nonparametric bootstrap confidence intervals
library(tidyr)    # data manipulation
library(plyr)     # data manipulation
library(dplyr)    # data manipulation
library(ggplot2)  # plotting
library(viridis)  # color palette
library(knitr)    # kable()
library(INLA)     # Bayesian Inference

## knitr options
opts_chunk$set(echo       = FALSE,
               message    = FALSE,
               warning    = FALSE,
               comment    = NA,
               fig.dpi    = 96,
               fig.width  = 5,
               fig.height = 3,
               fig.caption= FALSE,
               cache      = FALSE)

## ggplot2 options
theme_set(theme_bw())

## ----preprocessing, echo=FALSE, message=FALSE----------------------------
## Applies an square-root transform to CD (creating variables CDT)
## Buils object mdat in long format (with variables: Year; Variable; Value)
## Open with:
# file.show(file=system.file('vignettes/cleanup.R',
#                            package = 'RisingAshes'))
source('cleanup.R')

## ----setup-model-list, echo=FALSE, message=FALSE-------------------------
# Object for storing data, models and solutions
ml <- list()

## Base dataset for analyzing the CD variable (transformed)
ml$dat <- droplevels(
  transform(
    subset(mdat, variable == 'CDT'),
    variable = NULL,
    BLC = factor(BLC),
    PROVxBLC = factor(PROV:factor(BLC)),
    FAMxBLC = factor(FAM:factor(BLC))
  )
)

## ----CD-progress, fig.cap = 'Crown dieback (CD) progress over time', fig.width=4----
pdat <- mdat %>% 
  filter(variable == 'CD', !is.na(value)) %>% 
  droplevels() %>% 
  group_by(Year, value) %>% 
  summarise(N = n()) %>% 
  transform(Year = as.numeric(levels(Year))[Year],
            CD_status = as.factor(value))
levels(pdat$CD_status) <- sub("1", "1 (dead)", levels(pdat$CD_status))

ggplot(pdat, aes(Year, N, fill = CD_status)) +
  geom_area(na.rm = TRUE, colour = 'white') +
  scale_fill_grey(start = 0.9, end = 0) +
  theme(legend.title=element_blank()) +
  guides(fill = guide_legend(reverse = TRUE,
                             override.aes = list(colour = NULL))) +
  ylab('Number of trees')

## ----Pliura-study--------------------------------------------------------

ml$ref$res1 <- lmer(
  value ~ Year + PROV + BLC + PROVxBLC + (1|FAM) + (1|FAMxBLC),
  ml$dat
)

summary(ml$ref$res1)

## ----Pliura-profiling, eval = FALSE--------------------------------------
#  # for computation of confidence intervals
#  ml$ref$prof <- profile(ml$ref$res1)

## ----Pliura-profiling-cache-recover, echo=FALSE--------------------------

## Load into ml$ref$prof the cached result of the previous chunk
ml$ref$prof <- RisingAshes:::Pliura_profile

## ----Pliura-profiling-check, eval=FALSE----------------------------------
#  ## The fact that all plots show an almost linear shape
#  ## is a sign that a normal approximation would work fine
#  ## Douglas M. Bates (2014). lme4: Mixed-effects modeling with R. Sec. 1.5.
#  xyplot(ml$ref$prof)
#  densityplot(ml$ref$prof)
#  

## ----Pliura-variance-components------------------------------------------
confint(ml$ref$prof, 1:3)

## ----Pliura-residuals, fig.width=4, fig.height=3-------------------------
plot(ml$ref$res1,type=c("p","smooth"))
plot(ml$ref$res1,
     sqrt(abs(resid(.)))~fitted(.), ## scale-location
     type=c("p","smooth"))
qqmath(ml$ref$res1, id=0.05)
shapiro.test(resid(ml$ref$res1))

## ----Pliura-heritability-function----------------------------------------
## Function giving the point estimate of the (narrow sense) heritability
## given the model fit
h2.fun <- function(., exclude.interaction = FALSE) {
  v.hat.reml <- with(as.data.frame(VarCorr(.)),
                     structure(vcov,
                               names = grp))
  phe.var <- ifelse(exclude.interaction,
                    v.hat.reml["FAM"] + v.hat.reml["Residual"],
                    sum(v.hat.reml))
  ## Additive-genetic variance estimate
  sigma2_a <- 4*v.hat.reml["FAM"]
  structure(sigma2_a/phe.var,
            names = 'h^2')
  }

## ----Pliura-heritability, eval = FALSE-----------------------------------
#  
#  ml$ref$boot.h2 <- bootMer(ml$ref$res1,
#                            h2.fun,
#                            nsim = 500,
#                            parallel = 'multicore',
#                            ncpus = 8L)

## ----Pliura-heritability-cache-recover, echo=FALSE-----------------------

## Load into ml$ref$boot.h2 the cached result of the previous chunk
ml$ref$boot.h2 <- RisingAshes:::Pliura_heritability

## ----Pliura-population---------------------------------------------------
confint(ml$ref$prof, c('PROVc', 'PROVv'))

## ----Pliura-year---------------------------------------------------------
year.params <- paste0('Year', c('.L', '.Q', '.C', '^4'))
confint(ml$ref$prof, year.params)
year.effects <- contr.poly(5) %*% fixef(ml$ref$res1)[year.params]


## ----Pliura-year-bootstrap, eval = FALSE---------------------------------
#  year.eff.fun <- function(.) contr.poly(5) %*% fixef(.)[year.params]
#  ml$ref$boot1 <- bootMer(
#    ml$ref$res1,
#    year.eff.fun,
#    nsim = 500,
#    parallel = 'multicore',
#    ncpus = 8L
#  )

## ----Pliura-year-cache-recover, echo=FALSE-------------------------------

## Load into ml$ref$boot1 the cached result of the previous chunk
ml$ref$boot1 <- RisingAshes:::Pliura_year

## ----Pliura-year-plot, fig.width=4, fig.height=3-------------------------
ml$ref$bootCI <- t(sapply(
  1:5,
  function(ind) boot.ci(
    ml$ref$boot1,
    type = 'norm',
    index = ind)$normal[-1])
)
colnames(ml$ref$bootCI) <- c('ymin', 'ymax')

ggplot(data.frame(year = factor(2010:2014),
                  effect = year.effects,
                  ml$ref$bootCI),
       aes(year, effect)) + 
  geom_pointrange(aes(x = year, ymin = ymin, ymax = ymax))

## ----Pliura-best-submodel------------------------------------------------
ml$ref$res2 <- lmer(value ~ Year + BLC + (1|FAM) + (1|FAMxBLC),
                    ml$dat)
summary(ml$ref$res2)
confint(ml$ref$res2, method = 'Wald')
anova(ml$ref$res1, ml$ref$res2)

## ----Pliura-family, fig.width=5, fig.height=4----------------------------
dotplot(ranef(ml$ref$res2, condVar = TRUE, whichel = 'FAM'))

## ----Pliura-INLA, fig.width=4, fig.height=3------------------------------
ml$ref$res.inla <- 
  inla(
    formula = value ~ 
      Year + BLC + f(FAM, model = 'iid') + f(FAMxBLC, model = 'iid'),
    data = ml$dat,
    control.predictor = list(compute = TRUE),
    control.fixed = list(
      expand.factor.strategy = 'inla'), # account for NAs in factor
    control.compute = list(
      dic = TRUE,
      waic = TRUE,
      cpo = TRUE,
      po = TRUE,
      config = TRUE)
    #, inla.call = 'submit'
    )

summary(ml$ref$res.inla)

## Compare Family point estimates
ggplot(data.frame(lme4 = ranef(ml$ref$res2)$FAM$'(Intercept)',
                  INLA = ml$ref$res.inla$summary.random$FAM$mean),
       aes(lme4, INLA)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, col = 'darkgray')

## ----ref-model-CD-heritability, fig.width=3, fig.height=2, fig.cap = 'Posterior heritability'----
# summary(ml$add$res)
N.samples <- 5e3
hyperpar.sample <- inla.hyperpar.sample(N.samples, ml$ref$res.inla)

var.sample <- structure(1/hyperpar.sample,
                        dimnames = list(NULL, c('sigma2', 'sigma2_f', 'sigma2_fb')))

h2.fun <- function(x, exclude.interaction = FALSE) {
  sigma2_a <- 4*x["sigma2_f"]
  phe.var <- ifelse(exclude.interaction,
                    x["sigma2"] + x["sigma2_f"],
                    sum(x))
  ## Additive-genetic variance estimate
  structure(min(sigma2_a/phe.var, 1),
            names = 'h^2')
  }

ml$ref$post.h2 <- apply(var.sample, 1, h2.fun, exclude.interaction = FALSE)
quantile(ml$ref$post.h2, probs = c(.025, .975))
# ggplot(data.frame(h2 = ml$ref$post.h2), aes(x = h2)) + geom_density()


## ----genetic-structure---------------------------------------------------
## Creates objects:
##   gen_str
##     $Ainv
##     $recoding
##     $contraint.matrix
##   prec.A
## Open with:
# file.show(file=system.file('vignettes/genetic_structure',
#                            package = 'RisingAshes'))
source('genetic_structure.R')


## ----additive-genetic-fit------------------------------------------------

## Recoding individuals in data!
ml$add$dat <- transform(
  ml$dat,
  Y11 = as.numeric(Year == '2011'),
  Y12 = as.numeric(Year == '2012'),
  Y13 = as.numeric(Year == '2013'),
  Y14 = as.numeric(Year == '2014'),
  recode = gen_str$recoding
)

## Model formula
ml$add$fml <- 
  value ~ 0 + Year + BLC + 
  f(recode, model = "generic0", hyper = list(theta = prec.A),
    constr = TRUE, Cmatrix = gen_str$Ainv)

ml$add$res <- inla(
  formula = ml$add$fml,
  data = ml$add$dat,
  control.predictor = list(compute = TRUE),
  control.fixed = list(
    expand.factor.strategy = 'inla'  # account for NAs in factor
  ),
  control.compute = list(
    dic = TRUE,
    waic = TRUE,
    cpo = TRUE,
    po = TRUE,
    config = TRUE
  ),
  control.family = list(hyper = list(prec = prec.A))
  #, inla.call = 'submit'
)
summary(ml$add$res)

## ----additive-model-CD-heritability, fig.width=3, fig.height=2, fig.cap = 'Posterior heritability'----
# summary(ml$add$res)

N.samples <- 5e3
hyperpar.sample <- inla.hyperpar.sample(N.samples, ml$add$res)
var.sample <- structure(1/hyperpar.sample[, grep('Gaussian|recode',
                                                 colnames(hyperpar.sample))],
                        dimnames = list(NULL, c('sigma2', 'sigma2_a')))
ml$add$h2.sample <- apply(var.sample, 1, 
                          function(x) x['sigma2_a']/(x['sigma2_a'] + x['sigma2']))
ggplot(data.frame(h2 = ml$add$h2.sample), aes(x = h2)) + geom_density()


## ----spatial-structure---------------------------------------------------
## Generates an object
## sp_str
##   $loc
##   $mesh
##   $spde
sp_str <- Devecey_spde()


## ----spatiotemporal-structure--------------------------------------------
## mesh space (mesh$n * nYears = 942 * 5 = 4710 rows)
st.index <- inla.spde.make.index(
  name = 'theta',
  n.spde = sp_str$spde$n.spde,
  n.group = nlevels(ml$dat$Year)
)

## Maps from the observation to the prediction space
st.Amat <- inla.spde.make.A(
  sp_str$mesh, loc = sp_str$loc, 
  index = rep(1:nrow(sp_str$loc), nlevels(ml$dat$Year)),
  group = as.numeric(ml$dat$Year)
)


## ----spatiotemporal-model-CD-setup---------------------------------------
## I changed a bit the priors of the spatial effect
## See src/stmodel.R
ml$st$dat <- transform(
  ml$dat[!is.na(ml$dat$BF), c('Seq', 'FAM', 'PROV', 'BF98', 'Year', 'value')],
  mu = 1,
  recode = gen_str$recoding[Seq],
  FAM = as.numeric(FAM),
  BF  = BF98,
  Y11 = as.numeric(Year == '2011'),
  Y12 = as.numeric(Year == '2012'),
  Y13 = as.numeric(Year == '2013'),
  Y14 = as.numeric(Year == '2014')
)

ml$st$stack <-
  inla.stack(data    = list(y = ml$st$dat$value),
             A       = list(st.Amat[!is.na(ml$dat$BF),], 1),
             effects = list(st.index, ml$st$dat),
             tag     = 'est')

ml$st$stack.idx <- inla.stack.index(ml$st$stack, 'est')


## ----spatiotemporal-model-CD-fit, eval = FALSE---------------------------
#  ## Competing models
#  ml$st$var <- list(
#    wif = list(desc = 'Family integrated into blups',
#               fml  = y ~ 0 + Year +
#                 f(theta, model=sp_str$spde,
#                   group = theta.group,  # variable in the mesh space
#                   control.group = list(model = 'exchangeable')) +
#                 f(recode, model = "generic0", hyper = list(theta = prec.A),
#                   constr = TRUE, Cmatrix = gen_str$Ainv),
#               res = NULL),
#    wpe = list(desc = 'Model with Provenance fixed effect',
#               fml  = y ~ 0 + Year + PROV +
#                 f(theta, model=sp_str$spde,
#                   group = theta.group,  # variable in the mesh space
#                   control.group = list(model = 'exchangeable')) +
#                 f(recode, model = "generic0", hyper = list(theta = prec.A),
#                   constr = TRUE, Cmatrix = gen_str$Ainv),
#               res = NULL),
#    wbf = list(desc = 'Model with Bud-flush effect',
#               fml  = y ~ 0 + BF + Y11 + Y12 + Y13 + Y14 +
#                 f(theta, model=sp_str$spde,
#                   group = theta.group,  # variable in the mesh space
#                   control.group = list(model = 'exchangeable')) +
#                 f(recode, model = "generic0", hyper = list(theta = prec.A),
#                   constr = TRUE, Cmatrix = gen_str$Ainv),
#               res = NULL),
#    wef = list(desc = 'Explicit Family effect and centered blups',
#               fml  = y ~ 0 + BF + Y11 + Y12 + Y13 + Y14 +
#                 f(FAM, model = 'iid', constr = TRUE) +
#                 f(theta, model=sp_str$spde,
#                   group = theta.group,  # variable in the mesh space
#                   control.group = list(model = 'exchangeable')) +
#                 f(recode, model = "generic0", hyper = list(theta = prec.A), constr = FALSE,
#                   extraconstr = list(A=const.mat, e=rep(0, n.fam + 1)), Cmatrix = gen_str$Ainv),
#               res = NULL),
#    wof = list(desc = 'Without differences between families',
#               fml  = y ~ 0 + BF + Y11 + Y12 + Y13 + Y14 +
#                 f(theta, model=sp_str$spde,
#                   group = theta.group,  # variable in the mesh space
#                   control.group = list(model = 'exchangeable')) +
#                 f(recode, model = "generic0", hyper = list(theta = prec.A), constr = FALSE,
#                   extraconstr = list(A=const.mat, e=rep(0, n.fam + 1)), Cmatrix = gen_str$Ainv),
#               res = NULL))
#  
#  
#  
#  ## Send jobs
#  ## This takes a while. Commented lines allow sending jobs to a server.
#  ## control.comput$config = TRUE allows simulating from the posterior distributions
#  ## in order to compute the posterior density for the heritability
#  for(x in seq_along(ml$st$var)) {
#    ml$st$var[[x]]$res <- inla(
#      formula = ml$st$var[[x]]$fml,
#      data = inla.stack.data(ml$st$stack),
#      control.predictor = list(
#        A = inla.stack.A(ml$st$stack),
#        compute = TRUE),
#      control.fixed = list(
#        expand.factor.strategy = 'inla'), # account for NAs in factor
#      control.compute = list(
#        dic = TRUE,
#        waic = TRUE,
#        cpo = TRUE,
#        po = TRUE,
#        config = TRUE),
#      control.family = list(
#        hyper = list(prec = prec.A))
#      # , inla.call = 'submit'
#    )
#  }
#  

## ----spatiotemporal-model-CD-remote-retrieve, eval = FALSE---------------
#  
#  ## Retrieve jobs from server in due case
#  for(x in seq_along(ml$st$var)) {
#    ## Wait while the job is running
#    #   inla.qstat()
#    #   while(inla.qstat(ml$st$var[[x]]$res)[[1]]$status == "Running") {
#    #     Sys.sleep(60)
#    #     }
#    ml$st$var[[x]]$res <- inla.qget(ml$st$var[[x]]$res)
#    ml$st$var[[x]]$spde <- inla.spde2.result(ml$st$var[[x]]$res, 'theta', sp_str$spde)
#  }
#  

## ----spatiotemporal-model-CD-cache-save-recover, echo=FALSE--------------

## Load into ml$st$var the cached version of the two previous chunks
ml$st$var <- RisingAshes:::stmodels_CD


## ----test-spdexAR-model, eval = FALSE------------------------------------
#  
#  ## An AR1 model for the temporal evolution seems more natural
#  ## however, the model badly overfits the observations
#  ## Look at the effective number of paramters, and the high precisions
#  
#  fml  = y ~ 0 + Year +
#    f(theta, model=sp_str$spde,
#      group = theta.group,  # variable in the mesh space
#      control.group = list(model = 'ar1')) +
#    f(recode, model = "generic0", hyper = list(theta = prec.A),
#      constr = TRUE, Cmatrix = gen_str$Ainv)
#  tr <- inla(formula = fml,
#             data = inla.stack.data(ml$st$stack),
#             control.predictor = list(A = inla.stack.A(ml$st$stack),
#                                      compute = TRUE),
#             control.fixed = list(expand.factor.strategy = 'inla'), # account for NAs in factor
#             control.compute = list(dic = TRUE,
#                                    waic = TRUE,
#                                    cpo = TRUE,
#                                    po = TRUE,
#                                    config = TRUE),
#             control.family = list(hyper = list(prec = prec.A))
#             , inla.call = 'submit'
#             )
#  
#  summary(tr)

## ----spatiotemporal-model-CD-summary-------------------------------------
summary(ml$st$var$wbf$res)

## ----spatiotemporal-model-CD-variances, fig.width = 4, fig.height = 2.5, fig.cap = 'Posterior variances within and between families'----
# Marginals of the variances (inverse precisions)
ml$st$marginal.variances <- 
  with(
    ml$st$var,
    c(wef$res$marginals.hyperpar[c('Precision for recode',
                                   'Precision for FAM')],
      wbf$res$marginals.hyperpar[c('Precision for the Gaussian observations',
                                   'Precision for recode')])) %>% 
  lapply(function(y) inla.tmarginal(function(x) 1/x, y)) %>% 
  structure(names = c('Within Family',
                      'Between Families',
                      'Residual',
                      'Additive-genetic')) %>% 
  c(.,
    list('Spatio-temporal' =
         ml$st$var$wbf$spde$marginals.variance.nominal$variance.nominal.1))


ggplot(ldply(ml$st$marginal.variances[1:2], identity), aes(x, y)) + 
  geom_line(aes(col = .id))

## ----spatiotemporal-model-CD-summary-variances---------------------------
## Mode and 95% HPD Credible Intervals

ml$st$summary.variances <- 
  rbind(
    sapply(ml$st$marginal.variances, inla.mmarginal),
    sapply(ml$st$marginal.variances, function(x) inla.hpdmarginal(.95, x))
  ) %>% 
  structure(dimnames = list(c('mode', 'lo', 'hi'),
                            names(ml$st$marginal.variances)))

print_estimate <- function(x) 
  paste0(unname(round(ml$st$summary.variances[1, x], 3)),
        ' (', paste(round(ml$st$summary.variances[2:3, x], 3), collapse = ', '), ')')


## ----spatiotemporal-model-CD-variance-components, fig.width=4, fig.height=2.5, fig.cap = 'Posterior modes and 95% HPD Credible Intervals of the variance components.'----
## improved estimates of hyperparamters
# ml$st$wif$res <- inla.hyperpar(ml$st$var$wbf$res)

ml$st$summary.variances %>% 
  adply(2, .id = 'Component') %>% 
  rename(Variance = mode) %>% 
  arrange(desc(Variance)) %>% 
  mutate(Component = factor(Component, levels = Component)) %>% 
  ggplot(aes(Component, Variance)) +
  geom_pointrange(aes(ymin = lo, ymax = hi)) +
  coord_flip() +
  expand_limits(y = 0)


## ----spatiotemporal-model-CD-variance-components-table-------------------

kable(ml$st$summary.variances, digits = 3)


## ----spatiotemporal-model-CD-year, fig.width=4, fig.height=2, fig.cap = 'Posterior mean effect of the Year and 95% Credible Intervals. Transformed scale.'----
data.frame(
  Year = 2010:2014,
  rbind(ml$st$var$wbf$res$summary.fixed[3, ],
        ml$st$var$wbf$res$summary.fixed[3, 'mean'] +
          ml$st$var$wbf$res$summary.fixed[6:9, ])
) %>% 
  ggplot(aes(Year, mean, ymin = X0.025quant, ymax = X0.975quant)) + 
  geom_pointrange() + 
  geom_line() +
  labs(y=NULL)

## ----spatiotemporal-model-CD-BF, fig.width=4, fig.cap = 'Posterior mean effect of the Budflush and 95% Credible Intervals for Year 2010. Transformed scale.'----
ggplot(data.frame(BF = 1:5,
                  ml$st$var$wbf$res$summary.fixed[1:5,]),
       aes(BF, mean, ymin = X0.025quant, ymax = X0.975quant)) + 
  geom_pointrange() + 
  geom_line(col = 'darkgray') +
  labs(y=NULL)

## ----spatiotemporal-model-CD-BFxYear-oscale, fig.width=4, fig.height=2, fig.cap = 'Posterior mean effect of the Budflush and 95% Credible Intervals, for each year. Original scale.'----
inverse_transf <- function(x) 1/x**2 - 0.1


BFmar_latent <- 
  sapply(c(0, ml$st$var$wbf$res$summary.fixed[6:9, 'mean']),
         function(x) shift_marginal(
           ml$st$var$wbf$res$marginals.fixed[1:5],
           x
         )
  )

BFmar <- lapply(BFmar_latent,
                function(x) inla.tmarginal(inverse_transf, x))
BFhpd <- rbind(
  sapply(BFmar,
         function(x) inla.hpdmarginal(.95, x)),
  sapply(BFmar,
         function(x) inla.emarginal(identity, x))
)

plotdat <- 
  BFhpd %>% 
  t() %>% 
  as.data.frame() %>% 
  rename(ymin = V1, ymax = V2, mean = V3) %>% 
  mutate(BF = rep(1:5, 5),
         Year = factor(rep(2010:2014, each = 5)))

plotdat %>% 
  ggplot(aes(BF, mean, ymin = ymin, ymax = ymax)) + 
  geom_ribbon(aes(group = Year), alpha = 0.3) + 
  geom_line(aes(linetype = Year)) +
  labs(y='CD') +
  scale_linetype_discrete(guide = guide_legend(reverse = TRUE))

## ----spatiotemporal-model-CD-fixed-means-table---------------------------

kable(addmargins(xtabs(mean~Year+BF, data = plotdat), 
                 FUN = mean, 
                 quiet = TRUE),
      caption = 'mean predicted CD by Year and BF.',
      digits = 3)

## ----sptaiotemporal-model-CD-heritability, fig.width=3, fig.height=2, fig.cap = 'Posterior heritability'----
# summary(ml$st$var$wbf$res)

N.samples <- 5e3
hyperpar.sample <- inla.hyperpar.sample(N.samples, ml$st$var$wbf$res)
var.sample <- structure(
  cbind(1/hyperpar.sample[, grep('Gaussian|recode',
                                 colnames(hyperpar.sample))],
        inla.rmarginal(N.samples, ml$st$marginal.variances$`Spatio-temporal`)),
  dimnames = list(NULL, c('sigma2', 'sigma2_a', 'sigma2_st')))

sample_h2 <- function(fun)
  apply(var.sample, 1, fun)

ml$st$var$wbf$h2.sample <-
  c(with_ST = function(x)
    x['sigma2_a']/(x['sigma2_a']+x['sigma2']+x['sigma2_st']),
    wout_ST = function(x)
      x['sigma2_a']/(x['sigma2_a']+x['sigma2'])) %>% 
  sapply(sample_h2) %>% 
  as.data.frame()
  


ggplot(gather(ml$st$var$wbf$h2.sample, key, h2), aes(h2, colour = key)) +
  geom_density() +
  ylab(NULL) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())


## ----sptaiotemporal-model-CD-heritability-journal, fig.height=2----------

# Return the mean and the centered 95% credible interval
# for a sample x
mean_ci <- function(x) {
  c(mean = mean(x), quantile(x, probs = c(.025, .975)))
}

plotdat <- 
  sapply(ml$st$var$wbf$h2.sample, mean_ci) %>% 
  t() %>% 
  as.data.frame()

ggplot(mutate(plotdat, key = rownames(plotdat)), aes(key, y=mean)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`)) +
  coord_flip() +
  ylim(c(0, 1)) +
  ylab('Heritability') +
  xlab(NULL)


## ----sptaiotemporal-model-CD-heritability-table--------------------------

kable(plotdat, digits = 2)

## ----spatiotemporal-model-CD-fitvsobs, fig.width=4, fig.height=3, fig.cap = 'Fitted vs. observed values.Transformed scale.'----
plotdat <- data.frame(
  Observed = ml$st$dat$value,
  Predicted = ml$st$var$wbf$res$summary.fitted.values$mean[ml$st$stack.idx$data],
  Year = ml$st$dat$Year) 

ggplot(plotdat, aes(Observed, Predicted)) +
  geom_jitter(aes(color = Year), alpha = .5) +
  geom_abline(intercept = 0, slope = 1, colour = 'darkgray')

## ----spatiotemporal-model-CD-fitvsobs-oscale, fig.width=4, fig.height=3, fig.cap = 'Fitted vs. observed values. Original scale.'----

data.frame(Year = plotdat$Year,
           mutate_each(plotdat[,1:2], funs(inverse_transf))) %>% 
  ggplot(aes(Observed, Predicted)) +
  geom_jitter(aes(color = Year), alpha = .5) +
  geom_abline(intercept = 0, slope = 1, colour = 'darkgray')

## ----spatiotemporal-model-CD-fam, fig.width=4, fig.height=3, fig.cap = 'Posterior family effects'----
# # Identify the most distinguishable families at the extremes
# fam.idx <- order(ml$st$var$wef$res$summary.random$FAM$mean)[c(1:2,22:23)]
# message(paste('Two leftmost and two rightmost families: ',
#               paste(levels(dat$FAM)[fam.idx], collapse = ', ')))
# sizes <- rep(1, 23)
# sizes[fam.idx] <- 2
# ggplot(transform(ldply(ml$st$var$wef$res$marginals.random$FAM, data.frame),
#                  .id = factor(.id, levels = unique(.id), labels = levels(dat$FAM))),
#        aes(x, y, col = .id)) +
#   geom_line(aes(size = sizes[.id])) +
#   scale_size(range = c(.3, 1), guide = FALSE) +
#   guides(color = guide_legend(keyheight = 0.5))
# 
# dotplot(ranef(ml$ref$res2, condVar = TRUE, whichel = 'FAM'))

## Rankings from low to high from M1 and M3
## Not used
rank_change.idx <- which(order(ranef(ml$ref$res2,
                                     condVar = TRUE, whichel = 'FAM')$FAM) !=
                           order(ml$st$var$wef$res$summary.random$FAM$mean))
lwd <- rep(0.5, length(levels(dat$FAM)))
lwd[rank_change.idx] <- 1

ggplot(transform(ml$st$var$wef$res$summary.random$FAM,
                 Family = reorder(factor(levels(dat$FAM)), mean, 'mean')),
       aes(x = Family, y = mean)) +
  geom_pointrange(aes(ymin = X0.025quant, ymax = X0.975quant)) + 
  coord_flip() +
  labs(y=NULL)

## ----spatiotemporal-model-CD-genetic, fig.width = 7, fig.height = 6, fig.cap = 'Posterior mean and 95% HPD Credible Interval of individual additive-genetic effects'----
hpds <- sapply(ml$st$var$wbf$res$marginals.random$recode, 
               function(x) inla.hpdmarginal(.95, x))
recode.idx <- gen_str$recoding

ml$st$blups <- data.frame(
 Seq  = dat$Seq,
 Family = dat$FAM,
 CD14 = dat$CD_2014,
 mean = ml$st$var$wbf$res$summary.random$recode$mean[recode.idx],
 ymin = hpds[1, recode.idx],
 ymax = hpds[2, recode.idx])

ggplot(ml$st$blups, aes(x = mean, y = mean, ymin = ymin, ymax = ymax)) + 
  geom_pointrange() + 
  facet_wrap(~ Family) +
  coord_flip() +
  labs(x=NULL, y=NULL)

## ----check-genetic-vs-phenotype, eval = FALSE----------------------------
#  ggplot(ml$st$blups, aes(CD14, mean)) +
#    geom_point()

## ----spatiotemporal-model-CD-hyperpar, fig.width=4, fig.cap = 'Posterior hyperparameters of the spatio-temporal field'----

## Spatial-temporal field hyperparameters
stopifnot(exists('spde', ml$st$var$wbf))
res.spde.marginals <- ldply(
  tail(ml$st$var$wbf$spde, 2),
  function(x) data.frame(x[[1]])
)
res.spde.marginals <- transform(
  res.spde.marginals,
  .id=factor(res.spde.marginals$.id,
             labels = c('Range', 'Variance'))
)
res.spde.marginals <- rbind(
  cbind(type = 'Posterior', res.spde.marginals),
  cbind(type = 'Prior', .id = 'Range', as.data.frame(sp_str$priors$rho)),
  cbind(type = 'Prior', .id = 'Variance', as.data.frame(sp_str$priors$sigma))
)
# qplot(x, y, data=res.spde.marginals, geom='line', col = type) + facet_wrap(~.id, scales='free')


ggplot(res.spde.marginals, aes(x, y, lty = type)) +
  geom_line() +
  facet_wrap(~.id, scales='free', ncol=1) + 
  labs(x=NULL, y=NULL) +
  theme(legend.title = element_blank())


## ----spatiotemporal-model-CD-spatial, fig.width=7, fig.height=5, fig.cap = 'Posterior mean of the spatio-temporal field. Transformed scale.'----
proj.vec = inla.mesh.projector(sp_str$mesh, loc=sp_str$loc)
sp.f <- function(x) inla.mesh.project(proj.vec, field = x)

ml$st$post_spat <- ldply(1:5, function(x) 
  data.frame(
    X = dat$X,
    Y = dat$Y,
    Year = (2010:2014)[x],
    PostMean= sp.f(ml$st$var$wbf$spde$summary.values$mean[st.index$theta.group == x])))

ggplot(ml$st$post_spat, aes(X, Y)) +
  geom_tile(aes(fill = PostMean, color = PostMean)) +
  coord_fixed() +
  labs(x='Field column', y='Field row') +
  facet_wrap(~ Year) +
  scale_fill_gradient2(low='#832424FF', high='#3A3A98FF', name = 'Posterior\nmean') +
  scale_color_gradient2(low='#832424FF', high='#3A3A98FF', name = 'Posterior\nmean') 


# ggplot(transform(dat, theta.mean = ss.mean), aes(X, Y)) + geom_tile(aes(fill=theta.mean)) + coord_fixed()
# ggplot(transform(dat, theta.sd = ss.sd), aes(X, Y)) + geom_tile(aes(fill=theta.sd)) + coord_fixed()

# # sample some points in the trial
# ggplot(ldply(sample(nrow(dat), 30), function(x) {
#   idx <- which(ml$st$post_spat$X == dat[x, 'X'] & ml$st$post_spat$Y == dat[x, 'Y'])
#   data.frame(Tree = paste(dat[x, c('X', 'Y')], collapse = '_'),
#              ml$st$post_spat[idx, c('Year', 'PostMean')])
#   }
#   ),
#   aes(Year, PostMean, group = Tree)) + 
#   geom_point() + 
#   geom_line(color = 'darkgray')


## ----spatiotemporal-model-CD-spatial-oscale, fig.width=7, fig.height=5, fig.cap = 'Posterior mean of the spatio-temporal field. Original scale, for a fixed BF value of 3.'----

spmar <- ml$st$var$wbf$spde$marginals.values

## mean values of linear predictor for BF=3 and for each Year
latent_means <- with(ml$st$var$wbf$res, 
                     summary.fixed[3, 'mean'] + 
                       c(0, summary.fixed[6:9, 'mean']))
spmar.oscale <- 
  lapply(
    mapply(shift_marginal, spmar, latent_means[st.index$theta.group],
           SIMPLIFY = FALSE),
    function(x) inla.tmarginal(inverse_transf, x))

spmeans.oscale <- sapply(spmar.oscale, 
                         function(x) inla.emarginal(identity, x))

ml$st$post_spat_oscale <- ldply(
  1:5, function(x) 
    data.frame(
      X = dat$X,
      Y = dat$Y,
      Year = (2010:2014)[x],
      PostMean= sp.f(spmeans.oscale[st.index$theta.group == x])))

ggplot(ml$st$post_spat_oscale, aes(X, Y)) +
  geom_tile(aes(fill = PostMean, color = PostMean)) +
  coord_fixed() +
  labs(x='Field column', y='Field row') +
  facet_wrap(~ Year) +
  scale_fill_viridis(name = 'Posterior\nmean') +
  scale_color_viridis(name = 'Posterior\nmean') 



## ----spatiotemporal-model-CD-residuals, fig.width=4, fig.height=3--------
resid <- ml$st$dat$value - ml$st$var$wbf$res$summary.fitted.values$mean[ml$st$stack.idx$data]
qqmath(resid, xlab = "Standard normal quantiles",
       prepanel = prepanel.qqmathline,
       panel = function(x, subscripts, ...) {
         panel.qqmathline(x, ...)
         panel.qqmath(x, ...)
       })
shapiro.test(resid)

## ----summary-table, fig.height = 3, fig.width = 4, fig.cap = 'Model selection criteria.'----
st <- data.frame(Model = c(paste0('M', 1:5), paste(' M5', 1:2, sep = '.')),
                 Description = c("Pliūra et al. 2011",
                                 "Individual genetic effect",
                                 "Spatio-temporal effect",
                                 "Fixed provenance effect",
                                 "Fixed Budflush effect",
                                 "    explicit family effect",
                                 "    null variance between families"),
                 `Linear predictor` = c("Year + Block + _Fam_ + _Fam:Block_",
                                        "Year + Block + _BV_",
                                        "Year + _ST_ + _BV_",
                                        "Year + Prov + _ST_ + _BV_",
                                        "Year + BF + _ST_ + _BV_",
                                        "  Year + BF + _ST_ + _Fam_ + _BV_",
                                        "  Year + BF + _ST_ + (_Fam_=0) + _BV_"),
                 DIC = c(ml$ref$res.inla$dic$dic,
                         ml$add$res$dic$dic,
                         ml$st$var$wif$res$dic$dic,
                         ml$st$var$wpe$res$dic$dic,
                         ml$st$var$wbf$res$dic$dic,
                         ml$st$var$wef$res$dic$dic,
                         ml$st$var$wof$res$dic$dic),
                 p_D = c(ml$ref$res.inla$dic$p.eff,
                         ml$add$res$dic$p.eff,
                         ml$st$var$wif$res$dic$p.eff,
                         ml$st$var$wpe$res$dic$p.eff,
                         ml$st$var$wbf$res$dic$p.eff,
                         ml$st$var$wef$res$dic$p.eff,
                         ml$st$var$wof$res$dic$p.eff),
                 WAIC = c(ml$ref$res.inla$waic$waic,
                         ml$add$res$waic$waic,
                         ml$st$var$wif$res$waic$waic,
                         ml$st$var$wpe$res$waic$waic,
                         ml$st$var$wbf$res$waic$waic,
                         ml$st$var$wef$res$waic$waic,
                         ml$st$var$wof$res$waic$waic),
                 p_W = c(ml$ref$res.inla$waic$p.eff,
                         ml$add$res$waic$p.eff,
                         ml$st$var$wif$res$waic$p.eff,
                         ml$st$var$wpe$res$waic$p.eff,
                         ml$st$var$wbf$res$waic$p.eff,
                         ml$st$var$wef$res$waic$p.eff,
                         ml$st$var$wof$res$waic$p.eff),
                 `-2log-lik` = -2*c(ml$ref$res.inla$mlik[1],
                         ml$add$res$mlik[1],
                         ml$st$var$wif$res$mlik[1],
                         ml$st$var$wpe$res$mlik[1],
                         ml$st$var$wbf$res$mlik[1],
                         ml$st$var$wef$res$mlik[1],
                         ml$st$var$wof$res$mlik[1]))
kable(st, digits = 0)
ggplot(st, aes(DIC, X.2log.lik, label = Model)) +
  geom_text() +
  labs(y = "-2*Marginal likelihood")

## ----heritability-comparison, fig.height = 3, fig.width = 5, fig.cap = 'Posterior mean and 95% credible interval of heritability by model. Excluding ST variance.'----
dat.h2 <- rbind(data.frame(model = 'M1',
                           h2 = ml$ref$post.h2),
                data.frame(model = 'M2',
                           h2 = ml$add$h2.sample),
                data.frame(model = 'M3',
                           h2 = ml$st$var$wbf$h2.sample$wout_ST))

dat.h2 %<>%
  group_by(model) %>%
  summarise(mean.h2 = mean(h2),
            h.min = quantile(h2, .025),
            h.max = quantile(h2, .975))

ggplot(dat.h2, aes(x = model, y = mean.h2)) +
  geom_pointrange(aes(ymin = h.min, ymax = h.max)) +
  coord_flip() + 
  scale_x_discrete(limits = c('M3', 'M2', 'M1')) + 
  ylab('Posterior heritability') + 
  ylim(0, 1)


## ----family-variance-M1vsM3, fig.height = 2.5, fig.width = 4.5, fig.cap = 'Posterior variance between families.'----
ggplot(rbind(data.frame(Model = 'M1',
                        inla.tmarginal(function(x) 1/x,
                                       ml$ref$res.inla$marginals.hyperpar[['Precision for FAM']])),
             data.frame(Model = 'M3.1',
                        inla.tmarginal(function(x) 1/x,
                                       ml$st$var$wef$res$marginals.hyperpar[['Precision for FAM']]))),
       aes(x, y, col = Model)) + 
  geom_line()
       

## ----family-effect-M1vsM3, fig.height = 3, fig.width = 4.5, fig.cap = 'Posterior family effects by model.'----

ggplot(data.frame(M1 = ml$ref$res.inla$summary.random$FAM$mean,
                 M3.1 = ml$st$var$wef$res$summary.random$FAM$mean,
                 Family = levels(mdat$FAM)),
       aes(x = M1, y = M3.1, label = Family)) +
  geom_text(hjust = -1) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, col = 'darkgray')

