## ----setup, echo=FALSE, message=FALSE, cache=FALSE-----------------------
## Required packages
library(RisingAshes)
library(tidyr)    # data manipulation
library(plyr)     # data manipulation
library(dplyr)    # data manipulation
library(ggplot2)  # plotting
library(grid)     # plotting
library(gridExtra)# plotting
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

## ----preprocessing, echo=FALSE-------------------------------------------
## Applies an square-root transform to CD (creating variables CDT)
## Buils object mdat in long format (with variables: Year; Variable; Value)
## Open with:
# file.show(file=system.file('vignettes/cleanup.R',
#                            package = 'RisingAshes'))
source('cleanup.R')

## ----setup-model-list, echo=FALSE, message=FALSE-------------------------
# Object for storing models and solutions
ml <- list()

## Base dataset for analyzing the CL variable
ml$dat <- 
  mdat %>% 
  filter(
    variable %in% c('CL', 'BC')
    ,Year %in% c('2012', '2013')
  ) %>%
  spread(variable, value) %>% 
  transform(BLC = factor(BLC),
            PROVxBLC = factor(PROV:factor(BLC)),
            FAMxBLC = factor(FAM:factor(BLC))) %>%
  droplevels()


## ----data-description----------------------------------------------------
ggplot(ml$dat, aes(x = CL)) + geom_histogram(bins = 30)

# # Only the non-null data
# ggplot(data.frame(CL = ml$dat$CL[ml$dat$CL > 0]), aes(x = CL)) + geom_histogram()
# # Its log-transform (not very Gaussian...)
# ggplot(data.frame(CL = log(ml$dat$CL[ml$dat$CL > 0])), aes(x = CL)) + geom_histogram()
# tf <- function(x) log(x) - log(1 - x + 0.02)
# ggplot(data.frame(CL = tf(ml$dat$CL[ml$dat$CL > 0])), aes(x = CL)) + geom_histogram()
# # Its logit-transform (not very Gaussian...)
# ggplot(data.frame(CL = logit(ml$dat$CL[ml$dat$CL > 0])), aes(x = CL)) + geom_histogram()


## ----data-arrangement----------------------------------------------------
ggplot(ml$dat, aes(X, Y)) +
  facet_wrap(~ Year) +
  geom_tile(aes(fill = CL, color = CL)) +
  coord_fixed() +
  scale_fill_viridis() +
  scale_color_viridis()



## ----genetic-structure, include=FALSE------------------------------------
## Creates objects:
##   gen_str
##     $Ainv
##     $recoding
##     $contraint.matrix
##   prec.A
## Open with:
# file.show(file=system.file('vignettes/genetic_structure.R',
#                            package = 'RisingAshes'))
source('genetic_structure.R')


## ----spatial-structure---------------------------------------------------
## Generates an object
## sp_str
##   $loc
##   $mesh
##   $spde
##   $priors
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


## ----mst3-bin-BC-setup---------------------------------------------------

## Models for the binary component
## with Basal Circumference as explanatory variable

ml$bin_BC$dat <- 
  ml$dat[, c('Seq', 'FAM', 'PROV', 'BF98', 'Year', 'BC')] %>% 
  transform(
    Year2013 = as.numeric(Year == 2013),
    mu = 1,
    BC.cat = cut(BC, c(0, 30, 45, 60, 75, 90, 150), ordered = TRUE),
    BC.sc = scale(BC),
    recode = gen_str$recoding,
    FAM = as.numeric(FAM),
    CL0 = as.numeric(ml$dat$CL == 0)
  )

## Binomial stack
ml$bin_BC$stk <-
  inla.stack(
    data    = list(y = ml$bin_BC$dat$CL0),
    A       = list(st.Amat, 1),
    effects = list(
      suffix(st.index, 'bin'),
      suffix(ml$bin_BC$dat, 'bin')
    ),
    tag     = 'bin'
  )

## Base Formula: without BC
ml$bin_BC$fml <- 
  y ~ 0 + Year.bin + 
  f(theta.bin, model=sp_str$spde,
    group = theta.group.bin,  # variable in the mesh space
    control.group = list(model = 'exchangeable')) +
  f(recode.bin, model = "generic0", hyper = list(theta = prec.A),
    constr = TRUE, Cmatrix = Cmat)

## ----mst3-bin-BC-fit, eval = FALSE---------------------------------------
#  ## Fit models for comparison
#  ## Open with:
#  # file.show(file=system.file('vignettes/CL_bin_vs_BC_models.R',
#  #                            package = 'RisingAshes'))
#  source('CL_bin_vs_BC_models.R')

## ----mst3-bin-BC-remote-retrieve, eval = FALSE---------------------------
#  ## Retrieve results from server
#  ml$bin_BC$res <- c(lapply(ml$bin_BC$res, inla.qget))

## ----mst3-bin-BC-cache-save-recover--------------------------------------

## Load into ml$bin_BC the cached version of the two previous chunks
ml$bin_BC <- RisingAshes:::BCmodels_CLbin


## ----mst3-bin-BC-model-comparison, fig.cap = 'Relative performance of models. The reference model without genetic and spatio-temporal effects is off the scale.'----

lapply(ml$bin_BC$res, summary)

binBC_mcdat <- 
  data.frame(
    DIC = sapply(ml$bin_BC$res,
                 function(x) x$dic$dic),
    p_D = sapply(ml$bin_BC$res,
                 function(x) x$dic$p.eff),
    WAIC = sapply(ml$bin_BC$res,
                  function(x) x$waic$waic),
    p_W = sapply(ml$bin_BC$res,
                 function(x) x$waic$p.eff),
    marginal_loglik = sapply(ml$bin_BC$res,
                             function(x) x$mlik[1, 1]),
    LPML = sapply(ml$bin_BC$res,
                  function(x) sum(log(x$cpo$cpo), na.rm = TRUE))
  )

kable(
  binBC_mcdat,
  digits = 0
)

ggplot(binBC_mcdat[-1, ], aes(DIC, marginal_loglik)) +
  geom_point() + 
  geom_text(aes(label = rownames(binBC_mcdat)[-1]), hjust = -.1)


## ----mst3-bin-BC-ffect-estimates-latent-scale, fig.width = 6, fig.cap = 'Baseline effect of the Basal Circumference by Year and model in the latent scale.'----

BCdat <- local({
  ## Values of Basal Circumference to evaluate predictions
  x <- seq(min(ml$bin_BC$dat$BC, na.rm = TRUE),
           max(ml$bin_BC$dat$BC, na.rm = TRUE),
           length = 51)
  
  
  ## Predictions of the categorical model in the values of x
  x.cat <-  cut(x, c(0, 30, 45, 60, 75, 90, 150), ordered = TRUE)
  
  cat.values <- 
    with(
      ml$bin_BC$res$cat,
      c(summary.fixed[1:6, 'mean'],
        summary.fixed[7, 'mean'] + summary.fixed[1:6, 'mean'])
    )
  
  cat <- c(cat.values[x.cat], cat.values[7:12][x.cat])
  
  catxyear.values <- 
    with(
      ml$bin_BC$res$catxyear,
      as.vector(t(matrix(summary.fixed[, 'mean'], nrow = 2)))
    )
  
  catxyear <- c(catxyear.values[x.cat], catxyear.values[7:12][x.cat])
  
  cat2.values <-
    with(ml$bin_BC$res$cat2,
         as.vector(t(matrix(summary.fixed[, 'mean'], nrow = 2)))
    )
  
  cat2 <- c(cat2.values[cut(x, c(0, 30, 150))],
            cat2.values[3:4][cut(x, c(0, 30, 150))])
  
  ## Predictions of the linear and quadratic models in the values of x
  x.lin <- 
    with(
      attributes(scale(ml$bin_BC$dat$BC)),
      scale(x,
            center = `scaled:center`,
            scale = `scaled:scale`)
    )
  
  lin <- 
    with(
      ml$bin_BC$res$lin,
      c(summary.fixed['Year.bin2012', 'mean'] + 
          summary.fixed['BC.sc.bin', 'mean'] * x.lin,
        summary.fixed['Year.bin2013', 'mean'] + 
          summary.fixed['BC.sc.bin', 'mean'] * x.lin)
    )
  
  quad <- 
    with(
      ml$bin_BC$res$quad,
      c(summary.fixed['Year.bin2012', 'mean'] + 
          summary.fixed['BC.sc.bin', 'mean'] * x.lin +
          summary.fixed['I(BC.sc.bin^2)', 'mean'] * x.lin^2,
        summary.fixed['Year.bin2013', 'mean'] + 
          summary.fixed['BC.sc.bin', 'mean'] * x.lin +
          summary.fixed['I(BC.sc.bin^2)', 'mean'] * x.lin^2)
    )
  
  np  <-
    with(
      ml$bin_BC$res$np,
      c(summary.fixed['Year.bin2012', 'mean'] +
          summary.random$`cut(BC.sc.bin, 51)`[, 'mean'],
        summary.fixed['Year.bin2013', 'mean'] +
          summary.random$`cut(BC.sc.bin, 51)`[, 'mean'])
    )
  
  data.frame(
    Year = factor(rep(2012:2013, each = 51)),
    BC   = rep(x, 2),
    cat,
    cat2,
    catxyear,
    lin,
    quad,
    np) %>% 
    gather(model, y, cat:np)
})

ggplot(BCdat, aes(BC, y, colour = model)) +
  geom_line() + 
  facet_wrap(~ Year)


## ----mst3-bin-BC-effect-estimates-prob-scale, fig.width = 6, fig.cap = 'Baseline effect of the Basal Circumference by Year and model in the probability scale. The black dots are the observed proportions of zeros.'----
## Back-transform the posterior means, into the prob. scale
## Note this is not the posterior mean in that scale, but approx.

invlogit <- function(x) exp(x)/(1+exp(x))


tdat <- filter(mdat,
               variable %in% c('CD', 'CL', 'BC')
               ,Year %in% c('2012', '2013')
) %>%
  spread(variable, value) %>% 
  transform(CD = as.ordered(CD),
            BC = cut(BC, c(0, 30, 45, 60, 75, 90, 150), ordered = TRUE),
            CLbin = factor(CL > 0, labels = c('0', '>0'))) %>%
  filter(!is.na(CL), !is.na(CD)) %>% 
  droplevels() %>% 
  group_by(Year, BC) %>% 
  summarise(NCL0 = sum(CL == 0),
            TCL = n()) %>% 
  mutate(y = NCL0/TCL,
         PT.se = sqrt(y*(1-y)/TCL),
         BC = (15 * (.5 + 1:6))[BC])

ggplot(mutate(BCdat, y = invlogit(y)), aes(BC, y)) +
  geom_line(aes(colour = model)) + 
  geom_point(dat = tdat, stat = 'identity') +
  facet_wrap(~ Year) +
  ylab("Probability")



## ----mst3-bin-PROV-setup-------------------------------------------------

ml$bin_PROV <- ml$bin_BC[c('dat', 'stk', 'fml')]


## ----mst3-bin-PROV-fit, eval = FALSE-------------------------------------
#  ## Run models for comparison
#  source('CL_bin_vs_PROV_models.R')

## ----mst3-bin-PROV-remote-retrieve, eval = FALSE-------------------------
#  ## Retrieve results from server
#  ml$bin_PROV$res <- c(lapply(ml$bin_PROV$res, inla.qget))
#  
#  ## Summary of spatial results
#  ml$bin_PROV$spderes <-
#    lapply(ml$bin_PROV$res[-1], inla.spde2.result, 'theta.bin', sp_str$spde)

## ----mst3-bin-PROV-cache-save-recover------------------------------------

## Load into ml$bin_PROV the cached version of the previous chunks
ml$bin_PROV <- RisingAshes:::PROVmodels_CLbin


## ----mst3-bin-PROV-model-comparison, fig.cap = 'Relative performance of models. The reference model without genetic and spatio-temporal effects is off the scale.'----

lapply(ml$bin_PROV$res, summary)

binPROV_mcdat <- 
    data.frame(
    DIC = sapply(ml$bin_PROV$res,
                 function(x) x$dic$dic),
    p_D = sapply(ml$bin_PROV$res,
                 function(x) x$dic$p.eff),
    WAIC = sapply(ml$bin_PROV$res,
                  function(x) x$waic$waic),
    p_W = sapply(ml$bin_PROV$res,
                 function(x) x$waic$p.eff),
    marginal_loglik = sapply(ml$bin_PROV$res,
                      function(x) x$mlik[1, 1]),
    LPML = sapply(ml$bin_PROV$res,
                      function(x) sum(log(x$cpo$cpo), na.rm = TRUE))
  )

kable(
  binPROV_mcdat,
  digits = 0
)

ggplot(binPROV_mcdat[-1, ], aes(DIC, marginal_loglik)) +
  geom_point() + 
  geom_text(aes(label = rownames(binPROV_mcdat)[-1]), hjust = -.1)


## ----mst3-bin-PROV-ffect-estimates-latent-scale, fig.width = 6, fig.cap = 'Baseline effect of the Provenance by Year and model in the latent scale.'----

PROVdat <- local({

  ## Predictions of the categorical model
  
  cat <- with(
    ml$bin_PROV$res$cat,
    rbind(summary.fixed[1:3, c('mean', '0.025quant', '0.975quant')],
          summary.fixed[4, 'mean'] +
            summary.fixed[1:3, c('mean', '0.025quant', '0.975quant')])
  )
  

  catxyear <- with(
    ml$bin_PROV$res$catxyear,
        summary.fixed[t(matrix(1:6, 2)),
                      c('mean', '0.025quant', '0.975quant')]
  )
  

  ## mean of the sum of ST and genetic random effects
  ## for 2012 and 2013, in the cat model.
  center.ranef <-
    sapply(
    1:2,
    function(x) 
      mean(ml$bin_PROV$spderes$cat$summary.values$mean[st.index$theta.group == x],
           na.rm = TRUE)) +
    mean(ml$bin_PROV$res$cat$summary.random$recode.bin$mean[gen_str$recoding], 
         na.rm = TRUE)
  
  

  data.frame(
    Year = factor(rep(2012:2013, each = 3)),
    centering = rep(center.ranef, each = 3),
    PROV   = factor(rep(levels(ml$dat$PROV), 2)),
    cat = cat,
    cy = catxyear ) %>% 
    gather(model, y, cat.mean:cy.0.975quant) %>% 
    separate(model, c('model', 'q'), extra = 'merge') %>% 
    mutate(y = y + centering) %>% 
    spread(q, y)
})

ggplot(PROVdat, aes(PROV, mean, colour = model)) +
  geom_pointrange(aes(ymin = `0.025quant`, ymax = `0.975quant`),
                  position = position_dodge(width = 0.5)) +
  facet_wrap(~ Year)


## ----mst3-bin-PROV-effect-estimates-prob-scale, fig.width = 6, fig.cap = 'Baseline effect of the Provenance by Year and model in the probability scale. The black dots are the observed proportions of zeros.'----
## Back-transform the posterior means, into the prob. scale
## Note this is not the posterior mean in that scale, but approx.

invlogit <- function(x) exp(x)/(1+exp(x))


tdat <- filter(mdat,
               variable %in% c('CD', 'CL')
               ,Year %in% c('2012', '2013')
) %>%
  spread(variable, value) %>% 
  transform(CLbin = factor(CL > 0, labels = c('0', '>0'))) %>%
  filter(!is.na(CL), !is.na(PROV)) %>% 
  droplevels() %>% 
  group_by(Year, PROV) %>% 
  summarise(NCL0 = sum(CL == 0),
            TCL = n()) %>% 
  mutate(y = NCL0/TCL,
         PT.se = sqrt(y*(1-y)/TCL))

ggplot(mutate(PROVdat,
              y = invlogit(mean)), aes(PROV, y)) +
  geom_bar(aes(fill = model), stat = 'identity', position = 'dodge') + 
  geom_point(dat = tdat, stat = 'identity') +
  facet_wrap(~ Year) +
  ylab("Probability")



## ----mst3-bin-BF-setup---------------------------------------------------

ml$bin_BF <- ml$bin_BC[c('dat', 'stk', 'fml')]


## ----mst3-bin-BF-fit, eval = FALSE---------------------------------------
#  ## Run models for comparison
#  source('CL_bin_vs_BF_models.R')

## ----mst3-bin-BF-remote-retrieve, eval = FALSE---------------------------
#  
#  ml$bin_BF$res <- c(lapply(ml$bin_BF$res, inla.qget))
#  
#  ## Summary of spatial results
#  ml$bin_BF$spderes <- lapply(ml$bin_BF$res[-1], inla.spde2.result, 'theta.bin', sp_str$spde)

## ----mst3-bin-BF-cache-save-recover--------------------------------------

## Load into ml$bin_BF the cached version of the previous chunks
ml$bin_BF <- RisingAshes:::BFmodels_CLbin


## ----mst3-bin-BF-model-comparison, fig.cap = 'Relative performance of models. The reference model without genetic and spatio-temporal effects is off the scale.'----

lapply(ml$bin_BF$res, summary)

BF_mcdat <- 
    data.frame(
    DIC = sapply(ml$bin_BF$res,
                 function(x) x$dic$dic),
    WAIC = sapply(ml$bin_BF$res,
                  function(x) x$waic$waic),
    marginal_loglik = sapply(ml$bin_BF$res,
                      function(x) x$mlik[1, 1]),
    LPML = sapply(ml$bin_BF$res,
                      function(x) sum(log(x$cpo$cpo), na.rm = TRUE))
  )

kable(
  BF_mcdat,
  digits = 0
)

ggplot(BF_mcdat[-1, ], aes(DIC, marginal_loglik)) +
  geom_point() + 
  geom_text(aes(label = rownames(BF_mcdat)[-1]), hjust = -.1)


## ----mst3-bin-BF-ffect-estimates-latent-scale, fig.width = 6, fig.cap = 'Baseline effect of the Bud flush by Year and model in the latent scale.'----

BFdat <- local({

  nlevels <- 5
  ## Predictions of the categorical model
  
  cat <- with(
    ml$bin_BF$res$cat,
    rbind(summary.fixed[1:nlevels, c('mean', '0.025quant', '0.975quant')],
          summary.fixed[nlevels+1, 'mean'] +
            summary.fixed[1:nlevels, c('mean', '0.025quant', '0.975quant')])
  )
  

  catxyear <- with(
    ml$bin_BF$res$catxyear,
        summary.fixed[t(matrix(1:(2*nlevels), 2)),
                      c('mean', '0.025quant', '0.975quant')]
  )
  
  
  center.ranef <-
    sapply(
      1:2,
      function(x) 
        mean(ml$bin_BF$spderes$cat$summary.values$mean[st.index$theta.group == x],
             na.rm = TRUE)) +
    mean(ml$bin_BF$res$cat$summary.random$recode.bin$mean[gen_str$recoding],
         na.rm = TRUE)
  
  

  data.frame(
    Year = factor(rep(2012:2013, each = nlevels)),
    centering = rep(center.ranef, each = nlevels),
    BF   = factor(rep(levels(ml$dat$BF98), 2)),
    cat = cat,
    cy = catxyear ) %>% 
    gather(model, y, cat.mean:cy.0.975quant) %>% 
    separate(model, c('model', 'q'), extra = 'merge') %>% 
    mutate(y = y + centering) %>% 
    spread(q, y)
})

ggplot(BFdat, aes(BF, mean, colour = model)) +
  geom_pointrange(aes(ymin = `0.025quant`, ymax = `0.975quant`),
                  position = position_dodge(width = 0.5)) +
  facet_wrap(~ Year)


## ----mst3-bin-BF-effect-estimates-prob-scale, fig.width = 6, fig.cap = 'Baseline effect of the Bud flush by Year and model in the probability scale. The black dots are the observed proportions of zeros.'----
## Back-transform the posterior means, into the prob. scale
## Note this is not the posterior mean in that scale, but approx.

invlogit <- function(x) exp(x)/(1+exp(x))


tdat <- filter(mdat,
               variable %in% c('CD', 'CL')
               ,Year %in% c('2012', '2013')
) %>%
  spread(variable, value) %>% 
  transform(CLbin = factor(CL > 0, labels = c('0', '>0'))) %>%
  filter(!is.na(CL), !is.na(BF98)) %>% 
  droplevels() %>% 
  group_by(Year, BF98) %>% 
  summarise(NCL0 = sum(CL == 0),
            TCL = n()) %>% 
  mutate(y = NCL0/TCL,
         PT.se = sqrt(y*(1-y)/TCL),
         BF = BF98)

ggplot(mutate(BFdat,
              y = invlogit(mean)), aes(BF, y)) +
  geom_bar(aes(fill = model), stat = 'identity', position = 'dodge') + 
  geom_point(dat = tdat, stat = 'identity') +
  facet_wrap(~ Year) +
  ylab("Probability")



## ----mst3-bin-msel-setup-------------------------------------------------

## Models for the binary component

ml$bin_msel$dat <- 
  ml$dat[, c('Seq', 'FAM', 'PROV', 'BF98', 'Year', 'BC')] %>% 
  transform(
    Year2013 = as.numeric(Year == 2013),
    mu = 1,
    BC = cut(BC, c(0, 30, 45, 60, 75, 90, 150), ordered = TRUE),
    BF_2 = as.numeric(BF98 == 2),
    BF_3 = as.numeric(BF98 == 3),
    BF_4 = as.numeric(BF98 == 4),
    BF_5 = as.numeric(BF98 == 5),
    recode = gen_str$recoding,
    FAM = as.numeric(FAM),
    CL0 = as.numeric(ml$dat$CL == 0))

## Binomial stack
ml$bin_msel$stk <-
  inla.stack(data    = list(y = ml$bin_msel$dat$CL0),
             A       = list(st.Amat, 1),
             effects = list(suffix(st.index, 'bin'),
                            suffix(ml$bin_msel$dat, 'bin')),
             tag     = 'bin')

## Base Formula: without BC
ml$bin_msel$fml <- y ~ 0 + Year.bin + 
               f(theta.bin, model=sp_str$spde,
                 group = theta.group.bin,  # variable in the mesh space
                 control.group = list(model = 'exchangeable')) +
               f(recode.bin, model = "generic0", hyper = list(theta = prec.A),
                 constr = TRUE, Cmatrix = Cmat)



## ----mst3-bin-msel-fit, eval = FALSE-------------------------------------
#  
#  ## Run the full model
#  source('CL_bin_models.R')

## ----mst3-bin-msel-remote-retrieve, eval = FALSE-------------------------
#  
#  ml$bin_msel$res <- c(lapply(ml$bin_msel$res, inla.qget))
#  
#  

## ----mst3-bin-msel-cache-save-recover------------------------------------

## Load into ml$bin_msel$res$full the cached version of the two previous chunks
ml$bin_msel$res$full <- RisingAshes:::full_CLbin

## These were already fitted
ml$bin_msel$res <- c(list(
  null = ml$bin_BC$res$noBC,
  BC   = ml$bin_BC$res$catxyear,
  BF   = ml$bin_BC$res$cat
),
ml$bin_msel$res
)

## It would be good to fit also all the combinations with two of those


## ----mst3-bin-msel-model-comparison, fig.cap = 'Relative performance of models for the binary component.'----

lapply(ml$bin_msel$res, summary)

BF_mcdat <- 
    data.frame(
    DIC = sapply(ml$bin_msel$res,
                 function(x) x$dic$dic),
    WAIC = sapply(ml$bin_msel$res,
                  function(x) x$waic$waic),
    marginal_loglik = sapply(ml$bin_msel$res,
                      function(x) x$mlik[1, 1]),
    LPML = sapply(ml$bin_msel$res,
                      function(x) sum(log(x$cpo$cpo), na.rm = TRUE))
  )

kable(
  BF_mcdat,
  digits = 0
)

ggplot(BF_mcdat, aes(DIC, marginal_loglik)) +
  geom_point() + 
  geom_text(aes(label = rownames(BF_mcdat)), hjust = -.1)


## ----mst3-cont-data------------------------------------------------------

ml$cont$dat <- 
  ml$dat[, c('Seq', 'FAM', 'PROV', 'BF98', 'Year', 'BC')] %>% 
  transform(
    Year2013 = as.numeric(Year == 2013),
    BC.cat = cut(BC, c(0, 30, 45, 60, 75, 90, 150), ordered = TRUE),
    BC.sc = scale(BC),
    BF_2 = as.numeric(BF98 == 2),
    BF_3 = as.numeric(BF98 == 3),
    BF_4 = as.numeric(BF98 == 4),
    BF_5 = as.numeric(BF98 == 5),
    recode = gen_str$recoding,
    CL0 = as.numeric(ml$dat$CL == 0))

ml$cont$dat$CL_gt0   		<- ml$dat$CL
ml$cont$dat$CL_gt0[ml$dat$CL == 0] <- NA

ml$cont$dat$yt <- ml$cont$dat$CL_gt0
ml$cont$dat$yt[ml$dat$CL == 1] <- .98   # forbiden value for a beta


## ----mst3-cont-liksel----------------------------------------------------

## Use the same covariates as before

fit_lik <- function(lik) {
  fml <- yt ~ 0 + FAM + PROV + Year:BC
  inla(formula = fml,
       data = filter(ml$cont$dat,
                     !is.na(BF98),
                     !is.na(FAM)),
       family = lik,
       control.predictor = list(
         compute = TRUE, link = 1),
       control.fixed = list(
         expand.factor.strategy = 'inla'), # account for NAs in factor
       control.compute = list(
         dic = TRUE,
         waic = TRUE,
         cpo = TRUE,
         po = TRUE,
         config = TRUE)
       # , inla.call = 'remote'
  )
}

liklst <- c('beta', 'gamma')
cont_res <- lapply(liklst, fit_lik)
names(cont_res) <- liklst

cont_comp <- 
    data.frame(
    DIC = sapply(cont_res,
                 function(x) x$dic$dic),
    p_D = sapply(cont_res,
                 function(x) x$dic$p.eff),
    WAIC = sapply(cont_res,
                  function(x) x$waic$waic),
    p_W = sapply(cont_res,
                 function(x) x$waic$p.eff),
    marginal_loglik = sapply(cont_res,
                      function(x) x$mlik[1, 1]),
    LPML = sapply(cont_res,
                      function(x) sum(log(x$cpo$cpo), na.rm = TRUE))
  )

kable(
  cont_comp,
  digits = 0
)


## ----mst3-cont-setup-----------------------------------------------------

## Models for the continuous component
## with Basal Circumference as explanatory variable


## continuous stack

ml$cont$stk <-
  inla.stack(
    data    = list(y = ml$cont$dat$yt),
    A       = list(st.Amat, 1),
    effects = list(suffix(st.index, 'cont'),
                   suffix(ml$cont$dat, 'cont')),
    tag     = 'cont')


## Prior on the precision of the overdispersion term
## We use a Gamma with shape 0.5 and rate 0.2, which gives a
## prior expectation (in the log scale) of 2.5 and a reasonably ample
## variance of 12.5.
## Furthermore, we checked that sensitivity to this prior specification is very
## limited.
prec.od <- list(fixed = FALSE,
                initial = 1,
                param = c(.5, .2))


## Base Formula: without any covariate but Year
ml$cont$fml <- 
  y ~ 0 + Year.cont + 
  f(theta.cont, model=sp_str$spde,
    group = theta.group.cont,  # variable in the mesh space
    control.group = list(model = 'exchangeable')) +
  f(recode.cont, model = "generic0", hyper = list(theta = prec.A),
    constr = TRUE, Cmatrix = Cmat) #+
  # f(Seq.cont, model = 'iid', hyper = list(theta = prec.od))



## ----mst3-cont-BC-fit, eval = FALSE--------------------------------------
#  ## Run models for comparison
#  source('CL_cont_vs_BC_models.R')

## ----mst3-cont-BC-remote-retrieve, eval = FALSE--------------------------
#  
#  ml$cont_BC$res <- c(lapply(ml$cont_BC$res, inla.qget))
#  

## ----mst3-cont-BC-cache-save-recover-------------------------------------

## Load into ml$cont_BC the cached version of the previous chunks
ml$cont_BC <- RisingAshes:::BCmodels_CLcont


## ----mst3-cont-BC-model-comparison, fig.cap = 'Relative performance of models. The reference model without genetic and spatio-temporal effects is off the scale.'----

lapply(ml$cont_BC$res, summary)

contBC_mcdat <- 
    data.frame(
    DIC = sapply(ml$cont_BC$res,
                 function(x) x$dic$dic),
    p_D = sapply(ml$cont_BC$res,
                 function(x) x$dic$p.eff),
    WAIC = sapply(ml$cont_BC$res,
                  function(x) x$waic$waic),
    p_W = sapply(ml$cont_BC$res,
                 function(x) x$waic$p.eff),
    marginal_loglik = sapply(ml$cont_BC$res,
                      function(x) x$mlik[1, 1]),
    LPML = sapply(ml$cont_BC$res,
                      function(x) sum(log(x$cpo$cpo), na.rm = TRUE))
  )

kable(
  contBC_mcdat,
  digits = 0
)

ggplot(contBC_mcdat[-1, ], aes(DIC, marginal_loglik)) +
  geom_point() + 
  geom_text(aes(label = rownames(contBC_mcdat)[-1]), hjust = -.1)


## ----mst3-cont-BC-ffect-estimates-latent-scale, fig.width = 6, fig.cap = 'Baseline effect of the Basal Circumference by Year and model in the latent scale.'----

BCdat <- local({
  ## Values of Basal Circumference to evaluate predictions
  x <- seq(min(ml$cont$dat$BC, na.rm = TRUE),
           max(ml$cont$dat$BC, na.rm = TRUE), 
           length = 51)
  
  
  ## Predictions of the categorical model in the values of x
  x.cat <-  cut(x, c(0, 30, 45, 60, 75, 90, 150), ordered = TRUE)
  
  cat.values = with(ml$cont_BC$res$cat,
                 c(summary.fixed[1:6, 'mean'],
                   summary.fixed[7, 'mean'] + summary.fixed[1:6, 'mean'])
  )
  
  cat <- c(cat.values[x.cat], cat.values[7:12][x.cat])
  

  catxyear.values = with(ml$cont_BC$res$catxyear,
                    as.vector(t(matrix(summary.fixed[, 'mean'], nrow = 2)))
  )
  
  catxyear <- c(catxyear.values[x.cat], catxyear.values[7:12][x.cat])
  
  cat2.values = with(ml$cont_BC$res$cat2,
                         as.vector(t(matrix(summary.fixed[, 'mean'], nrow = 2)))
  )
  
  cat2 <- c(cat2.values[cut(x, c(0, 30, 150))], cat2.values[3:4][cut(x, c(0, 30, 150))])
  

  ## Predictions of the linear and quadratic models in the values of x
  x.lin <- with(
    attributes(scale(ml$cont$dat$BC)),
    scale(x,
          center = `scaled:center`,
          scale = `scaled:scale`)
  )
  
  lin = with(
    ml$cont_BC$res$lin,
    c(summary.fixed['Year.cont2012', 'mean'] + 
        summary.fixed['BC.sc.cont', 'mean'] * x.lin,
      summary.fixed['Year.cont2013', 'mean'] + 
        summary.fixed['BC.sc.cont', 'mean'] * x.lin)
  )
  
  quad = with(
    ml$cont_BC$res$quad,
    c(summary.fixed['Year.cont2012', 'mean'] + 
        summary.fixed['BC.sc.cont', 'mean'] * x.lin +
        summary.fixed['I(BC.sc.cont^2)', 'mean'] * x.lin^2,
      summary.fixed['Year.cont2013', 'mean'] + 
        summary.fixed['BC.sc.cont', 'mean'] * x.lin +
        summary.fixed['I(BC.sc.cont^2)', 'mean'] * x.lin^2)
  )
  
  
  np = with(
    ml$cont_BC$res$np,
    c(summary.fixed['Year.cont2012', 'mean'] +
        summary.random$`cut(BC.sc.cont, 51)`[, 'mean'],
      summary.fixed['Year.cont2013', 'mean'] +
        summary.random$`cut(BC.sc.cont, 51)`[, 'mean'])
  )
  
  
  data.frame(
    Year = factor(rep(2012:2013, each = 51)),
    BC   = rep(x, 2),
    cat,
    cat2,
    catxyear,
    lin,
    quad,
    np) %>% 
    gather(model, y, cat:np)
})

ggplot(BCdat, aes(BC, y, colour = model)) +
  geom_line() + 
  facet_wrap(~ Year)


## ----mst3-cont-BC-effect-estimates-prob-scale, fig.width = 6, fig.cap = 'Baseline effect of the Basal Circumference by Year and model in the observed scale. The black dots are the observed mean of CL|CL>0.'----
## Back-transform the posterior means, into the prob. scale
## Note this is not the posterior mean in that scale, but approx.

tdat <- 
  filter(mdat,
         variable %in% c('CD', 'CL', 'BC')
         ,Year %in% c('2012', '2013')
  ) %>%
  spread(variable, value) %>% 
  transform(CD = as.ordered(CD),
            BC = cut(BC, c(0, 30, 45, 60, 75, 90, 150), ordered = TRUE),
            CLcont = factor(CL > 0, labels = c('0', '>0'))) %>%
  filter(!is.na(CL), !is.na(CD)) %>% 
  droplevels() %>% 
  group_by(Year, BC) %>% 
  summarise(y = mean(CL[CL > 0]),
            n = n()) %>% 
  mutate(BC = (15 * (.5 + 1:6))[BC])

ggplot(mutate(BCdat, y = exp(y)), aes(BC, y)) +
  geom_line(aes(colour = model)) + 
  geom_point(aes(size = n), dat = tdat, stat = 'identity') +
  facet_wrap(~ Year) +
  ylab("CL | CL > 0")



## ----mst3-cont-PROV-fit, eval = FALSE------------------------------------
#  ## Run models for comparison
#  source('CL_cont_vs_PROV_models.R')

## ----mst3-cont-PROV-remote-retrieve, eval = FALSE------------------------
#  
#  ml$cont_PROV$res <- c(lapply(ml$cont_PROV$res, inla.qget))
#  
#  ## Summary of spatial results
#  ml$cont_PROV$spderes <- lapply(ml$cont_PROV$res[-1], inla.spde2.result, 'theta.cont', sp_str$spde)

## ----mst3-cont-PROV-cache-save-recover-----------------------------------

## Load into ml$cont_PROV the cached version of the previous chunks
ml$cont_PROV <- RisingAshes:::PROVmodels_CLcont


## ----mst3-cont-PROV-model-comparison, fig.cap = 'Relative performance of models. The reference model without genetic and spatio-temporal effects is off the scale.'----

lapply(ml$cont_PROV$res, summary)

contPROV_mcdat <- 
    data.frame(
    DIC = sapply(ml$cont_PROV$res,
                 function(x) x$dic$dic),
    p_D = sapply(ml$cont_PROV$res,
                 function(x) x$dic$p.eff),
    WAIC = sapply(ml$cont_PROV$res,
                  function(x) x$waic$waic),
    p_W = sapply(ml$cont_PROV$res,
                 function(x) x$waic$p.eff),
    marginal_loglik = sapply(ml$cont_PROV$res,
                      function(x) x$mlik[1, 1]),
    LPML = sapply(ml$cont_PROV$res,
                      function(x) sum(log(x$cpo$cpo), na.rm = TRUE))
  )

kable(
  contPROV_mcdat,
  digits = 0
)

ggplot(contPROV_mcdat[-1, ], aes(DIC, marginal_loglik)) +
  geom_point() + 
  geom_text(aes(label = rownames(contPROV_mcdat)[-1]), hjust = -.1)


## ----mst3-cont-PROV-ffect-estimates-latent-scale, fig.width = 6, fig.cap = 'Baseline effect of the Provenance by Year and model in the latent scale.'----

PROVdat <- local({

  ## Predictions of the categorical model
  
  cat <- with(
    ml$cont_PROV$res$cat,
    rbind(summary.fixed[1:3, c('mean', '0.025quant', '0.975quant')],
          summary.fixed[4, 'mean'] +
            summary.fixed[1:3, c('mean', '0.025quant', '0.975quant')])
  )
  

  catxyear <- with(
    ml$cont_PROV$res$catxyear,
        summary.fixed[t(matrix(1:6, 2)),
                      c('mean', '0.025quant', '0.975quant')]
  )
  

  center.ranef <-
    sapply(
    1:2,
    function(x) 
      mean(ml$cont_PROV$spderes$cat$summary.values$mean[st.index$theta.group == x],
           na.rm = TRUE)) +
    mean(ml$cont_PROV$res$cat$summary.random$recode.cont$mean[gen_str$recoding],
           na.rm = TRUE)
  
  

  data.frame(
    Year = factor(rep(2012:2013, each = 3)),
    centering = rep(center.ranef, each = 3),
    PROV   = factor(rep(levels(ml$dat$PROV), 2)),
    cat = cat,
    cy = catxyear ) %>% 
    gather(model, y, cat.mean:cy.0.975quant) %>% 
    separate(model, c('model', 'q'), extra = 'merge') %>% 
    mutate(y = y + centering) %>% 
    spread(q, y)
})

ggplot(PROVdat, aes(PROV, mean, colour = model)) +
  geom_pointrange(aes(ymin = `0.025quant`, ymax = `0.975quant`),
                  position = position_dodge(width = 0.5)) +
  facet_wrap(~ Year)


## ----mst3-cont-PROV-effect-estimates-prob-scale, fig.width = 6, fig.cap = 'Baseline effect of the Provenance by Year and model in the probability scale. The black dots are the observed proportions of zeros.'----
## Back-transform the posterior means, into the prob. scale
## Note this is not the posterior mean in that scale, but approx.


tdat <- filter(mdat,
               variable %in% c('CD', 'CL')
               ,Year %in% c('2012', '2013')
) %>%
  spread(variable, value) %>% 
  transform(CLcont = factor(CL > 0, labels = c('0', '>0'))) %>%
  filter(!is.na(CL), !is.na(PROV)) %>% 
  droplevels() %>% 
  group_by(Year, PROV) %>% 
  summarise(y = mean(CL[CL > 0]),
            n = n())

ggplot(mutate(PROVdat, y = exp(mean)), aes(PROV, y)) +
  geom_bar(aes(fill = model), stat = 'identity', position = 'dodge') + 
  geom_point(dat = tdat, stat = 'identity') +
  facet_wrap(~ Year) +
  ylab("CL | CL > 0")



## ----mst3-cont-BF-fit, eval = FALSE--------------------------------------
#  ## Run models for comparison
#  source('CL_cont_vs_BF_models.R')

## ----mst3-cont-BF-remote-retrieve, eval = FALSE--------------------------
#  
#  ml$cont_BF$res <- c(lapply(ml$cont_BF$res, inla.qget))
#  
#  ## Summary of spatial results
#  ml$cont_BF$spderes <- lapply(ml$cont_BF$res[-1], inla.spde2.result, 'theta.cont', sp_str$spde)

## ----mst3-cont-BF-cache-save-recover-------------------------------------

# saveRDS(ml$cont_BF, file.path(path, 'cache', 'CL.ml.cont_BF.rds'))
ml$cont_BF <- RisingAshes:::BFmodels_CLcont


## ----mst3-cont-BF-model-comparison, fig.cap = 'Relative performance of models. The reference model without genetic and spatio-temporal effects is off the scale.'----

lapply(ml$cont_BF$res, summary)

BF_mcdat <- 
    data.frame(
    DIC = sapply(ml$cont_BF$res,
                 function(x) x$dic$dic),
    WAIC = sapply(ml$cont_BF$res,
                  function(x) x$waic$waic),
    marginal_loglik = sapply(ml$cont_BF$res,
                      function(x) x$mlik[1, 1]),
    LPML = sapply(ml$cont_BF$res,
                      function(x) sum(log(x$cpo$cpo), na.rm = TRUE))
  )

kable(
  BF_mcdat,
  digits = 0
)

ggplot(BF_mcdat[-1, ], aes(DIC, marginal_loglik)) +
  geom_point() + 
  geom_text(aes(label = rownames(BF_mcdat)[-1]), hjust = -.1)


## ----mst3-cont-BF-ffect-estimates-latent-scale, fig.width = 6, fig.cap = 'Baseline effect of the Bud flush by Year and model in the latent scale.'----

BFdat <- local({

  nlevels <- 5
  ## Predictions of the categorical model
  
  cat <- with(
    ml$cont_BF$res$cat,
    rbind(summary.fixed[1:nlevels, c('mean', '0.025quant', '0.975quant')],
          summary.fixed[nlevels+1, 'mean'] +
            summary.fixed[1:nlevels, c('mean', '0.025quant', '0.975quant')])
  )
  

  catxyear <- with(
    ml$cont_BF$res$catxyear,
        summary.fixed[t(matrix(1:(2*nlevels), 2)),
                      c('mean', '0.025quant', '0.975quant')]
  )
  

  center.ranef <-
    sapply(
    1:2,
    function(x) 
      mean(ml$cont_BF$spderes$cat$summary.values$mean[st.index$theta.group == x],
           na.rm = TRUE)) +
    mean(ml$cont_BF$res$cat$summary.random$recode.cont$mean[gen_str$recoding],
           na.rm = TRUE)
  
  

  data.frame(
    Year = factor(rep(2012:2013, each = nlevels)),
    centering = rep(center.ranef, each = nlevels),
    BF   = factor(rep(levels(ml$dat$BF98), 2)),
    cat = cat,
    cy = catxyear ) %>% 
    gather(model, y, cat.mean:cy.0.975quant) %>% 
    separate(model, c('model', 'q'), extra = 'merge') %>% 
    mutate(y = y + centering) %>% 
    spread(q, y)
})

ggplot(BFdat, aes(BF, mean, colour = model)) +
  geom_pointrange(aes(ymin = `0.025quant`, ymax = `0.975quant`),
                  position = position_dodge(width = 0.5)) +
  facet_wrap(~ Year)


## ----mst3-cont-BF-effect-estimates-prob-scale, fig.width = 6, fig.cap = 'Baseline effect of the Bud flush by Year and model in the probability scale. The black dots are the observed proportions of zeros.'----
## Back-transform the posterior means, into the prob. scale
## Note this is not the posterior mean in that scale, but approx.

tdat <- filter(mdat,
               variable %in% c('CD', 'CL')
               ,Year %in% c('2012', '2013')
) %>%
  spread(variable, value) %>% 
  transform(CLcont = factor(CL > 0, labels = c('0', '>0'))) %>%
  filter(!is.na(CL), !is.na(BF98)) %>% 
  droplevels() %>% 
  group_by(Year, BF98) %>% 
  summarise(y = mean(CL[CL>0])) %>% 
  mutate(BF = BF98)

ggplot(mutate(BFdat,
              y = exp(mean)), aes(BF, y)) +
  geom_bar(aes(fill = model), stat = 'identity', position = 'dodge') + 
  geom_point(dat = tdat, stat = 'identity') +
  facet_wrap(~ Year) +
  ylab("CL | CL > 0")



## ----mst-setup-----------------------------------------------------------

## Dataset
ml$mst$dat <- 
  ml$dat[, c('Seq', 'FAM', 'PROV', 'BF98', 'Year', 'BC')] %>% 
  transform(
    Year2013 = as.numeric(Year == 2013),
    BC = cut(BC, c(0, 30, 45, 60, 75, 90, 150), ordered = TRUE),
    BF_2 = as.numeric(BF98 == 2),
    BF_3 = as.numeric(BF98 == 3),
    BF_4 = as.numeric(BF98 == 4),
    BF_5 = as.numeric(BF98 == 5),
    recode = gen_str$recoding,
    FAM = as.numeric(FAM),
    CL0 = as.numeric(ml$dat$CL == 0))

ml$mst$dat$CL_gt0   		<- ml$dat$CL
ml$mst$dat$CL_gt0[ml$dat$CL == 0] <- NA

ml$mst$dat$yt <- ml$mst$dat$CL_gt0
ml$mst$dat$yt[ml$dat$CL == 1] <- .98   # forbiden value for a beta


## Binomial stack
stk.bin <-
  inla.stack(
    data    = list(y = cbind(as.numeric(ml$dat$CL == 0), NA)),
    A       = list(st.Amat, 1, st.Amat, 1),
    effects = list(suffix(st.index, 'bin'),
                   suffix(ml$mst$dat, 'bin'),
                   suffix(st.index, 'cont', makeNA = TRUE),
                   suffix(ml$mst$dat, 'cont', makeNA = TRUE)),
    tag     = 'bin')


## continuous stack

stk.cont <-
  inla.stack(data    = list(y = cbind(NA, ml$mst$dat$yt)),
             A       = list(st.Amat, 1, st.Amat, 1),
             effects = list(suffix(st.index, 'bin', makeNA = TRUE),
                            suffix(ml$mst$dat, 'bin', makeNA = TRUE),
                            suffix(st.index, 'cont'),
                            suffix(ml$mst$dat, 'cont')),
             tag     = 'cont')

ml$mst$stk <- inla.stack(stk.bin, stk.cont)


## Prior on the precision of the overdispersion term
## We use a Gamma with shape 0.5 and rate 0.2, which gives a
## prior expectation (in the log scale) of 2.5 and a reasonably ample
## variance of 12.5.
## Furthermore, we checked that sensitivity to this prior specification is very
## limited.
prec.od <- list(fixed = FALSE,
                initial = 1,
                param = c(.5, .2))


## Formula
ml$mst$fml <- y ~ 0 +
  Year.bin +
  # Year.bin:BC.bin +BF_2.bin + BF_3.bin + BF_4.bin + BF_5.bin +
  f(theta.bin, model=sp_str$spde,
    group = theta.group.bin,  # variable in the mesh space
    control.group = list(model = 'exchangeable')) +
  f(recode.bin, model = "generic0", hyper = list(theta = prec.A),
    constr = TRUE, Cmatrix = Cmat) +
  # f(Seq.bin, model = "iid", hyper = list(theta = prec.od)) + # only in the continous part!!
  Year.cont +
  f(theta.cont, model=sp_str$spde,
    group = theta.group.cont,  # variable in the mesh space
    control.group = list(model = 'exchangeable')) +
  f(recode.cont, model = "generic0", hyper = list(theta = prec.A),
    constr = TRUE, Cmatrix = Cmat) 
# + f(Seq.cont, model = "iid", hyper = list(theta = prec.A))  # Causes overfitting


# Use the right link function for the prediction when the response is NA
link.NA <- rep(NA, nrow(ml$dat)*2L)
nas.idx.bin <- which(is.na(ml$dat$CL))
nas.idx.cont <- which(is.na(ml$mst$dat$yt))   # there were additional NAs
link.NA[nas.idx.bin] <- 1
link.NA[nrow(ml$dat) + nas.idx.cont] <- 2


## ----mst-fit, eval = FALSE-----------------------------------------------
#  ## Introduce the "residual" variation into a latent
#  ## overdispersion term and fix the likelihood precision
#  ## to a high value
#  ml$mst$res <- inla(
#    formula = ml$mst$fml,
#    data = inla.stack.data(ml$mst$stk),
#    family = list("binomial", "gamma"),
#    Ntrials = 1,
#    control.family = list(
#      list(),
#      list(
#        hyper = list(theta = prec.od)
#        # hyper = list(theta = list(fixed = TRUE, initial = 10))
#        # causes overfitting
#      )
#    ),
#    ## Correction of the Laplace Approximation
#    ## for Binomial data with many zeros
#    ## http://arxiv.org/abs/1503.07307
#    control.inla = list(
#      correct = TRUE,
#      correct.factor = 10),
#    control.predictor = list(
#      A = inla.stack.A(ml$mst$stk),
#      compute = TRUE, link = link.NA),
#    control.fixed = list(
#      expand.factor.strategy = 'inla'), # account for NAs in factor
#    control.compute = list(
#      dic = TRUE,
#      waic = TRUE,
#      cpo = TRUE,
#      po = TRUE,
#      config = TRUE)
#    #, control.family = list(hyper = list(prec = prec.A))
#    # , inla.call = 'submit'
#  )
#  
#  

## ----mst-remote-retrieve, eval = FALSE-----------------------------------
#  ml$mst$res <- inla.qget(ml$mst$res)
#  

## ----mst-cache-----------------------------------------------------------

ml$mst$res <- RisingAshes:::mst_CL

## Summary of spatial results
ml$mst$spde.bin <- 
  inla.spde2.result(ml$mst$res, 'theta.bin', sp_str$spde)
ml$mst$spde.cont <-
  inla.spde2.result(ml$mst$res, 'theta.cont', sp_str$spde)

## ----mst-CL-summary------------------------------------------------------
summary(ml$mst$res)

## ----mst-CL-fitvsobs-compute---------------------------------------------


bin.idx <- inla.stack.index(ml$mst$stk, tag = 'bin')$data
cont.idx <- inla.stack.index(ml$mst$stk, tag = 'cont')$data

## If I didn't use link = 2 in the INLA call:
fitted.bin <- ml$mst$res$summary.fitted.values$mean[bin.idx]
for( i in nas.idx.bin ) {
  fitted.bin[i] <- inla.emarginal(inv.logit, ml$mst$res$marginals.fitted.values[[bin.idx[i]]])
}
fitted.cont <- ml$mst$res$summary.fitted.values$mean[cont.idx]
for( i in nas.idx.cont ) {
  fitted.cont[i] <- inla.emarginal(inv.logit, ml$mst$res$marginals.fitted.values[[cont.idx[i]]])
}

fvso.dat <- rbind(
  data.frame(
    component = 'bin',
    Observed = as.numeric(ml$dat$CL == 0),
    Predicted = fitted.bin,
    Residuals = as.numeric(ml$dat$CL == 0) - fitted.bin,
    Year = ml$dat$Year
  ),
  data.frame(
    component = 'cont',
    Observed = ml$mst$dat$yt,
    Predicted = fitted.cont,
    Residuals = ml$mst$dat$yt - fitted.cont,
    Year = ml$dat$Year
  )
)

g <- ggplot(fvso.dat, aes(Observed, Predicted)) +
  geom_violin(aes(group = Observed),
              data = subset(fvso.dat, component == 'bin'),
              fill = 'grey80') + 
  geom_point(aes(x=Observed),
             data = subset(fvso.dat, component == 'cont')) + 
  geom_abline(intercept = 0, slope = 1, col = 'darkgray',
              data = subset(fvso.dat, component == 'cont')) +
  facet_wrap(~ component, scales = 'free_x')


## ----mst-CL-fitvsobs, fig.cap = 'Fitted vs. observed values'-------------
## Tweak the x-axis for the left facet
## ggplot under the hood:
## http://zevross.com/blog/2014/11/20/under-the-hood-of-ggplot2-graphics-in-r/
gTable <- ggplot_gtable(ggplot_build(g))
## Determine the right grob
# pdf("grobs.pdf")
# for(i in 1:length(gTable$grobs)){
#   grid.draw(gTable$grobs[[i]])
#   grid.text(i, x = unit(0.1, "npc"), y = unit(0.1, "npc"))
#   grid.newpage()
# }
# dev.off()
gTable$grobs[[8]]$children[2]$axis$grobs[[2]]$children[[1]]$label <- 
  c('', '> 0', '', 'Zero', '')
gTable$grobs[[8]]$children[2]$axis$grobs[[1]]$y <- 
  unit(0, units = 'cm')
grid.draw(gTable)

## ----resids-vs-BC, fig.cap = 'Residuals vs Basal Circumference.'---------
ggplot(cbind(fvso.dat, BC = ml$dat$BC), aes(BC, Residuals, colour = Year)) +
  geom_point() +
  facet_wrap(~ component)

## ----resids-vs-BC-bins, fig.cap = 'Mean of residuals by category of Basal Circumference.'----
BC = cut(ml$dat$BC, c(0, 30, 45, 60, 75, 90, 150), ordered = TRUE)

cbind(fvso.dat, BC) %>% 
  group_by(BC, Year, component) %>% 
  summarise(mresid = mean(Residuals, na.rm = TRUE),
            N = n()) %>% 
  ggplot(aes(BC, mresid, colour = Year)) +
  geom_point(aes(size = N)) +
  facet_wrap(~ component) +
  ylab("Mean residuals")

## ----resids-vs-PROV, fig.cap = 'Residuals vs Provenance.'----------------
ggplot(cbind(fvso.dat, BC = ml$dat$PROV),
       aes(BC, Residuals, colour = Year)) +
  geom_violin() +
  facet_wrap(~ component)

## ----resids-vs-BF98, fig.cap = 'Residuals vs Bud flushing.'--------------
ggplot(cbind(fvso.dat, BC = ml$dat$BF98), aes(BC, Residuals, colour = Year)) +
  geom_violin() +
  facet_wrap(~ component)

## ----mst-CL-genetic, fig.width = 7, fig.height = 6, fig.cap = 'Posterior individual additive-genetic effects. Binary component.'----
hpds.bin <- sapply(ml$mst$res$marginals.random$recode.bin,
                   function(x) inla.hpdmarginal(.95, x))

ml$mst$blups.bin <- data.frame(
  Seq = dat$Seq,
  Family = dat$FAM,
  mean = ml$mst$res$summary.random$recode.bin$mean[gen_str$recoding],
  ymin = hpds.bin[1, gen_str$recoding],
  ymax = hpds.bin[2, gen_str$recoding])

ggplot(ml$mst$blups, aes(x = mean, y = mean, ymin = ymin, ymax = ymax)) + 
  geom_pointrange() + 
  facet_wrap(~ Family) +
  coord_flip() +
  labs(x=NULL, y=NULL)

## ----mst-CL-hyperpar, fig.width  = 7.205, fig.height = 2.5, fig.cap = 'Prior and posterior densities for the spatio-temporal field characteristics.'----

stopifnot(exists('spde.bin', ml$mst))
res.spde.bin.marginals <- 
  ldply(tail(ml$mst$spde.bin, 2),
        function(x) data.frame(x[[1]]))
res.spde.bin.marginals <- 
  transform(res.spde.bin.marginals,
            .id=factor(res.spde.bin.marginals$.id,
                       labels = c('Range', 'Variance')))
res.spde.bin.marginals <- 
  rbind(
    cbind(type = 'Posterior',
          res.spde.bin.marginals),
    cbind(type = 'Prior', .id = 'Range',
          as.data.frame(sp_str$priors$rho)),
    cbind(type = 'Prior', .id = 'Variance',
          as.data.frame(sp_str$priors$sigma))
  )

# qplot(x, y, data=res.spde.bin.marginals, geom='line', col = type) +
#   facet_wrap(~.id, scales='free')

stopifnot(exists('spde.cont', ml$mst))
res.spde.cont.marginals <- 
  ldply(tail(ml$mst$spde.cont, 2),
        function(x) data.frame(x[[1]]))
res.spde.cont.marginals <- 
  transform(res.spde.cont.marginals,
            .id=factor(res.spde.cont.marginals$.id,
                       labels = c('Range', 'Variance')))
res.spde.cont.marginals <- 
  rbind(
    cbind(type = 'Posterior',
          res.spde.cont.marginals),
    cbind(type = 'Prior', .id = 'Range',
          as.data.frame(sp_str$priors$rho)),
    cbind(type = 'Prior', .id = 'Variance',
          as.data.frame(sp_str$priors$sigma)[sp_str$priors$sigma$x < 1, ]))
# qplot(x, y, data=res.spde.cont.marginals, geom='line', col = type) +
#   facet_wrap(~.id, scales='free')

hyper.dat <- 
  rbind(
    cbind(
      component = 'bin',
      res.spde.bin.marginals),
    cbind(component = 'cont',
          res.spde.cont.marginals)
  )

grid.arrange(
  ggplot(filter(hyper.dat, .id == 'Range'), aes(x, y, lty = type)) + 
    geom_line(show.legend = FALSE) + 
    facet_wrap(component ~ .id, scales = 'free_x', ncol = 1) +
    labs(x=NULL, y=NULL) +
    xlim(0, 150) +
    theme(legend.title = element_blank()),
  ggplot(filter(hyper.dat, .id == 'Variance'), aes(x, y, lty = type)) + 
    geom_line() + 
    facet_wrap(component ~ .id, scales = 'free_x', ncol = 1) +
    labs(x=NULL, y="") +
    # xlim(0, 150) +
    theme(legend.title = element_blank()),
  nrow = 1,
  widths = 3:4
)


## ----mst-CL-spatial_bin--------------------------------------------------
proj.vec = inla.mesh.projector(sp_str$mesh, loc = sp_str$loc)
sp.f <- function(x) inla.mesh.project(proj.vec, field = x)

ml$mst$post_spat.bin <- ldply(1:2, function(x) 
  data.frame(X = dat$X,
             Y = dat$Y,
             Year = levels(ml$dat$Year)[x],
             PostMean= sp.f(ml$mst$spde.bin$summary.values$mean[st.index$theta.group == x])))


sp.bin <- ggplot(ml$mst$post_spat.bin, aes(X, Y)) +
  geom_tile(aes(fill = PostMean, color = PostMean)) +
  coord_fixed() +
  labs(x='', y='Field row') +
  facet_wrap(~ Year) +
  scale_fill_gradient2(low='#832424FF', high='#3A3A98FF', name = 'Binary\npost.\nmean', space = "Lab") +
  scale_color_gradient2(low='#832424FF', high='#3A3A98FF', name = 'Binary\npost.\nmean', space = "Lab") 


## ----mst-CL-spatial_bet--------------------------------------------------
proj.vec = inla.mesh.projector(sp_str$mesh, loc = sp_str$loc)
sp.f <- function(x) inla.mesh.project(proj.vec, field = x)

ml$mst$post_spat.cont <- ldply(1:2, function(x) 
  data.frame(X = dat$X,
             Y = dat$Y,
             Year = levels(ml$dat$Year)[x],
             PostMean= sp.f(ml$mst$spde.cont$summary.values$mean[st.index$theta.group == x])))


sp.cont <- ggplot(ml$mst$post_spat.cont, aes(X, Y)) +
  geom_tile(aes(fill = PostMean, color = PostMean)) +
  coord_fixed() +
  labs(x='Field column', y='Field row') +
  facet_wrap(~ Year) +
  scale_fill_gradient2(low='#3A3A98FF', high='#832424FF', name = 'Gamma\npost.\nmean', space = "Lab") +
  scale_color_gradient2(low='#3A3A98FF', high='#832424FF', name = 'Gamma\npost.\nmean', space = "Lab") 


## ----mst-CL-spatial_joint, fig.width=10, fig.height=14, fig.cap = 'Posterior mean of the spatio-temporal field. Latent scale.'----
grid.arrange(sp.bin, sp.cont)

## ----mst-CL-spatial_joint-journal, include=FALSE-------------------------
## shades of gray version for the journal
comp_plot <- 
  grid.arrange(
    ggplot(ml$mst$post_spat.bin, aes(X, Y)) +
      geom_tile(aes(fill = PostMean, color = PostMean)) +
      coord_fixed() +
      labs(x='', y='Field row') +
      facet_wrap(~ Year) +
      scale_fill_gradient(
        low = 'black', high = 'white',
        name = 'Binary\npost.\nmean', space = "Lab") +
      scale_color_gradient(
        low = 'black', high = 'white',
        name = 'Binary\npost.\nmean', space = "Lab")
    ,
    ggplot(ml$mst$post_spat.cont, aes(X, Y)) +
      geom_tile(aes(fill = PostMean, color = PostMean)) +
      coord_fixed() +
      labs(x='Field column', y='Field row') +
      facet_wrap(~ Year) +
      scale_fill_gradient(
        low = 'white', high = 'black',
        name = 'Gamma\npost.\nmean', space = "Lab") +
      scale_color_gradient(
        low = 'white', high = 'black',
        name = 'Gamma\npost.\nmean', space = "Lab")
  )

## ----mst-CL-spatial-oscale-----------------------------------------------

## A list of 942 * 2 * 2 = 3768 marginals
## corresponding for the SPDE vertices of
## bin-2012; bin-2013; cont-2012; cont-2013
## in that specific order
spmar <- c(ml$mst$spde.bin$marginals.values,
           ml$mst$spde.cont$marginals.values)

## mean values of linear predictor for each Year
latent_means <- ml$mst$res$summary.fixed[, 'mean']

## ST marginals centered on the Year's mean
spmar.shifted <- mapply(
  shift_marginal, 
  spmar, 
  latent_means[c(st.index$theta.group,
                 2 + st.index$theta.group)],
  SIMPLIFY = FALSE)

## ST means on the original scale for each Year and component
spmeans.oscale <- c(
  sapply(
    spmar.shifted[seq_len(2*sp_str$mesh$n)],             # bin marginals
    function(x) inla.emarginal(inv.logit, x)
  ),
  sapply(
    spmar.shifted[2*sp_str$mesh$n + seq_len(2*sp_str$mesh$n)],  # cont marginals
    function(x) inla.emarginal(exp, x)
  )
)

ml$mst$post_spat_oscale <- ldply(
  1:2, function(x) 
    data.frame(
      X = dat$X,
      Y = dat$Y,
      Year = (2012:2013)[x],
      PostMean_bin= sp.f(spmeans.oscale[seq_len(2*sp_str$mesh$n)][st.index$theta.group == x]),
      PostMean_cont= sp.f(spmeans.oscale[2*sp_str$mesh$n+seq_len(2*sp_str$mesh$n)][st.index$theta.group == x])
    )
)


sp.bin <- ggplot(ml$mst$post_spat_oscale, aes(X, Y)) +
  geom_tile(aes(fill = PostMean_bin, color = PostMean_bin)) +
  coord_fixed() +
  labs(x='', y='Field row') +
  facet_wrap(~ Year) +
  scale_fill_viridis(name = 'Binary\npost.\nmean') +
  scale_color_viridis(name = 'Binary\npost.\nmean') 


sp.cont <- ggplot(ml$mst$post_spat_oscale, aes(X, Y)) +
  geom_tile(aes(fill = PostMean_cont, color = PostMean_cont)) +
  coord_fixed() +
  labs(x='Field column', y='Field row') +
  facet_wrap(~ Year) +
  scale_fill_viridis(name = 'Gamma\npost.\nmean') +
  scale_color_viridis(name = 'Gamma\npost.\nmean') 


## ----mst-CL-spatial_joint-oscale, fig.width=10, fig.height=14, fig.cap = 'Posterior mean of the spatio-temporal field. Original scale.'----
grid.arrange(sp.bin, sp.cont)

## ----distribution-PBVs-CLbin, fig.cap='Distribution of PBV for the Binomial component of CL, by observed infection status.'----

plotdat <- ml$dat %>%
  dplyr::select(Seq, Year, CL) %>%
  spread(Year, CL) %>%
  inner_join(ml$mst$blups.bin, by = 'Seq') %>%
  mutate(Infection = factor((`2012` > 0) + 2*(`2013` > 0),
                            labels = c('Not infec.',
                                       'Only 2012',
                                       'Only 2013',
                                       'Both')))

ggplot(plotdat, aes(x = mean, fill = Infection)) +
  geom_histogram(binwidth = diff(range(plotdat$mean))/30) +
  scale_fill_grey(na.value = 'white') +
  xlab('Predicted Breeding Value')



## ----distribution-PBVs-CLbin-journal, include=FALSE----------------------
## Simplified version for journal
ggplot(plotdat[!is.na(plotdat$Infection) & plotdat$Infection != 'Only 2012', ], aes(x = mean)) +
  geom_histogram(binwidth = diff(range(plotdat$mean))/30) +
  facet_grid(Infection~.) +
  xlab('Predicted Breeding Value')


## ----mst-CL-heritability-logit-scale, fig.cap = 'Posterior heritability by model component, including (a) or not (b) the spatio-temporal variance in the denominator.'----

N.samples <- 5e3

hyper.dat <- rbind(
  cbind(component = 'bin',
        res.spde.bin.marginals),
  cbind(component = 'cont',
        res.spde.cont.marginals)
)


sbin.mar <- filter(hyper.dat,
                   component == 'bin',
                   type == 'Posterior',
                   .id == 'Variance')
sbet.mar <- filter(hyper.dat,
                   component == 'cont',
                   type == 'Posterior',
                   .id == 'Variance')
hyperpar.sample <- inla.hyperpar.sample(N.samples, ml$mst$res)

var.sample <- data.frame(
  sbin = inla.rmarginal(N.samples, sbin.mar[, c('x', 'y')]),
  sbet = inla.rmarginal(N.samples, sbet.mar[, c('x', 'y')]),
  gbin  = 1/hyperpar.sample[, grep('recode.bin',
                                   colnames(hyperpar.sample))],
  gbet  = 1/hyperpar.sample[, grep('recode.cont',
                                   colnames(hyperpar.sample))],
  omega = 1/(1 + hyperpar.sample[, grep('Gamma observations',
                                        colnames(hyperpar.sample))])
)


## Heritability in the observed scale
## Approach following Villemereuil et al. 2015

# Samples from fitted values (i.e. exp(latent))

## Numerical problems with some very spliked marginals
## for the binary component. E.g.:
# tm <- ml$mst$res$marginals.fitted.values[inla.stack.index(ml$mst$stk,
#                                                       'bin')$data][[31]]
# plot(tm, type = 'l')
# inla.qmarginal(.1, tm)
# str(inla.smarginal(tm))
# plot(inla.smarginal(tm), type = 'l')

# The following discussion can be useful to overcome these numerical problems:
# https://groups.google.com/forum/embed/?parenturl=http%3A%2F%2Fwww.r-inla.org%2Fcomments-1&service=jotspot&ul=1&theme=default&place=forum%2Fr-inla-discussion-group&showsearch=true#!topic/r-inla-discussion-group/9B8CazRaik4)
        
## Other marginals give samples out of scale. E.g.: #16
# tm <- ml$mst$res$marginals.fitted.values[inla.stack.index(ml$mst$stk,
#                                                       'bin')$data][[31]]
# stopifnot(inla.rmarginal(N.samples, tm) > 0)

## Return a sample of N values from the posterior
## marginals of the fitted.values of the given component
sample_fit <- function(comp, mar, stk, N) {
  sapply(mar[inla.stack.index(stk, comp)$data],
         function(x) {
           tryCatch(inla.rmarginal(N, x),
                    error = function(e) 
                      rep(NA, N))
           }
  )
}

## Sample the fitted values for each component
fit.sample <- lapply(
  c(bin = 'bin', cont = 'cont'),
  sample_fit,
  mar = ml$mst$res$marginals.fitted.values,
  stk = ml$mst$stk,
  N   = N.samples
)

# Due to numerical problems in some marginals, I replace the 
# random sampling by the mode of the posterior
idx <- lapply(fit.sample, apply, 2, function(x) all(is.na(x)))
stopifnot(sapply(idx, sum) < sapply(idx, length) / .05) # < 5 %
fit.sample$bin[, idx$bin] <- 
  rep(
    ml$mst$res$summary.fitted.values$mean[
      inla.stack.index(ml$mst$stk, 'bin')$data][
        idx$bin],
    each = N.samples
  )

# Auxiliar diagnostic function
plot_sample <- function(x) {
  x %>% 
    as.data.frame() %>% 
    gather(component) %>% 
    ggplot(aes(value)) +
    geom_histogram() +
    facet_wrap(~component, scales = 'free')
}

# Average derivative of expected value wrt latent value (Eqs. 16-17)
Psi.sample <- lapply(fit.sample, apply, 1, mean)
# plot_sample(Psi.sample)

# Expected-scale phenotypic variance (Eq. 7)
V_P_exp.sample <- lapply(fit.sample, apply, 1, var)
# plot_sample(V_P_exp.sample)

# Distribution-specific variance:
# (MC integrated wrt latent)
# Gamma distribution: mu^2/phi
# Bernoulli distribution: p(1-p)
phi.sample <- hyperpar.sample[, grep('Gamma observations',
                                     colnames(hyperpar.sample))]
V_P_spec.sample <- list(
  bin = apply(fit.sample$bin*(1-fit.sample$bin), 1, mean),
  cont = apply(fit.sample$cont^2/phi.sample, 1, mean)
)
# plot_sample(V_P_spec.sample)

# Observed-scale phenotypic variance
V_P_obs <- mapply(`+`, V_P_exp.sample, V_P_spec.sample,
                  SIMPLIFY = FALSE)
# plot_sample(V_P_obs)

# Additive-genetic variance (latent scale)
V_A <- list(
  bin = 1/hyperpar.sample[, grep('recode.bin', colnames(hyperpar.sample))],
  cont = 1/hyperpar.sample[, grep('recode.cont', colnames(hyperpar.sample))]
)
# plot_sample(V_A)

# Additive-genetic variance (observed/expected phenotypic scales)
V_A_obs <- mapply(function(x, y) x^2*y,
                  Psi.sample,V_A,
                  SIMPLIFY = FALSE)
# plot_sample(V_A_obs)

## Method a: dividing by ST variance
## Method b: excluding ST variance
ml$mst <- c(ml$mst[-grep('h2', names(ml$mst))],
            transmute(var.sample,
                      h2.sample.bin.b = gbin/(gbin + pi^2/3),
                      # h2.sample.cont.b = gbet/(gbet + omega*pi^2/3),
                      # This formula applied for beta, but not for gamma
                      h2.sample.bin.a = gbin/(gbin + sbin + pi^2/3),
                      # h2.sample.bin.a = V_A_obs$bin/V_P_obs$bin,
                      # This give very high values. Numerical issues?
                      h2.sample.cont.a = V_A_obs$cont/V_P_obs$cont))


plotdat <- as.data.frame(ml$mst[grep('h2', names(ml$mst))]) %>% 
  # dplyr::rename(bin = h2.sample.bin, cont = h2.sample.cont) %>% 
  gather(component, x) %>% 
  mutate(component = factor(gsub('h2.sample.', '', component))) %>% 
  separate(component, c('component', 'spatial_variance'), convert = TRUE)


ggplot(plotdat, aes(x, colour = component)) +
  geom_density() +
  facet_grid(`spatial_variance`~.)


## ----mst-CL-heritability-logit-scale-journal,  opts.label = 'journal', include=TRUE----
(plotdat %>% 
  group_by(component, spatial_variance) %>% 
  dplyr::summarise(mean = mean(x),
                   ymin = quantile(x, probs = c(.025)),
                   ymax = quantile(x, probs = c(.975))) -> plotdatsum)
            
ggplot(plotdatsum, aes(component, mean, ymin = ymin, ymax = ymax)) +
  geom_pointrange() +
  facet_grid(`spatial_variance`~.) +
  coord_flip() +
  ylim(c(0,1)) +
  ylab('Posterior heritability')

## ----bin-summary-table---------------------------------------------------

## Gathering of previous results
bin_mcdat <- rbind(
  ## Null model and 5 models for BC
  binBC_mcdat[-(c(1,4)), -6],
  ## 2 models for PROV
  binPROV_mcdat[3:4, -6]
)

## Model numbering
rownames(bin_mcdat) <- NULL
bin_mcdat$Model <- paste0('M', 
                          c(0, paste(1, 1:5, sep = '.'),
                            paste(2, 1:2, sep = '.')))

## Model descriptions
bin_mcdat$Specificity <-   
  c('Reference',
    'BC categories',
    'BC cat. x Year',
    'BC linear',
    'BC quadratic',
    'BC non-parametric',
    'PROV',
    'PROV x Year')

## Model components
bin_mcdat$`Linear predictor` <- 
  c("Year + _IBV_ + _ST_",
    "Year + BC + _IBV_ + _ST_",
    "Year * BC + _IBV_ + _ST_",
    "Year + $\\beta$ BC + _IBV_ + _ST_",
    "Year + $\\beta$ BC + $\\beta$ BC^2 + _IBV_ + _ST_",
    "Year + f(BC) + _IBV_ + _ST_",
    "Year + PROV + _IBV_ + _ST_",
    "Year * PROV + _IBV_ + _ST_"
  )

## -2 log-likelihood
bin_mcdat <- transform(bin_mcdat,
                       'Dev.' = -2*marginal_loglik)


## Reorder and select columns
bin_mcdat <- bin_mcdat[, c(6:8, 1:4, 9)]

kable(bin_mcdat, digits = 0,
      caption="Model selection criteria for the binary component")


## ----cont-summary-table, fig.height = 3, fig.width = 4, fig.cap = 'Model selection criteria.'----

## Gathering of previous results
cont_mcdat <- rbind(
  cont_comp[, -6],
  ## Null model and 5 models for BC
  contBC_mcdat[-(c(1,4)), -6],
  ## 2 models for PROV
  contPROV_mcdat[3:4, -6]
)

## Model numbering
rownames(cont_mcdat) <- NULL
cont_mcdat$Model <- paste0('M', 
                           c(paste(0, 1:3, sep = '.'), 
                             paste(1, 1:5, sep = '.'),
                             paste(2, 1:2, sep = '.')))

## Model descriptions
cont_mcdat$Specificity <-   
  c('Beta likelihood',
    'Gamma likelihood',
    'Reference',
    'BC categories',
    'BC cat. x Year',
    'BC linear',
    'BC quadratic',
    'BC non-parametric',
    'PROV',
    'PROV x Year')

## Model components
cont_mcdat$`Linear predictor` <- 
  c("BF98 + PROV + Year*BC",
    "BF98 + PROV + Year*BC",
    "Year + _IBV_ + _ST_",
    "Year + BC + _IBV_ + _ST_",
    "Year * BC + _IBV_ + _ST_",
    "Year + $\\beta$ BC + _IBV_ + _ST_",
    "Year + $\\beta$ BC + $\\beta$ BC^2 + _IBV_ + _ST_",
    "Year + f(BC) + _IBV_ + _ST_",
    "Year + PROV + _IBV_ + _ST_",
    "Year * PROV + _IBV_ + _ST_"
  )

## -2 log-likelihood
cont_mcdat <- transform(cont_mcdat,
                       'Dev.' = -2*marginal_loglik)


## Reorder and select columns
cont_mcdat <- cont_mcdat[, c(6:8, 1:4, 9)]

kable(cont_mcdat, digits = 0,
      caption="Model selection criteria for the continuous component")

