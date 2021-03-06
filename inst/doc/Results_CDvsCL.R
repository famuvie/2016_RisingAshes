## ----setup-data, echo=FALSE, message=FALSE, cache=FALSE------------------
### Setup
# library(sp)
# library(lme4)
# library(lattice)  # some lme4 plotting functions
# library(boot)     # bootstrapping
# library(tidyr)
library(plyr)       # data manipulation
library(dplyr)      # data manipulation
library(ggplot2)    # plotting
library(gridExtra)  # plotting
library(GGally)     # ggpairs
library(knitr)      # kable()
# library(INLA)
# library(nnet)
# library(car)
# library(MASS)
# library(rpf)      # ordinal.gamma() association maasure

opts_chunk$set(echo       = FALSE,
               message    = FALSE,
               warning    = FALSE,
               comment    = NA,
               fig.dpi    = 96,
               fig.width  = 6,
               fig.height = 5,
               fig.caption= FALSE,
               cache      = FALSE)

## ggplot2 options
theme_set(theme_bw())

## ----preprocessing-------------------------------------------------------
## Applies an square-root transform to CD (creating variables CDT)
## Buils object mdat in long format (with variables: Year; Variable; Value)
## Open with:
# file.show(file=system.file('vignettes/cleanup.R',
#                            package = 'RisingAshes'))
source('cleanup.R')

## ----setup-dataset-------------------------------------------------------

bidat <- 
  filter(
    mdat,
    variable %in% c('CD', 'CL')
  ) %>%
  droplevels %>%
  spread(variable, value)


## ----CL-barplots, fig.height = 3, fig.cap = 'Relationship between proportion of trees with CL and BC (left) or CD (right)'----

tdat <- 
  filter(
    mdat,
    variable %in% c('CD', 'CL', 'BC'),
    Year %in% c('2012', '2013')
  ) %>%
  spread(variable, value) %>% 
  transform(
    CD = as.ordered(CD),
    BC = cut(BC, c(0, 30, 45, 60, 75, 90, 150), ordered = TRUE),
    CLbin = factor(CL > 0, labels = c('0', '>0'))) %>%
  filter(!is.na(CL), !is.na(CD)) %>% 
  droplevels()

levels(tdat$BC)[c(1, 6)] <- c('<= 30', '> 90')
levels(tdat$BC) <- gsub("\\(", "]", levels(tdat$BC))

BCdat <- 
  tdat %>% 
  group_by(BC, Year) %>% 
  summarise(
    NCL = sum(CL > 0),
    TCL = n()
  ) %>% 
  mutate(
    PT = NCL/TCL,
    PT.se = sqrt(PT*(1-PT)/TCL)
  )

p.bc <- 
  ggplot(BCdat, aes(BC, PT, fill = Year)) +
  geom_bar(stat = 'identity', position = 'dodge', show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = PT - 1.96*PT.se, ymax = PT + 1.96*PT.se),
    width = 0.2,
    position = position_dodge(.9)) +
  scale_fill_grey(start = 0.4) +
  coord_cartesian(ylim=c(0, 1.04)) +
  xlab('Basal Circumference (cm)') + 
  ylab('% of trees with CL') +
  theme(axis.text = element_text(size = rel(0.6)),
        axis.title = element_text(size = rel(0.8)))

CDdat <- 
  tdat %>% 
  group_by(CD, Year) %>% 
  summarise(
    NCL = sum(CL > 0),
    TCL = n()
  ) %>% 
  mutate(
    PT = NCL/TCL,
    PT.se = sqrt(PT*(1-PT)/TCL))

p.cd <- 
  ggplot(CDdat, aes(CD, PT, fill = Year)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  geom_errorbar(aes(ymin = PT - 1.96*PT.se, ymax = PT + 1.96*PT.se),
                width = 0.2,
                position = position_dodge(.9)) +
  #   geom_pointrange(aes(ymin = PT - 2*PT.se, ymax = PT + 2*PT.se),
  #                 position = position_dodge(.9)) +
  scale_fill_grey(start = 0.4) +
  coord_cartesian(ylim=c(0, 1.04)) +
  xlab('Crown Dieback (CD) categories') + 
  ylab('') +
  theme(axis.text = element_text(size = rel(0.6)),
        axis.title = element_text(size = rel(0.8)))

grid.arrange(p.bc, p.cd, ncol = 2, widths = 3:4)

## ----tests-independence--------------------------------------------------

table(tdat[, c('CLbin', 'BC', 'Year')])

## Chi-sqared tests of independence by year
summary(table(tdat[tdat$Year == '2012', c('CLbin', 'BC')]))
summary(table(tdat[tdat$Year == '2013', c('CLbin', 'BC')]))

table(tdat[, c('CLbin', 'CD', 'Year')])

## Chi-sqared tests of independence by year
summary(table(tdat[tdat$Year == '2012', c('CLbin', 'CD')]))
summary(table(tdat[tdat$Year == '2013', c('CLbin', 'CD')]))


## ----spearman-rho--------------------------------------------------------
with(filter(BCdat, Year == '2012'), cor.test(PT, as.numeric(BC), method = 'spearman'))
with(filter(BCdat, Year == '2013'), cor.test(PT, as.numeric(BC), method = 'spearman'))

with(filter(CDdat, Year == '2012'), cor.test(PT, as.numeric(CD), method = 'spearman'))
with(filter(CDdat, Year == '2013'), cor.test(PT, as.numeric(CD), method = 'spearman'))


## ----data-description, fig.height = 3, fig.cap='Phenotypic correlation between symptoms.'----

ggplot(filter(bidat, Year %in% c('2012', '2013')), aes(CD, CL)) +
  geom_violin(aes(group = round_any(CD, .01)), fill = 'grey80', adjust = 2)


# ggplot(bidat %>% filter(Year %in% c('2012', '2013')), aes(CD, CL)) +
#   geom_abline(int = 0, sl = 1, col = 'darkgray') +
# #   geom_jitter(aes(shape = Year)) +
#   geom_jitter(aes(col = Year), size = 1) +
#   theme_bw() +
#   theme(text=element_text(size=10)) +
#   scale_colour_grey()
# #   +scale_shape(solid = FALSE)

cormet <- function(x, y) 
  with(bidat[bidat$Year %in% y, ],
       cor(CD, CL,
           use = 'pairwise',
           method = tolower(x)))

cortmet <- function(x, y)
  formatC(with(bidat[bidat$Year %in% y, ],
               cor.test(CD, CL,
                        use = 'pairwise',
                        alternative = 'greater',
                        method = tolower(x)))$p.value,
          format = 'e',
          digits = 1)

# cortmet('pearson')
# cortmet('spearman')
# cortmet('kendall')
# p-values < numerical precision: but warnings about ties!!

cortab <- transform(data.frame(Method = c('Pearson', 'Spearman', 'Kendall')),
                    Year = c('2012', rep('', 2), '2013', rep('', 2)),
                    'Sample correlation' = c(sapply(Method, cormet, '2012'),
                                             sapply(Method, cormet, '2013')),
                    'p-value' = c(sapply(Method, cortmet, '2012'),
                                  sapply(Method, cortmet, '2013')))


cortab.pooled <- transform(data.frame(Method = c('Pearson', 'Spearman', 'Kendall')),
                           'Sample correlation' = sapply(Method, cormet, c('2012', '2013')),
                           'p-value' = sapply(Method, cortmet, c('2012', '2013')))

kable(cortab, digits = 2)
kable(cortab.pooled, digits = 2)

## ----data-PBVs-----------------------------------------------------------
ml <- list(
  CD = RisingAshes:::stmodels_CD,
  CL = RisingAshes:::mst_CL)

## compute the pedigree once more
n.fam <- nlevels(dat$FAM)
ped <- with(
  dat[, c('Seq', 'FAM')],
  data.frame(
    self = c(Seq, max(Seq) + 1:n.fam),
    mum  = c(max(Seq) + as.numeric(FAM), rep(0, n.fam)),
    dad  = 0))
ped$mum[is.na(ped$mum)] <- 0
Ainv <- compute_Ainverse(ped)

BV.dat <- data.frame(
  Seq = Ainv$map[, 1],
  BV_TCD = ml$CD$wif$res$summary.random$recode$mean,
  BV_CL.bin = ml$CL$summary.random$recode.bin$mean,
  BV_CL.cont = ml$CL$summary.random$recode.cont$mean,
  type = c(rep('Founder', 23), rep('Progeny', 788)))

plotdat <- 
  bidat %>%
  filter(Year %in% c('2012', '2013')) %>%
  dplyr::select(Seq, Year, CL) %>%
  droplevels %>%
  spread(Year, CL) %>%
  inner_join(BV.dat, by = 'Seq') %>%
  mutate(
    Infection = factor(
      (`2012` > 0) + 2*(`2013` > 0),
      labels = c('Not infec.',
                 'Only 2012',
                 'Only 2013',
                 'Both')))


## ----genetic-correlation, fig.height = 3, fig.cap='Genetic correlation between PBV and their posterior distribution.'----

ggpairs(
  BV.dat %>% filter(type == 'Progeny'),
  columns = 2:4,
  columnLabels = c('TCD', 'CL Binomial', 'CL Gamma')
) +
  theme_bw(base_size = 10)


## ----genetic-correlation-tests-------------------------------------------

cordat <- BV.dat %>% 
  filter(type == 'Progeny') %>% 
  dplyr::select(contains('BV')) 

cor.test(~ BV_TCD + BV_CL.bin,
         data = cordat, alternative = 'greater')
cor.test(~ BV_TCD + BV_CL.cont,
         data = cordat, alternative = 'less')
cor.test(~ BV_CL.bin + BV_CL.cont,
         data = cordat, alternative = 'less')


