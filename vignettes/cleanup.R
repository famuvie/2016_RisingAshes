### Cleanup

## This script leaves the data in four diferent formats
## Devecey : Original dataset. Should not be used!!!
## dat     : Original (wide) format with transformed variables CDT and CDL
## mdat    : long-format with a meta-factor 'variable' with values CD, CL, ...
## tdat    : tidy-format (one column per variable)
##
## The problem with the 'tidy' format is that I have NAs for all variables
## except BF for many years where only BF was measured.
## While mdat is easier to select variables while getting rid of all those
## year-long NAs at the same time.
## Recommended practice is to select variables and then 'spread' the variables
## using spread(droplevels(mdat[sel,]), variable, value)

library(tidyr)
library(dplyr)



### Data transforms

## Sqrt Transform
tr.f <- function(x) 1/sqrt(x + 0.1)
dat <- transform(
  Devecey,
  CDT10 = tr.f(CD10),
  CDT11 = tr.f(CD11),
  CDT12 = tr.f(CD12),
  CDT13 = tr.f(CD13),
  CDT14 = tr.f(CD14)
)

### Consider BF98 among the identifying variables
### rather than as a measure of a phenotype at some point in time.
### This is in order to use it as a covariate.
dat <- transform(dat, BF98 = factor(BF98))

### Variables spanning multiple columns
id.vars <- c('Seq', 'X', 'Y', 'BLC', 'SBLC', 'FAM', 'PROV', 'BF98')
var.basenames <- c('CD', 'CL', 'BC', 'CDT')
var.cols  <- sapply(paste0('^', var.basenames, '\\d{2}$'),
                    grep, names(dat), value = TRUE)
stopifnot(all.equal(names(dat), unname(c(id.vars, unlist(var.cols)))))


### Rename variables to separate variable name from year
separate_name_year <- function(x) {
  c(strtrim(x, nchar(x)-2),
    substr(x, nchar(x)-1, nchar(x)))
}
convert_year <- function(x) {
  paste(ifelse(x>format(Sys.time(), '%y'), '19', '20'),
        x, sep='')
}
names(dat) <-
  c(id.vars,
    sapply(lapply(unlist(var.cols), separate_name_year),
           function(x) paste(x[1], convert_year(x[2]), sep = '_')))

## Gather variables into one meta-variable mixed with the Year
mdat <- tidyr::gather(dat, variable, value, contains("_"))

## Separate the meta-variable from the Year
mdat <-
  tidyr::separate(mdat, variable, c('variable', 'Year')) %>%
  transform(
    variable = as.factor(variable),
    Year     = ordered(Year)
  )

### Split the meta-variable into separate variables
## This is in "tidy"" form
## https://github.com/hadley/tidyr
tdat <- spread(mdat, variable, value)

