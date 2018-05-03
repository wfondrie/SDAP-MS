BTAP-MS
================
William E Fondrie
5/2/2018

# Introduction

## Setup

``` r
library(tidyverse)
library(broom)
library(ggridges)
library(knitr) # for the kable() function

# load auxillary functions
source("R/ggplotTheme.R") # Plot themes

# Set ggplot2 theme
theme_set(coolTheme)
```

## Load MaxQuant Protein Results

Protein quantitation data was imported from the MaxQuant
`proteinGroups.txt` files. Reverse and contaminant proteins were
filtered out. Additionally, we required proteins to be quantified by a
minimum of 2 peptides.

``` r
prot <- read_tsv("data/combined/txt/proteinGroups.txt") %>%
  filter(!str_detect(`Protein IDs`, "^(REV|CON)__"),
         Peptides >= 2) 
```

The LFQ intensity columns for each protein were then extracted and
reshaped from wide to long format for downstream modeling. The column
names indicate the bait and concentration in the following form
`[bait]_[3^x nM]`. Thus, `LFQ intensity LRP1B_1` indicates that the
GST-LRP1B-ICD was used as bait at a concentration of 3 nM. Samples with
the GST suffix indicate that only GST was used.

``` r
lfq <- prot %>%
  select(`Protein IDs`, starts_with("LFQ ")) %>%
  gather(samp, lfq, -`Protein IDs`) %>%
  mutate(samp = str_match(samp, " (LRP1.*_.*)$")[ , 2],
         bait = str_match(samp, "(LRP1.*)_")[ , 2], # extract bait used
         conc = str_match(samp, "_(.*)$")[ , 2], # Create concentration column
         conc = ifelse(conc == "GST", 0, 3^as.numeric(conc))) 

# This results in "NAs introduced by coercion" warning, but the problem is
# handled in the final `ifelse()` statement.
```

The number of proteins passing this stage of analysis is:

``` r
lfq %>% 
  group_by(bait) %>%
  summarize(proteins = length(unique(`Protein IDs`))) %>%
  kable()
```

| bait  | proteins |
| :---- | -------: |
| LRP1  |     3025 |
| LRP1B |     3025 |

## Handling Missing Values

Zeros, `0`, indicate missing values in the MaxQuant results. With the
BTAP-MS experimental design, we expect a number of both values that are
missing at random (MAR) and left censored (meaning below the limit of
detection). In an effort to distinguish between these, missing values
were left as `0` if all of the intensities at greater bait
concentrations were also `0`. All other missing values were assumed to
be MAR, and were not considered in modeling.

``` r
lfqCensored <- lfq %>%
  group_by(`Protein IDs`, bait) %>%
  filter(sum(lfq) > 0) %>%
  mutate(maxConc = max(conc[lfq > 0]),
         lfq = ifelse(lfq == 0 & conc < maxConc, NA, lfq),
         numNA = sum(is.na(lfq)),
         numZero = sum(lfq == 0, na.rm = T))
```

In an effort to ensure robust modeling, only proteins with a minimum of
5 valid data points (out of the 9 total), were considered. Additionally,
no more than 3 of the data points could be left censored, and in total 4
data points needed to be non-zero.

``` r
lfqFiltered <- lfqCensored %>%
  group_by(`Protein IDs`, bait) %>%
  filter(numNA <= 4,
         numZero <= 3,
         numNA + numZero <= 5)
```

The number of proteins passing this stage of analysis is:

``` r
lfqFiltered %>% 
  group_by(bait) %>%
  summarize(proteins = length(unique(`Protein IDs`))) %>%
  kable()
```

| bait  | proteins |
| :---- | -------: |
| LRP1  |     2194 |
| LRP1B |     2077 |

## Transforming LFQ Intensities

The measured LFQ intensities indicate the amount of unbound protein in
each sample. To convert this to be proportional to the amount of bound
protein, which is needed obtain familiar binding curve shapes, the
maximum LFQ intensity for each protein was subtracted and the signed was
inversed.

``` r
lfqModInput <- lfqFiltered %>%
  group_by(`Protein IDs`, bait) %>%
  mutate(response = -(lfq - max(lfq)))
```

## Nonlinear Modeling

Each protein was fit to the 1:1 equilibrium binding isotherm using
nonlinear least-squares regression. The binding isotherm model takes the
form of:

  
![ R = \\frac{ \[B\]\_t \* R\_{max} }{ \[B\]\_t + K\_d }
](https://latex.codecogs.com/png.latex?%20R%20%3D%20%5Cfrac%7B%20%5BB%5D_t%20%2A%20R_%7Bmax%7D%20%7D%7B%20%5BB%5D_t%20%2B%20K_d%20%7D%20
" R = \\frac{ [B]_t * R_{max} }{ [B]_t + K_d } ")  

Where given ![R](https://latex.codecogs.com/png.latex?R "R") (the
response) and ![\[B\]\_t](https://latex.codecogs.com/png.latex?%5BB%5D_t
"[B]_t") (the total bait concentration), we fit the curves to estimate
![R\_{max}](https://latex.codecogs.com/png.latex?R_%7Bmax%7D "R_{max}")
(the estimated maximal response) and
![K\_d](https://latex.codecogs.com/png.latex?K_d "K_d") (the equilibrium
dissociation constant).

``` r
mods <- lfqModInput %>%
  rename(Protein = `Protein IDs`) %>% #glance() doesn't seem to like the long one
  group_by(Protein, bait) %>%
  filter(!is.na(response)) %>%
  do(models = nls(response ~ (Rmax * conc) / (Kd + conc),
                  data = .,
                  start = list(Kd = 100, Rmax = 4e+07),
                  control = list(maxiter = 200,
                                 warnOnly = T)))

# Retrieve information about model fits
fitInfo <- mods %>% glance(models)

kable(fitInfo[1:5, ])
```

| Protein | bait  |   sigma | isConv |  finTol |     logLik |      AIC |      BIC |     deviance | df.residual |
| :------ | :---- | ------: | :----- | ------: | ---------: | -------: | -------: | -----------: | ----------: |
| A0AVT1  | LRP1  | 1931140 | TRUE   | 9.1e-06 | \-141.9021 | 289.8042 | 290.3959 | 2.610510e+13 |           7 |
| A0AVT1  | LRP1B | 2878347 | TRUE   | 8.3e-06 | \-145.4941 | 296.9881 | 297.5798 | 5.799417e+13 |           7 |
| A1L020  | LRP1  | 1875136 | TRUE   | 1.9e-06 | \-141.6373 | 289.2745 | 289.8662 | 2.461294e+13 |           7 |
| A1L020  | LRP1B | 5120751 | TRUE   | 6.4e-06 | \-150.6788 | 307.3577 | 307.9493 | 1.835546e+14 |           7 |
| A1X283  | LRP1B | 1359508 | TRUE   | 7.1e-06 | \-138.7432 | 283.4865 | 284.0781 | 1.293783e+13 |           7 |

``` r
# Retrieve model parameter values
fitVals <- mods %>% 
  tidy(models) %>%
  mutate(CV = std.error / estimate * 100)

kable(fitVals[1:5, ])
```

| Protein | bait  | term |       estimate |    std.error |   statistic |   p.value |          CV |
| :------ | :---- | :--- | -------------: | -----------: | ----------: | --------: | ----------: |
| A0AVT1  | LRP1  | Kd   |   1.761082e+00 | 3.492715e+00 |   0.5042159 | 0.6295894 |   198.32775 |
| A0AVT1  | LRP1  | Rmax |   3.084011e+06 | 8.976176e+05 |   3.4357739 | 0.0108995 |    29.10552 |
| A0AVT1  | LRP1B | Kd   |   4.217124e-01 | 1.071854e+00 |   0.3934419 | 0.7057016 |   254.16716 |
| A0AVT1  | LRP1B | Rmax |   5.383047e+06 | 1.226457e+06 |   4.3891033 | 0.0031989 |    22.78370 |
| A1L020  | LRP1  | Kd   | \-3.835195e-01 | 8.987731e-01 | \-0.4267145 | 0.6824040 | \-234.34873 |

Because many of the proteins we measured do not actually interact with
the bait proteins, there are a large number of model fits that failed to
converge. These were removed from consideration.

``` r
fits <- fitInfo %>%
  filter(isConv) %>%
  select(Protein, bait) %>%
  left_join(fitVals)
```

And lastly, we reshape the parameter values to a wide format.

``` r
RmaxTbl <- fits %>%
  filter(term == "Rmax") %>%
  select(Protein, bait, estimate, CV) %>%
  rename(Rmax = estimate, Rmax_CV = CV)

fitTbl <- fits %>%
  filter(term == "Kd") %>%
  left_join(RmaxTbl)

kable(fitTbl[1:5, ])
```

| Protein | bait  | term |      estimate |  std.error |   statistic |   p.value |          CV |    Rmax |  Rmax\_CV |
| :------ | :---- | :--- | ------------: | ---------: | ----------: | --------: | ----------: | ------: | --------: |
| A0AVT1  | LRP1  | Kd   |     1.7610824 |  3.4927152 |   0.5042159 | 0.6295894 |   198.32775 | 3084011 |  29.10552 |
| A0AVT1  | LRP1B | Kd   |     0.4217124 |  1.0718544 |   0.3934419 | 0.7057016 |   254.16716 | 5383047 |  22.78370 |
| A1L020  | LRP1  | Kd   |   \-0.3835195 |  0.8987731 | \-0.4267145 | 0.6824040 | \-234.34873 |  962690 |  75.59953 |
| A1L020  | LRP1B | Kd   | \-202.8166979 | 85.9678233 | \-2.3592164 | 0.0503994 |  \-42.38696 | 1316586 | 204.29184 |
| A1X283  | LRP1B | Kd   |   \-0.0621288 |  0.6093989 | \-0.1019509 | 0.9216546 | \-980.86440 | 2166870 |  25.42602 |

The number of proteins passing this stage of analysis is:

``` r
fitTbl %>% 
  group_by(bait) %>%
  summarize(proteins = length(unique(Protein))) %>%
  kable()
```

| bait  | proteins |
| :---- | -------: |
| LRP1  |     1430 |
| LRP1B |     1280 |

## Filtering For Sufficient Model Fits

While some estimated ![K\_d](https://latex.codecogs.com/png.latex?K_d
"K_d") values are physically impossible, such as those below zero,
others are outside of the range that this experiment was designed to
measure. Because the bait concentrations used were between 1 and 2187
nM, we filtered for ![K\_d](https://latex.codecogs.com/png.latex?K_d
"K_d") between 3 nM and 1000 nM. Additionally, higher coefficients of
variation (CV) are indicative of poor model fits so we use it as an
additional filter.

``` r
interactors <- fitTbl %>%
  filter(estimate > 3, 
         estimate < 1000,
         CV < 100) %>%
  arrange(CV)
```

The number of proteins passing this stage of analysis is:

``` r
interactors %>% 
  group_by(bait) %>%
  summarize(proteins = length(unique(Protein))) %>%
  kable()
```

| bait  | proteins |
| :---- | -------: |
| LRP1  |       28 |
| LRP1B |       81 |

# Results
