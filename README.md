BTAP-MS
================
William E Fondrie
5/2/2018

# Introduction

# Methods

## Setup

``` r
library(tidyverse)
library(broom)
library(ggridges)
library(knitr) # for the kable() function
library(seqinr) # calculate molecular weights

set.seed(190875) # to ensure reproducibility

# load auxillary functions
source("R/ggplotTheme.R") # Plot themes
full <- 7
half <- 3.33

# Set ggplot2 theme
theme_set(coolTheme)
```

## Load Skyline Results

Precursor quantitation data was imported from the Skyline (exported
“Transition Results” report) `skyline_results.csv`. Proteins were
summarized as the sum of all precursor intensities for the protein. The
`Replicate` column indicates the bait and concentration in the following
form `...[bait]-Rep1-[3^x nM]`. Thus, `...LRP1B-Rep1-1` indicates that
the GST-LRP1B-ICD was used as bait at a concentration of 3 nM. Samples
with the GST suffix indicate that only GST was used.

``` r
prot <- read_csv("data/skyline_results.csv") %>%
    group_by(Protein, Replicate) %>%
    summarize(intensity = sum(Area, na.rm = T)) %>%
    ungroup() %>%
    mutate(bait = paste0("GST", str_match(Replicate, "BTAP(-.+-)Rep")[ , 2], "ICD"),
           bait = str_replace(bait, "LRP-", "LRP1-"),
           conc = str_match(Replicate, "Rep1-(.*?)(-|_|$)")[ , 2],
           conc = ifelse(conc == "GST", 0, 3^as.numeric(conc)))

# Retrieve protein info
fasta <- read.fasta("data/SwissProt-Human_v2018-02_GST.fasta", seqtype = "AA", 
                    as.string =T, set.attributes = T)

fastaKey <- tibble(Protein = names(fasta), sequence = unlist(fasta), 
                   annotation = map_chr(fasta, ~ attributes(.)$Annot))

protInfo <- prot %>%
    group_by(Protein) %>%
    summarize(accession = str_match(Protein, "\\|(.*)\\|") [1 , 2]) %>%
    left_join(fastaKey) %>%
    group_by(Protein) %>%
    mutate(MW = pmw(unlist(str_split(sequence, ""))) / 1000, #kDa
           gn = str_match(annotation, "GN=(.*?) ")[ , 2]) %>%
    select(Protein, accession, gn, MW) # only the cols we'll use

# Median normalization
globalMedian <- log2(median(prot$intensity))

lfq <- prot %>% 
    group_by(Replicate) %>% 
    mutate(medianIntensity = log2(median(intensity)),
           lfq = 2^(log2(intensity) - medianIntensity + globalMedian))
```

The number of proteins passing this stage of analysis is:

``` r
prot %>% 
  group_by(bait) %>%
  summarize(proteins = length(unique(Protein))) %>%
  kable()
```

| bait          | proteins |
| :------------ | -------: |
| GST-LRP1-ICD  |     2558 |
| GST-LRP1B-ICD |     2558 |

## Transforming Intensities

The measured intensities indicate the amount of unbound protein in each
sample. To convert this to be proportional to the amount of bound
protein, which is needed obtain familiar binding curve shapes, the
maximum LFQ intensity for each protein was subtracted and the signed was
inversed.

``` r
lfqModInput <- lfq %>%
  group_by(Protein, bait) %>%
  mutate(response =  max(lfq) - lfq)
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

| Protein                | bait          |     sigma | isConv |    finTol |     logLik |      AIC |      BIC |     deviance | df.residual |
| :--------------------- | :------------ | --------: | :----- | --------: | ---------: | -------: | -------: | -----------: | ----------: |
| sp|A0AVT1|UBA6\_HUMAN  | GST-LRP1-ICD  |  78129789 | TRUE   | 0.0000081 | \-175.2045 | 356.4089 | 357.0006 | 4.272985e+16 |           7 |
| sp|A0AVT1|UBA6\_HUMAN  | GST-LRP1B-ICD | 450942051 | TRUE   | 0.0000065 | \-190.9812 | 387.9624 | 388.5540 | 1.423441e+18 |           7 |
| sp|A0FGR8|ESYT2\_HUMAN | GST-LRP1-ICD  | 689192407 | FALSE  | 0.7538879 | \-194.7988 | 395.5976 | 396.1893 | 3.324903e+18 |           7 |
| sp|A0FGR8|ESYT2\_HUMAN | GST-LRP1B-ICD | 420933846 | TRUE   | 0.0000057 | \-190.3614 | 386.7228 | 387.3145 | 1.240297e+18 |           7 |
| sp|A0MZ66|SHOT1\_HUMAN | GST-LRP1-ICD  |   2112758 | TRUE   | 0.0000048 | \-142.7111 | 291.4222 | 292.0138 | 3.124623e+13 |           7 |

``` r
# Retrieve model parameter values
fitVals <- mods %>% 
  tidy(models) %>%
  mutate(CV = std.error / estimate * 100)

kable(fitVals[1:5, ])
```

| Protein                | bait          | term |       estimate |    std.error |   statistic |   p.value |             CV |
| :--------------------- | :------------ | :--- | -------------: | -----------: | ----------: | --------: | -------------: |
| sp|A0AVT1|UBA6\_HUMAN  | GST-LRP1-ICD  | Kd   |   6.740890e-01 | 1.065241e+00 |   0.6328041 | 0.5469690 |   1.580268e+02 |
| sp|A0AVT1|UBA6\_HUMAN  | GST-LRP1-ICD  | Rmax |   1.925117e+08 | 3.399920e+07 |   5.6622425 | 0.0007647 |   1.766085e+01 |
| sp|A0AVT1|UBA6\_HUMAN  | GST-LRP1B-ICD | Kd   | \-1.881341e+02 | 1.827913e+02 | \-1.0292288 | 0.3376165 | \-9.716012e+01 |
| sp|A0AVT1|UBA6\_HUMAN  | GST-LRP1B-ICD | Rmax |   7.371655e+07 | 2.266648e+08 |   0.3252228 | 0.7545171 |   3.074816e+02 |
| sp|A0FGR8|ESYT2\_HUMAN | GST-LRP1-ICD  | Kd   | \-5.526149e+09 | 8.514646e+16 | \-0.0000001 | 1.0000000 | \-1.540792e+09 |

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

| Protein                | bait          | term |      estimate |   std.error |   statistic |   p.value |         CV |       Rmax |    Rmax\_CV |
| :--------------------- | :------------ | :--- | ------------: | ----------: | ----------: | --------: | ---------: | ---------: | ----------: |
| sp|A0AVT1|UBA6\_HUMAN  | GST-LRP1-ICD  | Kd   |     0.6740890 |   1.0652413 |   0.6328041 | 0.5469690 |  158.02680 |  192511729 |    17.66085 |
| sp|A0AVT1|UBA6\_HUMAN  | GST-LRP1B-ICD | Kd   | \-188.1340940 | 182.7913198 | \-1.0292288 | 0.3376165 | \-97.16012 |   73716553 |   307.48155 |
| sp|A0FGR8|ESYT2\_HUMAN | GST-LRP1B-ICD | Kd   | \-791.8501715 | 316.0120269 | \-2.5057596 | 0.0406483 | \-39.90806 | \-50084258 | \-497.90444 |
| sp|A0MZ66|SHOT1\_HUMAN | GST-LRP1-ICD  | Kd   |     0.5174501 |   0.8461471 |   0.6115369 | 0.5601723 |  163.52243 |    5580542 |    16.26683 |
| sp|A0MZ66|SHOT1\_HUMAN | GST-LRP1B-ICD | Kd   |     0.1016078 |   0.8007250 |   0.1268947 | 0.9025919 |  788.05479 |   56966820 |    25.59808 |

The number of proteins passing this stage of analysis is:

``` r
fitTbl %>% 
  group_by(bait) %>%
  summarize(proteins = length(unique(Protein))) %>%
  kable()
```

| bait          | proteins |
| :------------ | -------: |
| GST-LRP1-ICD  |     2242 |
| GST-LRP1B-ICD |     2122 |

## Protein Concentration Estimation

For the binding isotherm equation above to be valid, we must make the
assumption that the conentration of prey protein,
![\[P\]\_t](https://latex.codecogs.com/png.latex?%5BP%5D_t "[P]_t") is
much less than the ![K\_d](https://latex.codecogs.com/png.latex?K_d
"K_d"). In an effort to verify this assumption, the prey protein
concentrations were crudely estimated using the “Total Protein
Approach.” This approach uses the following estimation to calculate
the relative protein concentration in a shotgun proteomics study:

  
![ \\frac{Protein~Mass}{Total~Protein~Mass} \\approx
\\frac{Protein~MS~Signal}{Total~MS~Signal}](https://latex.codecogs.com/png.latex?%20%5Cfrac%7BProtein~Mass%7D%7BTotal~Protein~Mass%7D%20%5Capprox%20%5Cfrac%7BProtein~MS~Signal%7D%7BTotal~MS~Signal%7D
" \\frac{Protein~Mass}{Total~Protein~Mass} \\approx \\frac{Protein~MS~Signal}{Total~MS~Signal}")  

Thus, given the mass spec signal of a protein, the total mass spec
signal (the sum of intensities for a run) we can estimate the relative
contributation of a single protein to the total protein mass analyzed.
Because the experiment was performed at a total protein conentration of
1 ug/uL, we can then calculate the individual protein concentrations
using their molecular weights.

``` r
tpaQuan <- lfqModInput %>%
  group_by(bait, conc) %>%
  left_join(protInfo) %>%
  mutate(tpa = lfq / sum(lfq, na.rm = T) / (MW * 1000),
         preyConc = tpa * 10^9,
         logPreyConc = log10(preyConc))

titer <- tpaQuan %>%
  group_by(Protein) %>%
  summarize(maxConc = max(preyConc, na.rm = T))
```

## Filtering For Sufficient Model Fits

While some estimated ![K\_d](https://latex.codecogs.com/png.latex?K_d
"K_d") values are physically impossible, such as those below zero,
others are outside of the range that this experiment was designed to
measure. Because the bait concentrations used were between 1 and 2187
nM, we filtered for ![K\_d](https://latex.codecogs.com/png.latex?K_d
"K_d") between 1 nM and 1000 nM. Additionally, higher coefficients of
variation (CV) are indicative of poor model fits so we use it as an
additional filter.

``` r
interactors <- fitTbl %>%
  left_join(titer) %>%
  filter(estimate > 1, 
         estimate < 1000,
         estimate > maxConc * 10) %>%
  arrange(CV) 
```

The number of proteins passing this stage of analysis is:

``` r
interactors %>% 
  group_by(bait) %>%
  summarize(proteins = length(unique(Protein))) %>%
  kable()
```

| bait          | proteins |
| :------------ | -------: |
| GST-LRP1-ICD  |      148 |
| GST-LRP1B-ICD |      185 |

## Final interactor list

``` r
highConfCut <- quantile(interactors$CV, 0.05)

highConf <- interactors %>%
    ungroup() %>%
    filter(CV <= highConfCut) %>%
    arrange(CV) 

highConf %>% 
  group_by(bait) %>%
  summarize(proteins = length(unique(Protein))) %>%
  kable()
```

| bait          | proteins |
| :------------ | -------: |
| GST-LRP1-ICD  |        9 |
| GST-LRP1B-ICD |        8 |

``` r
a <- filter(highConf, bait == "GST-LRP-ICD") %>%
    mutate(gn = str_match(Protein, "\\|.*\\|(.*)?_")[ , 2])
```

# Results

## TPA Estimation Plots

``` r
boxPlot <- tpaQuan %>%
  ggplot(aes(x = as.factor(conc), y = logPreyConc)) +
  geom_violin(fill = ggColors(3)[3]) +
  geom_boxplot(fill = "grey", width = 0.2, outlier.shape = NA) +
  facet_wrap(~ bait, ncol = 1) +
  ylab(expression("log"[10]*"[Protein (nM)]")) +
  xlab("Bait (nM)")

savePlot(boxPlot, w = full/2, h = 4)

overallDist <- titer %>%
  ggplot(aes(x = log10(maxConc))) +
  geom_density(fill = ggColors(3)[3]) +
  xlab(expression("log"[10]*"[Protein (nM)]")) +
  ylab("Density")
savePlot(overallDist, w = full/2, h = 2)


# some concentration stats
titer %>% 
    summarize(Mean = mean(log10(maxConc)),
              `Mean (nM)` = 10^Mean,
              `Std Dev` = sd(log10(maxConc)),
              total = length(Protein),
              `< 1` = sum(maxConc < 1) / total * 100,
              `< 10` = sum(maxConc < 10) / total * 100,
              `< 100` = sum(maxConc < 100)/ total * 100) %>%
    kable
```

|      Mean | Mean (nM) |   Std Dev | total |     \< 1 |    \< 10 |   \< 100 |
| --------: | --------: | --------: | ----: | -------: | -------: | -------: |
| 0.5807987 |  3.808892 | 0.6964215 |  2558 | 19.82017 | 74.47224 | 97.53714 |

``` r
# park stats
tpaQuan %>%
    left_join(protInfo) %>%
    filter(gn == "PARK7") %>%
    ungroup() %>%
    summarize(`mean (nM)` = mean(preyConc),
              `se` = sd(preyConc) / sqrt(length(preyConc))) %>%
    kable
```

| mean (nM) |       se |
| --------: | -------: |
|   34.1956 | 1.226328 |

``` r
barPlot <- titer %>%
  summarize(total = length(Protein),
            `< 1` = sum(maxConc < 1) / total * 100,
            `< 10` = sum(maxConc < 10) / total * 100,
            `< 100` = sum(maxConc < 100)/ total * 100) %>%
  gather(conc, percent, starts_with("<")) %>%
  ggplot(aes(x = conc, y = percent, fill = conc)) +
  geom_col(position = "dodge", color = "black") +
  #geom_errorbar(aes(ymax = avg + ci, ymin = avg - ci), width = 0.25) +
  xlab("Protein (nM)") +
  ylab("Number of Proteins (%)") +
  theme(legend.position = "none")

savePlot(barPlot, w = full/2, h = 2)
```

## Evaluating Model Fit Plots

``` r
fitScatter <- interactors %>% 
  ggplot(aes(x = log10(estimate), y = CV, color = bait)) +
  geom_point() +
  geom_hline(yintercept = highConfCut, linetype = "dashed") + 
  scale_color_discrete(name = "Bait") +
  theme(legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_rect(color = "black"),
        legend.title.align = 0.5,
        legend.key.height = unit(0.5, "lines")) +
  xlab(expression("Log"[10]*"[K"[d]~" (nM)]")) +
  ylab("Coefficient of Variation (%)")

savePlot(fitScatter, w = half, h = half)

cvDens <- interactors %>%
  ggplot(aes(x = CV, fill = bait)) +
  geom_histogram() +
  geom_vline(xintercept = highConfCut, linetype = "dashed") +
  facet_wrap(~bait, ncol = 1) +
  theme(legend.position = "none") +
  ylab("Number of Proteins") +
  xlab("Coefficient of Variation (%)")
  
savePlot(cvDens, w = half, h = half)
```

## Binding Curve Plots

``` r
plotIsotherm <- function(dat, responseDat = lfqModInput, logScale = F, ncol = 3) {
  
    # Function to create binding curves
    curveFun <- function(Protein, Rmax, estimate, CV, gn, bait,..., conc) {
        response <- (Rmax * conc) / (estimate + conc)
        tibble(response = response,
               conc = conc,
               Protein = Protein,
               CV = CV,
               gn = gn,
               estimate = estimate,
               bait = bait)
    }
    
    
    if(logScale) {
      x <- 10^(seq(0,log10(2187), length.out = 100))
    } else {
      x <- seq(0, 2187, length.out = 100)
    }
    
    # create curve data
    curveDat <- dat %>%
        select(Protein, estimate, Rmax, CV, gn, bait)%>%
        pmap_df(curveFun, conc = x) %>%
        ungroup() %>%
        mutate(gn = fct_reorder2(gn, as.numeric(bait), estimate))
    
    # Add response data
    datR <- dat %>% 
        select(Protein, bait, estimate, Rmax, CV, gn) %>%
        left_join(responseDat) %>%
        ungroup() %>%
        mutate(gn = fct_reorder2(gn, as.numeric(bait), estimate))

    # make plot
    p <- datR %>%
        ggplot(aes(x = conc, y = response * 1e-7)) +
        geom_line(aes(color = bait), data = curveDat) +#, color = ggColors(3)[3]) +
        geom_point(size = 0.5) +
        facet_wrap(~gn, ncol = ncol, scales = "free_y") +
        ylab(expression("Response (10"^7*")")) +
        scale_color_discrete(name = "Bait")
      
    
    if(logScale){ 
        p <- p + scale_x_log10(breaks = c(1, 10, 100, 1000)) 
    } else {
        p <- p + scale_x_continuous(breaks = c(0, 1000, 2000))
    } 
    
    if(length(unique(dat$bait)) == 2) {
        p <- p + xlab("Bait (nM)") 
        
    } else {
        p <- p + xlab(paste0(dat$bait[1], " (nM)"))
    }
    
    p
}


allCurves <- highConf %>%
    left_join(protInfo) %>%
    plotIsotherm(logScale = T, ncol = 3) +
    theme(legend.position = c(0.75, 0),
          legend.justification = c(0, 0),
          legend.key.size = unit(0.5, "lines"),
          legend.background = element_rect(color = "black")) 
savePlot(allCurves, w = 2 * full / 3, h = 5)
```

## Plot Calculated Affinities and Error Bars

``` r
exBar <- highConf %>%
    left_join(protInfo) %>%
    ungroup() %>%
    mutate(gn = fct_reorder2(gn, 
                             fct_rev(as.factor(bait)), 
                             estimate, 
                             .desc = F),
           affinity = ifelse(estimate >= 50, "low", "high")) %>% 
    ggplot(aes(x = gn, y = estimate, color = bait)) +
    geom_pointrange(aes(ymin = estimate - std.error, ymax = estimate + std.error),
                  width = 1, fatten = 2) +
    theme(legend.position = "none",
          axis.title.y = element_blank()) +
    facet_wrap(~bait, scales = "free_y", drop = T, ncol = 1) +
    ylab(expression("K"[d]~"(nM)")) +
    coord_flip()
savePlot(exBar, w = full / 3, h = 5)
```

## Info on specific proteins

``` r
# lrp1 biogrid
bioGridLrp1 <- read_tsv("data/bioGrid_lrp1.txt")
lrp1Genes <- unique(c(bioGridLrp1$`Official Symbol Interactor A`, 
                      bioGridLrp1$`Official Symbol Interactor B`))

interactors %>%
    left_join(protInfo) %>%
    filter(gn %in% lrp1Genes,
           bait == "GST-LRP1-ICD") %>%
    kable
```

| Protein                | bait         | term |  estimate | std.error | statistic |   p.value |        CV |      Rmax | Rmax\_CV |   maxConc | accession | gn     |        MW |
| :--------------------- | :----------- | :--- | --------: | --------: | --------: | --------: | --------: | --------: | -------: | --------: | :-------- | :----- | --------: |
| sp|Q9H7Z7|PGES2\_HUMAN | GST-LRP1-ICD | Kd   |  13.30562 |  11.80215 | 1.1273897 | 0.2967340 |  88.70047 |  40738541 | 16.74347 | 0.5189505 | Q9H7Z7    | PTGES2 |  41.94261 |
| sp|P04350|TBB4A\_HUMAN | GST-LRP1-ICD | Kd   | 409.46895 | 621.99500 | 0.6583155 | 0.5313825 | 151.90285 |  77633193 | 51.48883 | 1.3914422 | P04350    | TUBB4A |  49.58524 |
| sp|P29353|SHC1\_HUMAN  | GST-LRP1-ICD | Kd   | 125.36619 | 279.03546 | 0.4492841 | 0.6668074 | 222.57633 | 104916094 | 57.74150 | 0.4258674 | P29353    | SHC1   |  62.82140 |
| sp|O60610|DIAP1\_HUMAN | GST-LRP1-ICD | Kd   |  15.10822 |  51.35260 | 0.2942056 | 0.7771335 | 339.89837 |  24207452 | 65.07567 | 0.2227142 | O60610    | DIAPH1 | 141.34574 |

``` r
# lrp1b biogrid
bioGridLrp1b <- read_tsv("data/bioGrid_lrp1b.txt")
lrp1bGenes <- unique(c(bioGridLrp1b$`Official Symbol Interactor A`, 
                      bioGridLrp1b$`Official Symbol Interactor B`))

interactors %>%
    left_join(protInfo) %>%
    filter(gn %in% lrp1bGenes,
           bait == "GST-LRP1B-ICD") %>%
    kable
```

| Protein               | bait          | term | estimate | std.error | statistic |   p.value |       CV |      Rmax | Rmax\_CV |  maxConc | accession | gn     |       MW |
| :-------------------- | :------------ | :--- | -------: | --------: | --------: | --------: | -------: | --------: | -------: | -------: | :-------- | :----- | -------: |
| sp|P30533|AMRP\_HUMAN | GST-LRP1B-ICD | Kd   | 352.8549 |  520.5056 | 0.6779081 | 0.5196015 | 147.5126 | 816802262 | 48.08483 | 4.952549 | P30533    | LRPAP1 | 41.46546 |

``` r
# all high conf
highConf %>%
    left_join(protInfo) %>%
    kable
```

| Protein                | bait          | term |   estimate |  std.error | statistic |   p.value |       CV |       Rmax |  Rmax\_CV |   maxConc | accession | gn       |        MW |
| :--------------------- | :------------ | :--- | ---------: | ---------: | --------: | --------: | -------: | ---------: | --------: | --------: | :-------- | :------- | --------: |
| sp|Q96EV2|RBM33\_HUMAN | GST-LRP1-ICD  | Kd   |  11.947681 |   5.565052 |  2.146913 | 0.0689271 | 46.57851 |  115853488 |  8.689141 | 0.6834141 | Q96EV2    | RBM33    | 129.98434 |
| sp|Q2M2I8|AAK1\_HUMAN  | GST-LRP1-ICD  | Kd   |   4.955850 |   2.622861 |  1.889483 | 0.1007519 | 52.92453 |   27982828 |  8.979485 | 0.2092634 | Q2M2I8    | AAK1     | 103.88373 |
| sp|Q7L775|EPMIP\_HUMAN | GST-LRP1-ICD  | Kd   |   3.909935 |   2.087738 |  1.872809 | 0.1032510 | 53.39574 |   41188998 |  8.815560 | 0.3306080 | Q7L775    | EPM2AIP1 |  70.36872 |
| sp|P18583|SON\_HUMAN   | GST-LRP1B-ICD | Kd   | 264.772196 | 158.188330 |  1.673778 | 0.1380881 | 59.74507 | 4326930615 | 18.145617 | 3.4229295 | P18583    | SON      | 263.82720 |
| sp|Q9H222|ABCG5\_HUMAN | GST-LRP1B-ICD | Kd   |  37.904156 |  22.918367 |  1.653877 | 0.1421291 | 60.46399 |  531148294 | 12.962521 | 1.7893926 | Q9H222    | ABCG5    |  72.50287 |
| sp|A6NKG5|RTL1\_HUMAN  | GST-LRP1-ICD  | Kd   | 131.310677 |  83.654718 |  1.569675 | 0.1604814 | 63.70748 |  599677205 | 16.672600 | 4.0917805 | A6NKG5    | RTL1     | 155.04569 |
| sp|Q9NY93|DDX56\_HUMAN | GST-LRP1B-ICD | Kd   |   8.867674 |   5.674284 |  1.562783 | 0.1620769 | 63.98842 |  161237908 | 11.560072 | 0.8417132 | Q9NY93    | DDX56    |  61.58879 |
| sp|Q86U86|PB1\_HUMAN   | GST-LRP1-ICD  | Kd   |  75.012471 |  51.539346 |  1.455441 | 0.1888807 | 68.70770 |  147514543 | 16.297577 | 0.7606412 | Q86U86    | PBRM1    | 192.94556 |
| sp|Q9UKA9|PTBP2\_HUMAN | GST-LRP1B-ICD | Kd   |   9.115144 |   6.325854 |  1.440935 | 0.1927961 | 69.39938 |   57346725 | 12.574494 | 0.2588348 | Q9UKA9    | PTBP2    |  57.49006 |
| sp|Q08209|PP2BA\_HUMAN | GST-LRP1-ICD  | Kd   |  24.474112 |  17.059386 |  1.434642 | 0.1945173 | 69.70380 |  144907396 | 14.123416 | 0.8168949 | Q08209    | PPP3CA   |  58.68719 |
| sp|P06280|AGAL\_HUMAN  | GST-LRP1B-ICD | Kd   |  13.996039 |  10.087266 |  1.387496 | 0.2078559 | 72.07229 |   90571788 | 13.681156 | 0.5763919 | P06280    | GLA      |  48.76626 |
| sp|P22033|MUTA\_HUMAN  | GST-LRP1B-ICD | Kd   |  22.607308 |  16.373158 |  1.380754 | 0.2098283 | 72.42418 |   66115314 | 14.533222 | 0.2132138 | P22033    | MUT      |  83.13350 |
| sp|Q5JTZ9|SYAM\_HUMAN  | GST-LRP1-ICD  | Kd   |  30.749012 |  22.471784 |  1.368339 | 0.2135041 | 73.08132 |   35113641 | 15.239823 | 0.3182354 | Q5JTZ9    | AARS2    | 107.33917 |
| sp|P17612|KAPCA\_HUMAN | GST-LRP1B-ICD | Kd   |  10.878119 |   8.559185 |  1.270929 | 0.2443572 | 78.68258 |   99172491 | 14.529560 | 0.6038588 | P17612    | PRKACA   |  40.58915 |
| sp|P50995|ANX11\_HUMAN | GST-LRP1-ICD  | Kd   |  50.984417 |  40.442356 |  1.260669 | 0.2478214 | 79.32297 |   28603734 | 17.730047 | 0.1938662 | P50995    | ANXA11   |  54.38910 |
| sp|Q8TED0|UTP15\_HUMAN | GST-LRP1-ICD  | Kd   | 145.541122 | 116.402317 |  1.250328 | 0.2513550 | 79.97899 |   30684894 | 21.351607 | 0.9984038 | Q8TED0    | UTP15    |  58.41450 |
| sp|P51116|FXR2\_HUMAN  | GST-LRP1B-ICD | Kd   |  54.693310 |  44.520047 |  1.228510 | 0.2589523 | 81.39944 |   86545241 | 18.383586 | 0.2926980 | P51116    | FXR2     |  74.22244 |

## Write tables

``` r
# to be added
```

# Session Info

``` r
devtools::session_info()
```

    ##  setting  value                       
    ##  version  R version 3.5.0 (2018-04-23)
    ##  system   x86_64, darwin15.6.0        
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_US.UTF-8                 
    ##  tz       America/New_York            
    ##  date     2018-05-17                  
    ## 
    ##  package    * version date       source        
    ##  ade4         1.7-11  2018-04-05 CRAN (R 3.5.0)
    ##  assertthat   0.2.0   2017-04-11 CRAN (R 3.5.0)
    ##  backports    1.1.2   2017-12-13 CRAN (R 3.5.0)
    ##  base       * 3.5.0   2018-04-24 local         
    ##  bindr        0.1.1   2018-03-13 CRAN (R 3.5.0)
    ##  bindrcpp   * 0.2.2   2018-03-29 CRAN (R 3.5.0)
    ##  broom      * 0.4.4   2018-03-29 CRAN (R 3.5.0)
    ##  cellranger   1.1.0   2016-07-27 CRAN (R 3.5.0)
    ##  cli          1.0.0   2017-11-05 CRAN (R 3.5.0)
    ##  colorspace   1.3-2   2016-12-14 CRAN (R 3.5.0)
    ##  compiler     3.5.0   2018-04-24 local         
    ##  crayon       1.3.4   2017-09-16 CRAN (R 3.5.0)
    ##  datasets   * 3.5.0   2018-04-24 local         
    ##  devtools     1.13.5  2018-02-18 CRAN (R 3.5.0)
    ##  digest       0.6.15  2018-01-28 CRAN (R 3.5.0)
    ##  dplyr      * 0.7.4   2017-09-28 CRAN (R 3.5.0)
    ##  evaluate     0.10.1  2017-06-24 CRAN (R 3.5.0)
    ##  forcats    * 0.3.0   2018-02-19 CRAN (R 3.5.0)
    ##  foreign      0.8-70  2017-11-28 CRAN (R 3.5.0)
    ##  ggplot2    * 2.2.1   2016-12-30 CRAN (R 3.5.0)
    ##  ggridges   * 0.5.0   2018-04-05 CRAN (R 3.5.0)
    ##  glue         1.2.0   2017-10-29 CRAN (R 3.5.0)
    ##  graphics   * 3.5.0   2018-04-24 local         
    ##  grDevices  * 3.5.0   2018-04-24 local         
    ##  grid         3.5.0   2018-04-24 local         
    ##  gtable       0.2.0   2016-02-26 CRAN (R 3.5.0)
    ##  haven        1.1.1   2018-01-18 CRAN (R 3.5.0)
    ##  highr        0.6     2016-05-09 CRAN (R 3.5.0)
    ##  hms          0.4.2   2018-03-10 CRAN (R 3.5.0)
    ##  htmltools    0.3.6   2017-04-28 CRAN (R 3.5.0)
    ##  httr         1.3.1   2017-08-20 CRAN (R 3.5.0)
    ##  jsonlite     1.5     2017-06-01 CRAN (R 3.5.0)
    ##  knitr      * 1.20    2018-02-20 CRAN (R 3.5.0)
    ##  labeling     0.3     2014-08-23 CRAN (R 3.5.0)
    ##  lattice      0.20-35 2017-03-25 CRAN (R 3.5.0)
    ##  lazyeval     0.2.1   2017-10-29 CRAN (R 3.5.0)
    ##  lubridate    1.7.4   2018-04-11 CRAN (R 3.5.0)
    ##  magrittr     1.5     2014-11-22 CRAN (R 3.5.0)
    ##  MASS         7.3-49  2018-02-23 CRAN (R 3.5.0)
    ##  memoise      1.1.0   2017-04-21 CRAN (R 3.5.0)
    ##  methods    * 3.5.0   2018-04-24 local         
    ##  mnormt       1.5-5   2016-10-15 CRAN (R 3.5.0)
    ##  modelr       0.1.1   2017-07-24 CRAN (R 3.5.0)
    ##  munsell      0.4.3   2016-02-13 CRAN (R 3.5.0)
    ##  nlme         3.1-137 2018-04-07 CRAN (R 3.5.0)
    ##  parallel     3.5.0   2018-04-24 local         
    ##  pillar       1.2.2   2018-04-26 CRAN (R 3.5.0)
    ##  pkgconfig    2.0.1   2017-03-21 CRAN (R 3.5.0)
    ##  plyr         1.8.4   2016-06-08 CRAN (R 3.5.0)
    ##  psych        1.8.3.3 2018-03-30 CRAN (R 3.5.0)
    ##  purrr      * 0.2.4   2017-10-18 CRAN (R 3.5.0)
    ##  R6           2.2.2   2017-06-17 CRAN (R 3.5.0)
    ##  Rcpp         0.12.16 2018-03-13 CRAN (R 3.5.0)
    ##  readr      * 1.1.1   2017-05-16 CRAN (R 3.5.0)
    ##  readxl       1.1.0   2018-04-20 CRAN (R 3.5.0)
    ##  reshape2     1.4.3   2017-12-11 CRAN (R 3.5.0)
    ##  rlang        0.2.0   2018-02-20 CRAN (R 3.5.0)
    ##  rmarkdown    1.9     2018-03-01 CRAN (R 3.5.0)
    ##  rprojroot    1.3-2   2018-01-03 CRAN (R 3.5.0)
    ##  rstudioapi   0.7     2017-09-07 CRAN (R 3.5.0)
    ##  rvest        0.3.2   2016-06-17 CRAN (R 3.5.0)
    ##  scales       0.5.0   2017-08-24 CRAN (R 3.5.0)
    ##  seqinr     * 3.4-5   2017-08-01 CRAN (R 3.5.0)
    ##  stats      * 3.5.0   2018-04-24 local         
    ##  stringi      1.2.2   2018-05-02 CRAN (R 3.5.0)
    ##  stringr    * 1.3.0   2018-02-19 CRAN (R 3.5.0)
    ##  tibble     * 1.4.2   2018-01-22 CRAN (R 3.5.0)
    ##  tidyr      * 0.8.0   2018-01-29 CRAN (R 3.5.0)
    ##  tidyselect   0.2.4   2018-02-26 CRAN (R 3.5.0)
    ##  tidyverse  * 1.2.1   2017-11-14 CRAN (R 3.5.0)
    ##  tools        3.5.0   2018-04-24 local         
    ##  utf8         1.1.3   2018-01-03 CRAN (R 3.5.0)
    ##  utils      * 3.5.0   2018-04-24 local         
    ##  withr        2.1.2   2018-03-15 CRAN (R 3.5.0)
    ##  xml2         1.2.0   2018-01-24 CRAN (R 3.5.0)
    ##  yaml         2.1.19  2018-05-01 CRAN (R 3.5.0)
