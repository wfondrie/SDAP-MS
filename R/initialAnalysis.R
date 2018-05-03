library(tidyverse)
library(broom)
source("R/ggplotTheme.R")

theme_set(coolTheme)

# Import MaxQuant data --------------------------------------------------------
prot <- read_tsv("data/combined/txt/proteinGroups.txt") %>%
  filter(!str_detect(`Protein IDs`, "^(REV|CON)__"),
         Peptides >= 2) 

# Reshape data for modeling ---------------------------------------------------
lfq <- prot %>%
  select(`Protein IDs`, starts_with("LFQ ")) %>%
  gather(samp, lfq, -`Protein IDs`) %>%
  mutate(samp = str_match(samp, " (LRP1.*_.*)$")[ , 2],
         bait = str_match(samp, "(LRP1.*)_")[ , 2],
         conc = str_match(samp, "_(.*)$")[ , 2],
         conc = ifelse(conc == "GST", 0, 3^as.numeric(conc)))

lfq %>% 
  group_by(bait) %>%
  summarize(proteins = length(unique(`Protein IDs`)))

# Perform censoring ---------------------------------------------------------
lfqCensored <- lfq %>%
  group_by(`Protein IDs`, bait) %>%
  filter(sum(lfq) > 0) %>%
  mutate(maxConc = max(conc[lfq > 0]),
         lfq = ifelse(lfq == 0 & conc < maxConc, NA, lfq),
         numNA = sum(is.na(lfq)),
         numZero = sum(lfq == 0, na.rm = T))

lfqCensored %>% 
  group_by(bait) %>%
  summarize(proteins = length(unique(`Protein IDs`)))

# Filter before modeling ----------------------------------------------------
lfqFiltered <- lfqCensored %>%
  group_by(`Protein IDs`, bait) %>%
  filter(numNA <= 4,
         numZero <= 3,
         numNA + numZero <= 5)

lfqFiltered %>% 
  group_by(bait) %>%
  summarize(proteins = length(unique(`Protein IDs`)))


# Transform lfq to response -------------------------------------------------
lfqModInput <- lfqFiltered %>%
  group_by(`Protein IDs`, bait) %>%
  mutate(response = -(lfq - max(lfq)))

# Modeling ------------------------------------------------------------------
mods <- lfqModInput %>%
  rename(Protein = `Protein IDs`) %>% #glance() doesn't seem to like the long one
  group_by(Protein, bait) %>%
  filter(!is.na(response)) %>%
  do(models = nls(response ~ (Rmax * conc) / (Kd + conc),
                  data = .,
                  start = list(Kd = 100, Rmax = 4e+07),
                  control = list(maxiter = 200,
                                 warnOnly = T)))

fitInfo <- mods %>% glance(models)

fitVals <- mods %>% 
  tidy(models) %>%
  mutate(CV = std.error / estimate * 100)

fits <- fitInfo %>%
  filter(isConv) %>%
  select(Protein, bait) %>%
  left_join(fitVals)

tblKd <- fits %>%
  filter(term == "Kd")

tblRmax <- fits %>%
  filter(term == "Rmax") %>%
  select(Protein, bait, estimate, CV) %>%
  rename(Rmax = estimate, Rmax_CV = CV)

fitTbl <- tblKd %>%
  left_join(tblRmax)

# Possible Interactions -------------------------------------------------------
interactors <- fitTbl %>%
  filter(estimate > 3, 
         estimate < 2187,
         CV < 100) %>%
  arrange(CV)

interactors %>% 
    group_by(bait) %>%
    summarize(proteins = length(unique(Protein)))

bothInt <- interactors %>%
    group_by(Protein) %>%
    summarize(number = length(Protein)) %>%
    filter(number == 2) %>%
    left_join(interactors) %>%
    group_by(Protein) %>%
    filter(sum(CV) < 200)

bothInt

hist(interactors$estimate)
hist(interactors$CV)

interactors %>%
    ggplot(aes(x = log10(estimate), y = CV)) +
    geom_point() +
    geom_hline(yintercept = 100, color = "blue")

# more stuff
tbl <- interactors %>%
    rename(`Protein IDs` = Protein) %>%
    left_join(select(prot, `Protein IDs`, `Gene names`)) %>%
    filter(bait == "LRP1B")

# Interesting interactions -----------------------------------------------------
interactors %>%
    ungroup() %>%
    rename(`Protein IDs` = Protein) %>%
    top_n(-9, CV) %>%
    left_join(lfqModInput) %>%
    left_join(select(prot, `Protein IDs`, `Gene names`)) %>%
    ggplot(aes(x = conc, y = response)) +
    geom_smooth(method = "nls",
                formula = y ~ (x * Rmax)/ (Kd + x),
                method.args = list(start = list(Kd = 100, Rmax = 4e+07),
                                   control = list(maxiter = 200, warnOnly = T)),
                se = F,
                color = rgb(249, 38, 114, maxColorValue = 255),
                size = 2) +
    geom_point() +
    facet_wrap(~`Gene names` + bait, ncol = 3, scales = "free_y")


interactors %>%
    ungroup() %>%
    rename(`Protein IDs` = Protein) %>%
    top_n(-25, CV) %>%
    left_join(lfqModInput) %>%
    left_join(select(prot, `Protein IDs`, `Gene names`)) %>%
    group_by(`Protein IDs`, bait) %>%
    mutate(logConc = log10(conc),
           response = response/max(response)) %>%
    ggplot(aes(x = logConc, y = response)) +
    geom_smooth(method = "nls",
                formula = y ~ (10^x * Rmax)/ (Kd + 10^x),
                method.args = list(start = list(Kd = 100, Rmax = 4e+07),
                                   control = list(maxiter = 200, warnOnly = T)),
                se = F,
                color = rgb(249, 38, 114, maxColorValue = 255),
                size = 2) +
    geom_point() +
    facet_wrap(~`Gene names` + bait, ncol = 5, scales = "free_y")
