library(tidyverse)
library(broom)

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
         numZero <= 2)

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
  filter(estimate > 1, 
         estimate < 2187) %>%
  arrange(CV)

hist(interactors$estimate)
hist(interactors$CV)
