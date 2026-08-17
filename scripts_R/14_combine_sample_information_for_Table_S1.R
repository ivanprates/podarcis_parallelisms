
## Packages:
library(tidyverse)

## Directory:
setwd(dir = "~/Dropbox/Science/MYPAPERS_ongoing/Podarcis_parallelisms/")

## Trait data:
t <- read.csv("morphological_analyses/data/Podarcis_greenbrowns_selected_traits_2026-03-15.csv", header = TRUE)
t$species_abbpop

## Locality data:
l <- read.csv("sample_information/sampling_sites.csv", header = TRUE)
l$abbpop
l$species_abbpop <- paste(l$species, l$abbpop, sep = "_")
l <- l[setdiff(names(l), "GB")]

## Combine:
m <- merge(l, t)

## Check:
#View(m)

## Save:
write.csv(m, file = "drafts/Table_S1_sample_information.csv", row.names = FALSE)

## End of script.
