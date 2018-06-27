# script to make dataframe of taxonomic data

library(tidyverse)
library(rvest)
library(ape)

# obtain 1KP sample data from official project website

# make vector of 1KP sample codes to use
# (all eupolypod II ferns)
# codes <- c("PSKY", "KJZG", "URCP", "FCHS", "UFJN", "VITX", "LHLE", "XXHP", "YOWV", "RICC", "HNDZ", "HEGQ", "OCZL", "HTFH", "MROH", "YQEC", "YJJY")

# (subset of eupolypod II ferns incl Aspleniaceae, Athyriaceae, Blechnaceae, and Woodsiaceae)
# plus two eupolypod I ferns as outgroups (OQWW, FQGQ)
# note that AFPO no longer is a separate sample; see "RWYZ: combined assembly of AFPO+VITX"
codes <- c("FQGQ", "PSKY", "KJZG", "YJJY", "YQEC", "URCP", "FCHS", "UFJN")

# read in html
onekp_parsed <- read_html("http://www.onekp.com/samples/list.php")
# extract tables
tables <- html_table(onekp_parsed, fill = TRUE)
# the table we want is the third one
onekp_data <- tables[[3]]
rm(onekp_parsed, tables)
# only need code and species name
# inspect table, see these are the first and fourth columns
onekp_data <- onekp_data[,c(1, 4)]
colnames(onekp_data) <- c("code", "species")

# subset data to only the transcriptomes we are using
onekp_data <- onekp_data[onekp_data$code %in% codes, ]

# add genus and family-level data
onekp_data <- onekp_data %>%
  separate(species, c("genus", "specific_epithet"), remove=FALSE, extra="merge") %>%
  mutate(
    family = case_when (
      genus == "Asplenium" ~ "Aspleniaceae",
      genus == "Athyrium" ~ "Athyriaceae",
      genus == "Blechnum" ~ "Blechnaceae",
      genus == "Cystopteris" ~ "Cystopteridaceae",
      genus == "Deparia" ~ "Athyriaceae",
      genus == "Diplazium" ~ "Athyriaceae",
      genus == "Gymnocarpium" ~ "Cystopteridaceae",
      genus == "Homalosorus" ~ "Diplaziopsidaceae",
      genus == "Onoclea" ~ "Onocleaceae",
      genus == "Thelypteris" ~ "Thelypteridaceae",
      genus == "Woodsia" ~ "Woodsiaceae",
      genus == "Davallia" ~ "Davalliaceae",
      genus == "Polystichum" ~ "Dryopteridaceae"
    )
  ) %>%
  mutate(
    group = case_when (
      family == "Davalliaceae" | family == "Dryopteridaceae" ~ "out",
      TRUE ~ "in"
    )
  )

# export data -------------------------------------------------------------

usethis::use_data(onekp_data, overwrite = TRUE)
