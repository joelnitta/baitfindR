# script to make small transcriptome dataset for testing
# read in 1KP transcriptomes, randomly trim down to 10% of size

library(ape)

# vector of 1KP sample codes to use (all eupolypod II ferns)
# codes <- c("PSKY", "KJZG", "URCP", "AFPO", "FCHS", "UFJN", "VITX", "LHLE", "XXHP", "YOWV", "RICC", "HNDZ", "HEGQ", "OCZL", "HTFH", "MROH", "YQEC", "YJJY")
# vector of 1KP sample codes to use (subset of eupolypod II ferns incl Aspleniaceae, Athyriaceae, and Woodsiaceae)
codes <- c("PSKY", "KJZG", "YJJY", "YQEC", "AFPO", "URCP", "FCHS", "UFJN")

# function to read in data, randomly subset to 2% of sequences
# we want a very small dataset just for testing, not for complete analysis
trim_fasta <- function (code) {
  seq <- read.dna(paste0("/Volumes/Storage/thelypteris/data/one_kp/", code, ".fa"), format="fasta")
  seq <- seq[sample(1:length(seq), 0.02*length(seq))]
}

trimmed_transcriptomes <- lapply(codes, trim_fasta)
names(trimmed_transcriptomes) <- codes

usethis::use_data(trimmed_transcriptomes)
