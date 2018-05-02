# script to make small transcriptome dataset for testing
# download 1KP transcriptomes, randomly trim down to specified fraction of original size

library(ape)

# make vector of 1KP sample codes to use
# (all eupolypod II ferns)
# codes <- c("PSKY", "KJZG", "URCP", "FCHS", "UFJN", "VITX", "LHLE", "XXHP", "YOWV", "RICC", "HNDZ", "HEGQ", "OCZL", "HTFH", "MROH", "YQEC", "YJJY")
# (subset of eupolypod II ferns incl Aspleniaceae, Athyriaceae, and Woodsiaceae)
# note that AFPO no longer is a separate sample; see "RWYZ: combined assembly of AFPO+VITX"
codes <- c("PSKY", "KJZG", "YJJY", "YQEC", "URCP", "FCHS", "UFJN")

# function to download transcriptome from 1KP website, load the data, and randomly downsize
# arguments
#   sample_code: the 1KP sample code
#   keep_frac: the percentage of transcripts to keep, e.g. 0.10 = 10%
get_trimmed_transcriptome <- function (sample_code, keep_frac) {
  # set a name for the to-be downloaded compressed transcriptome
  filename <- paste0(sample_code, "-SOAPdenovo-Trans-assembly.fa.bz2")
  # download compressed transcriptome
  download.file(url = paste0("http://206.12.96.204/onekp/", sample_code, "-SOAPdenovo-Trans-assembly.fa.bz2"), destfile = paste0("data-raw/", filename))
  # read in file. need bzfile() because it's compressed
  seq <- read.dna(bzfile(paste0("data-raw/", filename)), format="fasta")
  # randomly downsize to specified fraction of transcripts
  seq <- seq[sample(1:length(seq), keep_frac*length(seq))]
}

# loop across all codes, keeping 2% of transcripts
trimmed_transcriptomes <- lapply(codes, get_trimmed_transcriptome, keep_frac = 0.02)
# name the transcriptomes
names(trimmed_transcriptomes) <- codes
# export data
usethis::use_data(trimmed_transcriptomes, overwrite = TRUE)
