% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transdecoder.R
\name{cd_hit_est}
\alias{cd_hit_est}
\title{Cluster DNA sequences.}
\usage{
cd_hit_est(
  input,
  output,
  wd = here::here(),
  other_args = NULL,
  echo = pkgconfig::get_config("baitfindR::echo", fallback = FALSE),
  ...
)
}
\arguments{
\item{input}{Character vector of length one; the path to the input file for
cd-hit-est. Should be DNA or AA sequences in fasta format.}

\item{output}{Character vector of length one; the name to assign to the output.
Can include a path, in which case the output will be written there.}

\item{wd}{Character vector of length one; the directory where the command
will be run.}

\item{other_args}{Character vector; other arguments to pass to cd-hit-est.
Each should be an element of the vector.}

\item{echo}{Logical; should the standard output and error be printed to the screen?}

\item{...}{Additional other arguments. Not used by this function, but meant
to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.}
}
\value{
Within the R environment, a list with components specified in
\code{\link[processx]{run}}.

Externally, two files will be written: according to the CD-HIT user's guide,
"The output are two files: a fasta file of representative sequences and a text
file of list of clusters."

The fasta file will be named with the value of \code{output}; the list of clusters
will be the same, with \code{.clstr} appended.
}
\description{
This is a wrapper for the CD-HIT-EST algorithm. According to the CD-HIT user's guide,
"CD-HIT-EST clusters a nucleotide dataset into clusters that meet a user-defined
similarity threshold, usually a sequence identity." cd-hit-est comes bundled with
transdecoder, so it is run from there.
}
\examples{
\dontrun{
library(ape)
library(baitfindR)

# Make temp dir for storing output
temp_dir <- fs::dir_create(fs::path(tempdir(), "baitfindR_example"))
data("PSKY")

# Write downsized transcriptome to temp dir
write.FASTA(PSKY, fs::path(temp_dir, "PSKY"))

# Get CDS
transdecoder_long_orfs(
  transcriptome_file = fs::path(temp_dir, "PSKY"),
  wd = temp_dir
  )

# Cluster similar genes in CDS
cd_hit_est(
  input = fs::path(temp_dir, "PSKY.transdecoder_dir", "longest_orfs.cds"),
  output = fs::path(temp_dir, "PSKY.cd-hit-est"),
  wd = temp_dir,
  echo = TRUE
)

# Check output
list.files(temp_dir)
head(readr::read_lines(fs::path(temp_dir, "PSKY.cd-hit-est")))
head(readr::read_lines(fs::path(temp_dir, "PSKY.cd-hit-est.clstr")))

# Cleanup
fs::file_delete(temp_dir)
}
}
\references{
\url{http://www.bioinformatics.org/cd-hit/}, \url{http://transdecoder.github.io}
}
\author{
Joel H Nitta, \email{joelnitta@gmail.com}
}
