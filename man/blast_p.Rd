% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blast.R
\name{blast_p}
\alias{blast_p}
\title{Run a blastp query.}
\usage{
blast_p(
  query,
  database,
  out_file = NULL,
  outfmt = "6",
  other_args = NULL,
  echo = TRUE,
  wd,
  ...
)
}
\arguments{
\item{query}{Character vector of length one; the path to the fasta file
to use as the query sequence(s).}

\item{database}{Character vector of length one; the name of the blast
database.}

\item{out_file}{Character vector of length one; the name to use for the
results file.}

\item{outfmt}{Character vector of length one; value to pass to
\code{blastp} \code{outfmt} argument. Default = "6".}

\item{other_args}{Character vector; other arguments to pass on to
\code{blastp}. For a list of options, run \code{blastp -help}.}

\item{echo}{Logical; should standard error and output be printed?}

\item{wd}{Character vector of length one; working directory. The blast
search will be conducted here.}

\item{...}{Additional other arguments. Not used by this function,
but meant to be used by \code{\link[drake]{drake_plan}} for tracking
during workflows.}
}
\value{
A tab-separated text file with the results of the blastp
search, named with the value of \code{out_file}.
}
\description{
This is a wrapper for blastp.
}
\examples{
library(ape)

# Make temp dir for storing files
temp_dir <- fs::dir_create(fs::path(tempdir(), "baitfindR_example"))

# Write out ape::woodmouse dataset as amino acids
data(woodmouse)
woodmouse_aa <- trans(woodmouse, 2)
ape::write.FASTA(woodmouse_aa, fs::path(temp_dir, "woodmouse.fasta"))

# Make protein blast database
build_blast_db(
  fs::path(temp_dir, "woodmouse.fasta"),
  db_type = "prot",
  out_name = "wood",
  parse_seqids = TRUE,
  wd = temp_dir)

# Blast the original sequences against the database
blast_p(
  fs::path(temp_dir, "woodmouse.fasta"),
  database = "wood",
  out_file = "blastp_results",
  wd = temp_dir,
  echo = TRUE
)

# Take a look at the results.
readr::read_tsv(
  fs::path(temp_dir, "blastp_results"),
  col_names = FALSE
  )

# Cleanup.
fs::file_delete(temp_dir)
}
\references{
\url{https://www.ncbi.nlm.nih.gov/books/NBK279690/}
}
\author{
Joel H Nitta, \email{joelnitta@gmail.com}
}
