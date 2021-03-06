% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/yang_and_smith.R
\name{blast_to_mcl}
\alias{blast_to_mcl}
\title{Prepare BLAST results for MCL.}
\usage{
blast_to_mcl(
  path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"),
  blast_results,
  hit_fraction_cutoff,
  echo = pkgconfig::get_config("baitfindR::echo", fallback = FALSE),
  ...
)
}
\arguments{
\item{path_to_ys}{Character vector of length one; the complete path to the folder containing Y&S python scripts, e.g., "/Users/me/apps/phylogenomic_dataset_construction/"}

\item{blast_results}{Character vector of length one; the complete path to the tab-separated text file containing the results from an all-by-all blast search. If blast searches were run separately (i.e., one for each sample), the results should be concatenated into a single file. For the blast search, the output format should specified as: -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore'}

\item{hit_fraction_cutoff}{Numeric between 0 and 1. Indicates the minimum percentage overlap between query and target in blast results to be retained in the output. According to Y&S, "A low hit-fraction cutoff will output clusters with more incomplete sequences and much larger and sparser alignments, whereas a high hit-fraction cutoff gives tighter clusters but ignores incomplete or divergent sequences."}

\item{echo}{Logical; should the standard output and error be printed to the screen?}

\item{...}{Other arguments. Not used by this function, but meant to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.}
}
\value{
A tab-separated text file with three columns: the first two are the matching query and target from the all-by-all blast, and the third is the negative log e-value for that match. This file is named \code{<blast_results>.hit-frac<hit_fraction_cutoff>.minusLogEvalue}, where \code{<blast_results>} and \code{<hit_fraction_cutoff>} correspond to the values of those arguments. If possible contaminants (i.e., identical sequences between different samples) were found, these are written to \code{<blast_results>.ident.hit-frac<hit_fraction_cutoff>}. Output files will be written to the same folder containing \code{blast_results}.
}
\description{
Converts the output of an all-by-all blast query into a format that can be parsed by mcl to find clusters.
}
\details{
Wrapper for Yang and Smith (2014) blast_to_mcl.py
}
\examples{
\dontrun{blast_to_mcl(blast_results = "some/folder/blastresults.tab", hit_fraction_cutoff = 0.5)}
}
\references{
Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
}
\author{
Joel H Nitta, \email{joelnitta@gmail.com}
}
