% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/yang_and_smith.R
\name{prune_paralogs_MI}
\alias{prune_paralogs_MI}
\title{Identify orthologs using the "MI" method.}
\usage{
prune_paralogs_MI(
  path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"),
  tree_folder,
  tree_file_ending,
  relative_cutoff,
  absolute_cutoff,
  minimal_taxa = 4,
  outdir,
  overwrite = FALSE,
  get_hash = TRUE,
  echo = pkgconfig::get_config("baitfindR::echo", fallback = FALSE),
  ...
)
}
\arguments{
\item{path_to_ys}{Character vector of length one; the path to the folder containing Y&S python scripts, e.g., "/Users/me/apps/phylogenomic_dataset_construction/"}

\item{tree_folder}{Character vector of length one; the path to the folder containing the trees to be used for pruning.}

\item{tree_file_ending}{Character vector of length one; only tree files with this file ending will be used.}

\item{relative_cutoff}{Numeric vector of length one; tips on a branch 10 times longer than their sister AND longer than this value will be trimmed.}

\item{absolute_cutoff}{Numeric vector of length one; tips on branches longer than this value will be trimmed.}

\item{minimal_taxa}{Numeric; minimal number of taxa required for tree to be included. Default 4, the minimum number of taxa needed for an un-rooted tree.}

\item{outdir}{Character vector of length one; the path to the folder where the pruned trees should be written.}

\item{overwrite}{Logical; should previous output of this command be erased so new output can be written? Once erased it cannot be restored, so use with caution!}

\item{get_hash}{Logical; should the 32-byte MD5 hash be computed for all pruned tree files concatenated together? Used for by \code{\link[drake]{drake_plan}} for tracking during workflows. If \code{TRUE}, this function will return the hash.}

\item{echo}{Logical; should the standard output and error be printed to the screen?}

\item{...}{Other arguments. Not used by this function, but meant to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.}
}
\value{
For each tree file ending in \code{tree_file_ending} in \code{tree_folder}, putative orthologs will be extracted from the tree using the MI method and written to \code{outdir} with the file ending \code{.MIortho1.tre}. If \code{get_hash} is \code{TRUE}, the 32-byte MD5 hash be computed for all extracted tree files concatenated together will be returned.
}
\description{
Given a folder containing homolog trees, prune paralogs from the trees
using the maximum inclusion (MI) method. This method iteratively
extracts subtrees containing the greatest number of samples (taxa)
without duplication. Long branches will also be trimmed (i.e., removed)
according to \code{relative_cutoff} and \code{absolute_cutoff} values. Trees
that consist solely of one-to-one orthologs (i.e., no duplications within a sample)
will also be retained. This function will overwrite any output files with
the same name in \code{outdir}.
}
\details{
Wrapper for Yang and Smith (2014) prune_paralogs_MI.py
}
\examples{
\dontrun{prune_paralogs_MI(tree_folder = "some/folder/containing/tree/files", tree_file_ending = ".tre", relative_cutoff = 0.2, absolute_cutoff = 0.4, outdir = "some/folder")}
}
\references{
Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
}
\author{
Joel H Nitta, \email{joelnitta@gmail.com}
}
