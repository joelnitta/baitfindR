% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/yang_and_smith.R
\name{prune_paralogs_RT}
\alias{prune_paralogs_RT}
\title{Identify orthologs using the "RT" method.}
\usage{
prune_paralogs_RT(
  path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"),
  tree_folder,
  tree_file_ending,
  ingroup,
  outgroup,
  min_ingroup_taxa = 2,
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

\item{ingroup}{Character vector; names of ingroup taxa/samples.}

\item{outgroup}{Character vector; names of outgroup taxa/samples.}

\item{min_ingroup_taxa}{Numeric; minimal number of taxa in the ingroup required for an ortholog to be written. Default 2.}

\item{outdir}{Character vector of length one; the path to the folder where the pruned trees should be written.}

\item{overwrite}{Logical; should previous output of this command be erased so new output can be written? Once erased it cannot be restored, so use with caution!}

\item{get_hash}{Logical; should the 32-byte MD5 hash be computed for all pruned tree files concatenated together? Used for by \code{\link[drake]{drake_plan}} for tracking during workflows. If \code{TRUE}, this function will return the hash.}

\item{echo}{Logical; should the standard output and error be printed to the screen?}

\item{...}{Other arguments. Not used by this function, but meant to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.}
}
\value{
For each tree file ending in \code{tree_file_ending} in \code{tree_folder}, the following outputs are possible depending on the presence of outgroups in the homolog tree:
\describe{
  \item{No outgroups in homolog tree}{Unrooted ingroup clades without duplications (files ending in \code{unrooted-ortho.tre})}
  \item{Outgroups present in homolog tree}{Rooted ingroup clades (files ending in \code{inclade}) and one or more paralogs (files ending in \code{inclade.ortho.tre})}
}

If \code{get_hash} is \code{TRUE}, the 32-byte MD5 hash be computed for all extracted tree files concatenated together will be returned.
}
\description{
Given a folder containing homolog trees, prune paralogs from the trees
using the rooted ingroups (RT) method. For trees that have outgroups, this
method iteratively extracts subtrees with the highest number of ingroup taxa/samples.
This function will overwrite any output files with the same
name in \code{outdir}.
}
\details{
Wrapper for Yang and Smith (2014) prune_paralogs_RT.py
}
\examples{
\dontrun{prune_paralogs_RT(tree_folder = "some/folder/containing/tree/files", tree_file_ending = ".tre", outgroup = c("ABC", "EFG"), ingroup = c("HIJ", "KLM"), outdir = "some/folder")}
}
\references{
Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
}
\author{
Joel H Nitta, \email{joelnitta@gmail.com}
}
