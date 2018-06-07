#' fasta_to_tree
#'
#' Wrapper for Yang and Smith (2014) fasta_to_tree.py
#'
#' Given a folder containing unaligned sequences in fasta format (i.e., clusters),
#' aligns each cluster with mafft (small alignments) or pasta (large alignments),
#' excludes poorly aligned sites with phyutility, and infers a maximum-likelihood
#' tree with RAxML (small alignments) or fasttree (large alignments). Requires all
#' of these programs to be installed and included in the user's `$PATH`.
#'
#' @param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., "/Users/me/apps/phylogenomic_dataset_construction/"
#' @param seq_folder Character vector of length one; the name of the folder containing the fasta files.
#' @param number_cores Numeric; number of threads to use for RAxML.
#' @param seq_type Character vector of length one indicating type of sequences. Should either be "dna" for DNA or "aa" for proteins.
#' @param bootstrap Logical; choose whether or not to run a bootstrap analysis for the trees.
#'
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#'
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#'
#' @examples
#' \dontrun{fasta_to_tree(seq_folder = "some/folder/containing/fasta/seqs", number_cores = 1, seq_type = "dna", bootstrap = FALSE)}
#'
#' @export
fasta_to_tree <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), seq_folder, number_cores, seq_type, bootstrap = FALSE) {
  bootstrap <- ifelse(bootstrap, "y", "n")
  arguments <- paste0(path_to_ys, "fasta_to_tree.py", " ", seq_folder, " ", number_cores, " ", seq_type, " ", bootstrap)
  system2("python", arguments)
}

#' write_fasta_files_from_mcl
#'
#' Wrapper for Yang and Smith (2014) write_fasta_files_from_mcl.py
#'
#' Given the output from the mcl clustering algorthim and a concatenated fasta file
#' including all sequences used for clustering, output one fasta file per cluster
#' including the sequences in that cluster.
#'
#' @param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., "/Users/me/apps/phylogenomic_dataset_construction/"
#' @param all_fasta Character vector of length one; the path to the fasta file
#' including all query sequences concatenated together, i.e., the fasta file used to
#' create the "all-by-all" blast database.
#' @param mcl_outfile Character vector of length one; the path to the output from running mcl on blast distances.
#' @param minimal_taxa Numeric; minimal number of taxa required to be present or the clus
#' @param outdir Character vector of length one; the path to the folder where the clusters should be written.
#'
#' @export
write_fasta_files_from_mcl <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), all_fasta, mcl_outfile, minimal_taxa, outdir) {
  arguments <- paste0(path_to_ys, "write_fasta_files_from_mcl.py", " ", all_fasta, " ", mcl_outfile, " ", minimal_taxa, " ", outdir)
  system2("python", arguments)
}
