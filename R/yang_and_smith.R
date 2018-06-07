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
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{fasta_to_tree(seq_folder = "some/folder/containing/fasta/seqs", number_cores = 1, seq_type = "dna", bootstrap = FALSE)}
#' @export
fasta_to_tree <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), seq_folder, number_cores, seq_type, bootstrap = FALSE) {
  bootstrap <- ifelse(bootstrap, "y", "n")
  path_to_ys <- jntools::add_slash(path_to_ys)
  arguments <- c(paste0(path_to_ys, "fasta_to_tree.py"), seq_folder, number_cores, seq_type, bootstrap)
  system2("python", arguments)
}

#' write_fasta_files_from_mcl
#'
#' Wrapper for Yang and Smith (2014) write_fasta_files_from_mcl.py
#'
#' Given the output from the mcl clustering algorthim and a concatenated fasta file
#' including all sequences used for clustering, outputs one fasta file per cluster
#' including the sequences in that cluster.
#'
#' @param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., "/Users/me/apps/phylogenomic_dataset_construction/"
#' @param all_fasta Character vector of length one; the path to the fasta file including all query sequences concatenated together, i.e., the fasta file used to create the "all-by-all" blast database.
#' @param mcl_outfile Character vector of length one; the path to the output from running mcl on blast distances.
#' @param minimal_taxa Numeric; minimal number of taxa required to be present or the clus. Default 4, the minimum number of taxa needed for an un-rooted tree.
#' @param outdir Character vector of length one; the path to the folder where the clusters should be written.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{write_fasta_files_from_mcl(all_fasta = "some/folder/all.fasta", mcl_outfile = "some/folder/hit-frac0.4_I1.4_e5", minimal_taxa = 5, outdir = "some/folder")}
#' @export
write_fasta_files_from_mcl <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), all_fasta, mcl_outfile, minimal_taxa = 4, outdir) {
  path_to_ys <- jntools::add_slash(path_to_ys)
  arguments <- c(paste0(path_to_ys, "write_fasta_files_from_mcl.py"), all_fasta, mcl_outfile, minimal_taxa, outdir)
  system2("python", arguments)
}

#' blast_to_mcl
#'
#'Wrapper for Yang and Smith (2014) blast_to_mcl.py
#'
#'Converts the output of an all-by-all blast query into a format that can be parsed by mcl to find clusters.
#'
#'@param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., "/Users/me/apps/phylogenomic_dataset_construction/"
#'@param blast_results Character vector of length one; the path to the tab-separated text file containing the results from an all-by-all blast search. If blast searches were run separately (i.e., one for each sample), the results should be concatenated into a single file. For the blast search, the output format should specified as: -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore'
#'@param hit_fraction_cutoff Numeric between 0 and 1. Indicates the minimum percentage overlap between query and target in blast results to be retained in the output. According to Y&S, "A low hit-fraction cutoff will output clusters with more incomplete sequences and much larger and sparser alignments, whereas a high hit-fraction cutoff gives tighter clusters but ignores incomplete or divergent sequences."
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{blast_to_mcl(blast_results = "some/folder/blastresults.tab", hit_fraction_cutoff = 0.5)}
#' @export
blast_to_mcl <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), blast_results, hit_fraction_cutoff, ...) {
  path_to_ys <- jntools::add_slash(path_to_ys)
  arguments <- c(paste0(path_to_ys, "blast_to_mcl.py"), blast_results, hit_fraction_cutoff)
  system2("python", arguments)
}
