#' Infer trees from fasta files.
#'
#' Given a folder containing unaligned sequences in fasta format (i.e., clusters),
#' aligns each cluster with \code{mafft} (small clusters) or \code{pasta} (large clusters),
#' excludes poorly aligned sites with \code{phyutility}, and infers a maximum-likelihood
#' tree with \code{RAxML} (small clusters) or \code{fasttree} (large clusters). Requires all
#' of these programs to be installed and included in the user's \code{$PATH}. Assumes clusters are named \code{"cluster1.fa"}, \code{"cluster2.fa"}, etc. Clusters with fewer than 1,000 sequences are considered "small," and those with more are considered "large."
#'
#' Wrapper for Yang and Smith (2014) \code{fasta_to_tree.py}
#'
#' @param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., \code{"/Users/me/apps/phylogenomic_dataset_construction/"}
#' @param seq_folder Character vector of length one; the path to the folder containing the fasta files.
#' @param number_cores Numeric; number of threads to use for \code{RAxML} and \code{mafft}.
#' @param seq_type Character vector of length one indicating type of sequences. Should either be \code{"dna"} for DNA or \code{"aa"} for proteins.
#' @param bootstrap Logical; should run a bootstrap analysis be run for the trees?
#' @param overwrite Logical; should previous output of this command be erased so new output can be written? Once erased it cannot be restored, so use with caution!
#' @param get_hash Logical; should the 32-byte MD5 hash be computed for all output tree files concatenated together? Used for by \code{\link[drake]{drake_plan}} for tracking during workflows. If \code{TRUE}, this function will return the hash.
#' @param echo Logical; should the standard output and error be printed to the screen?
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.
#' @return For each input cluster \code{cluster1.fa} in \code{seq_folder}, \code{cluster1.fa.mafft.aln} (small clusters) or \code{cluster1.pasta.aln} (large clusters), \code{cluster1.fa.mafft.aln-cln} (small clusters) or \code{cluster1.fa.pasta.aln-cln} (large clusters), and \code{cluster1.raxml.tre} (small clusters) or \code{cluster1.fasttree.tre} (large clusters) will be written to \code{seq_folder}. If \code{get_hash} is \code{TRUE}, the 32-byte MD5 hash be computed for all \code{.tre} files concatenated together will be returned.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{fasta_to_tree(seq_folder = "some/folder/containing/fasta/seqs", number_cores = 1, seq_type = "dna", bootstrap = FALSE)}
#' @export
fasta_to_tree <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), seq_folder, number_cores, seq_type, bootstrap = FALSE, overwrite = FALSE, get_hash = TRUE, echo = pkgconfig::get_config("baitfindR::echo", fallback = FALSE), ...) {

  # error checking
  if(is.null(path_to_ys)) {
    stop("Must provide 'path_to_ys' (path to Yang & Smith Phylogenomic Dataset Analysis folder)")
  }

  # modify arguments
  bootstrap <- ifelse(bootstrap, "y", "n")
  path_to_ys <- fs::path_abs(path_to_ys)
  seq_folder <- fs::path_abs(seq_folder)

  # define search terms for output files
  search_terms <- paste("cluster.*\\.fa\\.mafft\\.aln$",
                        "cluster.*\\.fa\\.mafft\\.aln-cln$",
                        "cluster.*\\.fa\\.mafft\\.aln-cln.reduced$",
                        "cluster.*\\.raxml\\.tre$",
                        sep = "|")

  # optional: delete all previous output written in this folder
  if (isTRUE(overwrite)) {
    files_to_delete <- list.files(seq_folder, pattern = search_terms)
    if (length(files_to_delete) > 0) {
      files_to_delete <- paste0(seq_folder, files_to_delete)
      file.remove(files_to_delete)
    }
  }

  # call command
  arguments <- c(fs::path(path_to_ys, "fasta_to_tree.py"), seq_folder, number_cores, seq_type, bootstrap)
  processx::run("python", arguments, wd = seq_folder, echo = echo)

  # optional: get MD5 hash of output
  if (isTRUE(get_hash)) {
    output <- list.files(seq_folder, pattern = search_terms, full.names = TRUE)
    output <- if (length(output) > 0) {unlist( lapply(output, readr::read_file) )} else {output}
    hash <- digest::digest(output)
    return(hash)
  }
}

#' Write fasta files from MCL results.
#'
#' Given the output from the mcl clustering algorthim and a concatenated fasta file
#' including all sequences used for clustering, outputs one fasta file per cluster
#' including the sequences in that cluster.
#'
#' Wrapper for Yang and Smith (2014) write_fasta_files_from_mcl.py
#'
#' @param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., "/Users/me/apps/phylogenomic_dataset_construction/"
#' @param all_fasta Character vector of length one; the path to the fasta file including all query sequences concatenated together, i.e., the fasta file used to create the "all-by-all" blast database.
#' @param mcl_outfile Character vector of length one; the path to the output from running mcl on blast distances.
#' @param minimal_taxa Numeric; minimal number of taxa required to be present for the cluster to be written. Default 4, the minimum number of taxa needed for an un-rooted tree.
#' @param outdir Character vector of length one; the path to the folder where the clusters should be written.
#' @param overwrite Logical; should previous output of this command be erased so new output can be written? Once erased it cannot be restored, so use with caution!
#' @param get_hash Logical; should the 32-byte MD5 hash be computed for all clusters concatenated together? Used for by \code{\link[drake]{drake_plan}} for tracking during workflows. If \code{TRUE}, this function will return the hash.
#' @param echo Logical; should the standard output and error be printed to the screen?
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.
#' @return One fasta file per cluster (\code{cluster1.fa}, \code{cluster2.fa}, etc.) will be written to \code{outdir}. If \code{get_hash} is \code{TRUE}, the 32-byte MD5 hash be computed for all \code{.fa} files concatenated together will be returned.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{write_fasta_files_from_mcl(all_fasta = "some/folder/all.fasta", mcl_outfile = "some/folder/hit-frac0.4_I1.4_e5", minimal_taxa = 5, outdir = "some/folder")}
#' @export
write_fasta_files_from_mcl <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), all_fasta, mcl_outfile, minimal_taxa = 4, outdir, overwrite = FALSE, get_hash = TRUE, echo = pkgconfig::get_config("baitfindR::echo", fallback = FALSE), ...) {

  # error checking
  if(is.null(path_to_ys)) {
    stop("Must provide 'path_to_ys' (path to Yang & Smith Phylogenomic Dataset Analysis folder)")
  }

  # modify arguments
  path_to_ys <- fs::path_abs(path_to_ys)
  outdir <- fs::path_abs(outdir)

  # define search terms for output files
  search_terms <- "cluster\\d*\\.fa$"

  # optional: delete all previous output written in this folder
  if (isTRUE(overwrite)) {
    files_to_delete <- list.files(outdir, pattern = search_terms, full.names = TRUE)
    if (length(files_to_delete) > 0) {
      file.remove(files_to_delete)
    }
  }

  # call command
  arguments <- c(fs::path(path_to_ys, "write_fasta_files_from_mcl.py"), all_fasta, mcl_outfile, minimal_taxa, outdir)
  processx::run("python", arguments, echo = echo)

  # optional: get MD5 hash of output
  if (isTRUE(get_hash)) {
    output <- list.files(outdir, pattern = search_terms, full.names = TRUE)
    output <- if (length(output) > 0) {unlist( lapply(output, readr::read_file) )} else {output}
    hash <- digest::digest(output)
    return(hash)
  }
}

#' Prepare BLAST results for MCL.
#'
#' Converts the output of an all-by-all blast query into a format that can be parsed by mcl to find clusters.
#'
#' Wrapper for Yang and Smith (2014) blast_to_mcl.py
#'
#' @param path_to_ys Character vector of length one; the complete path to the folder containing Y&S python scripts, e.g., "/Users/me/apps/phylogenomic_dataset_construction/"
#' @param blast_results Character vector of length one; the complete path to the tab-separated text file containing the results from an all-by-all blast search. If blast searches were run separately (i.e., one for each sample), the results should be concatenated into a single file. For the blast search, the output format should specified as: -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore'
#' @param hit_fraction_cutoff Numeric between 0 and 1. Indicates the minimum percentage overlap between query and target in blast results to be retained in the output. According to Y&S, "A low hit-fraction cutoff will output clusters with more incomplete sequences and much larger and sparser alignments, whereas a high hit-fraction cutoff gives tighter clusters but ignores incomplete or divergent sequences."
#' @param echo Logical; should the standard output and error be printed to the screen?
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.
#' @return A tab-separated text file with three columns: the first two are the matching query and target from the all-by-all blast, and the third is the negative log e-value for that match. This file is named \code{<blast_results>.hit-frac<hit_fraction_cutoff>.minusLogEvalue}, where \code{<blast_results>} and \code{<hit_fraction_cutoff>} correspond to the values of those arguments. If possible contaminants (i.e., identical sequences between different samples) were found, these are written to \code{<blast_results>.ident.hit-frac<hit_fraction_cutoff>}. Output files will be written to the same folder containing \code{blast_results}.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{blast_to_mcl(blast_results = "some/folder/blastresults.tab", hit_fraction_cutoff = 0.5)}
#' @export
blast_to_mcl <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), blast_results, hit_fraction_cutoff, echo = pkgconfig::get_config("baitfindR::echo", fallback = FALSE), ...) {

  # modify arguments
  path_to_ys <- fs::path_abs(path_to_ys)
  arguments <- c(fs::path(path_to_ys, "blast_to_mcl.py"), blast_results, hit_fraction_cutoff)

  # call command
  processx::run("python", arguments, echo = echo)

  # Normally, the warning file with sequences that are identical between
  # samples (possible contamination) is output as "blast_output.ident",
  # but this results in files with the same name from different hit_faction_cutoff
  # values. Append the hit_fraction_cutoff value so we can tell them apart.

  blast_results <- fs::path_abs(blast_results)
  ident_file <- fs::path(blast_results, ext = "ident")
  ident_file_rename <- fs::path(ident_file, ext = glue::glue("hit-frac.{hit_fraction_cutoff}"))
  if (file.exists(ident_file)) {
    file.rename(ident_file, ident_file_rename)
  }
}

#' Shorten names in fasta headers.
#'
#' Does the same thing as Yang and Smith (2014) fix_names_from_transdecoder.py, but works on one fasta file at at time.
#'
#' @param transdecoder_output Character vector of length one; the path to the .transdecoder.cds file produced by \code{\link{transdecoder_predict}}. It is assumed that the first part of the filename (immediately preceding .transdecoder.cds) is the sample code.
#' @param mol_type Character vector of length one; "dna" for DNA or "aa" for proteins.
#'
#' @return Object of class \code{DNAbin} or \code{AAbin} with names shortened to \code{sample_code`"@"`gene}
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{fix_names_from_transdecoder(transdecoder_output = "some/folder/CODE.transdecoder.cds")}
#' @export
fix_names_from_transdecoder <- function (transdecoder_output, mol_type = "dna") {

  transdecoder_output <- fs::path_abs(transdecoder_output)

  if (!(mol_type == "dna" | mol_type == "aa")) {
    stop("Must choose 'dna' or 'aa' for mol_type")
  } else if (mol_type == "dna") {
    seq <- ape::read.FASTA (transdecoder_output)
  } else if (mol_type == "aa") {
    seq <- ape::read.FASTA (transdecoder_output, type = "AA")
  }

  # get sequence names from fasta file
  seq_names <- names(seq)

  # get sample code from input string
  # break up by slashes, and grab last element
  code <- stringr::str_split(transdecoder_output, "/", simplify = TRUE)
  code <- code[,length(code)]

  # check to make sure input file name format is correct
  if (grepl(".transdecoder.cds", code) == FALSE) {
    stop("transdecoder_output does not end in '.transdecoder.cds' ")
  }

  # isolate sample code
  code <- gsub(".transdecoder.cds", "", code)

  # recreate Y&S code to isolate gene number:
  # 1. grab first substring separated by spaces,
  # 2. then grab the last substring of that separated by periods
  newid <- purrr::map_chr(seq_names, function (x) stringr::str_split(x, pattern = " ")[[1]][[1]] )
  newid <- purrr::map_chr(newid, function (x) stringr::str_split(x, pattern = "\\.")[[1]][[ length(stringr::str_split(x, pattern = "\\.")[[1]]) ]] )
  new_names <- paste(code, newid, sep="@")

  names(seq) <- new_names

  return(seq)

}

#' Trim tips.
#'
#' Given a folder containing phylogenetic trees, exclude (i.e., "trim"),
#' tips on unusually long branches. Tips on a branch 10 times longer
#' than their sister AND longer than \code{relative_cutoff}, OR tips
#' that are longer than \code{absolute_cutoff} will be trimmed. This
#' function will overwrite any output files with the same name in
#' \code{tree_folder}.
#'
#' Wrapper for Yang and Smith (2014) \code{trim_tips.py}
#'
#' @param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., \code{"/Users/me/apps/phylogenomic_dataset_construction/"}
#' @param tree_folder Character vector of length one; the path to the folder containing the trees to trim.
#' @param tree_file_ending Character vector of length one; only tree files with this file ending will be used.
#' @param relative_cutoff Numeric vector of length one; tips on a branch 10 times longer than their sister AND longer than this value will be trimmed.
#' @param absolute_cutoff Numeric vector of length one; tips on branches longer than this value will be trimmed.
#' @param overwrite Logical; should previous output of this command be erased so new output can be written? Once erased it cannot be restored, so use with caution!
#' @param get_hash Logical; should the 32-byte MD5 hash be computed for all output trimmed tree files concatenated together? Used for by \code{\link[drake]{drake_plan}} for tracking during workflows. If \code{TRUE}, this function will return the hash.
#' @param echo Logical; should the standard output and error be printed to the screen?
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.
#' @return For each input tree with a file ending matching \code{tree_file_ending} in \code{tree_folder}, a trimmed tree with a file ending in \code{.tt} will be written to \code{tree_folder}. If \code{get_hash} is \code{TRUE}, the 32-byte MD5 hash be computed for all trimmed tree files concatenated together will be returned.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{trim_tips(tree_folder = "some/folder/containing/tree/files", tree_file_ending = ".tre", relative_cutoff = 0.2, absolute_cutoff = 0.4)}
#' @export
trim_tips <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), tree_folder, tree_file_ending, relative_cutoff, absolute_cutoff, overwrite = FALSE, get_hash = TRUE, echo = pkgconfig::get_config("baitfindR::echo", fallback = FALSE), ...) {

  # error checking
  if(is.null(path_to_ys)) {
    stop("Must provide 'path_to_ys' (path to Yang & Smith Phylogenomic Dataset Analysis folder)")
  }

  # modify arguments
  path_to_ys <- fs::path_abs(path_to_ys)
  tree_folder <- fs::path_abs(tree_folder)

  # define search terms for output files
  search_terms <- "\\.tt$"

  # optional: delete all previous output written in this folder
  if (isTRUE(overwrite)) {
    files_to_delete <- list.files(tree_folder, pattern = search_terms)
    if (length(files_to_delete) > 0) {
      files_to_delete <- paste0(tree_folder, files_to_delete)
      file.remove(files_to_delete)
    }
  }

  # call command
  arguments <- c(fs::path(path_to_ys, "trim_tips.py"), tree_folder, tree_file_ending, relative_cutoff, absolute_cutoff)
  processx::run("python", arguments, echo = echo)

  # optional: get MD5 hash of output
  if (isTRUE(get_hash)) {
    output <- list.files(tree_folder, pattern = search_terms)
    output <- if (length(output) > 0) {unlist( lapply(paste0(tree_folder, output), readr::read_file) )} else {output}
    hash <- digest::digest(output)
    return(hash)
  }
}

#' Mask tips in tree.
#'
#' Given a folder containing phylogenetic trees and their alignments, mask monophyletic
#' (and optionally, paraphyletic) tips belonging to the same taxon (i.e., keep only a single tip
#' to represent clades consisting of a single taxon). Tree files are assumed to end in \code{.tt}
#' (the output of \code{trim_tips}), and only tree files with this ending will be included.
#' Alignment files are assumed to end in \code{.aln-cln} (the output of \code{fasta_to_tree}),
#' and only alignment files with this ending will be included.
#' The tip with the fewest ambiguous characters in the alignment will be kept. This
#' function will overwrite any output files with the same name in \code{tree_folder}.
#'
#' Wrapper for Yang and Smith (2014) \code{mask_tips_by_taxonID_transcripts.py}
#'
#' @param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., \code{"/Users/me/apps/phylogenomic_dataset_construction/"}
#' @param tree_folder Character vector of length one; the path to the folder containing the trees to mask.
#' @param aln_folder Character vector of length one; the path to the folder containing the alignments used to make the trees.
#' @param mask_paraphyletic Logical; should paraphyletic tips belonging to the same taxon be masked?
#' @param overwrite Logical; should previous output of this command be erased so new output can be written? Once erased it cannot be restored, so use with caution!
#' @param get_hash Logical; should the 32-byte MD5 hash be computed for all output masked tree files concatenated together? Used for by \code{\link[drake]{drake_plan}} for tracking during workflows. If \code{TRUE}, this function will return the hash.
#' @param echo Logical; should the standard output and error be printed to the screen?
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.
#' @return For each input tree with a file ending in \code{.tt} in \code{tree_folder}, a trimmed tree with a file ending in \code{.mm} will be written to \code{tree_folder}. If \code{get_hash} is \code{TRUE}, the 32-byte MD5 hash be computed for all masked tree files concatenated together will be returned.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{mask_tips_by_taxonID_transcripts(tree_folder = "some/folder/containing/tree/files", aln_folder = "some/folder/containing/alignment/files")}
#' @export
mask_tips_by_taxonID_transcripts <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), tree_folder, aln_folder, mask_paraphyletic = TRUE, overwrite = FALSE, get_hash = TRUE, echo = pkgconfig::get_config("baitfindR::echo", fallback = FALSE), ...) {

  # error checking
  if(is.null(path_to_ys)) {
    stop("Must provide 'path_to_ys' (path to Yang & Smith Phylogenomic Dataset Analysis folder)")
  }

  # modify arguments
  path_to_ys <- fs::path_abs(path_to_ys)
  tree_folder <- fs::path_abs(tree_folder)
  aln_folder <- fs::path_abs(aln_folder)
  mask_paraphyletic <- ifelse(isTRUE(mask_paraphyletic), "y", "n")

  # define search terms for output files
  search_terms <- "\\.mm$"

  # optional: delete all previous output written in this folder
  if (isTRUE(overwrite)) {
    files_to_delete <- list.files(tree_folder, pattern = search_terms)
    if (length(files_to_delete) > 0) {
      files_to_delete <- paste0(tree_folder, files_to_delete)
      file.remove(files_to_delete)
    }
  }

  # call command
  arguments <- c(fs::path(path_to_ys, "mask_tips_by_taxonID_transcripts.py"), tree_folder, aln_folder, mask_paraphyletic)
  processx::run("python", arguments, echo = echo)

  # optional: get MD5 hash of output
  if (isTRUE(get_hash)) {
    output <- list.files(tree_folder, pattern = search_terms)
    output <- if (length(output) > 0) {unlist( lapply(paste0(tree_folder, output), readr::read_file) )} else {output}
    hash <- digest::digest(output)
    return(hash)
  }
}

#' Cut long internal branches in tree.
#'
#' Given a folder containing phylogenetic trees, split the trees into multiple subtrees
#' for nodes that bifurcate deeper than \code{internal_branch_length_cutoff}.
#' \code{tree_folder} and \code{outdir} should be different to avoid writing over input trees.
#' This function will overwrite any output files with the same name in \code{outdir}.
#'
#' Wrapper for Yang and Smith (2014) \code{cut_long_internal_branches.py}
#'
#' @param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., \code{"/Users/me/apps/phylogenomic_dataset_construction/"}
#' @param tree_folder Character vector of length one; the path to the folder containing the trees to cut.
#' @param tree_file_ending Character vector of length one; only tree files with this file ending will be used.
#' @param internal_branch_length_cutoff Numeric vector of length one; the depth at which cuts should be made (smaller numbers indicate greater depth).
#' @param minimal_taxa Numeric; minimal number of taxa required for tree to be cut. Default 4, the minimum number of taxa needed for an un-rooted tree.
#' @param outdir Character vector of length one; the path to the folder where the subtrees should be written.
#' @param overwrite Logical; should previous output of this command be erased so new output can be written? Once erased it cannot be restored, so use with caution!
#' @param get_hash Logical; should the 32-byte MD5 hash be computed for all output subtree files concatenated together? Used for by \code{\link[drake]{drake_plan}} for tracking during workflows. If \code{TRUE}, this function will return the hash.
#' @param echo Logical; should the standard output and error be printed to the screen?
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.
#' @return For each input tree with a file ending in \code{tree_file_ending} in \code{tree_folder}, one or more subtrees with a file ending in \code{.subtree} will be written to \code{tree_folder}. If \code{get_hash} is \code{TRUE}, the 32-byte MD5 hash be computed for all subtree files concatenated together will be returned.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{cut_long_internal_branches(tree_folder = "some/folder/containing/tree/files", tree_file_ending = ".mm", internal_branch_length_cutoff = 0.3, outdir = "some/other/folder/")}
#' @export
cut_long_internal_branches <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), tree_folder, tree_file_ending, internal_branch_length_cutoff, minimal_taxa = 4, outdir, overwrite = FALSE, get_hash = TRUE, echo = pkgconfig::get_config("baitfindR::echo", fallback = FALSE), ...) {

  # error checking
  if(is.null(path_to_ys)) {
    stop("Must provide 'path_to_ys' (path to Yang & Smith Phylogenomic Dataset Analysis folder)")
  }

  # modify arguments
  path_to_ys <- fs::path_abs(path_to_ys)
  tree_folder <- fs::path_abs(tree_folder)
  outdir <- fs::path_abs(outdir)

  # more error checking
  if(tree_folder == outdir) {
    stop("Must provide provide different paths for input and output folders")
  }

  # define search terms for output files
  search_terms <- "\\.subtree$"

  # optional: delete all previous output written in this folder
  if (isTRUE(overwrite)) {
    files_to_delete <- list.files(outdir, pattern = search_terms)
    if (length(files_to_delete) > 0) {
      files_to_delete <- paste0(outdir, files_to_delete)
      file.remove(files_to_delete)
    }
  }

  # call command
  arguments <- c(fs::path(path_to_ys, "cut_long_internal_branches.py"), tree_folder, tree_file_ending, internal_branch_length_cutoff, minimal_taxa, outdir)
  processx::run("python", arguments, echo = echo)

  # optional: get MD5 hash of output
  if (isTRUE(get_hash)) {
    output <- list.files(outdir, pattern = search_terms)
    output <- if (length(output) > 0) {unlist( lapply(paste0(outdir, output), readr::read_file) )} else {output}
    hash <- digest::digest(output)
    return(hash)
  }
}

#' Write fasta files from trees.
#'
#' Given a folder containing phylogenetic trees and a single concatenated fasta file
#' including all the sequences used to build the trees, output one fasta file per tree
#' with the sequences in that tree. This function will overwrite any output files with
#' the same name in \code{outdir}.
#'
#' Wrapper for Yang and Smith (2014) write_fasta_files_from_trees.py
#'
#' @param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., "/Users/me/apps/phylogenomic_dataset_construction/"
#' @param all_fasta Character vector of length one; the path to the fasta file including all the sequences that were originally used to build the trees.
#' @param tree_folder Character vector of length one; the path to the folder containing the trees to be used for extracting fasta sequences.
#' @param tree_file_ending Character vector of length one; only tree files with this file ending will be used.
#' @param outdir Character vector of length one; the path to the folder where the fasta files should be written.
#' @param overwrite Logical; should previous output of this command be erased so new output can be written? Once erased it cannot be restored, so use with caution!
#' @param get_hash Logical; should the 32-byte MD5 hash be computed for all output fasta files concatenated together? Used for by \code{\link[drake]{drake_plan}} for tracking during workflows. If \code{TRUE}, this function will return the hash.
#' @param echo Logical; should the standard output and error be printed to the screen?
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.
#' @return One fasta file per tree file ending in \code{tree_file_ending} in \code{tree_folder} will be written to \code{outdir}. If \code{get_hash} is \code{TRUE}, the 32-byte MD5 hash be computed for all output fasta files concatenated together will be returned.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{write_fasta_files_from_trees(all_fasta = "some/folder/all.fasta", tree_file_ending = ".subtree", tree_folder = "some/folder/containing/tree/files", outdir = "some/folder")}
#' @export
write_fasta_files_from_trees <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), all_fasta, tree_folder, tree_file_ending, outdir, overwrite = FALSE, get_hash = TRUE, echo = pkgconfig::get_config("baitfindR::echo", fallback = FALSE), ...) {

  # error checking
  if(is.null(path_to_ys)) {
    stop("Must provide 'path_to_ys' (path to Yang & Smith Phylogenomic Dataset Analysis folder)")
  }

  # modify arguments
  path_to_ys <- fs::path_abs(path_to_ys)
  outdir <- fs::path_abs(outdir)
  tree_folder <- fs::path_abs(tree_folder)

  # define search terms for output files
  search_terms <- "rr\\.fa$"

  # optional: delete all previous output written in this folder
  if (isTRUE(overwrite)) {
    files_to_delete <- list.files(outdir, pattern = search_terms)
    if (length(files_to_delete) > 0) {
      files_to_delete <- paste0(outdir, files_to_delete)
      file.remove(files_to_delete)
    }
  }

  # call command
  arguments <- c(fs::path(path_to_ys, "write_fasta_files_from_trees.py"), all_fasta, tree_folder, tree_file_ending, outdir)
  processx::run("python", arguments, echo = echo)

  # optional: get MD5 hash of output
  if (isTRUE(get_hash)) {
    output <- list.files(outdir, pattern = search_terms)
    output <- if (length(output) > 0) {unlist( lapply(paste0(outdir, output), readr::read_file) )} else {output}
    hash <- digest::digest(output)
    return(hash)
  }
}

#' Identify orthologs using the "one-to-one" method.
#'
#' Given a folder containing homolog trees, filter the trees to only those
#' containing one-to-one orthologs (i.e., no duplications within a sample).
#' This function will overwrite any output files with the same name in
#' \code{outdir}.
#'
#' Wrapper for Yang and Smith (2014) filter_1to1_orthologs.py
#'
#' @param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., "/Users/me/apps/phylogenomic_dataset_construction/"
#' @param tree_folder Character vector of length one; the path to the folder containing the trees to be used for filtering.
#' @param tree_file_ending Character vector of length one; only tree files with this file ending will be used.
#' @param minimal_taxa Numeric; minimal number of taxa required for tree to be included. Default 4, the minimum number of taxa needed for an un-rooted tree.
#' @param outdir Character vector of length one; the path to the folder where the filtered trees should be written.
#' @param overwrite Logical; should previous output of this command be erased so new output can be written? Once erased it cannot be restored, so use with caution!
#' @param get_hash Logical; should the 32-byte MD5 hash be computed for all filtered tree files concatenated together? Used for by \code{\link[drake]{drake_plan}} for tracking during workflows. If \code{TRUE}, this function will return the hash.
#' @param echo Logical; should the standard output and error be printed to the screen?
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.
#' @return For each tree file ending in \code{tree_file_ending} in \code{tree_folder}, that tree will be written to \code{outdir} if it consists solely of one-to-one orthologs with the file ending \code{.1to1ortho.tre}. If \code{get_hash} is \code{TRUE}, the 32-byte MD5 hash be computed for all filtered tree files concatenated together will be returned.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{filter_1to1_orthologs(tree_folder = "some/folder/containing/tree/files", tree_file_ending = ".tre", tree_folder = "some/folder/containing/tree/files", outdir = "some/folder")}
#' @export
filter_1to1_orthologs <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), tree_folder, tree_file_ending, minimal_taxa = 4, outdir, overwrite = FALSE, get_hash = TRUE, echo = pkgconfig::get_config("baitfindR::echo", fallback = FALSE), ...) {

  # error checking
  if(is.null(path_to_ys)) {
    stop("Must provide 'path_to_ys' (path to Yang & Smith Phylogenomic Dataset Analysis folder)")
  }

  # modify arguments
  path_to_ys <- fs::path_abs(path_to_ys)
  outdir <- fs::path_abs(outdir)
  tree_folder <- fs::path_abs(tree_folder)

  # define search terms for output files
  search_terms <- "1to1ortho\\.tre$"

  # optional: delete all previous output written in this folder
  if (isTRUE(overwrite)) {
    files_to_delete <- list.files(outdir, pattern = search_terms)
    if (length(files_to_delete) > 0) {
      files_to_delete <- paste0(outdir, files_to_delete)
      file.remove(files_to_delete)
    }
  }

  # call command
  arguments <- c(fs::path(path_to_ys, "filter_1to1_orthologs.py"), tree_folder, tree_file_ending, minimal_taxa, outdir)
  processx::run("python", arguments, echo = echo)

  # optional: get MD5 hash of output
  if (isTRUE(get_hash)) {
    output <- list.files(outdir, pattern = search_terms)
    output <- if (length(output) > 0) {unlist( lapply(paste0(outdir, output), readr::read_file) )} else {output}
    hash <- digest::digest(output)
    return(hash)
  }
}

#' Identify orthologs using the "MI" method.
#'
#' Given a folder containing homolog trees, prune paralogs from the trees
#' using the maximum inclusion (MI) method. This method iteratively
#' extracts subtrees containing the greatest number of samples (taxa)
#' without duplication. Long branches will also be trimmed (i.e., removed)
#' according to \code{relative_cutoff} and \code{absolute_cutoff} values. Trees
#' that consist solely of one-to-one orthologs (i.e., no duplications within a sample)
#' will also be retained. This function will overwrite any output files with
#' the same name in \code{outdir}.
#'
#' Wrapper for Yang and Smith (2014) prune_paralogs_MI.py
#'
#' @param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., "/Users/me/apps/phylogenomic_dataset_construction/"
#' @param tree_folder Character vector of length one; the path to the folder containing the trees to be used for pruning.
#' @param tree_file_ending Character vector of length one; only tree files with this file ending will be used.
#' @param relative_cutoff Numeric vector of length one; tips on a branch 10 times longer than their sister AND longer than this value will be trimmed.
#' @param absolute_cutoff Numeric vector of length one; tips on branches longer than this value will be trimmed.
#' @param minimal_taxa Numeric; minimal number of taxa required for tree to be included. Default 4, the minimum number of taxa needed for an un-rooted tree.
#' @param outdir Character vector of length one; the path to the folder where the pruned trees should be written.
#' @param overwrite Logical; should previous output of this command be erased so new output can be written? Once erased it cannot be restored, so use with caution!
#' @param get_hash Logical; should the 32-byte MD5 hash be computed for all pruned tree files concatenated together? Used for by \code{\link[drake]{drake_plan}} for tracking during workflows. If \code{TRUE}, this function will return the hash.
#' @param echo Logical; should the standard output and error be printed to the screen?
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.
#' @return For each tree file ending in \code{tree_file_ending} in \code{tree_folder}, putative orthologs will be extracted from the tree using the MI method and written to \code{outdir} with the file ending \code{.MIortho1.tre}. If \code{get_hash} is \code{TRUE}, the 32-byte MD5 hash be computed for all extracted tree files concatenated together will be returned.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{prune_paralogs_MI(tree_folder = "some/folder/containing/tree/files", tree_file_ending = ".tre", relative_cutoff = 0.2, absolute_cutoff = 0.4, outdir = "some/folder")}
#' @export
prune_paralogs_MI <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), tree_folder, tree_file_ending, relative_cutoff, absolute_cutoff, minimal_taxa = 4, outdir, overwrite = FALSE, get_hash = TRUE, echo = pkgconfig::get_config("baitfindR::echo", fallback = FALSE), ...) {

  # error checking
  if(is.null(path_to_ys)) {
    stop("Must provide 'path_to_ys' (path to Yang & Smith Phylogenomic Dataset Analysis folder)")
  }

  # modify arguments
  path_to_ys <- fs::path_abs(path_to_ys)
  outdir <- fs::path_abs(outdir)
  tree_folder <- fs::path_abs(tree_folder)

  # define search terms for output files
  search_terms <- "1to1ortho\\.tre$|MIortho.*\\.tre$"

  # optional: delete all previous output written in this folder
  if (isTRUE(overwrite)) {
    files_to_delete <- list.files(outdir, pattern = search_terms)
    if (length(files_to_delete) > 0) {
      files_to_delete <- paste0(outdir, files_to_delete)
      file.remove(files_to_delete)
    }
  }

  # call command
  arguments <- c(fs::path(path_to_ys, "prune_paralogs_MI.py"), tree_folder, tree_file_ending, relative_cutoff, absolute_cutoff, minimal_taxa, outdir)
  processx::run("python", arguments, echo = echo)

  # optional: get MD5 hash of output
  if (isTRUE(get_hash)) {
    output <- list.files(outdir, pattern = search_terms)
    output <- if (length(output) > 0) {unlist( lapply(paste0(outdir, output), readr::read_file) )} else {output}
    hash <- digest::digest(output)
    return(hash)
  }
}

#' Identify orthologs using the "MO" method.
#'
#' Given a folder containing homolog trees, prune paralogs from the trees
#' using the monophyletic outgroups (MO) method. For trees that have non-repeating, monophyletic
#' outgroups, this method extracts the largest subtree containing the outgroup. Trees
#' that consist solely of one-to-one orthologs (i.e., no duplications within a sample/taxon)
#' will also be retained. This function will overwrite any output files with the same
#' name in \code{outdir}.
#'
#' Wrapper for Yang and Smith (2014) prune_paralogs_MO.py
#'
#' @param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., "/Users/me/apps/phylogenomic_dataset_construction/"
#' @param tree_folder Character vector of length one; the path to the folder containing the trees to be used for pruning.
#' @param tree_file_ending Character vector of length one; only tree files with this file ending will be used.
#' @param ingroup Character vector; names of ingroup taxa/samples.
#' @param outgroup Character vector; names of outgroup taxa/samples.
#' @param minimal_taxa Numeric; minimal number of taxa required for tree to be included. Default 4, the minimum number of taxa needed for an un-rooted tree.
#' @param outdir Character vector of length one; the path to the folder where the pruned trees should be written.
#' @param overwrite Logical; should previous output of this command be erased so new output can be written? Once erased it cannot be restored, so use with caution!
#' @param get_hash Logical; should the 32-byte MD5 hash be computed for all pruned tree files concatenated together? Used for by \code{\link[drake]{drake_plan}} for tracking during workflows. If \code{TRUE}, this function will return the hash.
#' @param echo Logical; should the standard output and error be printed to the screen?
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.
#' @return For each tree file ending in \code{tree_file_ending} in \code{tree_folder}, putative orthologs will be extracted from the tree using the MO method and written to \code{outdir} with the file ending \code{.ortho.tre}; re-rooted trees will also be written with the file ending \code{.reroot}. If \code{get_hash} is \code{TRUE}, the 32-byte MD5 hash be computed for all extracted tree files concatenated together will be returned.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{prune_paralogs_MO(tree_folder = "some/folder/containing/tree/files", tree_file_ending = ".tre", outgroup = c("ABC", "EFG"), ingroup = c("HIJ", "KLM"), outdir = "some/folder")}
#' @export
prune_paralogs_MO <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), tree_folder, tree_file_ending, ingroup, outgroup, minimal_taxa = 4, outdir, overwrite = FALSE, get_hash = TRUE, echo = pkgconfig::get_config("baitfindR::echo", fallback = FALSE), ...) {

  # error checking
  if(is.null(path_to_ys)) {
    stop("Must provide 'path_to_ys' (path to Yang & Smith Phylogenomic Dataset Analysis folder)")
  }

  # modify arguments
  path_to_ys <- fs::path_abs(path_to_ys)
  outdir <- fs::path_abs(outdir)
  tree_folder <- fs::path_abs(tree_folder)

  # define search terms for output files
  search_terms <- "1to1ortho\\.tre$|\\.reroot$|\\.ortho\\.tre$"

  # optional: delete all previous output written in this folder
  if (isTRUE(overwrite)) {
    files_to_delete <- list.files(outdir, pattern = search_terms)
    if (length(files_to_delete) > 0) {
      files_to_delete <- paste0(outdir, files_to_delete)
      file.remove(files_to_delete)
    }
  }

  # modify actual script to change ingroups and outgroups
  ys_script <- readr::read_lines(fs::path(path_to_ys, "prune_paralogs_MO.py"))

  # look for chunks of comments and single comment lines, delete these
  # (some comment chunks contained old outgroup / ingroup assignments)
  comment_block_lines <- grep('\"\"\"', ys_script)
  comment_block_start <- comment_block_lines[seq(1,length(comment_block_lines),by=2)]
  comment_block_end <- comment_block_lines[seq(2,length(comment_block_lines),by=2)]
  comment_blocks <- unlist(purrr::map2(comment_block_start, comment_block_end, ~ seq(.x, .y)))

  single_comment_lines <- grep("^#", ys_script )

  comments <- unique(c(comment_block_start, comment_block_end, comment_blocks, single_comment_lines))

  ys_script <- ys_script[-comments]

  # find lines specifying outgroups and ingroups
  outgroup_line <- grep("^OUTGROUPS = ", ys_script)
  ingroup_line <- grep("^INGROUPS = ", ys_script)

  # error if outgroup and ingroup lines can't be found or aren't unique
  if (length(outgroup_line) != 1) {stop ("Unable to assign single outgroup; check prune_paralogs_MO.py")}
  if (length(ingroup_line) != 1) {stop ("Unable to assign single ingroup; check prune_paralogs_MO.py")}

  # replace outgroup and ingroup lines
  outgroup_replacement <- paste0(outgroup, collapse= '\",\"')
  outgroup_replacement <- paste0('\"', outgroup_replacement, '\"')
  outgroup_replacement <- paste0("OUTGROUPS = [", outgroup_replacement, "]")

  ingroup_replacement <- paste0(ingroup, collapse= '\",\"')
  ingroup_replacement <- paste0('\"', ingroup_replacement, '\"')
  ingroup_replacement <- paste0("INGROUPS = [", ingroup_replacement, "]")

  ys_script[outgroup_line] <- outgroup_replacement
  ys_script[ingroup_line] <- ingroup_replacement

  # write out new temporary python script
  readr::write_lines(ys_script, path = fs::path(path_to_ys, "prune_paralogs_MO_temp.py"))

  # call command
  arguments <- c(fs::path(path_to_ys, "prune_paralogs_MO_temp.py"), tree_folder, tree_file_ending, minimal_taxa, outdir)
  processx::run("python", arguments, echo = echo)

  # delete temporary script
  file.remove(fs::path(path_to_ys, "prune_paralogs_MO_temp.py"))

  # optional: get MD5 hash of output
  if (isTRUE(get_hash)) {
    output <- list.files(outdir, pattern = search_terms)
    output <- if (length(output) > 0) {unlist( lapply(paste0(outdir, output), readr::read_file) )} else {output}
    hash <- digest::digest(output)
    return(hash)
  }
}

#' Identify orthologs using the "RT" method.
#'
#' Given a folder containing homolog trees, prune paralogs from the trees
#' using the rooted ingroups (RT) method. For trees that have outgroups, this
#' method iteratively extracts subtrees with the highest number of ingroup taxa/samples.
#' This function will overwrite any output files with the same
#' name in \code{outdir}.
#'
#' Wrapper for Yang and Smith (2014) prune_paralogs_RT.py
#'
#' @param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., "/Users/me/apps/phylogenomic_dataset_construction/"
#' @param tree_folder Character vector of length one; the path to the folder containing the trees to be used for pruning.
#' @param tree_file_ending Character vector of length one; only tree files with this file ending will be used.
#' @param ingroup Character vector; names of ingroup taxa/samples.
#' @param outgroup Character vector; names of outgroup taxa/samples.
#' @param min_ingroup_taxa Numeric; minimal number of taxa in the ingroup required for an ortholog to be written. Default 2.
#' @param outdir Character vector of length one; the path to the folder where the pruned trees should be written.
#' @param overwrite Logical; should previous output of this command be erased so new output can be written? Once erased it cannot be restored, so use with caution!
#' @param get_hash Logical; should the 32-byte MD5 hash be computed for all pruned tree files concatenated together? Used for by \code{\link[drake]{drake_plan}} for tracking during workflows. If \code{TRUE}, this function will return the hash.
#' @param echo Logical; should the standard output and error be printed to the screen?
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.
#' @return For each tree file ending in \code{tree_file_ending} in \code{tree_folder}, the following outputs are possible depending on the presence of outgroups in the homolog tree:
#' \describe{
#'   \item{No outgroups in homolog tree}{Unrooted ingroup clades without duplications (files ending in \code{unrooted-ortho.tre})}
#'   \item{Outgroups present in homolog tree}{Rooted ingroup clades (files ending in \code{inclade}) and one or more paralogs (files ending in \code{inclade.ortho.tre})}
#' }
#'
#' If \code{get_hash} is \code{TRUE}, the 32-byte MD5 hash be computed for all extracted tree files concatenated together will be returned.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{prune_paralogs_RT(tree_folder = "some/folder/containing/tree/files", tree_file_ending = ".tre", outgroup = c("ABC", "EFG"), ingroup = c("HIJ", "KLM"), outdir = "some/folder")}
#' @export
prune_paralogs_RT <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), tree_folder, tree_file_ending, ingroup, outgroup, min_ingroup_taxa = 2, outdir, overwrite = FALSE, get_hash = TRUE, echo = pkgconfig::get_config("baitfindR::echo", fallback = FALSE), ...) {

  # error checking
  if(is.null(path_to_ys)) {
    stop("Must provide 'path_to_ys' (path to Yang & Smith Phylogenomic Dataset Analysis folder)")
  }

  # modify arguments
  path_to_ys <- fs::path_abs(path_to_ys)
  outdir <- fs::path_abs(outdir)
  tree_folder <- fs::path_abs(tree_folder)

  # define search terms for output files
  search_terms <- "\\.inclade\\d*|\\.unrooted-ortho\\.tre$"

  # optional: delete all previous output written in this folder
  if (isTRUE(overwrite)) {
    files_to_delete <- list.files(outdir, pattern = search_terms)
    if (length(files_to_delete) > 0) {
      files_to_delete <- paste0(outdir, files_to_delete)
      file.remove(files_to_delete)
    }
  }

  # this script needs outgroups and ingroups provided as a text file
  in_out <- c(paste0("OUT\t", outgroup), paste0("IN\t", ingroup))

  # write out temporary text file with ingroups and outgroups
  readr::write_lines(in_out, path = here::here("in_out_temp"))

  # call command
  arguments <- c(fs::path(path_to_ys, "prune_paralogs_RT.py"), tree_folder, tree_file_ending, outdir, min_ingroup_taxa, here::here("in_out_temp"))
  processx::run("python", arguments, echo = echo)

  # delete temporary in_out file
  file.remove(here::here("in_out_temp"))

  # optional: get MD5 hash of output
  if (isTRUE(get_hash)) {
    output <- list.files(outdir, pattern = search_terms)
    output <- if (length(output) > 0) {unlist( lapply(paste0(outdir, output), readr::read_file) )} else {output}
    hash <- digest::digest(output)
    return(hash)
  }
}

#' Write fasta files from ortholog trees.
#'
#' Given a folder containing ortholog trees, write out the fasta files
#' that correspond to the sequences in the trees.
#'
#' Wrapper for Yang and Smith (2014) write_ortholog_fasta_files.py
#'
#' @param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., "/Users/me/apps/phylogenomic_dataset_construction/"
#' @param all_fasta Character vector of length one; the path to the fasta file including all the sequences that were originally used to build the trees.
#' @param tree_folder Character vector of length one; the path to the folder containing the trees to be used for extracting sequences.
#' @param outdir Character vector of length one; the path to the folder where the fasta files should be written.
#' @param minimal_taxa Numeric; minimal number of taxa required for output sequences to be written (regardless of ingroup/outgroup status).
#' @param overwrite Logical; should previous output of this command be erased so new output can be written? Once erased it cannot be restored, so use with caution!
#' @param get_hash Logical; should the 32-byte MD5 hash be computed for all fasta files concatenated together? Used for by \code{\link[drake]{drake_plan}} for tracking during workflows. If \code{TRUE}, this function will return the hash.
#' @param echo Logical; should the standard output and error be printed to the screen?
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.
#' @return One fasta file per tree will be written to \code{outdir}. If \code{get_hash} is \code{TRUE}, the 32-byte MD5 hash be computed for all \code{.fa} files concatenated together will be returned.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{write_ortholog_fasta_files(all_fasta = "some/folder/all.fasta", tree_folder = "some/folder/containing/tree/files", outdir = "some/folder", minimal_taxa = 5)}
#' @export
write_ortholog_fasta_files <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), all_fasta, tree_folder, outdir, minimal_taxa = 4, overwrite = FALSE, get_hash = TRUE, echo = pkgconfig::get_config("baitfindR::echo", fallback = FALSE), ...) {

  # error checking
  if(is.null(path_to_ys)) {
    stop("Must provide 'path_to_ys' (path to Yang & Smith Phylogenomic Dataset Analysis folder)")
  }

  # modify arguments
  path_to_ys <- fs::path_abs(path_to_ys)
  outdir <- fs::path_abs(outdir)
  tree_folder <- fs::path_abs(tree_folder)

  # define search terms for output files
  search_terms <- "\\.fa$"

  # optional: delete all previous output written in this folder
  if (isTRUE(overwrite)) {
    files_to_delete <- list.files(outdir, pattern = search_terms)
    if (length(files_to_delete) > 0) {
      files_to_delete <- paste0(outdir, files_to_delete)
      file.remove(files_to_delete)
    }
  }

  # call command
  arguments <- c(fs::path(path_to_ys, "write_ortholog_fasta_files.py"), all_fasta, tree_folder, outdir, minimal_taxa)
  processx::run("python", arguments, echo = echo)

  # optional: get MD5 hash of output
  if (isTRUE(get_hash)) {
    output <- list.files(outdir, pattern = search_terms)
    output <- if (length(output) > 0) {unlist( lapply(paste0(outdir, output), readr::read_file) )} else {output}
    hash <- digest::digest(output)
    return(hash)
  }
}

#' Align all fasta files in a directory.
#'
#' Given a directory containing unaligned fasta files, align all fasta files in
#' the directory. If there are > 1000 sequences in the directory, use the
#' mafft \code{--auto} algorithm. If less, use the \code{--genafpair}
#' algorithm.
#'
#' Wrapper for Yang and Smith (2014) mafft_wrapper.py
#'
#' @param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., "/Users/me/apps/phylogenomic_dataset_construction/"
#' @param fasta_folder Character vector of length one; the path to the folder containing the fasta files to be aligned.
#' @param infile_ending Character vector of length one; only files with this ending will be included.
#' @param number_cores Numeric; number of threads to use for and \code{mafft}.
#' @param seq_type Character vector of length one indicating type of sequences. Should either be \code{"dna"} for DNA or \code{"aa"} for proteins.
#' @param overwrite Logical; should previous output of this command be erased so new output can be written? Once erased it cannot be restored, so use with caution!
#' @param get_hash Logical; should the 32-byte MD5 hash be computed for all aligned fasta files concatenated together? Used for by \code{\link[drake]{drake_plan}} for tracking during workflows. If \code{TRUE}, this function will return the hash.
#' @param echo Logical; should the standard output and error be printed to the screen?
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.
#' @return Aligned fasta files will be written to \code{fasta_folder} with the file ending \code{.aln}. If \code{get_hash} is \code{TRUE}, the 32-byte MD5 hash be computed for all \code{.aln} files concatenated together will be returned.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{mafft_wrapper(fasta_folder = "some/folder/with/fasta/files", number_cores = 2, seq_type = "dna")}
#' @export
mafft_wrapper <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), fasta_folder, infile_ending = "fa", number_cores, seq_type = "dna", overwrite = FALSE, get_hash = TRUE, echo = pkgconfig::get_config("baitfindR::echo", fallback = FALSE), ...) {

  # error checking
  if(is.null(path_to_ys)) {
    stop("Must provide 'path_to_ys' (path to Yang & Smith Phylogenomic Dataset Analysis folder)")
  }

  # modify arguments
  path_to_ys <- fs::path_abs(path_to_ys)
  fasta_folder <- fs::path_abs(fasta_folder)

  # define search terms for output files
  search_terms <- "\\.mafft\\.aln$"

  # optional: delete all previous output written in this folder
  if (isTRUE(overwrite)) {
    files_to_delete <- list.files(fasta_folder, pattern = search_terms)
    if (length(files_to_delete) > 0) {
      files_to_delete <- paste0(fasta_folder, files_to_delete)
      file.remove(files_to_delete)
    }
  }

  # call command
  arguments <- c(fs::path(path_to_ys, "mafft_wrapper.py"), fasta_folder, infile_ending, number_cores, seq_type)
  processx::run("python", arguments, echo = echo)

  # optional: get MD5 hash of output
  if (isTRUE(get_hash)) {
    output <- list.files(fasta_folder, pattern = search_terms)
    output <- if (length(output) > 0) {unlist( lapply(paste0(fasta_folder, output), readr::read_file) )} else {output}
    hash <- digest::digest(output)
    return(hash)
  }
}

#' Clean all alignments in a directory.
#'
#' Given a directory containing aligned fasta files, clean the alignments
#' by removing columns below the specified occupancy cutoff.
#'
#' Wrapper for Yang and Smith (2014) phyutility_wrapper.py
#'
#' @param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., "/Users/me/apps/phylogenomic_dataset_construction/"
#' @param fasta_folder Character vector of length one; the path to the folder containing the alignments (fasta files) to be cleaned. Alignment files must end in \code{.aln}.
#' @param min_col_occup Numeric; characters (columns of the alignment) with less than this occupancy (as a decimal) will be removed from each alignment in the folder.
#' @param seq_type Character vector of length one indicating type of sequences. Should either be \code{"dna"} for DNA or \code{"aa"} for proteins.
#' @param overwrite Logical; should previous output of this command be erased so new output can be written? Once erased it cannot be restored, so use with caution!
#' @param get_hash Logical; should the 32-byte MD5 hash be computed for all result files concatenated together? Used for by \code{\link[drake]{drake_plan}} for tracking during workflows. If \code{TRUE}, this function will return the hash.
#' @param echo Logical; should the standard output and error be printed to the screen?
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.
#' @return Cleaned alignments will be written to \code{fasta_folder} with the file ending \code{.aln-cln}. If \code{get_hash} is \code{TRUE}, the 32-byte MD5 hash be computed for all \code{.aln-cln} files concatenated together will be returned.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{phyutility_wrapper(fasta_folder = "some/folder/with/alignments/", min_col_occup = 0.3, seq_type = "dna")}
#' @export
phyutility_wrapper <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), fasta_folder, min_col_occup, seq_type = "dna", overwrite = FALSE, get_hash = TRUE, echo = pkgconfig::get_config("baitfindR::echo", fallback = FALSE), ...) {

  # modify arguments
  path_to_ys <- fs::path_abs(path_to_ys)
  fasta_folder <- fs::path_abs(fasta_folder)

  # error checking
  if(is.null(path_to_ys)) {
    stop("Must provide 'path_to_ys' (path to Yang & Smith Phylogenomic Dataset Analysis folder)")
  }

  if (length(list.files(fasta_folder, pattern = "\\.aln$")) == 0) {
    stop(glue::glue('No file ending with .aln found in {fasta_folder}'))
  }

  # define search terms for output files
  search_terms <- "\\.aln-cln$"

  # optional: delete all previous output written in this folder
  if (isTRUE(overwrite)) {
    files_to_delete <- list.files(fasta_folder, pattern = search_terms)
    if (length(files_to_delete) > 0) {
      files_to_delete <- paste0(fasta_folder, files_to_delete)
      file.remove(files_to_delete)
    }
  }

  # call command
  arguments <- c(fs::path(path_to_ys, "phyutility_wrapper.py"), fasta_folder, min_col_occup, seq_type)
  processx::run("python", arguments, echo = echo)

  # optional: get MD5 hash of output
  if (isTRUE(get_hash)) {
    output <- list.files(fasta_folder, pattern = search_terms)
    output <- if (length(output) > 0) {unlist( lapply(paste0(fasta_folder, output), readr::read_file) )} else {output}
    hash <- digest::digest(output)
    return(hash)
  }
}
