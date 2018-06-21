#' fasta_to_tree
#'
#' Wrapper for Yang and Smith (2014) \code{fasta_to_tree.py}
#'
#' Given a folder containing unaligned sequences in fasta format (i.e., clusters),
#' aligns each cluster with \code{mafft} (small clusters) or \code{pasta} (large clusters),
#' excludes poorly aligned sites with \code{phyutility}, and infers a maximum-likelihood
#' tree with \code{RAxML} (small clusters) or \code{fasttree} (large clusters). Requires all
#' of these programs to be installed and included in the user's \code{$PATH}. Assumes clusters are named \code{"cluster1.fa"}, \code{"cluster2.fa"}, etc. Clusters with fewer than 1,000 sequences are considered "small," and those with more are considered "large."
#'
#' @param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., \code{"/Users/me/apps/phylogenomic_dataset_construction/"}
#' @param overwrite Logical; should previous output of this command be erased so new output can be written? Once erased it cannot be restored, so use with caution!
#' @param seq_folder Character vector of length one; the path to the folder containing the fasta files.
#' @param number_cores Numeric; number of threads to use for \code{RAxML} and \code{mafft}.
#' @param seq_type Character vector of length one indicating type of sequences. Should either be \code{"dna"} for DNA or \code{"aa"} for proteins.
#' @param bootstrap Logical; should run a bootstrap analysis be run for the trees?
#' @param get_hash Logical; should the 32-byte MD5 hash be computed for all output tree files concatenated together? Used for by \code{\link{drake}} for tracking during workflows. If \code{TRUE}, this function will return the hash.
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link{drake}} for tracking during workflows.
#' @return For each input cluster \code{cluster1.fa} in \code{seq_folder}, \code{cluster1.fa.mafft.aln} (small clusters) or \code{cluster1.pasta.aln} (large clusters), \code{cluster1.fa.mafft.aln-cln} (small clusters) or \code{cluster1.fa.pasta.aln-cln} (large clusters), and \code{cluster1.raxml.tre} (small clusters) or \code{cluster1.fasttree.tre} (large clusters) will be written to \code{seq_folder}. If \code{get_hash} is \code{TRUE}, the 32-byte MD5 hash be computed for all \code{.tre} files concatenated together will be returned.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{fasta_to_tree(seq_folder = "some/folder/containing/fasta/seqs", number_cores = 1, seq_type = "dna", bootstrap = FALSE)}
#' @export
fasta_to_tree <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), overwrite = FALSE, seq_folder, number_cores, seq_type, bootstrap = FALSE, get_hash = TRUE, ...) {

  # error checking
  if(is.null(path_to_ys)) {
    stop("Must provide 'path_to_ys' (path to Yang & Smith Phylogenomic Dataset Analysis folder)")
  }

  # modify arguments
  bootstrap <- ifelse(bootstrap, "y", "n")
  path_to_ys <- jntools::add_slash(path_to_ys)
  seq_folder <- jntools::add_slash(seq_folder)

  # optional: delete all previous output written by a previous call to fasta_to_tree.py in this folder
  if (overwrite) {
    search_terms <- paste("cluster.*\\.fa\\.mafft\\.aln$",
                           "cluster.*\\.fa\\.mafft\\.aln-cln$",
                           "cluster.*\\.fa\\.mafft\\.aln-cln.reduced$",
                           "cluster.*\\.raxml\\.tre$",
                           sep = "|")
    files_to_delete <- list.files(seq_folder, pattern = search_terms)
    if (length(files_to_delete) > 0) {
      files_to_delete <- paste0(seq_folder, files_to_delete)
      file.remove(files_to_delete)
    } else {
      print("No files to overwrite, continuing")
    }
  }

  # call fasta_to_tree.py
  arguments <- c(paste0(path_to_ys, "fasta_to_tree.py"), seq_folder, number_cores, seq_type, bootstrap)
  processx::run("python", arguments, wd = seq_folder)

  # optional: get MD5 hash of concatenated trees
  if (get_hash) {
    trees <- list.files(seq_folder)
    trees <- trees[grep("\\.raxml\\.tre", trees)]
    trees <- paste0(seq_folder, trees)
    trees <- unlist(lapply(trees, readr::read_file))
    hash <- digest::digest(trees)
    return(hash)
  }
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
#' @param overwrite Logical; should previous output of this command be erased so new output can be written? Once erased it cannot be restored, so use with caution!
#' @param all_fasta Character vector of length one; the path to the fasta file including all query sequences concatenated together, i.e., the fasta file used to create the "all-by-all" blast database.
#' @param mcl_outfile Character vector of length one; the path to the output from running mcl on blast distances.
#' @param minimal_taxa Numeric; minimal number of taxa required to be present for the cluster to be written. Default 4, the minimum number of taxa needed for an un-rooted tree.
#' @param outdir Character vector of length one; the path to the folder where the clusters should be written.
#' @param get_hash Logical; should the 32-byte MD5 hash be computed for all clusters concatenated together? Used for by \code{\link{drake}} for tracking during workflows. If \code{TRUE}, this function will return the hash.
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link{drake}} for tracking during workflows.
#' @return One fasta file per cluster (\code{cluster1.fa}, \code{cluster2.fa}, etc.) will be written to \code{outdir}. If \code{get_hash} is \code{TRUE}, the 32-byte MD5 hash be computed for all \code{.fa} files concatenated together will be returned.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{write_fasta_files_from_mcl(all_fasta = "some/folder/all.fasta", mcl_outfile = "some/folder/hit-frac0.4_I1.4_e5", minimal_taxa = 5, outdir = "some/folder")}
#' @export
write_fasta_files_from_mcl <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), overwrite, all_fasta, mcl_outfile, minimal_taxa = 4, outdir, get_hash, ...) {

  # error checking
  if(is.null(path_to_ys)) {
    stop("Must provide 'path_to_ys' (path to Yang & Smith Phylogenomic Dataset Analysis folder)")
  }

  # modify arguments
  path_to_ys <- jntools::add_slash(path_to_ys)
  outdir <- jntools::add_slash(outdir)


  # optional: delete all previous output written by a previous call to fasta_to_tree.py in this folder
  if (overwrite) {
    files_to_delete <- list.files(outdir)
    search_terms <- "cluster\\d*\\.fa$"
    files_to_delete <- files_to_delete[grep(search_terms, files_to_delete)]
    files_to_delete <- paste0(outdir, files_to_delete)
    file.remove(files_to_delete)
  }

  # call write_fasta_files_from_mcl.py
  arguments <- c(paste0(path_to_ys, "write_fasta_files_from_mcl.py"), all_fasta, mcl_outfile, minimal_taxa, outdir)
  processx::run("python", arguments)

  # optional: get MD5 hash of concatenated clusters
  if (get_hash) {
    clusters <- list.files(outdir)
    clusters <- clusters[grep("\\.fa$", clusters)]
    clusters <- paste0(outdir, clusters)
    clusters <- unlist(lapply(clusters, readr::read_file))
    hash <- digest::digest(clusters)
    return(hash)
  }

}

#' blast_to_mcl
#'
#'Wrapper for Yang and Smith (2014) blast_to_mcl.py
#'
#'Converts the output of an all-by-all blast query into a format that can be parsed by mcl to find clusters.
#'
#' @param path_to_ys Character vector of length one; the complete path to the folder containing Y&S python scripts, e.g., "/Users/me/apps/phylogenomic_dataset_construction/"
#' @param blast_results Character vector of length one; the complete path to the tab-separated text file containing the results from an all-by-all blast search. If blast searches were run separately (i.e., one for each sample), the results should be concatenated into a single file. For the blast search, the output format should specified as: -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore'
#' @param hit_fraction_cutoff Numeric between 0 and 1. Indicates the minimum percentage overlap between query and target in blast results to be retained in the output. According to Y&S, "A low hit-fraction cutoff will output clusters with more incomplete sequences and much larger and sparser alignments, whereas a high hit-fraction cutoff gives tighter clusters but ignores incomplete or divergent sequences."
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link{drake}} for tracking during workflows.
#' @return A tab-separated text file with three columns: the first two are the matching query and target from the all-by-all blast, and the third is the negative log e-value for that match. This file is named \code{<blast_results>.hit-frac<hit_fraction_cutoff>.minusLogEvalue}, where \code{<blast_results>} and \code{<hit_fraction_cutoff>} correspond to the values of those arguments. If possible contaminants (i.e., identical sequences between different samples) were found, these are written to \code{<blast_results>.ident.hit-frac<hit_fraction_cutoff>}. Output files will be written to the same folder containing \code{blast_results}.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{blast_to_mcl(blast_results = "some/folder/blastresults.tab", hit_fraction_cutoff = 0.5)}
#' @export
blast_to_mcl <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), blast_results, hit_fraction_cutoff, ...) {
  path_to_ys <- jntools::add_slash(path_to_ys)
  arguments <- c(paste0(path_to_ys, "blast_to_mcl.py"), blast_results, hit_fraction_cutoff)
  processx::run("python", arguments)

  # Normally, the warning file with sequences that are identical between
  # samples (possible contamination) is output as "blast_output.ident",
  # but this results in files with the same name from different hit_faction_cutoff
  # values. Append the hit_fraction_cutoff value so we can tell them apart.
  #
  ident_file <- paste0(blast_results, ".ident")
  ident_file_rename <- paste0(ident_file, ".hit-frac", hit_fraction_cutoff)
  if (file.exists(ident_file)) {
    file.rename(ident_file, ident_file_rename)
  }
}

#' fix_names_from_transdecoder
#'
#' Shortens names in fasta headers.
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

#' trim_tips
#'
#' Wrapper for Yang and Smith (2014) \code{trim_tips.py}
#'
#' Given a folder containing phylogenetic trees, exclude (i.e., "trim"),
#' tips on unusually long branches. Tips on a branch 10 times longer
#' than their sister AND longer than \code{relative_cutoff}, OR tips
#' that are longer than \code{absolute_cutoff} will be trimmed. This
#' function will overwrite any output files with the same name in
#' \code{tree_folder}.
#'
#' @param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., \code{"/Users/me/apps/phylogenomic_dataset_construction/"}
#' @param tree_folder Character vector of length one; the path to the folder containing the trees to trim.
#' @param tree_file_ending Character vector of length one; only tree files with this file ending will be used.
#' @param relative_cutoff Numeric vector of length one; tips on a branch 10 times longer than their sister AND longer than this value will be cut.
#' @param absolute_cutoff Numeric vector of length one; tips on branches longer than this value will be cut.
#' @param get_hash Logical; should the 32-byte MD5 hash be computed for all output trimmed tree files concatenated together? Used for by \code{\link{drake}} for tracking during workflows. If \code{TRUE}, this function will return the hash.
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link{drake}} for tracking during workflows.
#' @return For each input tree with a file ending matching \code{tree_file_ending} in \code{tree_folder}, a trimmed tree with a file ending in \code{.tt} will be written to \code{tree_folder}. If \code{get_hash} is \code{TRUE}, the 32-byte MD5 hash be computed for all trimmed tree files concatenated together will be returned.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{trim_tips(tree_folder = "some/folder/containing/tree/files", tree_file_ending = ".tre", relative_cutoff = 0.2, absolute_cutoff = 0.4)}
#' @export
trim_tips <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), tree_folder, tree_file_ending, relative_cutoff, absolute_cutoff, get_hash = TRUE, ...) {

  # error checking
  if(is.null(path_to_ys)) {
    stop("Must provide 'path_to_ys' (path to Yang & Smith Phylogenomic Dataset Analysis folder)")
  }

  # modify arguments
  path_to_ys <- jntools::add_slash(path_to_ys)
  tree_folder <- jntools::add_slash(tree_folder)

  # call trim_tips.py
  arguments <- c(paste0(path_to_ys, "trim_tips.py"), tree_folder, tree_file_ending, relative_cutoff, absolute_cutoff)
  processx::run("python", arguments)

  # optional: get MD5 hash of concatenated trees
  if (isTRUE(get_hash)) {
    trimmed_trees <- list.files(tree_folder, pattern = "\\.tt")
    trimmed_trees <- paste0(tree_folder, trimmed_trees)
    trimmed_trees <- unlist(lapply(trimmed_trees, readr::read_file))
    hash <- digest::digest(trimmed_trees)
    return(hash)
  }
}

#' mask_tips_by_taxonID_transcripts
#'
#' Wrapper for Yang and Smith (2014) \code{mask_tips_by_taxonID_transcripts.py}
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
#' @param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., \code{"/Users/me/apps/phylogenomic_dataset_construction/"}
#' @param tree_folder Character vector of length one; the path to the folder containing the trees to mask.
#' @param aln_folder Character vector of length one; the path to the folder containing the alignments used to make the trees.
#' @param mask_paraphyletic Logical; should paraphyletic tips belonging to the same taxon be masked?
#' @param get_hash Logical; should the 32-byte MD5 hash be computed for all output masked tree files concatenated together? Used for by \code{\link{drake}} for tracking during workflows. If \code{TRUE}, this function will return the hash.
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link{drake}} for tracking during workflows.
#' @return For each input tree with a file ending in \code{.tt} in \code{tree_folder}, a trimmed tree with a file ending in \code{.mm} will be written to \code{tree_folder}. If \code{get_hash} is \code{TRUE}, the 32-byte MD5 hash be computed for all masked tree files concatenated together will be returned.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{mask_tips_by_taxonID_transcripts(tree_folder = "some/folder/containing/tree/files", aln_folder = "some/folder/containing/alignment/files")}
#' @export
mask_tips_by_taxonID_transcripts <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), tree_folder, aln_folder, mask_paraphyletic = TRUE, get_hash = TRUE, ...) {

  # error checking
  if(is.null(path_to_ys)) {
    stop("Must provide 'path_to_ys' (path to Yang & Smith Phylogenomic Dataset Analysis folder)")
  }

  # modify arguments
  path_to_ys <- jntools::add_slash(path_to_ys)
  tree_folder <- jntools::add_slash(tree_folder)
  aln_folder <- jntools::add_slash(aln_folder)
  mask_paraphyletic <- ifelse(isTRUE(mask_paraphyletic), "y", "n")

  # call mask_tips_by_taxonID_transcripts.py
  arguments <- c(paste0(path_to_ys, "mask_tips_by_taxonID_transcripts.py"), tree_folder, aln_folder, mask_paraphyletic)
  processx::run("python", arguments)

  # optional: get MD5 hash of concatenated trees
  if (isTRUE(get_hash)) {
    masked_trees <- list.files(tree_folder, pattern = "\\.mm$")
    masked_trees <- paste0(tree_folder, masked_trees)
    masked_trees <- unlist(lapply(masked_trees, readr::read_file))
    hash <- digest::digest(masked_trees)
    return(hash)
  }
}

#' cut_long_internal_branches
#'
#' Wrapper for Yang and Smith (2014) \code{cut_long_internal_branches.py}
#'
#' Given a folder containing phylogenetic trees, split the trees into multiple subtrees
#' for nodes that bifurcate deeper than \code{internal_branch_length_cutoff}.
#' \code{tree_folder} and \code{outdir} should be different to avoid writing over input trees.
#' This function will overwrite any output files with the same name in \code{outdir}.
#'
#' @param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., \code{"/Users/me/apps/phylogenomic_dataset_construction/"}
#' @param tree_folder Character vector of length one; the path to the folder containing the trees to cut.
#' @param tree_file_ending Character vector of length one; only tree files with this file ending will be used.
#' @param internal_branch_length_cutoff Numeric vector of length one; the depth at which cuts should be made (smaller numbers indicate greater depth).
#' @param minimal_taxa Numeric; minimal number of taxa required for tree to be cut. Default 4, the minimum number of taxa needed for an un-rooted tree.
#' @param outdir Character vector of length one; the path to the folder where the subtrees should be written.
#' @param get_hash Logical; should the 32-byte MD5 hash be computed for all output subtree files concatenated together? Used for by \code{\link{drake}} for tracking during workflows. If \code{TRUE}, this function will return the hash.
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link{drake}} for tracking during workflows.
#' @return For each input tree with a file ending in \code{tree_file_ending} in \code{tree_folder}, one or more subtrees with a file ending in \code{.subtree} will be written to \code{tree_folder}. If \code{get_hash} is \code{TRUE}, the 32-byte MD5 hash be computed for all subtree files concatenated together will be returned.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{cut_long_internal_branches(tree_folder = "some/folder/containing/tree/files", tree_file_ending = ".mm", internal_branch_length_cutoff = 0.3, outdir = "some/other/folder/")}
#' @export
cut_long_internal_branches <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), tree_folder, tree_file_ending, internal_branch_length_cutoff, minimal_taxa = 4, outdir, get_hash = TRUE, ...) {

  # error checking
  if(is.null(path_to_ys)) {
    stop("Must provide 'path_to_ys' (path to Yang & Smith Phylogenomic Dataset Analysis folder)")
  }

  # modify arguments
  path_to_ys <- jntools::add_slash(path_to_ys)
  tree_folder <- jntools::add_slash(tree_folder)
  outdir <- jntools::add_slash(outdir)

  # more error checking
  if(tree_folder == outdir) {
    stop("Must provide provide different paths for input and output folders")
  }

  # call cut_long_internal_branches.py
  arguments <- c(paste0(path_to_ys, "cut_long_internal_branches.py"), tree_folder, tree_file_ending, internal_branch_length_cutoff, minimal_taxa, outdir)
  processx::run("python", arguments)

  # optional: get MD5 hash of concatenated trees
  if (isTRUE(get_hash)) {
    sub_trees <- list.files(outdir, pattern = "\\.subtree$")
    sub_trees <- paste0(outdir, sub_trees)
    sub_trees <- unlist(lapply(sub_trees, readr::read_file))
    hash <- digest::digest(sub_trees)
    return(hash)
  }
}

#' write_fasta_files_from_trees
#'
#' Wrapper for Yang and Smith (2014) write_fasta_files_from_trees.py
#'
#' Given a folder containing phylogenetic trees and a single concatenated fasta file
#' including all the sequences used to build the trees, output one fasta file per tree
#' with the sequences in that tree. This function will overwrite any output files with
#' the same name in \code{outdir}.
#'
#' @param path_to_ys Character vector of length one; the path to the folder containing Y&S python scripts, e.g., "/Users/me/apps/phylogenomic_dataset_construction/"
#' @param all_fasta Character vector of length one; the path to the fasta file including all the sequences that were originally used to build the trees.
#' @param tree_folder Character vector of length one; the path to the folder containing the trees to be used for extracting fasta sequences.
#' @param tree_file_ending Character vector of length one; only tree files with this file ending will be used.
#' @param outdir Character vector of length one; the path to the folder where the fasta files should be written.
#' @param get_hash Logical; should the 32-byte MD5 hash be computed for all output fasta files concatenated together? Used for by \code{\link{drake}} for tracking during workflows. If \code{TRUE}, this function will return the hash.
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link{drake}} for tracking during workflows.
#' @return One fasta file per tree file ending in \code{tree_file_ending} in \code{tree_folder} will be written to \code{outdir}. If \code{get_hash} is \code{TRUE}, the 32-byte MD5 hash be computed for all output fasta files concatenated together will be returned.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @examples
#' \dontrun{write_fasta_files_from_trees(all_fasta = "some/folder/all.fasta", tree_file_ending = ".subtree", tree_folder = "some/folder/containing/tree/files", outdir = "some/folder")}
#' @export
write_fasta_files_from_trees <- function (path_to_ys = pkgconfig::get_config("baitfindR::path_to_ys"), all_fasta, tree_folder, tree_file_ending, outdir, get_hash, ...) {

  # error checking
  if(is.null(path_to_ys)) {
    stop("Must provide 'path_to_ys' (path to Yang & Smith Phylogenomic Dataset Analysis folder)")
  }

  # modify arguments
  path_to_ys <- jntools::add_slash(path_to_ys)
  outdir <- jntools::add_slash(outdir)
  tree_folder <- jntools::add_slash(tree_folder)

  # call write_fasta_files_from_mcl.py
  arguments <- c(paste0(path_to_ys, "write_fasta_files_from_trees.py"), all_fasta, tree_folder, tree_file_ending, outdir)
  processx::run("python", arguments)

  # optional: get MD5 hash of concatenated clusters
  if (get_hash) {
    clusters <- list.files(outdir)
    clusters <- clusters[grep("rr\\.fa$", clusters)]
    clusters <- paste0(outdir, clusters)
    clusters <- unlist(lapply(clusters, readr::read_file))
    hash <- digest::digest(clusters)
    return(hash)
  }
}
