#' filter_fasta
#'
#' Additional options for filtering fasta files produced by Yang and Smith 2014 "Paralogy pruning" scripts
#'
#' The minimal_taxa setting in Y&S Step 6 "Paralogy pruning" scripts (\code{filter_1to1_orthologs.py}, \code{prune_paralogs_MI.py}, etc) filters paralog trees by a minimum number of taxa without considering ingroup / outgroup status (except for \code{prune_paralogs_RT.py}).
#'
#' Use \code{filter_fasta} to further filter the results from Y&S "Paralogy pruning" scripts by outgroup/ingroup status. \code{filter_fasta} can also filter by higher-level taxonomic groups provided by the user in the taxonomy_data argument. For example, if the dataset has multiple species per genus, we may wish to filter alignments such that we only keep those with at least one species per unique ingroup genus. To do this, include a column called "genus" in \code{taxonomy_data}, and set \code{filter_level} to "genus."
#'
#' @param MCL_settings Settings used for mcl (Markov Cluster Algorithm) step in Y&S pipeline. Should be in the format "hit-frac0.3_I1.4_e5", where "hit-frac" is hit_fraction_cutoff used by blast_to_mcl.py, "I" is the inflation value used by mcl, and "e" specifies minimal log transformed evalue passed to -tf 'gq()' in mcl (for details, see Yang and Smith 2014).
#' @param prune_method Strategy used to identify homologs in Y&S pipeline. Must be one of the following: "ortho_121", "ortho_MO", "ortho_MI", or "ortho_RT".
#' @param taxonomy_data Dataframe matching taxonID to ingroup/outgroup status and (optional) higher-level taxonomy for filtering. The columns must follow this format:
#' \describe{
#'   \item{taxonID}{Short unique identifier for the taxon (i.e., source of the transcriptome analyzed in Y&S pipeline). Must be exactly the same taxonID used in the Y&S pipeline.}
#'   \item{group_status}{Either "IN" or "OUT" depending if that taxon is in the ingroup or outgroup.}
#'   \item{(user-selected taxonomic rank)}{The user can provide any taxonomic rank they wish to filter by. For example, alignments can be filtered by having at least one representative of each genus (family, order, etc.) in the dataset.}
#' }
#' @param filter_level A single character matching the name of the column to be used for filtering in \code{taxonomy_data}.
#' @param min_taxa Minimum number of ingroup taxa required to pass the filter.
#' @param exclude_short Logical; should extremely short sequences be excluded from the alignment during filtering? If \code{TRUE}, the minimum length is set to be within 1 standard deviation of the mean sequence length for a given alignment.
#' @return A named list including:
#' \describe{
#'   \item{filtered_fasta_files}{A list of fasta files that passed the filter. These are not modified in any way; they simply met the requirements of the filter.}
#'   \item{filtered_fasta_names}{A character vector of the names of fasta files that passed the filter.}
#' }
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#'
#' @export
filter_fasta <- function (MCL_settings, prune_method, taxonomy_data, filter_level = NULL, min_taxa = NULL, exclude_short = FALSE) {

  # error checking ----------------------------------------------------------
  # check that at least filter_level or min_taxa are provided
  if (is.null(filter_level) & is.null(min_taxa)) {
    stop("Must provide value for either filter_level, min_taxa, or both")
  }

  # check that prune_method is one of the accepted values
  if (!prune_method %in% c("ortho_121", "ortho_MO", "ortho_MI", "ortho_RT")) {
    stop ("prune_method not one of 'ortho_121', 'ortho_MO', 'ortho_MI', or 'ortho_RT'")
  }

  # check if ./MCL_settings/prune_method/ exists
  if(!any(grepl(paste0(MCL_settings, "/", prune_method), list.dirs() ))) {
    stop ("Can't find directory containing fasta files in working directory. Check MCL_settings and prune_method.")
  }

  # check if ./MCL_settings/prune_method/fasta/ exists
  if(!any(grepl(paste0(MCL_settings, "/", prune_method, "/fasta"), list.dirs() ))) {
    stop ("Can't find 'fasta' directory within MCL_settings/prune_method")
  }

  # get names of all fasta files in ./MCL_settings/prune_method/fasta/
  fasta_names <- list.files(paste0(MCL_settings, "/", prune_method, "/fasta/"))
  fasta_names <- fasta_names[grep(".fa$", fasta_names)]

  # check that there are fasta files to read
  if (length(fasta_names) == 0) {
    stop("No files named .fa or .fasta in MCL_settings/prune_method/fasta/")
  }

  # read in fasta files, set up filters -------------------------------------

  # read in all fasta files in ./MCL_settings/prune_method/fasta/
  fasta_files <- lapply(paste(MCL_settings, "/", prune_method, "/fasta/", fasta_names, sep=""), ape::read.dna, format="fasta")

  # make list of unique user-provided higher-level taxa (e.g., genus, family, etc)
  if (!is.null(filter_level)) {
    if (!filter_level %in% colnames(taxonomy_data)) { stop ("'filter_level' not present in 'taxonomy_data' column names") }
    taxa_filter_list <- sort(unique(taxonomy_data[[filter_level]][taxonomy_data$group_status == "IN"]))
  }

  # make list of ingroup/outgroup status
  ingroup <- taxonomy_data$taxonID[taxonomy_data$group_status == "IN"]

  # filter function that works on a single alignment ------------------------
  # input: alignment = candidate bait fasta file
  #        ingroup = vector of taxonIDs that are in the ingroup
  #        exclude_short = described in documentation for the whole function at the top.
  #        taxonomy_data = described in documentation for the whole function at the top.
  #        taxa_filter_list = vector of higher-level ingroup taxa to use for filtering
  # output: single logical value for whether or not the candidate alignment passes the filter
  filter_alignment <- function (alignment, ingroup, exclude_short, taxonomy_data, taxa_filter_list, min_taxa) {

    # optionally remove extremely short sequences from consideration during filtering
    if (exclude_short) {
      # minimum cutoff seq length is mean overall seq length - 1 sd (mean taken from all species, not just ingroup)
      min_length <- mean(as.numeric(lapply(alignment, length))) - sd(as.numeric(lapply(alignment, length)))
      alignment <- alignment[sapply(alignment, length) > min_length ]
    }

    # subset alignment to only ingroup taxa
    alignment <- alignment[names(alignment) %in% ingroup]

    # optionally filter by higher-level taxa (passes filter by default)
    if (!is.null(filter_level)) {
      # get list of user-specified ingroup higher-level taxa present in subsetted alignment
      ingroup_higher_taxa_in_alignment <- taxonomy_data[[filter_level]][match(names(alignment), taxonomy_data$taxonID)]
      # check if all ingroup higher-level taxa are included at least once in the subsetted alignment
      taxa_filter_result <- all(taxa_filter_list %in% ingroup_higher_taxa_in_alignment)
    } else {
      taxa_filter_result <- TRUE
    }

    # optionally filter by minimum number of ingroup taxa (passes filter by default)
    if (!is.null(min_taxa)) {
      min_taxa_filter_result <- length(alignment) >= min_taxa
    } else {
      min_taxa_filter_result <- TRUE
    }

    # combine the two filters for the final result: did the alignment pass both filters?
    all(taxa_filter_result, min_taxa_filter_result)
  }

  # apply above function to all alignments (ie, fasta files)
  pass_filter <- sapply(fasta_files, filter_alignment, ingroup=ingroup, exclude_short=TRUE, taxonomy_data=taxonomy_data, taxa_filter_list=taxa_filter_list, min_taxa=NULL)

  # make lists of passing fasta files and their names
  filtered_fasta_files <- fasta_files[pass_filter]
  filtered_fasta_names <- fasta_names[pass_filter]

  results <- list (
    filtered_fasta_files = filtered_fasta_files,
    filtered_fasta_names = filtered_fasta_names
  )

  return(results)
}
