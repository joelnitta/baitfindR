#' filter_fasta
#'
#' Filter a directory of fasta files by ingroup/outgroup status and taxonomic rank.
#'
#' Given a folder containing DNA sequences in multi-fasta format (i.e., each fasta file contains more than one sequence) and a dataframe including taxonomic data and ingroup/outgroup status, \code{filter_fasta()} outputs a list of those fasta files that pass one of two filters, or a combination of both. One filter excludes fasta files that do not contain greater than the minimum number of ingroup sequences. The other filter excludes fasta files that do not contain at least one sequence per ingroup taxon at the specified taxonomic rank.
#'
#' For example, if the dataset includes multiple ingroup genera each with multiple samples per genus, we may wish to filter alignments such that we only keep those with at least one sequence per ingroup genus. To do this, include a column called \code{"genus"} in \code{taxonomy_data}, and set \code{filter_col = "genus"}.
#'
#' @param seq_folder Character vector of length one; the path to the folder containing the fasta files (ending in \code{.fa} or \code{.fasta}) to filter.
#' @param taxonomy_data Dataframe matching sequences to ingroup/outgroup status and (optionally) higher-level taxonomic ranks for filtering. The columns must follow this format:
#' \describe{
#'   \item{sample}{Unique identifier for the source of the sequence, such as transcriptome IDs or species names. All sequences names must include such an identifier.}
#'   \item{group}{Either "in" or "out" (case-insensitive) depending if that sample is in the ingroup or outgroup.}
#'   \item{(user-selected taxonomic rank)}{The user can provide any taxonomic rank they wish to filter by. For example, alignments can be filtered by having at least one representative of each ingroup genus (family, order, etc.) in the dataset.}
#' }
#' @param sample_col Optional character; user-provided column name for \code{sample} in \code{taxonomy_data}.
#' @param group_col Optional character; user-provided column name for \code{group} \code{taxonomy_data}.
#' @param filter_col Optional character; the name of the column to be used for filtering by taxonomic rank in \code{taxonomy_data}.
#' @param min_taxa Minimum number of ingroup samples required to pass the filter.
#' @param exclude_short Logical; should extremely short sequences be excluded from the alignment during filtering? If \code{TRUE}, the minimum length is set to be within 1 standard deviation of the mean sequence length for a given alignment.
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.
#' @return A named list of DNA sequences of class \code{DNAbin} that passed the filter. These are not modified in any way; they simply met the requirements of the filter.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @examples
#' \dontrun{filter_fasta(
#'   seq_folder = "some/folder/",
#'   taxonomy_data = onekp_data,
#'   filter_col = "genus",
#'   min_taxa = 2)}
#'
#' @export
filter_fasta <- function (seq_folder, taxonomy_data, filter_col = NULL, min_taxa = NULL, exclude_short = FALSE, sample_col = "sample", group_col = "group", ...) {

  ### error checking and editing taxonomy data

  if(!is.data.frame(taxonomy_data)) {stop ("taxonomy_data must be of class data.frame")}

  # check column names
  if(!sample_col %in% colnames(taxonomy_data)) {stop ("sample_col must be in colnames(taxonomy_data)")}
  if(!group_col %in% colnames(taxonomy_data)) {stop ("group_col must be in colnames(taxonomy_data)")}
  if(length(colnames(taxonomy_data)) != length(unique(colnames(taxonomy_data)))) {stop ("colnames(taxonomy_data) must be unique")}

  # convert column names
  colnames(taxonomy_data)[grep(sample_col, colnames(taxonomy_data))] <- "sample"
  colnames(taxonomy_data)[grep(group_col, colnames(taxonomy_data))] <- "group"

  # convert to lowercase
  taxonomy_data$group <- tolower(taxonomy_data$group)

  # other checks on taxonomy_data
  checkr::check_data(taxonomy_data, values = list(
      sample = "a",
      group = c("out", "in")),
      order = FALSE, nrow = TRUE, exclusive = FALSE, key = "sample", error = TRUE)

  ### check that at least filter_col or min_taxa are provided
  if (is.null(filter_col) & is.null(min_taxa)) {
    stop("Must provide value for either filter_col, min_taxa, or both")
  }

  ### read fasta files

  seq_folder <- jntools::add_slash(seq_folder)

  fasta_names <- list.files(seq_folder, pattern = "\\.fa$|\\.fasta$")

  # check that there are fasta files to read
  if (length(fasta_names) == 0) {
    stop(paste("No files ending in .fa or .fasta in", seq_folder))
  }

  fasta_files <- lapply(paste0(seq_folder, fasta_names), ape::read.FASTA)

  # make list of unique user-provided higher-level taxa (e.g., genus, family, etc)
  if (!is.null(filter_col)) {
    if (!filter_col %in% colnames(taxonomy_data)) { stop ("'filter_col' not present in 'taxonomy_data' column names") }
    taxa_filter_list <- sort(unique(taxonomy_data[[filter_col]][taxonomy_data$group == "in"]))
  }

  # make list of ingroup/outgroup status
  ingroup <- taxonomy_data$sample[taxonomy_data$group == "in"]

  # filter function that works on a single alignment ------------------------
  # input: alignment = candidate bait fasta file
  #        ingroup = vector of samples that are in the ingroup
  #        exclude_short = described in documentation for the whole function at the top.
  #        taxonomy_data = described in documentation for the whole function at the top.
  #        taxa_filter_list = vector of higher-level ingroup taxa to use for filtering
  # output: single logical value for whether or not the candidate alignment passes the filter
  filter_alignment <- function (alignment, ingroup, exclude_short, taxonomy_data, taxa_filter_list, min_taxa) {

    # first convert aligment to list if it isn't already
    alignment <- as.list(alignment)

    # optionally remove extremely short sequences from consideration during filtering
    if (exclude_short) {
      # minimum cutoff seq length is mean overall seq length - 1 sd (mean taken from all species, not just ingroup)
      min_length <- mean(as.numeric(lapply(alignment, length))) - sd(as.numeric(lapply(alignment, length)))
      alignment <- alignment[sapply(alignment, length) > min_length ]
    }

    # subset alignment to only ingroup taxa
    alignment <- alignment[names(alignment) %in% ingroup]

    # optionally filter by higher-level taxa (passes filter by default)
    if (!is.null(filter_col)) {
      # get list of user-specified ingroup higher-level taxa present in subsetted alignment
      ingroup_higher_taxa_in_alignment <- taxonomy_data[[filter_col]][match(names(alignment), taxonomy_data$sample)]
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
  pass_filter <- sapply(fasta_files, filter_alignment, ingroup = ingroup, exclude_short = exclude_short, taxonomy_data = taxonomy_data, taxa_filter_list = taxa_filter_list, min_taxa = min_taxa)

  # make lists of passing fasta files and their names
  filtered_fasta_files <- fasta_files[pass_filter]
  filtered_fasta_names <- fasta_names[pass_filter]

  names(filtered_fasta_files) <- filtered_fasta_names

  return(filtered_fasta_files)
}
