#' filter_fasta
#'
#' The minimal_taxa setting in Y&S Step 6 "Paralogy pruning" scripts filter paralog trees
#' by a minimal number of taxa without considering ingroup / outgroup status (except for
#' RT). Use this script to further filter the results from filter_1to1_orthologs.py,
#' prune_paralogs_MI.py etc by outgroup/ingroup status.
#'
#' @param MCL_settings Settings used for MCL (Markov Cluster Algorithm) step in Yang & Smith
#' pipeline. Should be in the format "hit-frac0.3_I1.4_e5".
#' @param prune_method Strategy used to identify homologs in Yang & Smith pipeline.
#' Must be one of the following: "ortho_121", "ortho_MO", "ortho_MI", or "ortho_RT".
#' @param taxaset
#' @param filter_type Type of filter to be applied to fasta files. Choose to filter by
#' ingroup family, subfamily, or genus (e.g., a given alignment must contain at least one
#' representative of each ingroup family, subfamily, etc), or by minimum number of ingroup
#' species. Must be one of the following: "family", "subfamily", "genus", "min_taxa"
#' @param min_taxa If filter_type is "min_taxa", this sets the minimum number of ingroup taxa required
#' to pass the filter.
#' @param taxonomy_data Dataframe matching taxonID to species, genus, subfamily, and family.
#' Must also include outgroup/ingroup status. The columns must follow this format:
#' \describe{
#'   \item{taxonID}{Short unique identifier for the taxon (usually species). Must be exactly the
#'   same taxonID used in the Y&S pipeline.}
#'   \item{group_status}{Either "IN" or "OUT" depending if that taxon is in the ingroup or outgroup.}
#'   \item{species}{Species in latin binomial form, with genus and species separated by space or underscore
#'   (optional, only needed if filtering by genus).}
#'   \item{family}{Taxonomic family (optional, only needed if filtering by family).
#'   \item{subfamily}{Taxonomic subfamily (optional, only needed if filtering by subfamily). Does not need to be
#'   specified for every taxon.}
#' }
#'
# given a set of fasta files with some minimum TOTAL taxa cutoff from an ortholog pruning step,
# filter the fasta files down to only those including a minimum number of INGROUP taxa only

# function to extract cluster name(s) from alignment, including "cluster" part and "1rr" part
extract_cluster_name <- function (align) {
  # grab either names or rownames of alignment depending if its in matrix format or a list
  if (length(dim(align)) == 0) {
    otu <- names(align)
  } else if (length(dim(align)) > 0) {
    otu <- rownames(align)
  }

  otu <- otu[grep("cluster", otu)]
  otu <- sapply(otu, function (x) strsplit(x, split="_"))
  otu_first <- sapply(otu, function (x) x[2])
  otu_second <- sapply(otu, function (x) x[3])
  otu <- paste(otu_first, otu_second, sep="_")
  cluster_name <- unique(otu)
  return(cluster_name)
}

filter_fasta <- function (MCL_settings, prune_method, taxaset, filter_type="min_taxa", min_taxa) {

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
  fasta.names <- list.files(paste0(MCL_settings, "/", prune_method, "/fasta/"))
  fasta.names <- fasta.names[grep(".fa$", fasta.names)]

  # check that there are fasta files to read
  if (length(fasta.names) == 0) {
    stop("No files named .fa or .fasta in MCL_settings/prune_method/fasta/")
  }

  # read in all fasta files in ./MCL_settings/prune_method/fasta/
  fasta.files <- lapply(paste(MCL_settings, "/", prune_method, "/fasta/", fasta.names, sep=""), ape::read.dna, format="fasta")

  # get list of unique ingroup families in taxonomy data
  family_list <- taxonomy_data$family[taxonomy_data$group_status == "IN"]
  ingroup <- taxonomy_data$taxonID[taxonomy_data$group_status == "IN"]

  # The following functions are different ways to filter the fasta files

  # initial_filter
  # first part of filtering step that is used in all filtering methods
  # input: alignment = candidate bait fasta file
  #        ingroup = vector of taxonIDs that are in the ingroup
  #        length_cutoff = logical; should the sequences in the alignment be filtered by a minimum length?
  #                        If TRUE, the minimum length is set to be within 1 standard deviation of the mean sequence length

  initial_filter <- function (alignment, ingroup, length_cutoff) {
    # minimum cutoff seq length is mean overall seq length - 1 sd (mean taken from all species, not just ingroup)
    min_length <- mean(as.numeric(lapply(alignment, length))) - sd(as.numeric(lapply(alignment, length)))
    # subsetted alignment with just ingroup taxa
    ingroup_alignment <- alignment[names(alignment) %in% ingroup]
    # optionally filter by minimum length
    if (length_cutoff) {
      ingroup_alignment <- ingroup_alignment[sapply(ingroup_alignment, length) > min_length ] }
    return(ingroup_alignment)
  }

  # filter_by_family
  # filter fasta files based on having at least one species per ingroup family and optionally by sequence length
  # input: alignment, ingroup, length_cutoff are as in initial_filter above.
  #        taxonomy_data = described in documentation for the whole function at the top.
  #        family_list = vector of families in the ingroup
  # output: single logical value for whether or not the candidate alignment passes the filter

  filter_by_family <- function (alignment, ingroup, length_cutoff, taxonomy_data, family_list) {
    # run initial filter to trim to only ingroup taxa, optionally dropping short sequences
    ingroup_alignment <- initial_filter(alignment, ingroup, length_cutoff)
    # get list of families for subsetted alignment
    ingroup_families_in_alignment <- taxonomy_data$family[match(names(ingroup_alignment), taxonomy_data$taxonID)]
    # check if all ingroup families are included at least once in the subsetted alignment
    all(family_list %in% ingroup_families_in_alignment)
  }

  # filter fasta files based on having at least one species per eupoly II subfamily with seq length within 1 sd of overall mean seq length
  # input: candidate bait fasta file
  # output: T/F whether the candidate passes the filter
  filter_eupoly2_subfamily <- function (alignment) {
    # minimum cutoff seq length is mean overall seq length - 1 sd
    min_length <- mean(as.numeric(lapply(alignment, length))) - sd(as.numeric(lapply(alignment, length)))
    # subsetted alignment with just EuII taxa
    # eu2_alignment <- alignment[full_subfamilies %in% eu2_subfamily_list]
    eu2_alignment <- alignment[names(alignment) %in% ingroup]
    # subsetted alignment with just EuII taxa having min seq length
    eu2_alignment <- eu2_alignment[sapply(eu2_alignment, length) > min_length ]
    # get list of families for subsetted alignment
    eu2_subfamilies <- onekp_data$Family_subfam[match(names(eu2_alignment), onekp_data$Code)]
    # check if all 12 eu2 subfamilies are included at least once in subsetted alignment
    all(eu2_subfamily_list %in% eu2_subfamilies)
  }

  # filter fasta files based on having at least one species per eupoly II subfamily with seq length within 1 sd of overall mean seq length
  # input: candidate bait fasta file
  # output: T/F whether the candidate passes the filter
  filter_eupoly2_genus <- function (alignment) {
    # get list of genera for each seq in alignment
    # full_genera <- onekp_data$genus[match(names(alignment), onekp_data$Code)]
    # minimum cutoff seq length is mean overall seq length - 1 sd
    min_length <- mean(as.numeric(lapply(alignment, length))) - sd(as.numeric(lapply(alignment, length)))
    # subsetted alignment with just EuII taxa
    #eu2_alignment <- alignment[full_genera %in% eu2_genus_list]
    eu2_alignment <- alignment[names(alignment) %in% ingroup]
    # subsetted alignment with just EuII taxa having min seq length
    eu2_alignment <- eu2_alignment[sapply(eu2_alignment, length) > min_length ]
    # get list of genera for subsetted alignment
    eu2_genera <- onekp_data$genus[match(names(eu2_alignment), onekp_data$Code)]
    # check if all 24 eu2 subfamilies are included at least once in subsetted alignment
    all(eu2_genus_list %in% eu2_genera)
    # length(which(eu2_genus_list %in% eu2_genera))
  }

  filter_xinping <- function (alignment_names, xinping_seqs) {
    clusters <- sapply(alignment_names, function (x) strsplit(x, split="_"))
    clusters_first <- sapply(clusters, function (x) x[1])
    clusters_second <- sapply(clusters, function (x) x[2])
    clusters <- paste(clusters_first, clusters_second, sep="_")
    clusters %in% xinping_seqs
  }

  if (filter_type == "min_taxa" && !is.na(min_taxa)) {
    # filter fasta files based on minimum number of taxa in ingroup
    count.in <- lapply(fasta.files, FUN = function (x) length(names(x)[names(x) %in% taxaset]))
    count.total <- lapply(fasta.files, FUN = function (x) length(names(x)))
    pass_filter <- count.in >= min_taxa
  } else if (filter_type == "min_taxa" && is.na(min_taxa)) {
    stop("need to provide minimum number of taxa for filtering")
  } else if (filter_type == "family") {
    pass_filter <- sapply(fasta.files, filter_eupoly2_family)
  } else if (filter_type == "subfamily") {
    pass_filter <- sapply(fasta.files, filter_eupoly2_subfamily)
  } else if (filter_type == "genus") {
    pass_filter <- sapply(fasta.files, filter_eupoly2_genus)
  } else if (filter_type == "xinping") {
    pass_filter <- filter_xinping (fasta.names, xinping_seqs)
  } else {
    stop ("need to choose valid filtering method")
  }

  fasta.filtered <- fasta.files[pass_filter]
  fasta.names.filtered <- fasta.names[pass_filter]

  results <- list (
    fasta.filtered = fasta.filtered,
    fasta.names.filtered = fasta.names.filtered,
    fasta.filtered.length = length(fasta.filtered),
    fasta.raw.length = length(fasta.files)
  )
  return(results)
}
