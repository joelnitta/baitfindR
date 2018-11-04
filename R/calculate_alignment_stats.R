#' Calculate summary statistics for an alignment.
#'
#' Including the original alignment in the output with \code{include_aln}
#' can be useful for mapping \code{calculate_alignment_stats} over a list
#' of alignments with \code{\link[purrr]{map_df}} to sort and filter
#' alignments by their summary statistics.
#'
#' @param alignment Input alignment; must be a matrix of class "DNAbin".
#' @param cutoff Numeric value indicating minimum exon length (optional);
#' flag this alignment if any/all exons are less than the cutoff length.
#' @param cutoff_any Logical; Should the alignment
#' be flagged if any exons are shorter than the cutoff? The default, FALSE,
#' means that the alignment will only be flagged if all exons are shorter
#' than the cutoff value.
#' @param include_aln Logical; Should the original alignment
#' be included in the output list?
#' @return A list including the following summary statistics: \describe{
#'   \item{intron_lengths}{List including vector of intron lengths}
#'   \item{exon_lengths}{List including vector of exon lengths}
#'   \item{num_introns}{Number of introns}
#'   \item{num_exons}{Number of introns}
#'   \item{mean_dist}{Mean genetic distance between sequences in alignment}
#'   \item{max_dist}{Maximum genetic distance between sequences in alignment}
#'   \item{GC_content}{Total \%GC content}
#'   \item{pars_inf}{Fraction of sites that are parsimony informative}
#'   \item{total_exon_length}{Total exon length}
#'   \item{less_than_cutoff}{Logical flag indicating whether alignment passed
#'   the minimum exon length cutoff or not}
#'   \item{alignment}{The original input alignment}
#' }
#' @export
calculate_alignment_stats <-
  function (alignment, cutoff = 120, cutoff_any = FALSE, include_aln = FALSE) {

    ### Error-checking
    assertthat::assert_that(any("DNAbin" %in% class(alignment)),
                            msg = "alignment must be of class 'DNAbin'")

    assertthat::assert_that(is.matrix(alignment),
                            msg = "alignment must be matrix")
    ### Setup
    # exon columns are columns that are NOT all 'n'
    exon_cols <- !(apply(as.character(alignment), 2, function (x) all(x == "n")))

    # intron columns are columns that ARE all 'n'
    intron_cols <- apply(as.character(alignment), 2, function (x) all(x == "n"))

    # make exon-only alignment to calculate various stats
    alignment_exons <- alignment[,exon_cols]

    ### Simple stats
    # Calculate mean DNA distance
    # (missing values ignored, so introns don't matter)
    mean_dist <- mean(ape::dist.dna(alignment, model="raw",
                                    pairwise.deletion=TRUE), na.rm=TRUE)
    max_dist <- max(ape::dist.dna(alignment, model="raw",
                                  pairwise.deletion=TRUE), na.rm=TRUE)

    # Calculate % pars. inf. chars
    # (introns DO matter, so use exon-only alignment)
    pars_inf <- ips::pis(alignment_exons, what="frac")

    # Calculate %GC content (missing values ignored, so introns don't matter)
    GC_content <- ape::GC.content(alignment)

    ### Count number and length of exons, introns
    # Convert from true/false vector to number of each column that is an exon
    exons <- which(exon_cols)
    introns <- which(intron_cols)

    # Get length of regions that consecutively increase by 1
    # https://stackoverflow.com/questions/16118050/how-to-check-if-a-vector-contains-n-consecutive-numbers
    exon_lengths <- rle(diff(exons))[["lengths"]]
    # also returns a bunch of 1s; get rid of these
    exon_lengths <- exon_lengths[exon_lengths > 1]

    # Do same for introns
    intron_lengths <- rle(diff(introns))[["lengths"]]
    intron_lengths <- intron_lengths[intron_lengths > 1]

    # Use this to calculate total number of exons and introns
    num_introns <- length(intron_lengths)
    num_exons <- length(exon_lengths)

    ### Calculate total length including only exons
    # (since don't know what intron length will be in target seq)
    total_exon_length <- sum(exon_lengths)

    ### Cutoff flag
    # Optionally flag gene if any (or all) exons are less than cutoff length
    if (is.null(cutoff) | is.null(cutoff_any)) less_than_cutoff <- NA

    if (cutoff_any == TRUE) {
      if ((any(exon_lengths < cutoff))) {
        less_than_cutoff <- TRUE
      } else {
        less_than_cutoff <- FALSE
      }
    } else if (cutoff_any == FALSE) {
      if (!(any(exon_lengths > cutoff))) {
        less_than_cutoff <- TRUE
      } else {
        less_than_cutoff <- FALSE
      }
    }

    results <- list(
      intron_lengths = list(intron_lengths),
      exon_lengths = list(exon_lengths),
      num_introns = num_introns,
      num_exons = num_exons,
      mean_dist = mean_dist,
      max_dist = max_dist,
      GC_content = GC_content,
      pars_inf = pars_inf,
      total_exon_length = total_exon_length,
      less_than_cutoff = less_than_cutoff,
      alignment = list(alignment)
    )

    if (!isTRUE(include_aln)) {
      results[["alignment"]] <- NULL
    }

    return(results)

  }
