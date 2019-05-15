
#' Trim columns of exclusively "n" characters from one end of alignment
#'
#' @param alignment matrix of class DNAbin
#'
#' @return matrix of class DNAbin
#' @examples
#' library(ape)
#' dna_mat <- matrix(rep(c("n", "t", "a", "n"), 2), 2, 4, byrow = TRUE)
#' dna_mat
#' dna_mat <- as.DNAbin(dna_mat)
#' dna_mat
#' trim_one_end(dna_mat)
#' @export
trim_one_end <- function (alignment) {

  assertthat::assert_that(
    isTRUE(inherits (alignment, "DNAbin")),
    msg = "alignment must be of class DNAbin")

  assertthat::assert_that(is.matrix(alignment))

  # convert to character
  alignment <- as.character(alignment)

  # prepare loop: to_trim is list of TRUE/FALSE for sequences to trim
  # at_ends starts as TRUE, becomes false at
  # first instance of non-"n" character
  to_trim <- list()
  at_ends <- TRUE

  # start at one end of alignment, assign TRUE to each column if
  # still at "end" and all columns are "n"
  for (i in 1:ncol(alignment)) {
    if (all(grepl("n", alignment[, i])) && at_ends == TRUE) {
      to_trim[i] <- TRUE
    } else {
      to_trim[i] <- FALSE
      at_ends <- FALSE
    }
  }
  to_trim <- unlist(to_trim)
  alignment <- alignment[, !(to_trim)]
  as.DNAbin(alignment)
}

#' Trim both ends of alignment
#'
#' @param alignment matrix of class DNAbin
#'
#' @return matrix of class DNAbin
#' @examples
#' library(ape)
#' dna_mat <- matrix(rep(c("n", "t", "a", "n"), 2), 2, 4, byrow = TRUE)
#' dna_mat <- as.DNAbin(dna_mat)
#' dna_mat
#' trim_both_ends(dna_mat)
#' @export
trim_both_ends <- function (alignment) {
  alignment %>%
  trim_one_end() %>%
  ape::complement() %>%
  trim_one_end() %>%
  ape::complement()
}

#' Fill-in introns in an alignment
#'
#' The alignment must be produced by aligning sequence to a reference with
#' introns masked by 'n's, which results in gaps in the rest of the sequences.
#' This function replaces those gaps with 'n's and removes the reference
#' sequence.
#'
#' @param alignment matrix of class DNAbin
#' @param outgroup Character vector; names of outgroup sequences.
#' @param ref_pattern Pattern used for matching with grep to identify reference
#' sequences.
#' @param trim_outgroup Logical; should the outgroups be trimmed from the
#' alignment?
#'
#' @return Matrix of class DNAbin
#'
#' @examples
#' library(ape)
#' data(woodmouse)
#'
#' # Make reference sequence with a 50bp intron (a string of 'n's)
#' # in the middle.
#' woodmouse_ref <- as.character(woodmouse[1,])
#' woodmouse_ref <- c(woodmouse_ref[1:400],
#'   rep("n", 50),
#'   woodmouse_ref[401:length(woodmouse_ref)])
#' woodmouse_ref <- as.DNAbin(woodmouse_ref)
#'
#' # Align with other sequences that don't include the intron.
#' # (need to convert to list first)
#' woodmouse_ref <- as.list(woodmouse_ref)
#' names(woodmouse_ref) <- "ref"
#' woodmouse <- as.list(woodmouse)
#' woodmouse_with_introns <- ips::mafft(
#'   c(woodmouse, woodmouse_ref),
#'   path = "/usr/bin/mafft")
#'
#' # Image of the alignment shows that 'ref' has 'n's at positions 400-450,
#' # while other sequences have gaps ('-').
#' image(woodmouse_with_introns)
#'
#' # Fill-in introns
#' woodmouse_masked <- fill_introns(
#'   woodmouse_with_introns,
#'   ref_pattern = "ref"
#' )
#'
#' # After filling-in, the reference sequence is gone and
#' # all introns are 'n's.
#' image(woodmouse_masked)
#'
#' @export
fill_introns <- function (alignment, ref_pattern, outgroup = NULL, trim_outgroup = FALSE) {

  # Check input
  assertthat::assert_that(
    isTRUE(inherits (alignment, "DNAbin")),
    msg = "alignment must be of class DNAbin")

  assertthat::assert_that(is.character(outgroup) | is.null(outgroup))
  assertthat::assert_that(assertthat::is.string(ref_pattern))
  assertthat::assert_that(is.logical(trim_outgroup))
  assertthat::assert_that(is.matrix(alignment))

  # Extract the "reference genome" sequence, and remove from alignment
  # all transcriptome sequences are 4-letter codes.
  # Sacu = Salvinia, Azfi = Azolla, numbers 1-5 are Arabidopsis (numbered by chromosome)
  assertthat::assert_that(
    any(stringr::str_detect(rownames(alignment), ref_pattern)),
    msg = "Reference sequence not detected in alignment")

  ref_seq <- alignment[grep(ref_pattern, rownames(alignment)), ]

  assertthat::assert_that(
    nrow(ref_seq) == 1,
    msg = "Multiple reference sequences detected"
  )

  alignment <- alignment[grep(ref_pattern, rownames(alignment), invert=TRUE), ]

  # Rev-comp alignments if MAFFT flipped all of them
  if (all(grepl("_R_", rownames(alignment)))) {
    alignment <- ape::complement(alignment)
    rownames(alignment) <- gsub("_R_", "", rownames(alignment))
  }

  # Classify each character in ref-genome sequence as hard-masked intron ("n") or not
  introns <- apply(as.character(ref_seq), 2, function (x) x=="n"|x=="N")

  # Also make list of all transcript characters that are gap-only
  gaps <- apply(as.character(alignment), 2, function (x) all(x == "-"))

  # Only convert introns if all other characters are gaps
  introns[which(gaps==FALSE)] <- FALSE

  # Mask introns (convert to character first)
  alignment <- as.character(alignment)
  alignment[,introns] <- "n"

  # Remove outgroup taxa
  if(isTRUE(trim_outgroup)) alignment <- alignment[!(rownames(alignment) %in% outgroup), ]

  # Convert back to DNAbin for deleteEmptyCells
  alignment <- as.DNAbin(alignment)

  # Delete columns that are only empty cells
  alignment <- ips::deleteEmptyCells(alignment, nset="-")

  # Trim ends if they are only introns (all columns "n")
  trim_both_ends(alignment)

}
