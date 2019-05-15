
# trim columns of exclusively "n" characters from one end of alignment
trim_one_end <- function (alignment) {
  # prepare loop: to.trim is list of TRUE/FALSE for sequences to trim
  # at_ends starts as TRUE, becomes false at first instance of non-"n" character
  to.trim <- list()
  at_ends <- TRUE

  # start at one end of alignment, assign TRUE to each column if still at "end" and all columns are "n"
  for (i in 1:ncol(alignment)) {
    if (all(grepl("n", alignment[, i])) && at_ends == TRUE) {
      to.trim[i] <- TRUE
    } else {
      to.trim[i] <- FALSE
      at_ends <- FALSE
    }
  }
  to.trim <- unlist(to.trim)
  alignment <- alignment[, !(to.trim)]
  return(alignment)
}

# custom func that uses trim_one_end on both ends of alignment
trim_both_ends <- function (alignment) {
  require(ape)
  alignment <- trim_one_end(alignment)
  alignment <- ape::complement(alignment)
  alignment <- trim_one_end(alignment)
  alignment <- ape::complement(alignment)
  return(alignment)
}

# fill introns in alignment based on reference sequence,
# then trim that reference sequence and the outgroup
# from the alignment
fill_introns <- function (alignment, outgroup) {

  # extract the "reference genome" sequence, and remove from alignment
  # all transcriptome sequences are 4-letter codes.
  # Sacu = Salvinia, Azfi = Azolla, numbers 1-5 are Arabidopsis (numbered by chromosome)
  ref_seq <- alignment[grep("Sacu|Azfi|1|2|3|4|5", rownames(alignment)), ]
  alignment <- alignment[grep("Sacu|Azfi|1|2|3|4|5", rownames(alignment), invert=TRUE), ]

  # rev-comp alignments if MAFFT flipped all of them
  if (all(grepl("_R_", rownames(alignment)))) {
    alignment <- ape::complement(alignment)
    rownames(alignment) <- gsub("_R_", "", rownames(alignment))
  }

  # classify each character in ref-genome sequence as hard-masked intron ("n") or not
  introns <- apply(as.character(ref_seq), 2, function (x) x=="n"|x=="N")

  # also make list of all transcript characters that are gap-only
  gaps <- apply(as.character(alignment), 2, function (x) all(x == "-"))

  # only convert introns if all other characters are gaps
  introns[which(gaps==FALSE)] <- FALSE

  # mask introns (convert to character first)
  alignment <- as.character(alignment)
  alignment[,introns] <- "n"

  # remove outgroup taxa
  alignment <- alignment[!(rownames(alignment) %in% outgroup), ]

  # convert back to DNAbin for deleteEmptyCells
  alignment <- as.DNAbin(alignment)

  # delete columns that are only empty cells
  alignment <- ips::deleteEmptyCells(alignment, nset="-")

  # trim ends if they are only introns (all columns "n")
  alignment <- trim_both_ends(alignment)

  return(alignment)
}

# loop-ize fill_introns for drake workflow
fill_introns_loop <- function (alignment_list, outgroup, ...) {
  map(alignment_list, fill_introns, outgroup = outgroup)
}
