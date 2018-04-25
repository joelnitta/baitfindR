# pad_exons_carl25_genes.R

# Script to take Carl's alignments including both masked and unmasked introns,
# and add "padding" around short exons

# Input is "genomic" alignments including transcriptomes and masked and unmasked introns
# with ends that I have already trimmed manually
# (but intron masking hasn't been applied to other seqs yet)

library(ape)
library(ips)
library(xlsx)

setwd("/Volumes/Storage/thelypteris/working/yang_smith/rothfels_25_genes")

# get associated data for 1kp samples
onekp_data <- read.xlsx("/Volumes/Transcend/Documents/Projects/JSPS/one_kp/1kP-Sample-List.xlsx",
                        sheetIndex = 1,
                        colIndex = 1:18,
                        startRow=1,
                        endRow=84,
                        stringsAsFactors=FALSE)

outgroups <- onekp_data$Code[onekp_data$out_in == "OUT"]
outgroups <- outgroups[!(is.na(outgroups))]
ingroups <- onekp_data$Code[onekp_data$out_in == "IN"]
ingroups <- ingroups[!(is.na(ingroups))]

# function to only count length and position of exons/introns in final padded alignment
# choose to return warning about presence and position of short exons ("report"),
# actual lengths ("length"),
# or position ("position")
# can choose to trim results to just short exons or not
count_exon_intron <- function (alignment_file, return_type="length", trim=TRUE) {

  # read in manually trimmed alignment as DNA.bin matrix
  alignment <- read.dna(alignment_file, format="fasta", as.matrix=TRUE)

  ### get exon lengths
  # exons are columns in ingroup that are NOT all 'n'
  exons <- !(apply(as.character(alignment), 2, function (x) all(x == "n")))
  # convert from true/false vector to number of each column that is an exon
  exons_num <- which(exons)
  # handy trick to get length of regions that consecutively increase by 1
  # https://stackoverflow.com/questions/16118050/how-to-check-if-a-vector-contains-n-consecutive-numbers
  # turns out that in our case, the first item in the list is the exon lengths, and the second is intron lengths
  # (with some tweaking)
  exon_lengths <- unlist(rle(diff(exons_num))[1])
  intron_lengths <- unlist(rle(diff(exons_num))[2])
  exon_lengths <- exon_lengths[exon_lengths > 1]
  intron_lengths <- intron_lengths[intron_lengths > 1]
  exon_lengths <- exon_lengths+1
  intron_lengths <- intron_lengths-1
  names(exon_lengths) <- gsub("lengths", "exon", names(exon_lengths))
  names(intron_lengths) <- gsub("values", "intron", names(intron_lengths))

  # combine exon and intron lengths
  both_lengths <- c(exon_lengths, intron_lengths)
  # sort by number
  both_lengths <- both_lengths[order(as.numeric(gsub("exon|intron", "", names(both_lengths))))]

  # check for short exons
  short_exons <- list()
  for (i in 1:length(both_lengths)) {
    # check if it's an exon shorter than 120bp
    if (grepl("exon", names(both_lengths[i])) && both_lengths[i] < 120) {
      short_exons[[i]] <- names(both_lengths[i])
    } else {
      short_exons[[i]] <- NULL
    }
  }
  short_exons <- unlist(short_exons)

  positions <- cumsum(both_lengths)

  if (trim == TRUE) {
    if (length(short_exons) > 0) {
      positions <- positions[names(positions) %in% short_exons]
      both_lengths <- both_lengths[names(both_lengths) %in% short_exons]
    }
  }

  if (return_type == "length") {
    result <- both_lengths
  } else if (return_type == "report") {
    result <- short_exons
  } else if (return_type == "position") {
    result <- positions
  }

  return(result)
}

# function to pad short exons OR optionally return named vector of exon/intron lengths
pad_exons <- function (alignment_file, align_name="", return_lengths = FALSE, ref_masked_type="Azfi", ref_unmasked_type="WF01", extend_right_only="", extend_left_only="") {

  # convert extend right/left to empty if NA
  if(is.na(extend_right_only)) {extend_right_only <- ""}
  if(is.na(extend_left_only)) {extend_left_only <- ""}

  # read in manually trimmed alignment as DNA.bin matrix
  alignment <- read.dna(alignment_file, format="fasta", as.matrix=TRUE)

  ### subsetting
  # remove non-eupoly II ferns OTHER than ref seqs
  alignment <- alignment[grep(paste(outgroups,collapse="|"), rownames(alignment), invert=TRUE), ]

  # delete gap-only columns
  alignment <- deleteEmptyCells(alignment, nset="-")

  # make sure ref seqs present
  if ( any(grepl("Sacu|Azfi|WF01", rownames(alignment))) ) {

    # subset alignment without ref seqs
    alignment_noref <- alignment[grep("Sacu|Azfi|WF01", rownames(alignment), invert=TRUE), ]

    # extract ref seqs, one with and one without masked introns
    ref_seqs <- alignment[grep("Sacu|Azfi|WF01", rownames(alignment)), ]

    # extract masked or unmasked ref seq with search method depending on naming scheme of ref seq
    if (ref_masked_type == "Azfi" | ref_masked_type == "Sacu") {
      ref_masked <- ref_seqs[grep("Sacu|Azfi", rownames(ref_seqs)), ]
      ref_masked <- ref_masked[grep("unmasked", rownames(ref_masked), invert=TRUE), ]
    } else if (ref_masked_type == "WF01") {
      ref_masked <- ref_seqs[rownames(ref_seqs) == "WF01", ]
    }

    if (ref_unmasked_type == "Azfi" | ref_unmasked_type == "Sacu") {
      ref_unmasked <- ref_seqs[grep("unmasked", rownames(ref_seqs)), ]
    } else if (ref_unmasked_type == "WF01") {
      ref_unmasked <- ref_seqs[grep("WF01-", rownames(ref_seqs)), ]
    }

    rm(ref_seqs)

    ### mask introns
    # classify each character in ref-genome sequence as hard-masked intron ("n") or not
    introns <- apply(as.character(ref_masked), 2, function (x) x=="n")

    # also make list of all transcript characters that are 90% or greater gap-only
    percent_gap <- function (x) {length(which(as.character(x) == "-"))/length(x)}
    gaps <- apply(as.character(alignment_noref), 2, function (x) all(percent_gap(x) > 0.90))

    # only convert introns if all other characters are 95% gaps
    introns[which(gaps==FALSE)] <- FALSE

    # mask ingroup taxa introns (convert to character first)
    alignment <- as.character(alignment)
    # alignment[grep("unmasked", rownames(alignment), invert=TRUE), introns] <- "n"
    alignment[grep(paste(ingroups,collapse="|"), rownames(alignment)), introns] <- "n"

    # subset alignment without unmasked ref seq
    # alignment_no_unmasked <- alignment[grep("unmasked", rownames(alignment), invert=TRUE), ]
    alignment_no_unmasked <- alignment[grep(paste(ingroups,collapse="|"), rownames(alignment)), ]
    alignment_no_unmasked <- as.DNAbin(alignment_no_unmasked)

    # for troubleshooting, write out alignment
    # write.dna(alignment, file="rothfels_test_align.phy")

    ### get exon lengths
    # exons are columns in ingroup that are NOT all 'n'
    exons <- !(apply(as.character(alignment_no_unmasked), 2, function (x) all(x == "n")))
    # convert from true/false vector to number of each column that is an exon
    exons_num <- which(exons)
    # handy trick to get length of regions that consecutively increase by 1
    # https://stackoverflow.com/questions/16118050/how-to-check-if-a-vector-contains-n-consecutive-numbers
    # turns out that in our case, the first item in the list is the exon lengths, and the second is intron lengths
    # (with some tweaking)
    exon_lengths <- unlist(rle(diff(exons_num))[1])
    intron_lengths <- unlist(rle(diff(exons_num))[2])
    exon_lengths <- exon_lengths[exon_lengths > 1]
    intron_lengths <- intron_lengths[intron_lengths > 1]
    exon_lengths <- exon_lengths+1
    intron_lengths <- intron_lengths-1
    names(exon_lengths) <- gsub("lengths", "exon", names(exon_lengths))
    names(intron_lengths) <- gsub("values", "intron", names(intron_lengths))

    # combine exon and intron lengths
    both_lengths <- c(exon_lengths, intron_lengths)
    # sort by number
    both_lengths <- both_lengths[order(as.numeric(gsub("exon|intron", "", names(both_lengths))))]

    # define end positions of each exon and intron
    end_pos <- cumsum(both_lengths)

    # check if padding is needed
    if (any(exon_lengths < 120)) {

      padding_all <- list()

      for (i in 1:nrow(alignment_no_unmasked)) {
        ### pad introns: make vector of bases to pad
        # set up variables for loop to define sequences for padding
        padding <- list()
        pad_right <- list()
        pad_left <- list()

        start_right <- NULL
        end_right <- NULL
        start_left <- NULL
        end_left <- NULL

        pad_seqs_right <- NULL
        pad_seqs_left <- NULL
        gaps <- NULL
        pad_seqs <- list()

        # loop over all exons and introns
        for (j in 1:length(both_lengths)) {

          # check if it's an exon shorter than 120bp
          if (grepl("exon", names(both_lengths[j])) && both_lengths[j] < 120) {

            # check that there are no gaps in original exon sequence (DON'T pad if there are gaps)
            if (!any(grepl("-", as.character(alignment_no_unmasked[i,ifelse(j==1,1,end_pos[j-1]+1):end_pos[j]])))) {

              # check if want to extend on only one side or the other (extend_right_only and extend_left_only are manually set vector of exons to only extend L or R)
              if (j %in% extend_right_only) {
                expand <- 124-both_lengths[j]
                start_left <- NULL
                end_left <- NULL
                start_right <- ifelse ((end_pos[j]+1) < ncol(ref_masked), (end_pos[j]+1), ncol(ref_masked))
                end_right <- ifelse ((end_pos[j] + expand) > ncol(ref_masked), ncol(ref_masked), (end_pos[j] + expand))
              } else if (j %in% extend_left_only) {
                expand <- 124-both_lengths[j]
                start_right <- NULL
                end_right <- NULL
                if (j == 1) {
                  start_left <- NULL
                  end_left <- NULL
                } else if (j > 1) {
                  start_left <- ifelse ((end_pos[j-1]-expand) < 1, 1, (end_pos[j-1]-expand))
                  end_left <- end_pos[j-1]
                }

              # otherwise, extend on both sides
              } else {
                # number of bases to pad left and right: make it a little longer than 120 total just to be safe
                expand <- ceiling((124-both_lengths[j])/2)

                # define starting and ending points
                # on left side, starting position uses the previous end position, so must be greater than 1
                if (j == 1) {
                  start_left <- NULL
                  end_left <- NULL
                } else if (j > 1) {
                  start_left <- ifelse ((end_pos[j-1]-expand) < 1, 1, (end_pos[j-1]-expand))
                  end_left <- end_pos[j-1]
                }
                # on right side take into account max length of ref seq (e.g., if end point is longer than ref seq, use max ref seq length )
                start_right <- ifelse ((end_pos[j]+1) < ncol(ref_masked), (end_pos[j]+1), ncol(ref_masked))
                end_right <- ifelse ((end_pos[j] + expand) > ncol(ref_masked), ncol(ref_masked), (end_pos[j] + expand))
              }

              # vectors of bases to pad left and right (check for null values first)
              if (length(start_right) == 0) {
                pad_right <- NULL
              } else {
                  pad_right <- start_right:end_right
              }

              if (length(start_left) == 0) {
                pad_left <- NULL
              } else {
                pad_left <- start_left:end_left
              }

              # combine these into a unique, sorted vector: this is the padding for this sample
              padding[[j]] <- as.numeric(sort(unique(c(pad_right, pad_left))))
            }
            padding_all[[i]] <- padding
          }
        }
      }

      padding_all <- lapply(padding_all, unlist)

      ### fill-in padding
      # convert to character
      alignment_no_unmasked <- as.character(alignment_no_unmasked)
      ref_unmasked <- as.character(ref_unmasked)

      # loop through rows of alignment, if a base is in the padding group, replace with ref seq intron
      for (i in 1:nrow(alignment_no_unmasked)) {
        # check that there are sequences to pad
        if (length(padding_all[[i]]) > 0) {
          # for each column, check if that column is in padding
          for (j in 1:ncol(alignment_no_unmasked)) {
            if (j %in% padding_all[[i]]) {
              # replace that base the with ref seq
              # alignment_no_unmasked[i,j] <- ref_unmasked[1,j]
              alignment_no_unmasked[i,j] <- ref_unmasked[j]
            }
          }
        }
      }
      # test
      # write.dna(alignment_no_unmasked, colsep="", indent=0, blocksep=0, format = "fasta", file="testpad2.aln")
    }
  }

  ### final clean up
  # final "alignment" is now the padded alignment with ref seqs removed
  alignment <- alignment_no_unmasked

  # delete gap-only columns
  if (class(alignment) == "matrix") {
    alignment <- as.DNAbin(alignment)
  }
  alignment <- deleteEmptyCells(alignment, nset="-")

  # add name of cluster to taxa names in alignment (will fix any weird characters later)
  rownames(alignment) <- paste(rownames(alignment), align_name, sep="_")

  # choose to return list of intron/exons lengths or actual alignment
  if (return_lengths == TRUE) {
    result <- both_lengths
  } else {
    result <- alignment
  }

  return(result)
}

# NEW VERSION
# read in sampling info for baits
bait_files_data <- read.csv("baits_new_sampling.csv", stringsAsFactors=FALSE)

# loop through files, pad each alignment, and write it out
padded_alignment <- NULL
for (i in 1:nrow(bait_files_data)) {
  padded_alignment <- pad_exons(
    alignment_file = paste0("baits_new/", bait_files_data$file[i]),
    align_name = bait_files_data$region[i],
    ref_masked_type = bait_files_data$ref_masked[i],
    ref_unmasked_type = bait_files_data$ref_unmasked[i],
    extend_right_only=bait_files_data$pad_right_only[i],
    extend_left_only=bait_files_data$pad_left_only[i])

  write.dna(padded_alignment, colsep="", indent=0, blocksep=0, format = "fasta", file=paste0("baits_new/", bait_files_data$region[i], ".padded.aln") )
}
#
# # for testing only
# alignment_file = paste0("baits_new/", bait_files_data$file[i])
# align_name = bait_files_data$region[i]
# ref_masked_type = bait_files_data$ref_masked[i]
# ref_unmasked_type = bait_files_data$ref_unmasked[i]
# extend_right_only=bait_files_data$pad_right_only[i]
# extend_left_only=bait_files_data$pad_left_only[i]

# check for short exons (optionally, in final REVISED bait set)
setwd("/Volumes/Storage/thelypteris/working/yang_smith/")
short_exon_check <- list()
for (i in 1:nrow(bait_files_data)) {
  #short_exon_check[[i]] <- count_exon_intron(alignment_file = paste0("baits_new/", bait_files_data$region[i], ".padded.aln"), return_type="report")
  short_exon_check[[i]] <- count_exon_intron(alignment_file = paste0("final_baits_REVISED/", bait_files_data$region[i], ".aln"), return_type="report")
  }
names(short_exon_check) <- bait_files_data$region
short_exon_check_keep <- short_exon_check [which(lapply(short_exon_check, length) > 0)]

# check position of short exons
short_exon_position <- list()
for (i in 1:nrow(bait_files_data)) {
  # short_exon_position[[i]] <- count_exon_intron(alignment_file = paste0("baits_new/", bait_files_data$region[i], ".padded.aln"), return_type="position")
  short_exon_position[[i]] <- count_exon_intron(alignment_file = paste0("final_baits_REVISED/", bait_files_data$region[i], ".aln"), return_type="position")
  }
names(short_exon_position) <- bait_files_data$region
short_exon_position <- short_exon_position [which(lapply(short_exon_check, length) > 0)]

short_exon_check_keep
short_exon_position

# function to trim fasta seqs to just those in ingroup and align with MAFFT

align_eu2 <- function (alignment_file, ingroup_select) {

  # read in manually trimmed alignment as DNA.bin matrix
  alignment <- read.dna(alignment_file, format="fasta", as.matrix=FALSE)

  ### subsetting
  # remove non-eupoly II ferns OTHER than ref seqs
  alignment <- alignment[grep(paste(ingroup_select,collapse="|"), names(alignment)) ]
  alignment <- mafft(alignment)

  return(alignment)

}

# list all files in baits folder
bait_files <- list.files("baits_new")
# then select only unaligned sequences
bait_files <- bait_files[(grep("supercontig.fasta", bait_files, fixed=TRUE))]

# add hyb-piper data to ingroups
# all four eu-2 seqs, or just Thelytperis
#ingroups <- c(ingroups, "WF01", "WF20", "WF29", "WF30")
ingroups <- c(ingroups, "WF01", "Azfi", "Sacu")
# or try one at a time
ingroups.list <- list(
  c(ingroups, "WF01"),
  c(ingroups, "WF20"),
  c(ingroups, "WF29"),
  c(ingroups, "WF30")
)

# loop through files, make alignment, and write it out
new_alignment <- NULL
for (i in 1:length(bait_files)) {
  new_alignment <- align_eu2(paste("baits_new/", bait_files[i], sep=""), ingroups)
  write.dna(new_alignment, colsep="", indent=0, blocksep=0, format = "fasta", file=paste0("baits_new/", gsub("fasta", "ref.aln", bait_files[i])))
}

# or, just choose one region but align for each eupoly II hyb piper sequence
bait_file <- bait_files[grep("TPLATE", bait_files)]
new_name <- gsub("fasta", "", bait_file)

for (i in 1:length(ingroups.list)) {
  new_alignment <- align_eu2(paste("baits_new/", bait_file, sep=""), ingroups.list[[i]])
  write.dna(new_alignment, colsep="", indent=0, blocksep=0, format = "fasta", file=paste0("baits_new/", new_name, "euII.", i, ".aln") )
}



# # OLD VERSION
# # list all files in baits folder
# bait_files <- list.files("baits")
#
# # keep those that I manually trimmed and those that didn't need trimming (keep)
# keep <- c("NDUFS6.genome.unmasked.aln", "Hemera.genome.unmasked.aln")
# bait_files <- bait_files[c(which(bait_files %in% keep), grep("trim", bait_files))]
# bait_files <- bait_files[-(grep("meta", bait_files))]
#
# # names is just name of each gene
# names <- gsub(".genome.unmasked.trim.aln", "", bait_files, fixed=TRUE)
# names <- gsub(".genome.unmasked.normal_mafft.trim.aln", "", names, fixed=TRUE)
# names <- gsub(".genome.unmasked.aln", "", names, fixed=TRUE)
#
# # loop through files, pad each alignment, and write it out
# padded_alignment <- NULL
# for (i in 1:length(bait_files)) {
#     padded_alignment <- pad_exons(paste("baits/", bait_files[i], sep=""), names[i])
#     write.dna(padded_alignment, colsep="", indent=0, blocksep=0, format = "fasta", file=paste0("baits/", names[i], ".padded.aln") )
# }
#
# # optional: check on exon/intron lengths in each alignment
# # list all files in original baits folder
# bait_files <- list.files("baits")
# keep <- c("NDUFS6.genome.unmasked.aln", "Hemera.genome.unmasked.aln")
# bait_files <- bait_files[c(which(bait_files %in% keep), grep("trim", bait_files))]
# bait_files <- bait_files[-(grep("meta", bait_files))]
# bait_files <- paste0("baits/", bait_files)
# # names is just name of each gene
# bait_names <- gsub(".genome.unmasked.trim.aln", "", bait_files, fixed=TRUE)
# bait_names <- gsub(".genome.unmasked.normal_mafft.trim.aln", "", bait_names, fixed=TRUE)
# bait_names <- gsub(".genome.unmasked.aln", "", bait_names, fixed=TRUE)
# bait_names <- gsub("baits/", "", bait_names, fixed=TRUE)
# old_bait_lengths <- lapply(bait_files, pad_exons, return_lengths=TRUE)
# names(old_bait_lengths) <- bait_names
#
# # list all files in new baits folder
# bait_files <- list.files("baits_new")
# bait_files <- bait_files[grep("supercontig.masked.clean.aln", bait_files, fixed=TRUE)]
# bait_files <- paste0("baits_new/", bait_files)
#
# # names is just name of each gene
# bait_names <- gsub(".genome.supercontig.masked.clean.aln", "", bait_files, fixed=TRUE)
# bait_names <- gsub("baits_new/", "", bait_names, fixed=TRUE)
#
# new_bait_lengths <- lapply(bait_files, pad_exons, return_lengths=TRUE, ref_seq = "WF01")
# names(new_bait_lengths) <- bait_names

