# bedtools functions ------------------------------------------------------

#' Extract regions from a fasta file
#'
#' Wrapper for bedtools getfasta.
#'
#' @param bed_file Path to bed file with locations of regions to extract.
#' bed file is a tab-separated file with columns for chromosome (e.g., chr1),
#' start position (e.g., 1), and end position (e.g., 10), in that order.
#' No column headers are used.
#' @param fasta_file Path to file in fasta format to extract regions from.
#' @param out_fasta_file Path to write extracted regions (in fasta format).
#' @param ...
#' @return List; output of processx::run(). Externally, a fasta file will be
#' written to the path specified by `out_fasta_file`.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @examples
#' extract_regions_from_fasta(
#'   bed_file = "temp_dir/test_genes"
#'   fasta_file = "data_raw/Arabidopsis_thaliana.TAIR10.dna.toplevel.renamed.fasta"
#'   out_fasta_file = "temp_dir/test_masked_genes"
#' )
#' @export
extract_regions_from_fasta <- function (bed_file, fasta_file, out_fasta_file, ...) {

  # Check input
  assertthat::assert_that(assertthat::is.string(bed_file))
  bed_file <- fs::path_abs(bed_file)
  assertthat::assert_that(assertthat::is.readable(bed_file))

  assertthat::assert_that(assertthat::is.string(fasta_file))
  fasta_file <- fs::path_abs(fasta_file)
  assertthat::assert_that(assertthat::is.readable(fasta_file))

  assertthat::assert_that(assertthat::is.string(out_fasta_file))
  out_fasta_file <- fs::path_abs(out_fasta_file)
  assertthat::assert_that(assertthat::is.dir(fs::path_dir(out_fasta_file)))

  bed <- readr::read_tsv(bed_file, col_names = c("chr", "start", "end"), col_types = "cdd")

  checkr::check_data(
    bed,
    values = list(
      chr = "a",
      start = 1,
      end = 1
    ),
    order = TRUE)

  if (!(check_bed_genome_names(fasta_file, bed))) {
    stop ("Names don't match between bed file and fasta file headers")
  }

  # Run bedtools getfasta
  processx::run(
    command = "bedtools",
    args = c(
      "getfasta",
      "-fi", fasta_file,
      "-bed", bed_file,
      "-fo", out_fasta_file
    ),
    echo = TRUE
  )

}

#' Mask regions in a fasta file.
#'
#' Wrapper for bedtools maskfasta.
#'
#' All regions of the `fasta_file` specified by the `bed_file` will be
#' replaced ("hard-masked") with 'N's.
#'
#' The bed file is a tab-separated file with columns for chromosome (e.g., chr1),
#' start position (e.g., 1), and end position (e.g., 10), in that order.
#' No column headers are used.
#'
#' @param bed_file Path to bed file with locations of regions to mask.
#' @param fasta_file Path to unmasked fasta file.
#' @param out_fasta_file Path to write masked fasta file.
#' @param ...
#' @return List; output of processx::run(). Externally, a fasta file will be
#' written to the path specified by `out_fasta_file`.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @examples
#' # First write genes, introns, and exons out as tsv files
#'
#' dir.create("temp_dir")
#' find_bed_regions(
#'   gff3_file = "data_raw/Arabidopsis_thaliana.TAIR10.40.gff3.gz",
#'   source_select = "araport11",
#'   out_type = "write_all",
#'   out_dir = "temp_dir",
#'   prefix = "test"
#' )
#'
#' # Now mask the genome, using the bed file and genome fasta file.
#' mask_genome(
#'   bed_file = "temp_dir/test_introns",
#'   fasta_file = "data_raw/Arabidopsis_thaliana.TAIR10.dna.toplevel.renamed.fasta",
#'   out_fasta_file = "temp_dir/test_masked"
#' )
#' @export
mask_regions_in_fasta <- function (bed_file, fasta_file, out_fasta_file, ...) {

  # Check input
  assertthat::assert_that(assertthat::is.string(bed_file))
  bed_file <- fs::path_abs(bed_file)
  assertthat::assert_that(assertthat::is.readable(bed_file))

  assertthat::assert_that(assertthat::is.string(fasta_file))
  fasta_file <- fs::path_abs(fasta_file)
  assertthat::assert_that(assertthat::is.readable(fasta_file))

  assertthat::assert_that(assertthat::is.string(out_fasta_file))
  out_fasta_file <- fs::path_abs(out_fasta_file)
  assertthat::assert_that(assertthat::is.dir(fs::path_dir(out_fasta_file)))

  bed <- readr::read_tsv(bed_file, col_names = c("chr", "start", "end"), col_types = "cdd")

  checkr::check_data(
    bed,
    values = list(
      chr = "a",
      start = 1,
      end = 1
    ),
    order = TRUE)

  if (!(check_bed_genome_names(fasta_file, bed))) {
    stop ("Names don't match between bed file and fasta file headers")
  }

  # Run bedtools maskfasta
  processx::run(
    command = "bedtools",
    args = c(
      "maskfasta",
      "-fi", fasta_file,
      "-bed", bed_file,
      "-fo", out_fasta_file
    ),
    echo = TRUE
  )

}

#' Clean up data from a gff file and
#' convert to bed format
#'
#' Helper function for `find_bed_regions`. Merges overlapping regions and sorts
#' regions.
#'
#' @param region Dataframe; list of gene regions in "bed" format. Must include
#' the following columns in order: `chr` ('chromosome', character), `start`
#' (start position, numeric), and `end` (end position, numeric).
#' @param check.chr Logical; should coordinates be checked for chromosomal
#' format with "chr" prefix?
#' @param verbose Logical; should `bedr` functions output all messages?
#'
#' @return Dataframe in "bed" format.
#'
clean_gff <- function (region, check.chr = FALSE, verbose = FALSE) {

  # Check input format
  assertthat::assert_that(
    bedr:::determine.input(region, verbose = verbose) == 1,
    msg = "Input must be dataframe in 'bed' format")

  region %>%
    # Check if region is valid
    dplyr::filter(bedr::is.valid.region(., check.chr = check.chr, verbose = verbose)) %>%
    # Collapse overlapping regions
    bedr::bedr.merge.region(check.chr = check.chr, verbose = verbose) %>%
    # Sort regions
    bedr::bedr.sort.region(check.chr = check.chr, verbose = verbose) %>%
    # Convert to bed format
    bedr::convert2bed(check.chr = check.chr, verbose = verbose)
}

#' Find genes, exons, and introns in a gff3 file
#'
#' If tsv files are written out by selecting "write_all" for `out_type`,
#' they will overwrite any existing files with the same name in `out_dir`.
#'
#' @param gff3_file Path to input file in `gff3` format.
#' @param source_select Character vector; only use regions from these
#' sources. Must match values in `source` column of gff3 file. Optional.
#' @param gene_label String; value used to indicate genes in gff3 file.
#' Must match at least one value in `type` column of gff3 file. Default "gene".
#' @param exon_label String; value used to indicate exons in gff3 file.
#' Must match at least one value in `type` column of gff3 file. Default "exon".
#' @param verbose Logical; should `bedr` functions output all messages?
#' @param out_type Type of output to return:
#' "genes": dataframe in "bed" format of genes.
#' "introns": dataframe in "bed" format of introns.
#' "exons": dataframe in "bed" format of exons.
#' "write_all": write tab-separated files for each of `genes`, `introns`, and
#' `exons` to `out_dir`. The hash digest of the combined genes, introns, and
#' exons will be returned.
#' @param prefix String; prefix to attach to tsv files if `out_type` is
#' "write_all".
#' @param out_dir Directory to write tsv files if `out_type` is "write_all".
#' @param  ... Other arguments. Not used by this function, but meant to
#' be used by \code{\link[drake]{drake_plan}} for tracking during workflows.
#' @return Dataframe or character.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @examples
#' # Find genes
#' genes <- find_bed_regions(
#'   gff3_file = "/home/rstudio/baitfindR_simple/data_raw/Arabidopsis_thaliana.TAIR10.40.gff3.gz",
#'   source_select = "araport11",
#'   out_type = "genes"
#' )
#' tibble::as_tibble(genes)
#'
#' # Find introns
#' introns <- find_bed_regions(
#'   gff3_file = "/home/rstudio/baitfindR_simple/data_raw/Arabidopsis_thaliana.TAIR10.40.gff3.gz",
#'   source_select = "araport11",
#'   out_type = "introns"
#' )
#' tibble::as_tibble(introns)
#'
#' # Find exons
#' exons <- find_bed_regions(
#'   gff3_file = "/home/rstudio/baitfindR_simple/data_raw/Arabidopsis_thaliana.TAIR10.40.gff3.gz",
#'   source_select = "araport11",
#'   out_type = "exons"
#' )
#' tibble::as_tibble(exons)
#'
#' # Write genes, introns, and exons out as tsv files
#' dir.create("temp_dir")
#' find_bed_regions(
#'   gff3_file = "/home/rstudio/baitfindR_simple/data_raw/Arabidopsis_thaliana.TAIR10.40.gff3.gz",
#'   source_select = "araport11",
#'   out_type = "write_all",
#'   out_dir = "temp_dir",
#'   prefix = "test"
#' )
#' @export
find_bed_regions <- function (gff3_file,
                              source_select = NULL,
                              gene_label = "gene", exon_label = "exon",
                              verbose = FALSE,
                              prefix = NULL, out_dir = NULL,
                              out_type = c("genes", "introns", "exons", "write_all"),
                              ...) {

  # Check input
  assertthat::assert_that(assertthat::is.readable(gff3_file))
  assertthat::assert_that(assertthat::is.string(gene_label))
  assertthat::assert_that(assertthat::is.string(exon_label))
  assertthat::assert_that(assertthat::is.string(out_type))
  assertthat::assert_that(is.logical(verbose))
  assertthat::assert_that(out_type %in% c("genes", "introns", "exons", "write_all"),
                          msg = "'out_type' must be one of 'genes', 'introns', 'exons', or 'write_all'")


  # Read in gff3 file as dataframe
  gff3 <- ape::read.gff(gff3_file) %>%
    dplyr::mutate(chr = as.character(seqid))

  # Keep only annotations from selected source
  if (!is.null(source_select)) {
    assertthat::assert_that(is.character(source_select))
    assertthat::assert_that(all(source_select %in% gff3$source))
    gff3 <- dplyr::filter(gff3, source %in% source_select)
  }

  # Extract and clean up genes
  genes <- gff3 %>% dplyr::filter(type == gene_label) %>%
    dplyr::select(chr, start, end) %>%
    clean_gff(verbose = verbose)

  # Extract and clean up exons
  exons <- gff3 %>% dplyr::filter(type == exon_label) %>%
    dplyr::select(chr, start, end) %>%
    clean_gff(verbose = verbose)

  # Introns are genes - exons
  introns <- bedr::bedr.subtract.region(
    genes,
    exons,
    remove.whole.feature = FALSE,
    check.chr = FALSE,
    verbose = verbose)

  # Write out all regions and return hash of genes + exons + introns
  if (out_type == "write_all") {
    out_dir <- fs::path_abs(out_dir)
    assertthat::assert_that(assertthat::is.writeable(out_dir))
    assertthat::assert_that(assertthat::is.string(prefix))
    all_regions <- list(genes = genes,
                        exons = exons,
                        introns = introns)
    all_regions %>%
      purrr::set_names(fs::path(out_dir, paste0(prefix, "_", names(.)))) %>%
      purrr::iwalk(readr::write_tsv, col_names = FALSE)
    return(digest::digest(all_regions))
  }

  # Or, return a particular result type.
  return(switch(
    out_type,
    genes = genes,
    exons = exons,
    introns = introns
  ))

}

# Check that genome fasta headers and
# "chromosome" names in bed file match.
# (This must be true for
# bedtools maskfasta to work).
check_bed_genome_names <- function (fasta_file, bed) {
  # find all sequence headers in fasta file
  seq_names <-
    readr::read_lines(fasta_file) %>%
    purrr::keep(~ grepl(">", .x)) %>%
    purrr::map(~ gsub(">", "", .x))
  # make sure "chr" (chromosome) names of bed file are all in
  # fasta sequence headers
  chr <- unique(bed$chr)
  all(chr %in% seq_names)
}
