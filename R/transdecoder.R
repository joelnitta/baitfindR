#' Find long open reading frames.
#'
#' Extracts long open reading frames (ORFs) from a fasta file containing transcript
#' sequences (i.e., the transcriptome). This is a wrapper that calls
#' \code{TransDecoder.LongOrfs}.
#'
#' @param transcriptome_file Character vector of length one; the path to the
#' fasta file containing transcript sequences (i.e., the transcriptome).
#' @param wd Character vector of length one; the directory where the command
#' will be run, and the output folder created.
#' @param other_args Character vector; other arguments to pass to TransDecoder.
#' Each should be an element of the vector.
#' @param echo Logical; should the standard output and error be printed to the screen?
#' @param ... Additional other arguments. Not used by this function, but meant to
#' be used by \code{\link[drake]{drake_plan}} for tracking during workflows.
#'
#' @return
#' Within the R environment, a list with components specified in
#' \code{\link[processx]{run}}.
#'
#' Externally, an output folder in the working directory containing base frequencies
#' and .cds, .gff3, and .pep files for the recovered ORFs. The folder will be named
#' \code{<transcriptome>.transdecoder_dir}, where \code{<transcriptome>} is the value
#' of that argument.
#'
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references \url{http://transdecoder.github.io}
#' @examples
#' \dontrun{
#' library(ape)
#' temp_dir <- tempdir()
#' data("example_transcriptomes")
#' PSKY <- example_transcriptomes$PSKY
#' write.FASTA(PSKY, fs::path(temp_dir, "PSKY.fasta"))
#' transdecoder_long_orfs(
#'   transcriptome_file = fs::path(temp_dir, "PSKY.fasta"),
#'   wd = temp_dir
#'   )
#' list.files(temp_dir)
#' }
#' @export
transdecoder_long_orfs <- function (
  transcriptome_file,
  wd = here::here(),
  other_args = NULL,
  echo = pkgconfig::get_config("baitfindR::echo", fallback = FALSE),
  ...) {

  # check input
  wd <- fs::path_abs(wd)
  assertthat::assert_that(assertthat::is.dir(wd))

  transcriptome_file <- fs::path_abs(transcriptome_file)
  assertthat::assert_that(assertthat::is.readable(transcriptome_file))

  assertthat::assert_that(is.character(other_args) | is.null(other_args))

  # modify arguments
  arguments <- c("-t", transcriptome_file, other_args)

  # run command
  processx::run("TransDecoder.LongOrfs", arguments, wd = wd, echo = echo)
}

#' Predict coding regions.
#'
#' This is a wrapper for transdecoder.predict. It must be run in the same directory
#' containing the .transdecoder_dir output folders from
#' \code{\link{transdecoder_long_orfs}}. Optionally include the results of a
#' blastp search to make sure that peptides with a blastp hit against the
#' reference database are retained in the TransDecoder output.
#'
#' For a more detailed example, see vignettes.
#'
#' @param transcriptome_file Character vector of length one; the path to the fasta
#' file containing transcript sequences (i.e., the transcriptome).
#' @param blast_result Character vector of length one; the path to the tab-separated
#' text file containing the results from a blastp search of the transcriptome against
#' a reference blast protein database. For the blast search, the output format should
#' specified as: -outfmt 6.
#' @param wd Character vector of length one; the directory where the command will be
#' run. Must contain .transdecoder_dir folder with results from
#' \code{\link{transdecoder_long_orfs}}.
#' @param echo Logical; should the standard output and error be printed to the screen?
#' @param other_args Character vector; other arguments to pass to TransDecoder. Each
#' should be an element of the vector.
#' @param ... Additional other arguments. Not used by this function, but meant to be
#' used by \code{\link[drake]{drake_plan}} for tracking during workflows.
#'
#' @return
#' Within the R environment, a list with components specified in \code{\link[processx]{run}}.
#'
#' Externally, four output files (.cds, .pep, .bed, and .gff3) for each input
#' transcriptome file will be written to the working directory.
#'
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references \url{http://transdecoder.github.io}
#' @examples
#' \dontrun{
#' transdecoder_predict("some/transcriptome_file.fa", "some/blast_result.txt")
#' }
#'
#' @export
transdecoder_predict <- function (
  transcriptome_file,
  blast_result = NULL,
  wd = here::here(),
  other_args = NULL,
  echo = pkgconfig::get_config("baitfindR::echo", fallback = FALSE),
  ...) {

  # check input
  wd <- fs::path_abs(wd)
  assertthat::assert_that(assertthat::is.dir(wd))

  transcriptome_file <- fs::path_abs(transcriptome_file)
  assertthat::assert_that(assertthat::is.readable(transcriptome_file))

  assertthat::assert_that(is.character(other_args) | is.null(other_args))

  # modify arguments
  blast_argument <- if(is.null(blast_result)) { NULL } else {
    assertthat::assert_that(assertthat::is.string(blast_result))
    blast_result <- fs::path_abs(blast_result)
    assertthat::assert_that(assertthat::is.readable(blast_result))
    c("--retain_blastp_hits", blast_result)
  }

  arguments <- c("-t", transcriptome_file, blast_argument, other_args)

  # run command
  processx::run("TransDecoder.Predict", arguments, wd = wd, echo = echo)

}

#' Cluster DNA sequences.
#'
#' This is a wrapper for the CD-HIT-EST algorithm. According to the CD-HIT user's guide,
#' "CD-HIT-EST clusters a nucleotide dataset into clusters that meet a user-defined
#' similarity threshold, usually a sequence identity." cd-hit-est comes bundled with
#' transdecoder, so it is run from there.
#'
#' @param input Character vector of length one; the path to the input file for
#' cd-hit-est. Should be DNA or AA sequences in fasta format.
#' @param output Character vector of length one; the name to assign to the output.
#' Can include a path, in which case the output will be written there.
#' @param wd Character vector of length one; the directory where the command
#' will be run.
#' @param echo Logical; should the standard output and error be printed to the screen?
#' @param other_args Character vector; other arguments to pass to cd-hit-est.
#' Each should be an element of the vector.
#' @param ... Additional other arguments. Not used by this function, but meant
#' to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.
#'
#' @return
#' Within the R environment, a list with components specified in
#' \code{\link[processx]{run}}.
#'
#' Externally, two files will be written: according to the CD-HIT user's guide,
#' "The output are two files: a fasta file of representative sequences and a text
#' file of list of clusters."
#'
#' The fasta file will be named with the value of \code{output}; the list of clusters
#' will be the same, with \code{.clstr} appended.
#'
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references \url{http://www.bioinformatics.org/cd-hit/}, \url{http://transdecoder.github.io}
#' @examples
#' \dontrun{
#' library(ape)
#' library(baitfindR)
#'
#' # Make temp dir for storing output
#' temp_dir <- fs::dir_create(fs::path(tempdir(), "baitfindR_example"))
#' data("PSKY")
#'
#' # Write downsized transcriptome to temp dir
#' write.FASTA(PSKY, fs::path(temp_dir, "PSKY"))
#'
#' # Get CDS
#' transdecoder_long_orfs(
#'   transcriptome_file = fs::path(temp_dir, "PSKY"),
#'   wd = temp_dir
#'   )
#'
#' # Cluster similar genes in CDS
#' cd_hit_est(
#'   input = fs::path(temp_dir, "PSKY.transdecoder_dir", "longest_orfs.cds"),
#'   output = fs::path(temp_dir, "PSKY.cd-hit-est"),
#'   wd = temp_dir,
#'   echo = TRUE
#' )
#'
#' # Check output
#' list.files(temp_dir)
#' head(readr::read_lines(fs::path(temp_dir, "PSKY.cd-hit-est")))
#' head(readr::read_lines(fs::path(temp_dir, "PSKY.cd-hit-est.clstr")))
#'
#' # Cleanup
#' fs::file_delete(temp_dir)
#' }
#' @export
cd_hit_est <- function (
  input,
  output,
  wd = here::here(),
  other_args = NULL,
  echo = pkgconfig::get_config("baitfindR::echo", fallback = FALSE),
  ...) {

  # Check input
  assertthat::assert_that(assertthat::is.string(input))
  assertthat::assert_that(assertthat::is.string(output))
  assertthat::assert_that(is.character(other_args) | is.null(other_args))
  assertthat::assert_that(is.logical(echo))
  assertthat::assert_that(assertthat::is.string(wd))

  wd <- fs::path_abs(wd)
  assertthat::assert_that(assertthat::is.dir(wd))

  input <- fs::path_abs(input)
  assertthat::assert_that(assertthat::is.readable(input))

  # modify arguments
  arguments <- c("-i", input, "-o", output, other_args)

  # run command
  processx::run("cd-hit-est", arguments, wd = wd, echo = echo)

}
