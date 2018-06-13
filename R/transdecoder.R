#' transdecoder_long_orfs
#'
#' Wrapper to call TransDecoder.LongOrfs
#'
#' Extracts long open reading frames (ORFs) from a fasta file containing transcript sequences (i.e., the transcriptome).
#'
#' @param path_to_transdecoder Character vector of length one; the path to the folder containing TransDecoder scripts, e.g., \code{"/Users/me/apps/TransDecoder/"}
#' @param transcriptome_file Character vector of length one; the path to the fasta file containing transcript sequences (i.e., the transcriptome).
#' @param wd Character vector of length one; the directory where the command will be run, and the output folder created.
#' @param other_args Character vector; other arguments to pass to TransDecoder. Each should be an element of the vector.
#' @param ... Additional other arguments. Not used by this function, but meant to be used by \code{\link{drake}} for tracking during workflows.
#'
#' @return An output folder in the working directory containing base frequencies and .cds, .gff3, and .pep files for the recovered ORFs. The folder will be named \code{<transcriptome>.transdecoder_dir}, where \code{<transcriptome>} is the value of that argument.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references \url{http://transdecoder.github.io}
#' @examples
#' \dontrun{transdecoder_long_orfs("some/transcriptome_file.fa")}
#'
#' @export
transdecoder_long_orfs <- function (path_to_transdecoder = pkgconfig::get_config("baitfindR::path_to_transdecoder"), transcriptome_file, wd = here::here(), other_args = NULL, ...) {

  # error checking
  if(is.null(path_to_transdecoder)) {
    stop("Must provide 'path_to_transdecoder' (path to TransDecoder folder)")
  }

  # modify arguments
  path_to_transdecoder <- jntools::add_slash(path_to_transdecoder)
  arguments <- c("-t", transcriptome_file, other_args)

  # modify command
  command <- paste0(path_to_transdecoder, "TransDecoder.LongOrfs")

  # run command
  processx::run(command, arguments, wd = wd)
}

#' transdecoder_predict
#'
#' Wrapper to call transdecoder.predict
#'
#' This has to be run in the same directory containing the .transdecoder_dir output folders from \code{\link{transdecoder_long_orfs}}. Optionally include the results of a blastp search to make sure that peptides with a blastp hit against the reference database are retained in the TransDecoder output.
#'
#' @param path_to_transdecoder Character vector of length one; the path to the folder containing TransDecoder scripts, e.g., \code{"/Users/me/apps/TransDecoder/"}
#' @param transcriptome_file Character vector of length one; the path to the fasta file containing transcript sequences (i.e., the transcriptome).
#' @param blast_result Character vector of length one; the path to the tab-separated text file containing the results from a blastp search of the transcriptome against a reference blast protein database. For the blast search, the output format should specified as: -outfmt 6.
#' @param wd Character vector of length one; the directory where the command will be run. Must contain .transdecoder_dir folder with results from \code{\link{transdecoder_long_orfs}}.
#' @param other_args Character vector; other arguments to pass to TransDecoder. Each should be an element of the vector.
#' @param ... Additional other arguments. Not used by this function, but meant to be used by \code{\link{drake}} for tracking during workflows.
#'
#' @return Four output files (.cds, .pep, .bed, and .gff3) for each input transcriptome file will be written to the working directory.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references \url{http://transdecoder.github.io}
#' @examples
#' \dontrun{transdecoder_predict_with_blast("some/transcriptome_file.fa", "some/blast_result.txt")}
#'
#' @export
transdecoder_predict <- function (path_to_transdecoder = pkgconfig::get_config("baitfindR::path_to_transdecoder"), transcriptome_file, blast_result = NULL, wd = here::here(), other_args = NULL, ...) {

  # error checking
  if(is.null(path_to_transdecoder)) {
    stop("Must provide 'path_to_transdecoder' (path to TransDecoder folder)")
  }

  # modify arguments
  path_to_transdecoder <- jntools::add_slash(path_to_transdecoder)
  blast_argument <- if(is.null(blast_result)) {NULL} else {c("--retain_blastp_hits", blast_result)}

  arguments <- c("-t", transcriptome_file, blast_argument, other_args)

  # modify command
  command <- paste0(path_to_transdecoder, "TransDecoder.Predict")

  # run command
  processx::run(command, arguments, wd = wd)

}
