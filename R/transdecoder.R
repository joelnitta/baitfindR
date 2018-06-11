#' transdecoder_long_orfs
#'
#' Wrapper to call TransDecoder.LongOrfs
#'
#' Extracts long open reading frames (ORFs) from a fasta file containing transcript sequences (i.e., the transcriptome).
#'
#' @param path_to_transdecoder Character vector of length one; the path to the folder containing TransDecoder scripts, e.g., \code{"/Users/me/apps/TransDecoder/"}
#' @param transcriptome Character vector of length one; the path to the fasta file containing transcript sequences (i.e., the transcriptome).
#' @param wd Character vector of length one; the directory where the command will be run, and the output folder created.
#' @param other_args Character vector; other arguments to pass to TransDecoder. Each should be an element of the vector.
#' @param ... Additional other arguments. Not used by this function, but meant to be used by \code{\link{drake}} for tracking during workflows.
#' @return An output folder in the working directory containing base frequencies and .cds, .gff3, and .pep files for the recovered ORFs. The folder will be named \code{<transcriptome>.transdecoder_dir}, where \code{<transcriptome>} is the value of that argument.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references \url{http://transdecoder.github.io}
#' @examples
#' \dontrun{transdecoder_long_orfs("some/transcriptome_file.fa")}
transdecoder_long_orfs <- function (path_to_transdecoder = pkgconfig::get_config("baitfindR::path_to_transdecoder"), transcriptome, wd, other_args = NULL, ...) {

  # error checking
  if(is.null(path_to_transdecoder)) {
    stop("Must provide 'path_to_transdecoder' (path to TransDecoder folder)")
  }

  # modify arguments
  path_to_transdecoder <- jntools::add_slash(path_to_transdecoder)
  arguments <- c("-t", transcriptome, other_args)

  # modify command
  command <- paste0(path_to_transdecoder, "TransDecoder.LongOrfs")

  # run command
  processx::run(command, arguments, wd = wd)
}
