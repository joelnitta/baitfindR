#' tcs
#'
#' Wrapper to call t_coffee to calculate Transitive Consistency Score (TCS)
#'
#' @param alignment Character vector of length one; the path to the input alignment
#' @param number_cores Number of cores to use by t_coffee
#' @param other_args Character vector; other arguments to pass to t_coffee. Each should be an element of the vector.
#' @param wd Character vector of length one; the directory where the command will be run, and the external output written.
#' @param ... Additional other arguments. Not used by this function, but meant to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.
#'
#' @return
#' A numeric value; the overall TCS score of the input alignment.
#'
#' Externally, a plain text file with TCS scores for each sequence, position, and the whole alignment. The file will be named \code{<alignment>.score_ascii}, where \code{<alignment>} is the filename of the input alignment.
#'
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Chang JM, Di Tommaso P, Notredame C. TCS: a new multiple sequence alignment reliability measure to estimate alignment accuracy and improve phylogenetic tree reconstruction. Mol Biol Evol. 2014 Jun. 31(6):1625-37. \url{http://www.tcoffee.org/Projects/tcs/}
#' @examples
#' \dontrun{tcs("some/alignment.fasta")}
#'
#' @export
tcs <- function (alignment, number_cores = 1, other_args = NULL, wd = here::here(), ...) {

  # assemble arguments
  arguments <- c("-infile", alignment, "-n_core", number_cores, "-output", "score_ascii", "-score", other_args)

  # execute command
  output <- processx::run("t_coffee", arguments, wd = wd)

  # isolate overall score from stdout and return
  stdout <- output$stdout
  stdout <- unlist(stringr::str_split(stdout, "\n"))
  score_line <- stdout[grep("SCORE", stdout)]
  score_line <- unlist(stringr::str_split(score_line, ","))
  score <- score_line[grep("SCORE", score_line)]
  score <- stringr::str_extract(score, "\\d*$")
  score <- as.numeric(score)

  return(score)
}

#' tcs_loop
#'
#' Call t_coffee to calculate Transitive Consistency Scores (TCS) for all alignments in a folder
#'
#' @param folder Character vector of length one; the path to the folder containing the alignments.
#' @param pattern An optional regular expression. Only alignment files with names that match the regular expression will be included.
#' @param number_cores Number of cores for t-coffee to use.
#' @param other_args Character vector; other arguments to pass to t_coffee. Each should be an element of the vector.
#' @param ... Additional other arguments. Not used by this function, but meant to be used by \code{\link[drake]{drake_plan}} for tracking during workflows.
#'
#' @return A numeric vector; the overall TCS scores for all alignments in the folder.
#'
#' Externally, plain text files that include TCS scores for each sequence, position, and the whole alignment overall, for each alignment. The files will be written to \code{folder}, and named \code{<alignment>.score_ascii}, where \code{<alignment>} is the filename of the input alignment.
#'
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Chang JM, Di Tommaso P, Notredame C. TCS: a new multiple sequence alignment reliability measure to estimate alignment accuracy and improve phylogenetic tree reconstruction. Mol Biol Evol. 2014 Jun. 31(6):1625-37. \url{http://www.tcoffee.org/Projects/tcs/}
#' @examples
#' \dontrun{tcs_loop(folder = "some/folder/with/alignments/", pattern = "fa.mafft.aln-cln$")}
#'
#' @export
tcs_loop <- function (folder, pattern = NULL, number_cores = 1, other_args = NULL, ...) {

  # add trailing slash
  folder <- jntools::add_slash(folder)

  # run loop
  alignments <- list.files(folder, pattern = pattern)
  scores <- purrr::map_dbl(alignments, baitfindR::tcs, wd = folder, number_cores = number_cores, other_args = other_args)

  return(scores)
}
