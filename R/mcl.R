#' mcl
#'
#' Wrapper to call mcl
#'
#' Calls the Markov Cluster Algorithm (mcl), a clustering algorithm for graphs. Here, it is meant to be used on genetic distances from BLAST.
#'
#' @param mcl_input Character vector of length one; the path to the input file for mcl clustering.
#' @param mcl_output Character vector of length one; the path to the output file produced by the mcl algorithm.
#' @param i_value Numeric or character vector of length one; the inflation value.
#' @param e_value Numeric or character vector of length one; the minimal -log transformed evalue to be considered by the algorithm.
#' @param other_args Character vector; other arguments to pass to mcl. Each should be an element of the vector. For example, to pass "-abc" to specify the input file format and "--te" to specify number of threads, use \code{c("--abc", "-te", "2")}.
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link{drake}} for tracking during workflows.
#'
#' @return A plain text file of tab-separated values, where each value on a line belongs to the same cluster.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Stijn van Dongen, A cluster algorithm for graphs. Technical Report INS-R0010, National Research Institute for Mathematics and Computer Science in the Netherlands, Amsterdam, May 2000. \url{https://micans.org/mcl/}
#' @examples
#' \dontrun{mcl(mcl_input = "some/folder/distance.file", mcl_output = "some/folder/mcl_output.txt", i_value = 1.4, evalue = 5)}
#' @export
mcl <- function (mcl_input, mcl_output, i_value, e_value, other_args = NULL, ...) {

  # modify arguments
  arguments <- c(mcl_input, "-I", i_value, "-tf", paste0("'gq(", e_value, ")'"), "-o", mcl_output, other_args )

  # run command
  system2("mcl", arguments)
}
