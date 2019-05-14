#' Set the default path to Y&S scripts.
#'
#' This will overwrite the default value of the \code{path_to_ys} argument
#' for all \code{baitfindR} functions.
#'
#' @param path Character vector of length one; the complete path to the
#' folder containing Y&S python scripts, e.g.,
#' \code{"/Users/me/apps/phylogenomic_dataset_construction/"}
#' @examples
#' set_ys_path("/Users/me/apps/phylogenomic_dataset_construction/")
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @export
set_ys_path <- function (path) {
  pkgconfig::set_config("baitfindR::path_to_ys" = path)
}

#' Make a directory.
#'
#' Creates a directory in the working directory, and adds a hidden \code{.keep}
#' file. The purpose of the \code{.keep} file is to allow for tracking by
#'  \code{\link[drake]{drake_plan}} during workflows, because
#'  \code{\link[drake]{drake_plan}} can only track files, not folders.
#'
#' @param dir_name Name of the directory to be created.
#' @param ... Other arguments. Not used by this function, but meant to be used
#' by \code{\link[drake]{drake_plan}} for tracking during workflows.
#'
#' @return \code{NULL} in the R environment; externally, creates a
#' directory \code{dir_name}.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @examples
#' \dontrun{make_dir("new_dir")}
#' @export
make_dir <- function (dir_name, ...) {
  dir_name <- fs::path_abs(dir_name)
  if(!dir.exists(dir_name)) {
    dir.create(dir_name)
    sink(file = fs::path(dir_name, ".keep"), type = "output")
    cat(dir_name)
    sink()
  }
}


#' Concatenate files.
#'
#' Concatenate a list and write out the result as single file. Equivalent
#' of \code{cat file1 file2}.
#'
#' @param input_file_list A list. Could be any list, but it's meant to used
#' for lists of files that have been read into R as character vectors, e.g.,
#' using \code{\link[readr]{read_file}}.
#' @param output_file Path to write output file.
#'
#' @return A character vector of length one in the R environment; externally,
#' the concatenated output file.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @examples
#' \dontrun{cat_files(list("a", "b", "c"), "my_list.txt")}
#' @export
cat_files <- function (input_file_list, output_file) {
  unlisted_files <- unlist(input_file_list)
  collapsed_files <- paste(unlisted_files, collapse="")
  readr::write_file(collapsed_files, output_file)
}

#' Delete old output in a folder
#'
#' All file names matching `terms` will be deleted
#'
#' @param folder Path to folder
#' @param terms A regular expression: search terms to find files to delete.
#'
#' @return The deleted paths (invisibly).
#'
#' @examples
#' \dontrun{
#' make_dir("test")
#' fs::file_create("test/a")
#' fs::file_create("test/b")
#' fs::file_create("test/1")
#' list.files("test")
#' delete_old_output("test", "c")
#' list.files("test")
#' delete_old_output("test", "a|b")
#' list.files("test")
#' fs::file_delete("test")
#' }
delete_old_output <- function (folder, terms) {
  files_to_delete <- list.files(folder, pattern = terms, full.names = TRUE)
  fs::file_delete(files_to_delete)
}

#' Get MD5 hash of all files in a folder
#'
#' All file names matching `terms` will included
#'
#' @param folder Path to folder
#' @param terms A regular expression: search terms to find files to digest.
#'
#' @return Character vector of length 1; the MD5 hash.
#'
#' @examples
#' \dontrun{
#' make_dir("test")
#' fs::file_create("test/a")
#' fs::file_create("test/b")
#' fs::file_create("test/1")
#' list.files("test")
#' get_out_hash("test", "c")
#' get_out_hash("test", "a|b")
#' get_out_hash("test", "d")
#' get_out_hash("test", "e")
#' fs::file_delete("test")
#' }
get_out_hash <- function(folder, terms) {
  output <- list.files(folder, pattern = terms, full.names = TRUE)
  output <- if (length(output) > 0) unlist(lapply(output, readr::read_file)) else output
  return(digest::digest(output))
}

#' Write out a list of fasta files to a directory
#'
#' Optionally assign names to the files if `fasta_list` isn't already named
#' or if you want to over-write the original names.
#'
#' @param fasta_list List of DNA sequences to write out. Each item in list
#' must of class \code{\link[ape]{DNAbin}}.
#' @param out_dir Path to directory to write out DNA sequences.
#' @param fasta_names Optional character vector of file names to use when
#' writing out DNA sequences.
#' @param ext Optional character vector of length one; file extension to
#' append to DNA sequences (e.g., "fasta").
#' @param get_hash Logical; should the MD5 hash of `fasta_list` be returned?
#' @param ... Additional other arguments. Not used by this function,
#' but meant to be used by \code{\link[drake]{drake_plan}} for tracking
#' during workflows.
#'
#' @return None (invisible ‘NULL’) or character vector if `hash` is `TRUE`.
#' Externally, fasta files will be written to `out_dir`.
#'
#' @examples
#' # Load some example DNA sequences.
#' library(ape)
#' data(woodmouse)
#'
#' # Make a temporary working directory to write out files.
#' temp_dir <- fs::dir_create(fs::path(tempdir(), "baitfindR_example"))
#'
#' # Make list of DNA samples.
#' dna_list <- list(a = woodmouse, b = woodmouse)
#'
#' # Write out list:
#' # names as-is
#' write_fasta_files(dna_list, temp_dir)
#' list.files(temp_dir)
#'
#' # add extension
#' write_fasta_files(dna_list, temp_dir, ext = "fasta")
#' list.files(temp_dir)
#'
#' # new names and extension
#' write_fasta_files(
#'   dna_list,
#'   temp_dir,
#'   fasta_names = c("ho", "ge"),
#'   ext = "fasta")
#' list.files(temp_dir)
#'
#' # Cleanup
#' fs::file_delete(temp_dir)
#' @export
write_fasta_files <- function (fasta_list, out_dir,
                               fasta_names = NULL, ext = NULL,
                               get_hash = TRUE,
                               ...) {

  # Check input
  assertthat::assert_that(is.logical(get_hash))

  assertthat::assert_that(is.list(fasta_list))
  assertthat::assert_that(
    unique(purrr::map_chr(fasta_list, class)) == "DNAbin",
    msg = "All items in fasta_list must of class 'DNAbin'"
  )

  assertthat::assert_that(assertthat::is.string(ext) | is.null(ext))
  if(is.null(ext)) ext <- ""

  assertthat::assert_that(assertthat::is.string(out_dir))
  out_dir <- fs::path_abs(out_dir)
  assertthat::assert_that(assertthat::is.dir(out_dir))

  assertthat::assert_that(!is.null(fasta_names) | !is.null(names(fasta_list)),
                          msg = "fasta list must have names or fasta_names must
                          be provided")

  # Set new names
  if (!is.null(fasta_names)) {
    assertthat::assert_that(
      length(fasta_names) == length(fasta_list),
      msg = "fasta_names and fasta_list must be of same length"
    )
    names(fasta_list) <- fasta_names
  }

  # Write out fasta seqs
  fasta_list %>%
    purrr::set_names(fs::path(out_dir, names(.), ext = ext)) %>%
    purrr::iwalk(ape::write.FASTA)

  # optional: get MD5 hash of output
  if (isTRUE(get_hash)) unlist(fasta_list) %>% digest::digest()

}
