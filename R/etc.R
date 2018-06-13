#' set_ys_path
#'
#' Set the default path to the folder containing Yang and Smith (2014) python scripts.
#'
#' This will overwrite the default value of the \code{path_to_ys} argument for all \code{baitfindR} functions.
#'
#' @param path Character vector of length one; the complete path to the folder containing Y&S python scripts, e.g., \code{"/Users/me/apps/phylogenomic_dataset_construction/"}
#' @examples
#' set_ys_path("/Users/me/apps/phylogenomic_dataset_construction/")
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @export
set_ys_path <- function (path) {
  pkgconfig::set_config("baitfindR::path_to_ys" = path)
}

#' set_transdecoder_path
#'
#' Set the default path to the folder containing \code{transdecoder} scripts.
#'
#' This will overwrite the default value of the \code{path_to_transdecoder} argument for all \code{baitfindR} functions.
#'
#' @param path Character vector of length one; the complete path to the folder containing \code{transdecoder}, e.g., \code{"/Users/me/apps/transdecoder/"}
#' @examples
#' set_ys_path("/Users/me/apps/phylogenomic_dataset_construction/")
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references Yang, Y. and S.A. Smith. 2014. Orthology inference in non-model organisms using transcriptomes and low-coverage genomes: improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31:3081-3092. \url{https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview}
#' @export
set_transdecoder_path <- function (path) {
  pkgconfig::set_config("baitfindR::path_to_transdecoder" = path)
}

#' make_dir
#'
#' Make a directory
#'
#' Creates a directory in the working directory, and adds a hidden \code{.name} file, which is a plain text file containing the name of the directory. The purpose of the \code{.name} file is to allow for tracking by \code{\link{drake}} during workflows, because \code{\link{drake}} can only track files, not folders.
#'
#' @param dir_name Name of the directory to be created.
#' @param ... Other arguments. Not used by this function, but meant to be used by \code{\link{drake}} for tracking during workflows.
#'
#' @return \code{NULL} in the R environment; externally, creates a directory \code{dir_name}.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @examples
#' \dontrun{make_dir("new_dir")}
#' @export
make_dir <- function (dir_name, ...) {
  if(!dir.exists(dir_name)) {
    dir.create(dir_name)
    sink(file = paste0(here::here(dir_name), "/.name"), type = "output")
    cat(dir_name)
    sink()
  }
}
