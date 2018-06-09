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
  pkgconfig::set_config("baitfindR::set_transdecoder_path" = path)
}
