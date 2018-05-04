#' Fern transcriptomes
#'
#' Selected, downsized fern transcriptomes from the 1KP dataset to use for testing code.
#'
#' @format A list of eight fern transcriptomes stored in binary DNA format (class ape::DNAbin).
#' All transcriptomes have been randomly subsampled to 2\% of their original size.
#' @source \url{http://www.onekp.com/public_data.html}
"onekp_ferns"

#' Taxonomic data for fern transcriptomes
#'
#' Taxonomic data associated with the \code{\link{onekp_ferns}} dataset. Family-level taxonomy follows PPG I (2016).
#'
#' @format A data frame with 7 rows and 5 variables:
#' \describe{
#'   \item{code}{4 letter code used to identify the sample}
#'   \item{species}{Binomial species name}
#'   \item{genus}{Genus}
#'   \item{specific_epithet}{Specific epithet}
#'   \item{family}{Family}
#' }
#' @source \url{http://www.onekp.com/public_data.html}, Pteridophyte Phylogeny Group I. 2016. A community-derived classification for extant lycophytes and ferns. Journal of Systematics and Evolution 54:563–603.
"onekp_data"

#' Lygodium japonicum proteome
#'
#' Complete proteome (i.e., all coding genes) of Lygodium japonicum. Duplicate sequence IDs were identified and removed.
#'
#' @format A single text file in fasta format.
#' @source \url{http://bioinf.mind.meiji.ac.jp/kanikusa/}
"lygodium_proteome"

#' Arabidopsis thaliana proteome
#'
#' Complete proteome (i.e., all coding genes) of Arabidopsis thaliana.
#'
#' @format A single text file in fasta format.
#' @source \url{https://www.arabidopsis.org/}
"arabidopsis_proteome"
