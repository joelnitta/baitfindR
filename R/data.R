#' Taxonomic data for fern transcriptomes
#'
#' Taxonomic data for seven fern transcripomes downloaded from the onekp project. Family-level taxonomy follows PPG I (2016).
#'
#' @format A data frame with 8 rows and 6 variables:
#' \describe{
#'   \item{code}{4 letter code used to identify the sample}
#'   \item{species}{Binomial species name}
#'   \item{genus}{Genus}
#'   \item{specific_epithet}{Specific epithet}
#'   \item{family}{Family}
#'   \item{group}{Ingroup or outgroup status}
#' }
#' @source \url{http://www.onekp.com/public_data.html}
#'
#' Pteridophyte Phylogeny Group I. 2016. A community-derived classification for extant lycophytes and ferns. Journal of Systematics and Evolution 54:563–603.
"onekp_data"

#' Fern transcriptome
#'
#' List of fern transcriptomes from the oneKp project that have been randomly
#' downsized to 5 percent of original size. Named by 1kP code.
#'
#' @format List of lists, each of class 'DNAbin".
#' @source \url{http://www.onekp.com/public_data.html}
"example_transcriptomes"
