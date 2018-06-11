#' build_blast_db
#'
#' Wrapper to call makeblastdb
#'
#' Takes a list of input sequences and makes them into a BLAST database. Assumes that all BLAST commands are in the user's \code{$PATH}.
#'
#' @param in_seqs Character vector of length one; the path to the fasta file containing the sequences to be included in the database.
#' @param out_name Character vector of length one; name of BLAST database to be created.
#' @param db_type Character vector of length one; "nucl" for DNA or "prot" for amino acids (proteins).
#' @param other_args Character vector; other arguments to pass on to \code{makeblastdb}. For a list of options, run \code{makeblastdb -help}.
#' @param ... Additional other arguments. Not used by this function, but meant to be used by \code{\link{drake}} for tracking during workflows.
#' @return A series of files starting with \code{out_name} and ending in .phr, .pin, .pog, .psd, .psi, psq (for proteins) or .nhr, .nin, .nog, .nsd, .nsi, and .nsq (for DNA) that constitute the BLAST database in the directory specified by \code{wd}.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references \url{https://www.ncbi.nlm.nih.gov/books/NBK279690/}
#' @export
build_blast_db <- function (in_seqs, out_name, db_type = "nucl", other_args = NULL, ...) {

  # modify arguments
  arguments <- c(paste("-in", in_seqs), paste("-out", out_name), paste("-dbtype", db_type), other_args)

  # run command
  system2("makeblastdb", arguments)

}

#' blastp
#'
#' Wrapper to call blastp
#'
#' @param query Character vector of length one; the path to the fasta file to use as the query sequence(s).
#' @param database Character vector of length one; the name of the blast database.
#' @param out_file Character vector of length one; the name to use for the results file.
#' @param other_args Character vector; other arguments to pass on to \code{makeblastdb}. For a list of options, run \code{blastp -help}.
#' @param ... Additional other arguments. Not used by this function, but meant to be used by \code{\link{drake}} for tracking during workflows.
#' @return A tab-separated text file with the results of the blastp search, named with the value of \code{out_file}.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references \url{https://www.ncbi.nlm.nih.gov/books/NBK279690/}
#' @export
blastp <- function (query, database, out_file, other_args = NULL, ...) {

  # modify arguments
  arguments <- c(paste("-query", query), paste("-db", database), paste("-out", out_file), other_args)

  # run command
  system2("makeblastdb", arguments)

}

# # blast n single file, save as outformat 6
# blast_cdhitest <- function (query, blast_result_file, database, other_args = "-evalue 10 -num_threads 1 -max_target_seqs 1000 -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore'", ...) {
#   arguments <- paste0("-query ", query, " -db ", database, " ", other_args, " > ", blast_result_file)
#   system2("blastn", arguments)
# }
