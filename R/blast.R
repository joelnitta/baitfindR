# blast -------------------------------------------------------------------

#' Build a BLAST database.
#'
#' This is a wrapper for makeblastdb.
#'
#' @param in_seqs Character vector of length one; the path to the fasta
#' file containing the sequences to be included in the database.
#' @param db_type Character vector of length one; "nucl" for DNA or "prot"
#' for amino acids (proteins).
#' @param out_name Character vector of length one; name of BLAST database
#' to be created (optional). This will be used to name all database files;
#' if omitted, the name of the `in_seqs` file will be used instead.
#' @param title Character vector of length one; title of BLAST database
#' to be created (optional).
#' @param parse_seqids Logical; should the makeblastdb flag
#' "parse_seqids" be used?
#' @param wd Character vector of length one; working directory. The blast
#' database will be made here.
#' @param echo Logical; should standard error and output be printed?
#' @param ... Additional other arguments. Not used by this function, but
#' meant to be used by \code{\link[drake]{drake_plan}} for tracking
#' during workflows.
#' @return A series of files starting with \code{out_name} and ending in
#' .phr, .pin, .pog, .psd, .psi, psq (for proteins) or .nhr, .nin, .nog,
#' .nsd, .nsi, and .nsq (for DNA) that constitute the BLAST database in
#' the working directory.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references \url{https://www.ncbi.nlm.nih.gov/books/NBK279690/}
#' @examples
#' library(ape)
#' data(woodmouse)
#' temp_dir <- fs::dir_create(fs::path(tempdir(), "baitfindR_example"))
#' ape::write.FASTA(woodmouse, fs::path(temp_dir, "woodmouse.fasta"))
#' list.files(temp_dir)
#' build_blast_db(
#'   fs::path(temp_dir, "woodmouse.fasta"),
#'   title = "test db",
#'   out_name = "wood",
#'   parse_seqids = TRUE,
#'   wd = temp_dir)
#' list.files(temp_dir)
#' fs::file_delete(temp_dir)
#' @export
build_blast_db <- function (in_seqs,
                            db_type = "nucl",
                            out_name = NULL,
                            title = NULL,
                            parse_seqids = FALSE,
                            echo = TRUE,
                            wd, ...) {

  # Check input
  assertthat::assert_that(assertthat::is.string(in_seqs))
  assertthat::assert_that(assertthat::is.string(db_type))
  assertthat::assert_that(assertthat::is.string(wd))
  assertthat::assert_that(is.logical(parse_seqids))
  assertthat::assert_that(assertthat::is.string(title) | is.null(title))
  assertthat::assert_that(assertthat::is.string(out_name) | is.null(out_name))

  assertthat::assert_that(db_type %in% c("nucl", "prot"))

  wd <- fs::path_abs(wd)
  assertthat::assert_that(assertthat::is.dir(wd))

  in_seqs <- fs::path_abs(in_seqs)
  assertthat::assert_that(assertthat::is.readable(in_seqs))

  # Prepare arguments
  parse_seqids <- if(isTRUE(parse_seqids)) "-parse_seqids" else NULL
  title <- if(!is.null(title)) c("-title", title) else NULL
  out_name <- if(!is.null(out_name)) c("-out", out_name) else NULL

  arguments <- c("-in", in_seqs,
                 "-dbtype", db_type,
                 parse_seqids,
                 out_name,
                 title)

  # run command
  processx::run("makeblastdb", arguments, wd = wd, echo = echo)

}

#' Run a blastp query.
#'
#' This is a wrapper for blastp.
#'
#' @param query Character vector of length one; the path to the fasta file
#' to use as the query sequence(s).
#' @param database Character vector of length one; the name of the blast
#' database.
#' @param out_file Character vector of length one; the name to use for the
#' results file.
#' @param outfmt Character vector of length one; value to pass to
#' \code{blastp} \code{outfmt} argument. Default = "6".
#' @param other_args Character vector; other arguments to pass on to
#' \code{blastp}. For a list of options, run \code{blastp -help}.
#' @param wd Character vector of length one; working directory. The blast
#' search will be conducted here.
#' @param echo Logical; should standard error and output be printed?
#' @param ... Additional other arguments. Not used by this function,
#' but meant to be used by \code{\link[drake]{drake_plan}} for tracking
#' during workflows.
#' @return A tab-separated text file with the results of the blastp
#' search, named with the value of \code{out_file}.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references \url{https://www.ncbi.nlm.nih.gov/books/NBK279690/}
#' @examples
#' library(ape)
#'
#' # Make temp dir for storing files
#' temp_dir <- fs::dir_create(fs::path(tempdir(), "baitfindR_example"))
#'
#' # Write out ape::woodmouse dataset as amino acids
#' data(woodmouse)
#' woodmouse_aa <- trans(woodmouse, 2)
#' ape::write.FASTA(woodmouse_aa, fs::path(temp_dir, "woodmouse.fasta"))
#'
#' # Make protein blast database
#' build_blast_db(
#'   fs::path(temp_dir, "woodmouse.fasta"),
#'   db_type = "prot",
#'   out_name = "wood",
#'   parse_seqids = TRUE,
#'   wd = temp_dir)
#'
#' # Blast the original sequences against the database
#' blast_p(
#'   fs::path(temp_dir, "woodmouse.fasta"),
#'   database = "wood",
#'   out_file = "blastp_results",
#'   wd = temp_dir,
#'   echo = TRUE
#' )
#'
#' # Take a look at the results.
#' readr::read_tsv(
#'   fs::path(temp_dir, "blastp_results"),
#'   col_names = FALSE
#'   )
#'
#' # Cleanup.
#' fs::file_delete(temp_dir)
#' @export
blast_p <- function (query,
                     database,
                     out_file = NULL,
                     outfmt = "6",
                     other_args = NULL,
                     echo = TRUE,
                     wd,
                     ...) {

  # Check input

  assertthat::assert_that(assertthat::is.string(query))
  assertthat::assert_that(assertthat::is.string(database))
  assertthat::assert_that(assertthat::is.string(out_file) | is.null(out_file))
  assertthat::assert_that(assertthat::is.string(outfmt))
  assertthat::assert_that(is.character(other_args) | is.null(other_args))
  assertthat::assert_that(is.logical(echo))
  assertthat::assert_that(assertthat::is.string(wd))

  wd <- fs::path_abs(wd)
  assertthat::assert_that(assertthat::is.dir(wd))

  query <- fs::path_abs(query)
  assertthat::assert_that(assertthat::is.readable(query))

  # modify arguments
  if(!is.null(out_file)) out_file <- c("-out", out_file)

  arguments <- c("-query", query,
                 "-db", database,
                 "-outfmt", outfmt,
                 out_file,
                 other_args)

  # run command
  processx::run("blastp", arguments, wd = wd, echo = echo)

}

#' Run a blastn query.
#'
#' This is a wrapper for blastn.
#'
#' @param query Character vector of length one; the path to the fasta
#' file to use as the query sequence(s).
#' @param database Character vector of length one; the name of the blast
#' database.
#' @param out_file Character vector of length one; the name to use for
#' the results file.
#' @param outfmt Character vector of length one; value to pass to
#' \code{blastn} \code{outfmt} argument. Default = "6".
#' @param other_args Character vector; other arguments to pass on to
#' \code{blastn}.
#' Must be formatted so that each argument name and its value are
#' separate, consecutive elements of the vector, e.g.,
#' \code{c("-evalue", 10, "-num_threads", 1)}.
#' The argument name must be preceded by a hyphen.
#' For a list of options, run \code{blastn -help}.
#' @param wd Character vector of length one; working directory. The blast
#' search will be conducted here.
#' @param echo Logical; should standard error and output be printed?
#' @param ... Additional other arguments. Not used by this function,
#' but meant to be used by \code{\link[drake]{drake_plan}} for tracking
#' during workflows.
#' @return A tab-separated text file with the results of the blastn
#' search, named with the value of \code{out_file}.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references \url{https://www.ncbi.nlm.nih.gov/books/NBK279690/}
#' @examples
#' library(ape)
#'
#' # Make temp dir for storing files
#' temp_dir <- fs::dir_create(fs::path(tempdir(), "baitfindR_example"))
#'
#' # Write out ape::woodmouse dataset as DNA
#' data(woodmouse)
#' ape::write.FASTA(woodmouse, fs::path(temp_dir, "woodmouse.fasta"))
#'
#' # Make blast database
#' build_blast_db(
#'   fs::path(temp_dir, "woodmouse.fasta"),
#'   db_type = "nucl",
#'   out_name = "wood",
#'   parse_seqids = TRUE,
#'   wd = temp_dir)
#'
#' # Blast the original sequences against the database
#' blast_n(
#'   fs::path(temp_dir, "woodmouse.fasta"),
#'   database = "wood",
#'   out_file = "blastn_results",
#'   wd = temp_dir,
#'   echo = TRUE
#' )
#'
#' # Take a look at the results.
#' readr::read_tsv(
#'   fs::path(temp_dir, "blastn_results"),
#'   col_names = FALSE
#'   )
#'
#' # Cleanup.
#' fs::file_delete(temp_dir)
#' @export
blast_n <- function (query,
                     database,
                     out_file = NULL,
                     outfmt = "6",
                     other_args = NULL,
                     echo = TRUE,
                     wd,
                     ...) {

  # Check input

  assertthat::assert_that(assertthat::is.string(query))
  assertthat::assert_that(assertthat::is.string(database))
  assertthat::assert_that(assertthat::is.string(out_file) | is.null(out_file))
  assertthat::assert_that(assertthat::is.string(outfmt))
  assertthat::assert_that(is.character(other_args) | is.null(other_args))
  assertthat::assert_that(is.logical(echo))
  assertthat::assert_that(assertthat::is.string(wd))
  assertthat::assert_that(
    length(other_args) > 1 | is.null(other_args),
    msg = "other_args not formatted correctly. Check help file by running ?baitfindr::blast_n")

  wd <- fs::path_abs(wd)
  assertthat::assert_that(assertthat::is.dir(wd))

  query <- fs::path_abs(query)
  assertthat::assert_that(assertthat::is.readable(query))

  # modify arguments
  if(!is.null(out_file)) out_file <- c("-out", out_file)

  arguments <- c("-query", query,
                 "-db", database,
                 "-outfmt", outfmt,
                 out_file,
                 other_args)

  # run command
  processx::run("blastn", arguments, wd = wd, echo = echo)

}

#' Run blastn on all fasta files in a folder.
#'
#' Output is written to the same folder containing the input files.
#'
#' @param fasta_folder Path to the folder containing fasta files to BLAST.
#' @param fasta_pattern Optional; pattern used for matching with grep. Only
#' files with names matching the pattern will be included in the
#' BLAST search.
#' @param database_path Path to the BLAST database, including the database
#' name.
#' @param out_ext File extension used for BLAST results files. The result of
#' each BLAST search will be a file with the same name as the input fasta files,
#' but with this extension appended.
#' @param outfmt String; format to use for BLAST output. See
#' https://www.ncbi.nlm.nih.gov/books/NBK279684/ (Table C1) for details.
#' @param other_args Character vector; other arguments to pass on to
#' \code{blastn}. For a list of options, run \code{blastn -help}.
#' @param overwrite Logical: should old output be erased before running this
#' function? "Old output" will be determined by matching any file names with
#' `out_ext`.
#' @param echo Logical; should standard error and output be printed?
#' @param get_hash Logical; if TRUE, the MD5 hash of the output will be
#' returned.
#' @param ... Additional other arguments. Not used by this function,
#' but meant to be used by \code{\link[drake]{drake_plan}} for tracking
#' during workflows.
#' @return NULL or character vector if `get_hash` is TRUE.
#' Externally, a text file file with the results of the blastn
#' search, named by adding `out_ext` to each input fasta file name.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references \url{https://www.ncbi.nlm.nih.gov/books/NBK279690/}
#'
#' @examples
#' library(ape)
#'
#' # Make temp dir for storing files
#' temp_dir <- fs::dir_create(fs::path(tempdir(), "baitfindR_example"))
#'
#' # Write out ape::woodmouse dataset as DNA
#' data(woodmouse)
#' ape::write.FASTA(woodmouse, fs::path(temp_dir, "woodmouse.fasta"))
#' ape::write.FASTA(woodmouse, fs::path(temp_dir, "woodmouse2.fasta"))
#'
#' # Make blast database
#' build_blast_db(
#'   fs::path(temp_dir, "woodmouse.fasta"),
#'   db_type = "nucl",
#'   out_name = "wood",
#'   parse_seqids = TRUE,
#'   wd = temp_dir)
#'
#' # Blast the original sequences against the database
#' blast_n_list(
#'   fasta_folder = temp_dir,
#'   fasta_pattern = "fasta",
#'   database_path = fs::path(temp_dir, "wood")
#' )
#'
#' # Take a look at the results.
#' readr::read_tsv(
#'   fs::path(temp_dir, "woodmouse.tsv"),
#'   col_names = FALSE
#'   )
#'
#' readr::read_tsv(
#'   fs::path(temp_dir, "woodmouse2.tsv"),
#'   col_names = FALSE
#'   )
#'
#' # Cleanup.
#' fs::file_delete(temp_dir)
#' @export
blast_n_list <- function (fasta_folder,
                          fasta_pattern,
                          database_path,
                          out_ext = "tsv",
                          outfmt = "6",
                          other_args = NULL,
                          overwrite = FALSE,
                          echo = FALSE,
                          get_hash = TRUE, ...) {

  # Check input
  assertthat::assert_that(assertthat::is.dir(fasta_folder))
  assertthat::assert_that(assertthat::is.string(fasta_pattern))
  assertthat::assert_that(assertthat::is.string(database_path))
  assertthat::assert_that(assertthat::is.string(out_ext))
  assertthat::assert_that(assertthat::is.string(outfmt))
  assertthat::assert_that(is.character(other_args) | is.null(other_args))
  assertthat::assert_that(is.logical(overwrite))
  assertthat::assert_that(is.logical(echo))
  assertthat::assert_that(is.logical(get_hash))

  # Setup arguments
  fasta_folder <- fs::path_abs(fasta_folder)
  database_path <- fs::path_abs(database_path)
  search_terms <- glue::glue("\\.{out_ext}$")

  # optional: delete all previous output written in this folder based on output extension
  if (isTRUE(overwrite)) delete_old_output(fasta_folder, search_terms)

  # List all fasta files in folder
  fasta_files <- list.files(fasta_folder, pattern = fasta_pattern, full.names = TRUE)

  # Make list of result files to write
  # (Just append out_ext to each input fasta file)
  results_list <- purrr::map_chr(fasta_files, ~fs::path_ext_set(., ext = out_ext))

  # Loop blast search over lists
  purrr::walk2(
    fasta_files,
    results_list,
    ~ blast_n(
      query = .x,
      out_file = .y,
      wd = fasta_folder,
      database = database_path,
      outfmt = outfmt,
      other_args = other_args,
      echo = echo)
  )

  # Optional: get MD5 hash of output
  if (isTRUE(get_hash)) get_out_hash(fasta_folder, search_terms)

}

#' Extract top blast hits from multiple blast output files
#'
#' BLAST results files must be in tabular format (e.g., outfmt 6). For each
#' BLAST results file, the single top hit (sorted by ascending evalue, then
#' descending bitscore) will be extracted from the BLAST database and written
#' to `out_dir`.
#'
#' @param blast_results_dir Path to folder containing BLAST results files.
#' @param blast_results_pattern Optional; pattern used for matching with grep.
#' Only files with names matching the pattern will be used to extract top blast
#' hits.
#' @param blast_cols Character vector; column names of BLAST results. See
#' https://www.ncbi.nlm.nih.gov/books/NBK279684/ (Table C1) for details. Must
#' include at least 'sseqid', 'evalue', and 'bitscore'. Defaults to standard
#' columns for output format 6.
#' @param database_path Path to the BLAST database, including the database
#' name.
#' @param out_dir Path to folder to write results.
#' @param out_ext File extension used when writing out top BLAST hit.
#' @param ... Additional other arguments. Not used by this function,
#' but meant to be used by \code{\link[drake]{drake_plan}} for tracking
#' during workflows.
#'
#' @return TRUE when runs successfully; externally, the top blast hit for each
#' fasta result will be written to `out_dir`.
#'
#' @examples
#' library(ape)
#'
#' # Make temp dir for storing files
#' temp_dir <- fs::dir_create(fs::path(tempdir(), "baitfindR_example"))
#'
#' # Write out ape::woodmouse dataset as DNA
#' data(woodmouse)
#' ape::write.FASTA(woodmouse, fs::path(temp_dir, "woodmouse.fasta"))
#' ape::write.FASTA(woodmouse, fs::path(temp_dir, "woodmouse2.fasta"))
#'
#' # Make blast database
#' build_blast_db(
#'   fs::path(temp_dir, "woodmouse.fasta"),
#'   db_type = "nucl",
#'   out_name = "wood",
#'   parse_seqids = TRUE,
#'   wd = temp_dir)
#'
#' # Blast the original sequences against the database
#' blast_n_list(
#'   fasta_folder = temp_dir,
#'   fasta_pattern = "fasta",
#'   database_path = fs::path(temp_dir, "wood")
#' )
#'
#' # Extract the top BLAST hit for each fasta file.
#' extract_blast_hits(
#'   blast_results_dir = temp_dir,
#'   blast_results_pattern = "\\.tsv$",
#'   database_path = fs::path(temp_dir, "wood"),
#'   out_dir = temp_dir
#' )
#'
#' list.files(temp_dir)
#' ape::read.FASTA(fs::path(temp_dir, "woodmouse.tsv.bestmatch.fasta"))
#'
#' # Cleanup.
#' fs::file_delete(temp_dir)
#' @export
extract_blast_hits <- function (
  blast_results_dir,
  blast_results_pattern,
  blast_cols = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                 "qstart", "qend", "sstart", "send", "evalue", "bitscore"),
  database_path,
  out_dir,
  out_ext = "bestmatch.fasta", ...) {

  # Check input
  assertthat::assert_that(assertthat::is.string(blast_results_pattern))
  assertthat::assert_that(assertthat::is.string(database_path))
  assertthat::assert_that(assertthat::is.string(out_ext))
  assertthat::assert_that(assertthat::is.dir(out_dir))
  assertthat::assert_that(assertthat::is.dir(blast_results_dir))
  assertthat::assert_that(is.character(blast_cols))


  out_dir <- fs::path_abs(out_dir)
  blast_results_dir <- fs::path_abs(blast_results_dir)

  # Make list of blastdbcmd calls to extract top blast hits
  command <-
    list.files(blast_results_dir, pattern = blast_results_pattern, full.names = TRUE) %>%
    purrr::set_names(.) %>%
    purrr::map_df(readr::read_tsv, col_names = blast_cols, .id = "file_name") %>%
    dplyr::mutate(
      file_name = fs::path_file(file_name),
      output = fs::path(out_dir, file_name, ext = out_ext)
    ) %>%
    dplyr::group_by(file_name) %>%
    dplyr::arrange(evalue, desc(bitscore)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    # format the command to extract the blast hit as a column
    dplyr::mutate(command = glue::glue('blastdbcmd -db {fs::path_file(database_path)} -entry "{sseqid}" > {output}')) %>%
    dplyr::pull(command)

  # Write out bash script
  # Name file after hash of commands vector, should be unique.
  temp_bash_file <- digest::digest(command)
  readr::write_lines(
    path = fs::path(fs::path_dir(database_path), temp_bash_file),
    x = c("#!/bin/bash", command)
  )

  # Run bash script
  processx::run("bash", temp_bash_file, echo = TRUE, wd = fs::path_dir(database_path))

  # Cleanup
  file.remove(fs::path(fs::path_dir(database_path), temp_bash_file))
}

#' Realign top blast hit of multi-fasta file with that fasta file.
#'
#' Top hits and original fasta files are matched based on the first part
#' of the filename separated by periods (i.e., the filename without any
#' extension).
#'
#' @param best_hits_dir Path to directory containing top blast hits.
#' @param best_hits_pattern Pattern used for matching with grep. Only
#' files with names matching the pattern will be included as the top
#' blast hit.
#' @param fasta_dir Path to directory containing fasta files for realignment.
#' @param fasta_pattern Pattern used for matching with grep. Only
#' files with names matching the pattern will be included for realignment.
#' @param ... Additional other arguments. Not used by this function,
#' but meant to be used by \code{\link[drake]{drake_plan}} for tracking
#' during workflows.
#'
#' @return List of lists, each of which is of class `DNAbin`.
#'
#' @examples
#'
#' library(ape)
#'
#' # Make temp dir for storing files
#' temp_dir <- fs::dir_create(fs::path(tempdir(), "baitfindR_example"))
#'
#' # Write out ape::woodmouse dataset as DNA
#' data(woodmouse)
#' ape::write.FASTA(woodmouse, fs::path(temp_dir, "woodmouse.fasta"))
#' ape::write.FASTA(woodmouse, fs::path(temp_dir, "woodmouse2.fasta"))
#'
#' # Make blast database
#' build_blast_db(
#'   fs::path(temp_dir, "woodmouse.fasta"),
#'   db_type = "nucl",
#'   out_name = "wood",
#'   parse_seqids = TRUE,
#'   wd = temp_dir)
#'
#' # Blast the original sequences against the database
#' blast_n_list(
#'   fasta_folder = temp_dir,
#'   fasta_pattern = "fasta",
#'   database_path = fs::path(temp_dir, "wood")
#' )
#'
#' # Extract the top BLAST hit for each fasta file.
#' extract_blast_hits(
#'   blast_results_dir = temp_dir,
#'   blast_results_pattern = "\\.tsv$",
#'   database_path = fs::path(temp_dir, "wood"),
#'   out_dir = temp_dir,
#'   out_ext = "bestmatch"
#' )
#'
#' realign_with_best_hits(
#'   best_hits_dir = temp_dir,
#'   best_hits_pattern = "bestmatch",
#'   fasta_dir = temp_dir,
#'   fasta_pattern = "fasta"
#' )
#'
#' # Cleanup.
#' fs::file_delete(temp_dir)
#' @export
realign_with_best_hits <- function (best_hits_dir,
                                    best_hits_pattern = "bestmatch",
                                    fasta_dir,
                                    fasta_pattern = "\\.fa$", ...) {

  # Check input
  assertthat::assert_that(assertthat::is.string(best_hits_dir))
  assertthat::assert_that(assertthat::is.string(best_hits_pattern))
  assertthat::assert_that(assertthat::is.string(fasta_dir))
  assertthat::assert_that(assertthat::is.string(fasta_pattern))

  best_hits_dir <- fs::path_abs(best_hits_dir)
  assertthat::assert_that(assertthat::is.dir(best_hits_dir))
  best_hits_dir <- fs::path_abs(best_hits_dir)

  fasta_dir <- fs::path_abs(fasta_dir)
  assertthat::assert_that(assertthat::is.dir(fasta_dir))
  fasta_dir <- fs::path_abs(fasta_dir)

  # Get file names of top blast hits
  best_hits_files <- list.files(best_hits_dir, best_hits_pattern, full.names = TRUE)

  # Read in seqs of top blast hits
  blast_top_matches <- purrr::map(best_hits_files, ape::read.FASTA)

  # Get file names of alignments to which to add best hits
  fasta_files <- list.files(fasta_dir, fasta_pattern, full.names = TRUE)

  # Read in alignments to which to add best hits
  fasta_to_add <- purrr::map(fasta_files, ape::read.FASTA)

  # Match sequence to add to alignment by first part of filename
  # (without extension).
  fasta_to_add_names <-
    fasta_files %>%
    fs::path_file() %>%
    stringr::str_split("\\.") %>%
    purrr::map_chr(1)

  best_hits_names <-
    best_hits_files %>%
    fs::path_file() %>%
    stringr::str_split("\\.") %>%
    purrr::map_chr(1)

  select <- fasta_to_add_names == best_hits_names

  assertthat::assert_that(
    any(select),
    msg = "No names match between top blast hits and fasta sequences to realign"
  )

  # combine and re-align blast-filtered alignments with their top matches
  purrr::map2(fasta_to_add[select], blast_top_matches[select], c) %>%
    purrr::map(ips::mafft, path = "/usr/bin/mafft", options = "--adjustdirection") %>%
    purrr::set_names(best_hits_names[select])
}
