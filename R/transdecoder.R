# Wrapper to call transdecoder for a particular sequence
# This produces one .transdecoder_dir folder as output per transcriptome in the working directory
transdecoder_long_orfs <- function (transcriptome, transdecoder_result_folder, path_to_transdecoder = "~/apps/TransDecoder/", other_args = NULL, outfile) {

  arguments <- paste0("-t ", transcriptome, " ", other_args)

  # Normally I would use system2() but for some reason it doesn't like the path in the command.
  # system() seems to have no issues with this
  system(paste0(path_to_transdecoder, "TransDecoder.LongOrfs ", arguments))

  # Copy the output (transdecoder_result_folder) to translation subfolder,
  # using rsync to overwrite any existing files
  system2("rsync", c("-a", here::here(transdecoder_result_folder), here::here("translation")))

  # Delete the copy remaining in the working folder
  system2("rm", c("-r ", here::here(transdecoder_result_folder)))
}
