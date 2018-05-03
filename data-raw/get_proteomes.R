# download arabidopsis and lygodium proteomes
# these will be used to make blastp database in Yang and Smith (2014) workflow

library(readr)

# download arabidopsis proteome
download.file(url = "https://www.arabidopsis.org/download_files/Sequences/TAIR10_blastsets/TAIR10_pep_20110103_representative_gene_model_updated", destfile = "data-raw/TAIR10_pep_20110103_representative_gene_model_updated.fasta")
# can verify this is (probably) the same as used in Y&S (2014) Atha.fa because
# grep ">" TAIR10_pep_20110103_representative_gene_model_updated.txt -c
# and
# grep ">" Atha.fa -c
# both yield 27416

# download lygodium proteome (note, but don't correct, 'potein' typo)
download.file(url = "http://bioinf.mind.meiji.ac.jp/kanikusa/data/download/lygodium_predicted_potein_ver1.0RC.fasta.tar.gz", destfile = "data-raw/lygodium_predicted_potein_ver1.0RC.fasta.tar.gz")

# lygodium proteome has some issues: repeated sequence names. unzip and remove these before importing.
# unzip
system("tar -zxf data-raw/lygodium_predicted_potein_ver1.0RC.fasta.tar.gz -C data-raw/")

# remove repeat names with awk
system("awk '/^>/{f=!d[$1];d[$1]=1}f' data-raw/lygodium_predicted_potein_ver1.0RC.fasta > data-raw/lygodium_predicted_potein_ver1.0RC.clean.fasta")

# read in cleaned lygodium proteome as character vector
# (use seqinr::read.fasta instead to store within in R as a peptide sequence)
lygodium_proteome <- read_file("data-raw/lygodium_predicted_potein_ver1.0RC.clean.fasta")

# note that when we want to use this, it is a plain text, so use write_file
# write_file(lygodium_proteome_txt, "lygodium_proteome.fasta")

# read in arabidopsis proteome as character vector like we did for lygodium
arabidopsis_proteome <- read_file("data-raw/TAIR10_pep_20110103_representative_gene_model_updated.fasta")

# export data
usethis::use_data(lygodium_proteome, overwrite = TRUE, compress="xz")
usethis::use_data(arabidopsis_proteome, overwrite = TRUE, compress="xz")
