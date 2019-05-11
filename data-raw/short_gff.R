# Make small version of Arabidopsis GFF file for testing bed functions.
system("cat data-raw/Arabidopsis_thaliana.TAIR10.40.gff3 | head -n 1000 > inst/extdata/Arabidopsis_thaliana_TAIR10_40_small.gff3")
