make\_small\_transcriptome
================
Joel Nitta
5/12/2019

Make small version of a transcriptome file for testing functions.

The transriptome will be randomly downsized to the fraction set by `keep_frac`.

``` r
# Set fraction of transcriptome to keep
keep_frac <- 0.05

# Set seed
set.seed(5394)

# Clear out old files, if present
fs::file_delete(here::here("data-raw/PSKY-SOAPdenovo-Trans-assembly.fa.bz2"))

# Download transcriptome
download.file(
  url = "https://drive.google.com/uc?export=download&id=1Avi1oAeW9RPExzzreWZMht7EplvKEX0J",
  destfile = here::here("data-raw/PSKY-SOAPdenovo-Trans-assembly.fa.bz2")
)

# Read in file with ape. Use bzfile() because it's compressed
PSKY <- ape::read.FASTA(
  bzfile(here::here("data-raw/PSKY-SOAPdenovo-Trans-assembly.fa.bz2"))
    )

# Randomly downsize to specified fraction of transcripts
PSKY <- PSKY[sample(1:length(PSKY), keep_frac*length(PSKY))]

# Write out to data/
usethis::use_data(PSKY, overwrite = TRUE)
```

    ## ✔ Setting active project to '/home/rstudio/baitfindR'
    ## ✔ Saving 'PSKY' to 'data/PSKY.rda'

``` r
sessionInfo()
```

    ## R version 3.5.2 (2018-12-20)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Debian GNU/Linux 9 (stretch)
    ## 
    ## Matrix products: default
    ## BLAS: /usr/lib/openblas-base/libblas.so.3
    ## LAPACK: /usr/lib/libopenblasp-r0.2.19.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C             
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.0       lattice_0.20-38  ape_5.2          here_0.1        
    ##  [5] clisymbols_1.2.0 crayon_1.3.4     digest_0.6.18    rprojroot_1.3-2 
    ##  [9] grid_3.5.2       nlme_3.1-137     backports_1.1.3  magrittr_1.5    
    ## [13] evaluate_0.12    stringi_1.3.1    fs_1.2.6         rmarkdown_1.11  
    ## [17] tools_3.5.2      stringr_1.4.0    glue_1.3.0       parallel_3.5.2  
    ## [21] xfun_0.4         yaml_2.2.0       compiler_3.5.2   htmltools_0.3.6 
    ## [25] knitr_1.21       usethis_1.4.0
