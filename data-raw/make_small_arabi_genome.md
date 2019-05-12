make\_small\_arabi\_genome
================
Joel Nitta
5/11/2019

Make small version of Arabidopsis genome file for testing functions.

This only includes the first ca. 100,000 bases of chromosome 1.

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.1.0       ✔ purrr   0.3.2  
    ## ✔ tibble  2.1.1       ✔ dplyr   0.8.0.1
    ## ✔ tidyr   0.8.1       ✔ stringr 1.4.0  
    ## ✔ readr   1.3.1       ✔ forcats 0.3.0

    ## ── Conflicts ────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
# Download genome
download.file(
  url = "ftp://ftp.ensemblgenomes.org/pub/release-40/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz",
  destfile = here::here("data-raw/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz")
  )

# Unzip
R.utils::gunzip(
  filename = here::here("data-raw/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"),
  destname = here::here("data-raw/Arabidopsis_thaliana.TAIR10.dna.toplevel.fasta")
)

# Read in fasta
arabidopsis_genome <- ape::read.FASTA(here::here("data-raw/Arabidopsis_thaliana.TAIR10.dna.toplevel.fasta"))

# Fix names (only take first part of name separated by space)
names(arabidopsis_genome) <-
  names(arabidopsis_genome) %>% str_split(" ") %>% map_chr(1)

# Keep only first 100,000 nucleotides of chromosome 1
arabidopsis_genome_small <- arabidopsis_genome["1"]
arabidopsis_genome_small[["1"]] <- arabidopsis_genome_small[["1"]][1:100000]

# Write out FASTA file to inst/extdata
ape::write.FASTA(
  arabidopsis_genome,
  here::here("inst/extdata/Arabidopsis_thaliana_TAIR10_40_small.fasta"))

sessionInfo()
```

    ## R version 3.5.2 (2018-12-20)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS High Sierra 10.13.6
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] forcats_0.3.0   stringr_1.4.0   dplyr_0.8.0.1   purrr_0.3.2    
    ## [5] readr_1.3.1     tidyr_0.8.1     tibble_2.1.1    ggplot2_3.1.0  
    ## [9] tidyverse_1.2.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_0.2.5  xfun_0.5          haven_1.1.2      
    ##  [4] lattice_0.20-38   colorspace_1.3-2  htmltools_0.3.6  
    ##  [7] yaml_2.2.0        rlang_0.3.1       R.oo_1.22.0      
    ## [10] pillar_1.3.1      glue_1.3.1        withr_2.1.2      
    ## [13] R.utils_2.7.0     modelr_0.1.2      readxl_1.1.0     
    ## [16] plyr_1.8.4        munsell_0.5.0     gtable_0.2.0     
    ## [19] cellranger_1.1.0  rvest_0.3.2       R.methodsS3_1.7.1
    ## [22] evaluate_0.13     knitr_1.21        parallel_3.5.2   
    ## [25] broom_0.5.0       Rcpp_1.0.1        scales_1.0.0     
    ## [28] backports_1.1.3   jsonlite_1.6      hms_0.4.2        
    ## [31] digest_0.6.18     stringi_1.4.3     grid_3.5.2       
    ## [34] rprojroot_1.3-2   here_0.1          cli_1.1.0        
    ## [37] tools_3.5.2       magrittr_1.5      lazyeval_0.2.1   
    ## [40] crayon_1.3.4      ape_5.3           pkgconfig_2.0.2  
    ## [43] xml2_1.2.0        lubridate_1.7.4   assertthat_0.2.1 
    ## [46] rmarkdown_1.11    httr_1.4.0        rstudioapi_0.9.0 
    ## [49] R6_2.4.0          nlme_3.1-137      compiler_3.5.2
