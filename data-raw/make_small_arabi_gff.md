make\_small\_arabi\_gff
================
Joel Nitta
5/11/2019

Make small version of Arabidopsis GFF file for testing bed functions.

This only includes the first ca. 1000 features in the GFF file, which all happen to be on chromosome 1.

``` r
# Download gff3
download.file(
  url = "ftp://ftp.ensemblgenomes.org/pub/release-40/plants/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.40.gff3.gz",
  destfile = here::here("data-raw/Arabidopsis_thaliana.TAIR10.40.gff3.gz")
)

# Unzip
R.utils::gunzip(
  filename = here::here("data-raw/Arabidopsis_thaliana.TAIR10.40.gff3.gz"),
  destname = here::here("data-raw/Arabidopsis_thaliana.TAIR10.40.gff3")
)

# Copy first 1000 lines of file to inst/extdata/
arabidopsis_gff <- readr::read_lines(here::here("data-raw/Arabidopsis_thaliana.TAIR10.40.gff3"))

arabidopsis_gff_small <- arabidopsis_gff[1:1000]

readr::write_lines(
  arabidopsis_gff_small,
  here::here("inst/extdata/Arabidopsis_thaliana_TAIR10_40_small.gff3")
)

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
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.1        here_0.1          crayon_1.3.4     
    ##  [4] digest_0.6.18     rprojroot_1.3-2   R.methodsS3_1.7.1
    ##  [7] R6_2.4.0          backports_1.1.3   magrittr_1.5     
    ## [10] evaluate_0.13     pillar_1.3.1      rlang_0.3.1      
    ## [13] stringi_1.4.3     R.oo_1.22.0       R.utils_2.7.0    
    ## [16] rmarkdown_1.11    tools_3.5.2       readr_1.3.1      
    ## [19] stringr_1.4.0     hms_0.4.2         xfun_0.5         
    ## [22] yaml_2.2.0        compiler_3.5.2    pkgconfig_2.0.2  
    ## [25] htmltools_0.3.6   knitr_1.21        tibble_2.1.1
