make\_small\_transcriptomes
================
Joel Nitta
5/12/2019

Make small versions of transcriptome files for testing functions.

The transcriptomes will be randomly downsized to the fraction set by `keep_frac`.

Load packages.

``` r
library(rvest)
```

    ## Loading required package: xml2

``` r
library(tidyverse)
```

    ## ── Attaching packages ──────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.1.0       ✔ purrr   0.3.2  
    ## ✔ tibble  2.1.1       ✔ dplyr   0.8.0.1
    ## ✔ tidyr   0.8.1       ✔ stringr 1.4.0  
    ## ✔ readr   1.3.1       ✔ forcats 0.3.0

    ## ── Conflicts ─────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter()         masks stats::filter()
    ## ✖ readr::guess_encoding() masks rvest::guess_encoding()
    ## ✖ dplyr::lag()            masks stats::lag()
    ## ✖ purrr::pluck()          masks rvest::pluck()

``` r
library(XML)
```

    ## 
    ## Attaching package: 'XML'

    ## The following object is masked from 'package:rvest':
    ## 
    ##     xml

``` r
library(baitfindR)
library(ape)
```

Define helper function for downsizing transcriptomes.

``` r
downsize_transcriptome <- function (file, keep_frac) {
  # Read in file. Use bzfile() because it's compressed
  seq <- ape::read.FASTA(bzfile(file))
  # Randomly downsize to specified fraction of transcripts
  seq <- seq[sample(1:length(seq), keep_frac*length(seq))]
}
```

Scrape sample codes and URLs to download transcriptome assemblies from [oneKp website](http://www.onekp.com/public_data.html).

``` r
# First, scrape links contained within oneKp data website.

# Set URL for oneKp data website
url <- "http://www.onekp.com/public_data.html"

# Scrape links for assemblies.
# Use Selectorgadget to identify the CSS component with links.
# Note: run vignette(“selectorgadget”) for a refresher on how to do this.
links <- read_html(url) %>% html_nodes("td:nth-child(5) a")

# Next, read in the main table (just plain text, without links).
tbls_xml <- readHTMLTable(url)

# Extract the table and add links to assemblies,
# keeping only transcriptomes in baitfindR example dataset .
onekp_download_data <- 
  tbls_xml[[1]] %>% 
  as_tibble %>%
  mutate_all(as.character) %>%
  mutate(assembly_link = map_chr(links, xml_attrs)) %>%
  select(code = `1kP_Code`, assembly_link) %>%
  filter(code %in% baitfindR::onekp_data$code)
```

Download transcriptomes.

``` r
purrr::walk2(
  .x = onekp_download_data$assembly_link, 
  .y = fs::path(
    "data-raw", 
    glue::glue("{onekp_download_data$code}-SOAPdenovo-Trans-assembly.fa.bz2")) %>% 
    here::here(),
  ~ download.file(url = .x, destfile = .y)
)
```

Downsize transcriptomes.

``` r
# Set fraction of transcriptome to keep
keep_frac <- 0.05

# Set seed
set.seed(5394)

# Downsize and write out transcriptomes as list.
example_transcriptomes <-
  purrr::map(
    fs::path(
      "data-raw", 
      glue::glue("{onekp_download_data$code}-SOAPdenovo-Trans-assembly.fa.bz2")) %>% 
      here::here() %>%
      set_names(onekp_download_data$code),
    ~ downsize_transcriptome(., keep_frac = keep_frac)
  )
```

Write out list of downsized transcriptomes to `data/`.

``` r
usethis::use_data(example_transcriptomes, overwrite = TRUE)
```

    ## ✔ Setting active project to '/Users/joelnitta/R/baitfindR'
    ## ✔ Saving 'example_transcriptomes' to 'data/example_transcriptomes.rda'

``` r
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
    ##  [1] ape_5.3         baitfindR_0.1.0 XML_3.98-1.16   forcats_0.3.0  
    ##  [5] stringr_1.4.0   dplyr_0.8.0.1   purrr_0.3.2     readr_1.3.1    
    ##  [9] tidyr_0.8.1     tibble_2.1.1    ggplot2_3.1.0   tidyverse_1.2.1
    ## [13] rvest_0.3.2     xml2_1.2.0     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] clisymbols_1.2.0 tidyselect_0.2.5 xfun_0.5         haven_1.1.2     
    ##  [5] lattice_0.20-38  colorspace_1.3-2 usethis_1.4.0    htmltools_0.3.6 
    ##  [9] yaml_2.2.0       rlang_0.3.1      pillar_1.3.1     glue_1.3.1      
    ## [13] withr_2.1.2      selectr_0.4-1    modelr_0.1.2     readxl_1.1.0    
    ## [17] plyr_1.8.4       munsell_0.5.0    gtable_0.2.0     cellranger_1.1.0
    ## [21] evaluate_0.13    knitr_1.21       parallel_3.5.2   curl_3.3        
    ## [25] broom_0.5.0      Rcpp_1.0.1       scales_1.0.0     backports_1.1.3 
    ## [29] jsonlite_1.6     fs_1.2.6         hms_0.4.2        digest_0.6.18   
    ## [33] stringi_1.4.3    rprojroot_1.3-2  grid_3.5.2       here_0.1        
    ## [37] cli_1.1.0        tools_3.5.2      magrittr_1.5     lazyeval_0.2.1  
    ## [41] crayon_1.3.4     pkgconfig_2.0.2  lubridate_1.7.4  assertthat_0.2.1
    ## [45] rmarkdown_1.11   httr_1.4.0       rstudioapi_0.9.0 R6_2.4.0        
    ## [49] nlme_3.1-137     compiler_3.5.2
