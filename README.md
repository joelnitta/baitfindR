# baitfindR

`baitfindR` is an [R package](https://www.r-project.org/about.html) that extends the [Yang and Smith orthology inference pipeline](https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview) (2014 Mol. Biol. Evol. 31:3081-3092) to find appropriate loci for use as baits in sequence capture.

## Installation

You can install `baitfindR` as follows:

``` r
install.packages("devtools")
devtools::install_github("joelnitta/baitfindR")
```

## Dependencies

Most functions are wrappers for other external (to R) programs, so these must be installed and on the users' `PATH`. `baitfindR` requires all of the [dependencies for the  Yang & Smith workflow](https://bitbucket.org/yangya/phylogenomic_dataset_construction/src/master/tutorials/part1_dependencies.md).

## Docker image

[A docker image](https://hub.docker.com/r/joelnitta/baitfindr) is available that includes the latest version of `baitfindR` and all required, properly versioned dependencies.

## Examples

`baitfindR` works best in conjunction with [drake](https://ropensci.github.io/drake/) to manage a bait-finding workflow. Please see [the example project](https://github.com/joelnitta/baitfindR_example).
