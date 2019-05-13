# baitfindR

`baitfindR` is an R package that extends the [Yang and Smith orthology inference pipeline](https://bitbucket.org/yangya/phylogenomic_dataset_construction/overview) (2014 Mol. Biol. Evol. 31:3081-3092) to find appropriate loci for use as baits in sequence capture.

## Installation

You can install `baitfindR` as follows:

``` r
install.packages("devtools")
devtools::install_github("joelnitta/baitfindR")
```

## Examples

`baitfindR` works best in conjunction with [drake](https://ropensci.github.io/drake/) to manage a bait-finding workflow. Please see [the example repo](https://github.com/joelnitta/baitfindR_example).
