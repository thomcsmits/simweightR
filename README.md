# simweightR: similarity-based counts weighting of T cell receptor data for differential abundance analysis

simweightR is a lightweight R implementation of the data pre-processing pipeline described in [Buytenhuijs et al., "Differential T cell receptor clonotype abundance analysis with similarity-based counts weighting"](https://www.biorxiv.org/content/10.1101/2025.03.28.645951).

The package can be used to pre-process immunological sequencing data from T cell receptor sequencing data. The pre-processing involves an adjustment of count values based on the expression levels of highly similar TCRs, thereby improving detection of differential abundance analysis via various methods, including the Wilcoxon rank-sum test and DESeq2. 

## Installation
```r
remotes::install_github("thomcsmits/simweightR")
library(simweightR)
```

## Usage

simweightR follows [AIRR Standards 1.6](https://docs.airr-community.org/en/stable/datarep/rearrangements.html). Adjusting the counts as follows:

```r
library(simweightR)
x <- simweight(mouse_pbs_srbc_mini)
```

The weighted counts are now in `x$weighted`.

```r
head(x$sequences)
head(x$counts[, 1:4])
head(x$weighted[, 1:4])
```

See our accompanying vignette for more examples.
