# simweightR: similarity-based counts weighting of T-cell receptor data for differential abundance analysis

simweightR is a lightweight R implementation of the data pre-processing pipeline described in [Buytenhuijs et al. (2025). "Differential T cell receptor gene expression analysis using the Wilcoxon test with similarity-based weighting"](https://www.biorxiv.org/content/10.1101/2025.03.28.645951v1).

The package can be used to pre-process immunological sequencing data from T-cell receptor sequencing data. The pre-processing involves an adjustment of count values based on the expression levels of highly similar TCRs, thereby improving detection of differential abundance analysis via various methods, including the Wilcoxon signed-rank test and DESeq2. 

## Installation
```r
remotes::install_github("thomcsmits/simweightR")
library(simweightR)
```

## Usage

simweightR follows [AIRR Standards 1.6](https://docs.airr-community.org/en/stable/datarep/rearrangements.html). Adjusting the counts as follows:

```r
library(simweightR)
results <- adjust_counts(data)
```

The weighted counts are now in `results$wrc`. For integration with downstream methods that require a matrix, simply convert to the appropriate format using:

```r
count_matrix <- as_counts_matrix(results)
```

See our accompanying vignette for more examples.
