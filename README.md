# TCRsimilift

## Name
TCRsimilift - TCR sequencing data preprocessing via similarity-based count adjustment

## Description
The TCRsimilift package is a lightweight R implementation of the data pre-processing pipeline described in [Buytenhuijs et al. (2024). "Differential T cell receptor gene expression analysis using the Wilcoxon test with similarity-based weighting"](https://www.biorxiv.org/content/10.1101/2025.03.28.645951v1).

The package can be used to pre-process immunological sequencing data from T Cell Receptor samples. The pre-processing involves an adjustment of count values based on the expression levels of highly similar TCRs, thereby improving detection of differential gene expression via various methods, including edgeR, deseq2 and the Wilcoxon test. 

In essence, the approach involves the calculation of similarity scores between each TCR sequence and all other similar sequences in a given dataset. The counts of each TCR are then adjusted based on a weighted average of their neighbours in the similarity graph. For details regarding the theory behind the method, please see the pre-print.

## Installation
```
library(TCRsimilift)
```

## Usage

The easiest way to use the package is by running the `TCRsimilift_calculate()` function on a dataframe containing TCR sequencing data, where the column namings respect the [AIRR Standards 1.6](https://docs.airr-community.org/en/stable/datarep/rearrangements.html). This function will output a dataframe containing aggregated counts, with an additional column `wrc` that contains the adjusted counts based on similarity scores. Please see the help file of the `TCRsimilift_calculate()` for details about all relevant parameters and their defaults.

Here is an example of running this function on the example dataset provided with the package:

```
results <- TCRsimilift_calculate(mouse_PBSvTCZ_data, sim_method="HAMMING", cutoff = 0.77, export_results=TRUE, output_directory="my_outputs")
```

### Alternatively: running the functions in the pipeline separately

It is also possible to run the functions of this pipeline manually to see the outputs created along the way. These are as follows:
 
* `datacheck()` : will throw an error if the data is formatted incorrectly. 
* `dataprep()` : creates a dataframe with the appropriate columns and contents to process downstream.
* `net_update_data()` : adjust counts based on similarity, sample by sample, and return results.

If a different method for count adjustment is desired than the weighted average we use, it is also possible to call the functions that output the similarity matrices directly via `blosum_similarity()` or `hamming_similarity()`. See respective help files for details. Bear in mind that all cdr3 amino acid sequences fed to these two functions must have the same length in order to be valid inputs.

In the package code, the similarity is calculated for all TCRs with the same CDR3 amino acid sequence length separately, via a loop like so:
```
for (seq_length in all_sequence_lengths) {
  # calculate similarity matrices and adjust counts
}
```

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.

## License
For open source projects, say how it is licensed.
