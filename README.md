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
```
results <- TCRsimilift_calculate(mouse_PBSvTCZ_data, sim_method="HAMMING", cutoff = 0.77, export_results=TRUE, output_directory="my_outputs")
```

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.

## License
For open source projects, say how it is licensed.
