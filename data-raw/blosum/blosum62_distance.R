# BLOSUM62 distance score matrix
# we normalized original target probabilities (Henikoff & Henikoff, 1992) per amino acid to get the probability vector
# for each combination of amino acid, we computed the squared distance between the probability vectors, and normalized to [0,1]

blosum62_distance <- as.matrix(read.csv("data-raw/blosum/blosum_derived_score_hen.csv", row.names = 1))
usethis::use_data(blosum62_distance, overwrite = TRUE)
