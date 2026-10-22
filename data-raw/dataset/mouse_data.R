# Creates example datasets mouse_pbs_srbc and mouse_pbs_srbc_mini using mouse TCR-seq from da manuscript
# Run as: Rscript data-raw/dataset/mouse_data.R

raw_dir <- "data-raw/dataset"

# Use TCZ, first 4 mice
zone <- "TCZ"
n_mice <- 4
# use one V gene for whole, and 2 J genes for mini
v_gene <- "TRBV2"
mini_j_genes <- c("TRBJ1-1", "TRBJ2-1")
mini_length <- 9

cols <- c("CDR3.amino.acid.sequence", "Read.count", "V.segments", "J.segments", "treatment", "zone", "mouse")
raw <- rbind(read.csv(file.path(raw_dir, "pbs.csv"))[, cols],
             read.csv(file.path(raw_dir, "srbc.csv"))[, cols])
raw <- raw[raw$zone == zone, ]

mice <- unlist(lapply(c("PBS", "SRBC_4d"), function(tr) {
  head(sort(unique(raw$mouse[raw$treatment == tr])), n_mice)
}))
raw <- raw[raw$mouse %in% mice & raw$V.segments == v_gene, ]

# AIRR column names
mouse_pbs_srbc <- data.frame(
  cdr3_aa = raw$CDR3.amino.acid.sequence,
  v_call = raw$V.segments,
  j_call = raw$J.segments,
  consensus_count = raw$Read.count,
  sample_processing_id = paste0(raw$mouse, raw$zone, raw$treatment),
  treatment = raw$treatment
)
rownames(mouse_pbs_srbc) <- NULL

mini <- mouse_pbs_srbc$j_call %in% mini_j_genes & nchar(mouse_pbs_srbc$cdr3_aa) == mini_length
mouse_pbs_srbc_mini <- mouse_pbs_srbc[mini, ]
rownames(mouse_pbs_srbc_mini) <- NULL

usethis::use_data(mouse_pbs_srbc, mouse_pbs_srbc_mini, overwrite = TRUE)
