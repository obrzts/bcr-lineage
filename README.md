# bcr-lineage
B cell clonal group and tree inference

*parse_mixcr.R* - reads MiXCR output and extract mutations (either from alignment or from mutationsDetailed* if available)
*clusterize.R* - takes clonotype table as input, run MIR and defines cluster (cluster_id added as a new column to clonotype table)
