## DADA2

This is a workflow to run the DADA2[^1] (v1.26.0) software for amplicon data.
It currently runs the following DADA2 functions under the hood:

- `filterAndTrim` for filtering and trimming sequence files
- `derepFastq` for dereplication
- `learnErrors` to learn the error rates from the data
- `dada` the main DADA2 function, to infer sample composition
- `mergePairs` to merge paired-end reads
- `makeSequenceTable` to build the per-sample sequence table
- `removeBimeraDenovo` to remove sequencing chimeras
- `assignTaxonomy` and `addSpecies` for taxonomic classification of ASVs

The workflow outputs two data files, one corresponding to the per-sequence table after chimera removal
(`asv_table.csv`) and one for the taxonomic classification (`species_table.csv`)

[^1]:
    Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016).
    “DADA2: High-resolution sample inference from Illumina amplicon data.”
    Nature Methods, 13, 581-583. https://doi.org/10.1038/nmeth.3869.
