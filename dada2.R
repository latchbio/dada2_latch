library(dada2)

args <- commandArgs(trailingOnly = TRUE)

sample_name <- args[1]

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files("/root/data/MiSeq_SOP/", pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files("/root/data/MiSeq_SOP/", pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

##### Filter and Trim

filtFs <- file.path("/root/data", "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("/root/data", "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)

##### Dereplicate

derepF1 <- derepFastq(filtFs, verbose=TRUE)
derepR1 <- derepFastq(filtRs, verbose=TRUE)

##### Learn Errors

errF <- learnErrors(derepF1, multithread=TRUE)
errR <- learnErrors(derepR1, multithread=TRUE)

##### Infer ASVs

dadaFs <- dada(derepF1, err=errF, multithread=TRUE)
dadaRs <- dada(derepR1, err=errR, multithread=TRUE)

##### Remove Chimeras

mergers <- mergePairs(dadaFs, derepF1, dadaRs, derepR1, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

write.table(seqtab.nochim, paste0(sample_name, "_asv_table.tsv"))

###### Assign taxonomy

taxa <- assignTaxonomy(seqtab.nochim, "/root/data/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

taxa <- addSpecies(taxa, "/root/data/silva_species_assignment_v138.1.fa.gz")

write.table(taxa, paste0(sample_name, "_species_table.tsv"))