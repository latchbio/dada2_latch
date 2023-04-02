library(dada2)

args <- commandArgs(trailingOnly = TRUE)

read_dir <- args[1]
output_dir <- args[2]
taxonomy_ref_fasta <- args[3]
minLen <- as.numeric(args[4])
maxN <- as.numeric(args[5])
minQ <- as.numeric(args[6])
maxEE <- as.numeric(args[7])
truncQ <- as.numeric(args[8])
trimLeft <- as.numeric(args[9])
trimRight <- as.numeric(args[10])
species_assignment_fasta <- args[11]

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <-
  sort(list.files(read_dir, pattern = "1.fastq", full.names = TRUE))
fnRs <-
  sort(list.files(read_dir, pattern = "2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

##### Filter and Trim

filtFs <-
  file.path("/root", "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <-
  file.path("/root", "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(
  fnFs,
  filtFs,
  fnRs,
  filtRs,
  minLen = minLen,
  maxN = maxN,
  minQ = minQ,
  maxEE = maxEE,
  truncQ = truncQ,
  trimLeft = trimLeft,
  trimRight = trimRight,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

##### Dereplicate

derepF1 <- derepFastq(filtFs, verbose = TRUE)
derepR1 <- derepFastq(filtRs, verbose = TRUE)

##### Learn Errors

errF <- learnErrors(derepF1, multithread = TRUE)
errR <- learnErrors(derepR1, multithread = TRUE)

errF_plot <- plotErrors(errF, nominalQ=TRUE)
errR_plot <- plotErrors(errR, nominalQ=TRUE)

ggplot2::ggsave(paste0(output_dir, "/error_rates_forward.pdf"), errF_plot)
ggplot2::ggsave(paste0(output_dir, "/error_rates_reverse.pdf"), errR_plot)

##### Infer ASVs

dadaFs <- dada(derepF1, err = errF, multithread = TRUE)
dadaRs <- dada(derepR1, err = errR, multithread = TRUE)

##### Merge Pairs and Remove Chimeras

mergers <-
  mergePairs(dadaFs, derepF1, dadaRs, derepR1, verbose = TRUE)

seqtab <- makeSequenceTable(mergers)

seqtab.nochim <-
  removeBimeraDenovo(seqtab,
                     method = "consensus",
                     multithread = TRUE,
                     verbose = TRUE)

write.csv(as.data.frame(seqtab.nochim),
          paste0(output_dir, "/asv_table.csv"),
          col.names = NA)

###### Assign taxonomy

taxa <-
  assignTaxonomy(seqtab.nochim, taxonomy_ref_fasta, multithread = TRUE)

if (is.na(species_assignment_fasta) == FALSE) {
  taxa <- addSpecies(taxa, species_assignment_fasta)
}

write.csv(taxa, paste0(output_dir, "/species_table.csv"), col.names = NA)