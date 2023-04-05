library(dada2)
library(phyloseq)

args <- commandArgs(trailingOnly = TRUE)

read_dir1 <- args[1]
read_dir2 <- args[2]
output_dir <- args[3]
taxonomy_ref_fasta <- args[4]
minLen <- as.numeric(args[5])
maxN <- as.numeric(args[6])
minQ <- as.numeric(args[7])
maxEE <- as.numeric(args[8])
truncQ <- as.numeric(args[9])
trimLeft <- as.numeric(args[10])
trimRight <- as.numeric(args[11])
species_assignment_fasta <- args[12]

fnFs <-
  sort(list.files(read_dir1, full.names = TRUE))
fnRs <-
  sort(list.files(read_dir2, full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

##### Filter and Trim

filtFs <-
  file.path("/root", "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <-
  file.path("/root", "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

message("--- Filtering and Trimming Data ---")

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

message("--- Dereplicating FASTQ data ---")

derepF1 <- derepFastq(filtFs, verbose = TRUE)
derepR1 <- derepFastq(filtRs, verbose = TRUE)

##### Learn Errors

message("--- Learning Error Rates from the Dataset ---")

errF <- learnErrors(derepF1, multithread = TRUE)
errR <- learnErrors(derepR1, multithread = TRUE)

errF_plot <- plotErrors(errF, nominalQ=TRUE)
errR_plot <- plotErrors(errR, nominalQ=TRUE)

ggplot2::ggsave(paste0(output_dir, "/error_rates_forward.pdf"), errF_plot)
ggplot2::ggsave(paste0(output_dir, "/error_rates_reverse.pdf"), errR_plot)

##### Infer ASVs

message("--- Inferring sample composition ---")

dadaFs <- dada(derepF1, err = errF, multithread = TRUE)
dadaRs <- dada(derepR1, err = errR, multithread = TRUE)

##### Merge Pairs and Remove Chimeras

message("--- Merging read pairs and removing chimeras ---")

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

message("--- Assigning taxonomy to ASVs ---")

taxa <-
  assignTaxonomy(seqtab.nochim, taxonomy_ref_fasta, multithread = TRUE)

if (is.na(species_assignment_fasta) == FALSE) {
  taxa <- addSpecies(taxa, species_assignment_fasta)
}

write.csv(taxa, paste0(output_dir, "/species_table.csv"), col.names = NA)


ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               tax_table(as.matrix(taxa)))

save(ps, file = paste0(output_dir, "/phyloseq_object.RData"))
