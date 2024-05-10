library(dada2)
library(phyloseq)

args <- commandArgs(trailingOnly = TRUE)

# read_dir1 <- args[1]
# read_dir2 <- args[2]
# output_dir <- args[3]
# taxonomy_ref_fasta <- args[4]
# minLen <- as.numeric(args[5])
# maxN <- as.numeric(args[6])
# minQ <- as.numeric(args[7])
# maxEE <- as.numeric(args[8])
# truncQ <- as.numeric(args[9])
# trimLeft <- as.numeric(args[10])
# trimRight <- as.numeric(args[11])
# filtering_setting <- args[12]
# maxMismatch <- as.numeric(args[13])
# minOverlap <- as.numeric(args[14])
# derep <- args[15]
# chimeric_setting <- args[16]
# species_assignment_fasta <- args[17]

reads_dir <- args[1]
output_dir <- args[2]
taxonomy_ref_fasta <- args[3]
truncLen <- as.numeric(args[4])
OMEGA_A <- as.numeric(args[5])
singletons_setting <- args[6]
derep <- args[7]
chimeric_setting <- args[8]
pooling_param <- args[9]
merge_reads <- args[10]
maxMismatch <- as.numeric(args[11])
minOverlap <- as.numeric(args[12])
minLen <- as.numeric(args[13])
species_assignment_fasta <- args[14]


print("Read files found:")
list.files(reads_dir)

merge_setting = FALSE

if (merge_reads == "merge") {
  merge_setting = TRUE

}


all_files <- list.files(reads_dir, full.names = TRUE, recursive = TRUE)

# Check if "_R1_001.fastq" is present first in any of the files
has_R1_001_first <- any(grepl("_R1_001\\.fastq", all_files) & grepl("^.*_R1_001\\.fastq", all_files))

# Print the result
if (has_R1_001_first) {
  message("At least one file has '_R1_001.fastq' present first.\n Checking for reads that contain _R1_001.fastq and _R2_001.fastq.")
  fnFs <- sort(list.files(reads_dir, pattern="_R1_001.fastq", full.names = TRUE, recursive=TRUE))
  fnRs <- sort(list.files(reads_dir, pattern="_R2_001.fastq", full.names = TRUE,recursive=TRUE))

  order_check <- all(gsub("_R1_001.fastq", "", fnFs) == gsub("_R2_001.fastq", "", fnRs))

  # Take action if order_check is FALSE
  if (!order_check) {
    # Code to handle the case where the order is not correct
    message("Error: The forward and reverse files are not in the same order.\n")
    fnFs <- sort(list.files(reads_dir, pattern="_R1_001.fastq", full.names = TRUE, recursive = TRUE))
    fnRs <- sort(list.files(reads_dir, pattern="_R2_001.fastq", full.names = TRUE, recursive = TRUE))
    # You can add additional code here to handle the error, such as reordering the files or exiting the script.
  }
} else {
  message("No file has '_R1_001.fastq' present.\n Checking for reads that contain _1.fastq and _2.fastq.")
  fnFs <- sort(list.files(reads_dir, pattern="_1.fastq", full.names = TRUE, recursive=TRUE))
  fnRs <- sort(list.files(reads_dir, pattern="_2.fastq", full.names = TRUE,recursive=TRUE))

  order_check <- all(gsub("_1.fastq", "", fnFs) == gsub("_2.fastq", "", fnRs))

  # Take action if order_check is FALSE
  if (!order_check) {
    # Code to handle the case where the order is not correct
    message("Error: The forward and reverse files are not in the same order.\n")
    fnFs <- sort(list.files(reads_dir, pattern="_1.fastq", full.names = TRUE, recursive = TRUE))
    fnRs <- sort(list.files(reads_dir, pattern="_2.fastq", full.names = TRUE, recursive = TRUE))
    # You can add additional code here to handle the error, such as reordering the files or exiting the script.
  }
}

single_sample = FALSE

if (length(fnFs) == 1) {
  single_sample = TRUE
}

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
print(sample.names)
# sample.names <- tools::file_path_sans_ext(basename(fnFs))

qual_plot_f <- plotQualityProfile(fnFs[1])
qual_plot_r <- plotQualityProfile(fnRs[1])
ggplot2::ggsave(paste0(output_dir, "/qual_plot_f.pdf"), qual_plot_f)
ggplot2::ggsave(paste0(output_dir, "/qual_plot_r.pdf"), qual_plot_r)

##### Filter and Trim

message("--- Filtering and Trimming Data ---")

# maxEE = c(200,200)

# out <- filterAndTrim(
#   fnFs,
#   filtFs,
#   fnRs,
#   filtRs,
#   minLen = 10,
#   maxN = 3,
#   # minQ = minQ,
#   # # maxEE = maxEE,
#   # truncQ = 1,
#   # truncLen = c(145,145),
#   trimLeft = c(41, 45),
#   trimRight = c(10, 10),
#   # trimRight = trimRight,
#   rm.phix = FALSE,
#   compress = TRUE,
#   multithread = TRUE,
# )

filtFs <- file.path("/root", "filtered", basename(fnFs))
filtRs <- file.path("/root", "filtered", basename(fnRs))

# names(filtFs) <- sample.names
# names(filtRs) <- sample.names


out <- filterAndTrim(
  fnFs,
  filtFs,
  fnRs,
  filtRs,
  truncLen = truncLen,
  compress = TRUE,
  multithread = TRUE,
  minLen = minLen,
)
print(head(out))



##### Dereplicate


if (derep == "derep") {
  message("--- Dereplicating FASTQ data ---")

  derepF1 <- derepFastq(filtFs, verbose = TRUE)
  derepR1 <- derepFastq(filtRs, verbose = TRUE)
} else {
  message("--- Not dereplicating FASTQ data ---")
  derepF1 <- filtFs
  derepR1 <- filtRs
}


##### Learn Errors

message("--- Learning Error Rates from the Dataset ---")

errF <- learnErrors(derepF1, multithread = TRUE)
errR <- learnErrors(derepR1, multithread = TRUE)

errF_plot <- plotErrors(errF, nominalQ=TRUE)
errR_plot <- plotErrors(errR, nominalQ=TRUE)

ggplot2::ggsave(paste0(output_dir, "/error_rates_forward.pdf"), errF_plot)
ggplot2::ggsave(paste0(output_dir, "/error_rates_reverse.pdf"), errR_plot)

##### Infer ASVs

if (pooling_param == "TRUE") {
  pooling_param_val = TRUE
} else if (pooling_param == "FALSE") {
  pooling_param_val = FALSE
} else {
  pooling_param_val = "pseudo"
}


message("--- Inferring sample composition ---")
if (singletons_setting == "detect_singletons") {
  dadaFs <- dada(derepF1, err = errF, OMEGA_A=OMEGA_A, pool=pooling_param_val, DETECT_SINGLETONS=TRUE, multithread = TRUE)
  dadaRs <- dada(derepR1, err = errR, OMEGA_A=OMEGA_A, pool=pooling_param_val, DETECT_SINGLETONS=TRUE, multithread = TRUE)
} else {
  dadaFs <- dada(derepF1, err = errF, OMEGA_A=OMEGA_A, pool=pooling_param_val, DETECT_SINGLETONS=FALSE, multithread = TRUE)
  dadaRs <- dada(derepR1, err = errR, OMEGA_A=OMEGA_A, pool=pooling_param_val, DETECT_SINGLETONS=FALSE, multithread = TRUE)
}

##### Merge Pairs and Remove Chimeras

counter = 1

if (single_sample == TRUE) {

  if (merge_setting == TRUE) {
    message("--- Merging reads ---")
    mergers <-
      mergePairs(dadaFs, derepF1, dadaRs, derepR1, verbose = TRUE, maxMismatch=maxMismatch, minOverlap=minOverlap)
    seqtab <- makeSequenceTable(mergers)
  } else {
    seqtab <- makeSequenceTable(dadaFs)
  }

  sample_name = sample.names[1]

  if (chimeric_setting == "remove") {
      message("--- Removing chimeras ---")

      seqtab.nochim <-
        removeBimeraDenovo(seqtab,
                          method = "consensus",
                          multithread = TRUE,
                          verbose = TRUE)

      print(dim(seqtab.nochim))
      print(sum(seqtab.nochim)/sum(seqtab))
    } else {
      message("--- Not removing chimeras ---")
      seqtab.nochim = seqtab
    }

    subDir <- sample_name
    folder_name <- strsplit(subDir, "_")[[1]][1]

    dir.create(file.path(output_dir, folder_name), showWarnings = FALSE)

    write.csv(as.data.frame(seqtab.nochim),
          paste0(output_dir, "/", folder_name, "/asv_table.csv"),
          col.names = NA)

    message("--- Assigning taxonomy to ASVs ---")

    taxa <-
      assignTaxonomy(seqtab.nochim, taxonomy_ref_fasta, multithread = TRUE, minBoot=0, tryRC=TRUE)

    if (is.na(species_assignment_fasta) == FALSE) {
      message("--- Assigning species to ASVs ---")
      taxa <- addSpecies(taxa, species_assignment_fasta)
    }

    message(paste0(output_dir, "/", folder_name, "/species_table.csv"))
    write.csv(taxa, paste0(output_dir, "/", folder_name, "/species_table.csv"), col.names = NA)


    ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
                  tax_table(as.matrix(taxa)))

    save(ps, file = paste0(output_dir, "/", folder_name, "/phyloseq_object.RData"))


} else {

  if (merge_setting == TRUE) {
    message("--- Merging reads ---")
    mergers <-
        mergePairs(dadaFs, derepF1, dadaRs, derepR1, verbose = TRUE, maxMismatch=maxMismatch, minOverlap=minOverlap)
    all_samples = mergers
  } else {
    all_samples = dadaFs
  }



  for (sample in all_samples) {



      sample_name = names(dadaFs)[counter]
      print(sample_name)
      counter = counter +1
      seqtab <- makeSequenceTable(sample)

    if (chimeric_setting == "remove") {
      message("--- Removing chimeras ---")

      seqtab.nochim <-
        removeBimeraDenovo(seqtab,
                          method = "consensus",
                          multithread = TRUE,
                          verbose = TRUE)

      print(dim(seqtab.nochim))
      print(sum(seqtab.nochim)/sum(seqtab))
    } else {
      message("--- Not removing chimeras ---")
      seqtab.nochim = seqtab
    }

    subDir <- sample_name
    folder_name <- strsplit(subDir, "_")[[1]][1]

    dir.create(file.path(output_dir, folder_name), showWarnings = FALSE)

    write.csv(as.data.frame(seqtab.nochim),
          paste0(output_dir, "/", folder_name, "/asv_table.csv"),
          col.names = NA)

    message("--- Assigning taxonomy to ASVs ---")

    taxa <-
      assignTaxonomy(seqtab.nochim, taxonomy_ref_fasta, multithread = TRUE, minBoot=0, tryRC=TRUE)

    if (is.na(species_assignment_fasta) == FALSE) {
      message("--- Assigning species to ASVs ---")
      taxa <- addSpecies(taxa, species_assignment_fasta)
    }

    message(paste0(output_dir, "/", folder_name, "/species_table.csv"))
    write.csv(taxa, paste0(output_dir, "/", folder_name, "/species_table.csv"), col.names = NA)


    ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
                  tax_table(as.matrix(taxa)))

    save(ps, file = paste0(output_dir, "/", folder_name, "/phyloseq_object.RData"))

}

}
