import subprocess
from pathlib import Path
from datetime import datetime
import shutil

from latch import medium_task, workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import LatchDir, LatchFile
from typing import List, Optional

from wf.docs import wf_docs
from wf.types import Sample


@medium_task
def run_dada2(
    samples: List[Sample],
    taxonomy_ref_fasta: LatchFile,
    species_assignment_fasta: Optional[LatchFile],
    minLen: int,
    maxN: int,
    minQ: int,
    maxEE: int,
    truncQ: int,
    trimLeft: int,
    trimRight: int,
) -> LatchDir:
    """Task to run dada2"""

    # Move all reads into same directory
    read_dir = "dada2_reads"
    read_dirpath = Path(read_dir).resolve()
    read_dirpath.mkdir(parents=True, exist_ok=True)

    for sample in samples:
        shutil.move(sample.read1.local_path, str(read_dirpath))
        shutil.move(sample.read2.local_path, str(read_dirpath))

    # Run dada2
    dt_string = datetime.now().strftime("%d.%m.%Y_%H.%M.%S")
    output_dir = f"dada2_{dt_string}"
    output_dirpath = Path(output_dir).resolve()
    output_dirpath.mkdir(parents=True, exist_ok=True)

    _run_cmd = [
        "Rscript",
        "/root/dada2.R",
        str(read_dirpath),
        str(output_dirpath),
        taxonomy_ref_fasta.local_path,
        str(minLen),
        str(maxN),
        str(minQ),
        str(maxEE),
        str(truncQ),
        str(trimLeft),
        str(trimRight),
    ]

    if species_assignment_fasta:
        _run_cmd.append(species_assignment_fasta.local_path)

    subprocess.run(_run_cmd)

    return LatchDir(str(output_dirpath), f"latch:///dada2_results/{output_dir}")


@workflow(wf_docs)
def dada2(
    samples: List[Sample],
    taxonomy_ref_fasta: LatchFile,
    species_assignment_fasta: Optional[LatchFile] = None,
    minLen: int = 50,
    maxN: int = 0,
    minQ: int = 0,
    maxEE: int = 2,
    truncQ: int = 2,
    trimLeft: int = 0,
    trimRight: int = 0,
) -> LatchDir:
    """A workflow for fast and accurate sample inference from amplicon data

    DADA2
    ------

    This is a workflow to run the DADA2[^1] software for amplicon data.
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

    [^1]: Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016).
    “DADA2: High-resolution sample inference from Illumina amplicon data.”
    Nature Methods, 13, 581-583. https://doi.org/10.1038/nmeth.3869.
    """
    return run_dada2(
        samples=samples,
        taxonomy_ref_fasta=taxonomy_ref_fasta,
        species_assignment_fasta=species_assignment_fasta,
        minLen=minLen,
        maxN=maxN,
        minQ=minQ,
        maxEE=maxEE,
        truncQ=truncQ,
        trimLeft=trimLeft,
        trimRight=trimRight,
    )


LaunchPlan(
    dada2,
    "Fungus metagenome",
    {
        "samples": [
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/PRJNA377530_trimmed/SRR5314314_trim_1.fastq.gz"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/PRJNA377530_trimmed/SRR5314314_trim_2.fastq.gz"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/PRJNA377530_trimmed/SRR5314336_trim_1.fastq.gz"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/PRJNA377530_trimmed/SRR5314336_trim_2.fastq.gz"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/PRJNA377530_trimmed/SRR5838532_trim_1.fastq.gz"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/PRJNA377530_trimmed/SRR5838532_trim_2.fastq.gz"
                ),
            ),
        ],
        "taxonomy_ref_fasta": LatchFile(
            "s3://latch-public/test-data/4318/sh_general_release_dynamic_29.11.2022.fasta"
        ),
    },
)

LaunchPlan(
    dada2,
    "Human Gut Microbiome",
    {
        "samples": [
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/F3D0_S188_L001_1.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/F3D0_S188_L001_2.fastq"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/F3D149_S215_L001_1.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/F3D149_S215_L001_2.fastq"
                ),
            ),
        ],
        "taxonomy_ref_fasta": LatchFile(
            "s3://latch-public/test-data/4318/silva_nr99_v138.1_train_set.fa.gz"
        ),
        "species_assignment_fasta": LatchFile(
            "s3://latch-public/test-data/4318/silva_species_assignment_v138.1.fa.gz"
        ),
    },
)
