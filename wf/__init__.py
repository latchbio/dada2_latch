import subprocess
from pathlib import Path
from datetime import datetime
import shutil
from urllib.parse import urlparse
from latch import medium_task, workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import LatchDir, LatchFile
from typing import List, Optional

from wf.docs import wf_docs
from wf.types import (
    Sample,
    TaxonomyReference,
    SpeciesAssignmentReference,
    TAXONOMY_REF_DICT,
    SPECIES_ASSIGN_DICT,
)


def aws_cp(src: str, dst: str, show_progress: bool = False):
    cmd = ["aws", "s3", "cp", src, dst]
    if not show_progress:
        cmd.append("--no-progress")
    subprocess.run(cmd, check=True)


@medium_task
def run_dada2(
    samples: List[Sample],
    taxonomy_reference: TaxonomyReference,
    species_assignment: Optional[SpeciesAssignmentReference],
    taxonomy_ref_fasta: Optional[LatchFile],
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

    if taxonomy_ref_fasta:
        taxonomy_loc = taxonomy_ref_fasta.local_path
    else:
        remote_path = TAXONOMY_REF_DICT[taxonomy_reference.value]
        local_path = Path.cwd() / Path(urlparse(remote_path).path).name
        taxonomy_loc = str(local_path.resolve())
        aws_cp(remote_path, taxonomy_loc)

    if species_assignment:
        remote_path = SPECIES_ASSIGN_DICT[species_assignment.value]
        local_path = Path.cwd() / Path(urlparse(remote_path).path).name
        species_assign_loc = str(local_path.resolve())
        aws_cp(remote_path, species_assign_loc)
    elif species_assignment_fasta:
        species_assign_loc = species_assignment_fasta.local_path
    else:
        species_assign_loc = None

    # Move all reads into respective directory
    read_dirpath1 = Path("dada2_forwardreads").resolve()
    read_dirpath2 = Path("dada2_reversereads").resolve()
    read_dirpath1.mkdir(parents=True, exist_ok=True)
    read_dirpath2.mkdir(parents=True, exist_ok=True)

    for sample in samples:
        shutil.move(sample.read1.local_path, str(read_dirpath1))
        shutil.move(sample.read2.local_path, str(read_dirpath2))

    # Run dada2
    dt_string = datetime.now().strftime("%d.%m.%Y_%H.%M.%S")
    output_dir = f"dada2_{dt_string}"
    output_dirpath = Path(output_dir).resolve()
    output_dirpath.mkdir(parents=True, exist_ok=True)

    _run_cmd = [
        "Rscript",
        "/root/dada2.R",
        str(read_dirpath1),
        str(read_dirpath2),
        str(output_dirpath),
        taxonomy_loc,
        str(minLen),
        str(maxN),
        str(minQ),
        str(maxEE),
        str(truncQ),
        str(trimLeft),
        str(trimRight),
    ]

    if species_assign_loc:
        _run_cmd.append(species_assign_loc)

    subprocess.run(_run_cmd)

    return LatchDir(str(output_dirpath), f"latch:///dada2_results/{output_dir}")


@workflow(wf_docs)
def dada2(
    samples: List[Sample],
    tax_ref_fork: str,
    taxonomy_reference: TaxonomyReference,
    species_assign_fork: str = "none",
    species_assignment: Optional[SpeciesAssignmentReference] = None,
    taxonomy_ref_fasta: Optional[LatchFile] = None,
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

    [^1]: Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016).
    “DADA2: High-resolution sample inference from Illumina amplicon data.”
    Nature Methods, 13, 581-583. https://doi.org/10.1038/nmeth.3869.
    """
    return run_dada2(
        samples=samples,
        taxonomy_reference=taxonomy_reference,
        species_assignment=species_assignment,
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
        "taxonomy_reference": TaxonomyReference.UNITE,
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
                    "s3://latch-public/test-data/4318/F3D8_S196_L001_1.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/F3D8_S196_L001_2.fastq"
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
        "taxonomy_reference": TaxonomyReference.SILVA,
        "species_assignment": SpeciesAssignmentReference.SILVA,
    },
)
