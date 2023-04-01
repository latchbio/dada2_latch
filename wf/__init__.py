import subprocess
from pathlib import Path
from datetime import datetime
import shutil

from latch import medium_task, workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import LatchDir, LatchFile
from typing import List

from wf.docs import wf_docs
from wf.types import Sample


@medium_task
def run_dada2(
    samples: List[Sample],
    taxonomy_ref_fasta: LatchFile,
    species_assignment_fasta: LatchFile,
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
    dt_string = datetime.now().strftime("%d/%m/%Y_%H.%M.%S")
    output_dir = f"dada2_{dt_string}"
    output_dirpath = Path(output_dir).resolve()
    output_dirpath.mkdir(parents=True, exist_ok=True)

    _run_cmd = [
        "Rscript",
        "/root/dada2.R",
        str(read_dirpath),
        str(output_dirpath),
        taxonomy_ref_fasta.local_path,
        species_assignment_fasta.local_path,
    ]

    subprocess.run(_run_cmd)

    return LatchDir(str(output_dirpath), f"latch:///dada2_results/{output_dirpath}")


@workflow(wf_docs)
def dada2(
    samples: List[Sample],
    taxonomy_ref_fasta: LatchFile,
    species_assignment_fasta: LatchFile,
) -> LatchDir:
    """Workflow to do X

    Header
    ------

    This is a workflow that does X.

    """
    return run_dada2(
        samples=samples,
        taxonomy_ref_fasta=taxonomy_ref_fasta,
        species_assignment_fasta=species_assignment_fasta,
    )


LaunchPlan(
    dada2,
    "Test Data",
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
