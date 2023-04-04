from wf import dada2
from wf.types import Sample, TaxonomyReference, SpeciesAssignmentReference
from latch.types import LatchFile


dada2(
    samples=[
        Sample(
            read1=LatchFile("s3://latch-public/test-data/4318/F3D0_S188_L001_1.fastq"),
            read2=LatchFile("s3://latch-public/test-data/4318/F3D0_S188_L001_2.fastq"),
        ),
        Sample(
            read1=LatchFile("s3://latch-public/test-data/4318/F3D8_S196_L001_1.fastq"),
            read2=LatchFile("s3://latch-public/test-data/4318/F3D8_S196_L001_2.fastq"),
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
    taxonomy_reference=TaxonomyReference.SILVA_138_1,
    species_assignment=SpeciesAssignmentReference.SILVA_138_1,
)

dada2(
    samples=[
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
    taxonomy_ref_fasta=LatchFile(
        "s3://latch-public/test-data/4318/sh_general_release_dynamic_29.11.2022.fasta"
    ),
)
