from wf import dada2
from wf.types import Sample
from latch.types import LatchFile


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

# dada2(
#     samples=[
#         Sample(
#             read1=LatchFile("s3://latch-public/test-data/4318/F3D0_S188_L001_1.fastq"),
#             read2=LatchFile("s3://latch-public/test-data/4318/F3D0_S188_L001_2.fastq"),
#         ),
#         Sample(
#             read1=LatchFile("s3://latch-public/test-data/4318/F3D8_S196_L001_1.fastq"),
#             read2=LatchFile("s3://latch-public/test-data/4318/F3D8_S196_L001_2.fastq"),
#         ),
#         Sample(
#             read1=LatchFile(
#                 "s3://latch-public/test-data/4318/F3D149_S215_L001_1.fastq"
#             ),
#             read2=LatchFile(
#                 "s3://latch-public/test-data/4318/F3D149_S215_L001_2.fastq"
#             ),
#         ),
#     ],
#     taxonomy_ref_fasta=LatchFile(
#         "s3://latch-public/test-data/4318/silva_nr99_v138.1_train_set.fa.gz"
#     ),
#     species_assignment_fasta=LatchFile(
#         "s3://latch-public/test-data/4318/silva_species_assignment_v138.1.fa.gz"
#     ),
# )
