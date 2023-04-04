from dataclasses import dataclass
from enum import Enum

from dataclasses_json import dataclass_json
from latch.types import LatchFile


class TaxonomyReference(Enum):
    SILVA_138_1 = "s3://latch-public/test-data/4318/silva_nr99_v138.1_train_set.fa.gz"
    RDP = "s3://latch-public/test-data/4318/rdp_train_set_18.fa.gz"
    UNITE = (
        "s3://latch-public/test-data/4318/sh_general_release_dynamic_29.11.2022.fasta"
    )


class SpeciesAssignmentReference(Enum):
    SILVA_138_1 = (
        "s3://latch-public/test-data/4318/silva_species_assignment_v138.1.fa.gz"
    )
    RDP = "s3://latch-public/test-data/4318/rdp_species_assignment_18.fa.gz"


@dataclass_json
@dataclass
class Sample:
    read1: LatchFile
    read2: LatchFile
