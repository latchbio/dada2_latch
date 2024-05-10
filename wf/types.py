from dataclasses import dataclass
from enum import Enum

from latch.types import LatchFile

TAXONOMY_REF_DICT = {
    "SILVA": "s3://latch-public/test-data/4318/silva_nr99_v138.1_train_set.fa.gz",
    "RDP": "s3://latch-public/test-data/4318/rdp_train_set_18.fa.gz",
    "UNITE": "s3://latch-public/test-data/4318/sh_general_release_dynamic_29.11.2022.fasta",
    "GREENGENES": "s3://latch-public/test-data/23882/gg_13_8_train_set_97.fa.gz",
}

SPECIES_ASSIGN_DICT = {
    "SILVA": "s3://latch-public/test-data/4318/silva_species_assignment_v138.1.fa.gz",
    "RDP": "s3://latch-public/test-data/4318/rdp_species_assignment_18.fa.gz",
}


class TaxonomyReference(Enum):
    SILVA = "SILVA"
    RDP = "RDP"
    UNITE = "UNITE"
    GREENGENES = "GREENGENES"


class SpeciesAssignmentReference(Enum):
    SILVA = "SILVA"
    RDP = "RDP"


@dataclass
class Sample:
    read1: LatchFile
    read2: LatchFile
