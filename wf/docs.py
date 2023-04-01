from latch.types.metadata import (
    LatchAuthor,
    LatchMetadata,
    LatchParameter,
    Params,
    Section,
    Text,
    Spoiler,
)

PARAMS = {
    "samples": LatchParameter(
        display_name="Samples",
    ),
    "taxonomy_ref_fasta": LatchParameter(display_name="Taxonomy Reference FASTA"),
    "species_assignment_fasta": LatchParameter(display_name="Species Assignment FASTA"),
    "minLen": LatchParameter(
        display_name="minLen",
        description="Remove reads with length less than this value",
    ),
    "maxN": LatchParameter(
        display_name="maxN",
        description="After truncation, sequences with more than this number of Ns will be discarded",
    ),
    "minQ": LatchParameter(
        display_name="minQ",
        description="Reads that contain a quality score less than this will be discarded",
    ),
    "maxEE": LatchParameter(
        display_name="maxEE",
        description="After truncation, reads with higher expected errors than this will be discarded.",
    ),
    "truncQ": LatchParameter(
        display_name="truncQ",
        description="Truncate reads at the first instance of a quality score less than or equal to this",
    ),
    "trimLeft": LatchParameter(
        display_name="trimLeft",
        description="The number of nucleotides to remove from the start of each read.",
    ),
    "trimRight": LatchParameter(
        display_name="trimRight",
        description="The number of nucleotides to remove from the end of each read.",
    ),
}

FLOW = [
    Section(
        "Samples",
        Text(
            "Samples provided have to include"
            " two files corresponding to the reads (paired-end),"
            " which follow the pattern '*1.fastq' and '*2.fastq'."
        ),
        Params("samples"),
    ),
    Section(
        "Reference Files",
        Text(
            "The taxonomy reference FASTA should be a training set of reference sequences with known taxonomy"
            " which will be used by the `assignTaxonomy` function. The Species Assignment FASTA will be used"
            " as input to the `addSpecies` function. See further explanation on the"
            " [DADA2 documentation](http://benjjneb.github.io/dada2/tutorial.html#assign-taxonomy)."
        ),
        Params("taxonomy_ref_fasta", "species_assignment_fasta"),
    ),
    Spoiler(
        "filterAndTrim parameters",
        Text(
            "Parameters provided to the `filterAndTrim` function"
            " in the DADA2 workflow"
        ),
        Params(
            "minLen",
            "maxN",
            "minQ",
            "maxEE",
            "truncQ",
            "trimLeft",
            "trimRight",
        ),
    ),
]

WORKFLOW_NAME = "DADA2"

wf_docs = LatchMetadata(
    display_name=WORKFLOW_NAME,
    documentation=f"https://github.com/jvfe/{WORKFLOW_NAME}_latch/blob/main/README.md",
    author=LatchAuthor(
        name="jvfe",
        github="https://github.com/jvfe",
    ),
    repository=f"https://github.com/jvfe/{WORKFLOW_NAME}_latch",
    license="MIT",
    parameters=PARAMS,
    tags=["NGS", "16S", "Amplicon", "ASV"],
    flow=FLOW,
)
