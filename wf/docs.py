from latch.types.metadata import (
    LatchAuthor,
    LatchMetadata,
    LatchParameter,
    Params,
    Section,
    Text,
)

PARAMS = {
    "samples": LatchParameter(
        display_name="Samples",
    ),
    "taxonomy_ref_fasta": LatchParameter(display_name="Taxonomy Reference FASTA"),
    "species_assignment_fasta": LatchParameter(display_name="Species Assignment FASTA"),
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
