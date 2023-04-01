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
        batch_table_column=True,
    ),
    "taxonomy_ref_fasta": LatchParameter(display_name="Taxonomy Reference FASTA"),
    "species_assignment_fasta": LatchParameter(
        display_name="Species Assignment Reference FASTA"
    ),
}

FLOW = [
    Section(
        "Samples",
        Text(
            "Sample provided has to include an identifier for the sample (Sample name)"
            " and two files corresponding to the reads (paired-end)"
        ),
        Params("sample"),
    )
]

WORKFLOW_NAME = "dada2"

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
    tags=["NGS"],
    flow=FLOW,
)
