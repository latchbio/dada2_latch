from latch.types.metadata import (
    Fork,
    ForkBranch,
    LatchAuthor,
    LatchMetadata,
    LatchParameter,
    Params,
    Section,
    Text,
    Spoiler,
    LatchRule,
)

PARAMS = {
    "samples": LatchParameter(display_name="Samples", samplesheet=True),
    "tax_ref_fork": LatchParameter(),
    "species_assign_fork": LatchParameter(),
    "taxonomy_reference": LatchParameter(display_name="Taxonomy Reference Database"),
    "species_assignment": LatchParameter(display_name="Species Assignment Database"),
    "taxonomy_ref_fasta": LatchParameter(
        display_name="Taxonomy Reference FASTA",
        detail="(.fa, .fasta, .fa.gz, .fasta.gz)",
        rules=[
            LatchRule(
                regex="(.fa|.fasta|.fa.gz|.fasta.gz)$",
                message="Must be a valid FASTA file",
            )
        ],
    ),
    "species_assignment_fasta": LatchParameter(
        display_name="Species Assignment FASTA",
        detail="(.fa, .fasta, .fa.gz, .fasta.gz)",
        rules=[
            LatchRule(
                regex="(.fa|.fasta|.fa.gz|.fasta.gz)$",
                message="Must be a valid FASTA file",
            )
        ],
    ),
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
            " two files corresponding to the reads (paired-end)."
        ),
        Params("samples"),
    ),
    Section(
        "Reference Files",
        Text(
            "The taxonomy reference should be a training set of reference sequences with known taxonomy"
            " which will be used by the `assignTaxonomy` function. The species assignment reference will be used"
            " as input to the `addSpecies` function. See further explanation on the"
            " [DADA2 documentation](http://benjjneb.github.io/dada2/tutorial.html#assign-taxonomy)."
        ),
        Text(
            "We already provide [SILVA v138.1](https://doi.org/10.5281/zenodo.4587954),"
            " [RDP 11.5](https://zenodo.org/record/4310151) and [UNITE 9.0](https://unite.ut.ee/repository.php)."
        ),
        Fork(
            "tax_ref_fork",
            "Choose a reference database",
            taxonomy_reference=ForkBranch("Database", Params("taxonomy_reference")),
            taxonomy_ref_fasta=ForkBranch(
                "FASTA Reference", Params("taxonomy_ref_fasta")
            ),
        ),
        Fork(
            "species_assign_fork",
            "Choose a species assignment database",
            species_assignment=ForkBranch("Database", Params("species_assignment")),
            species_assignment_fasta=ForkBranch(
                "FASTA Reference", Params("species_assignment_fasta")
            ),
        ),
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
    tags=["16S", "Amplicon", "ASV"],
    flow=FLOW,
)
