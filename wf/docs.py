from latch.types.metadata import (
    Fork,
    ForkBranch,
    LatchAuthor,
    LatchMetadata,
    LatchParameter,
    LatchRule,
    Params,
    Section,
    Spoiler,
    Text,
)

PARAMS = {
    # "input_dir": LatchParameter(
    #     display_name="Input Folder", description="Directory containing FASTQ files"
    # ),
    # "samples": LatchParameter(
    #     display_name="Samples",
    #     description="Input FASTQ files",
    #     # batch_table_column=True,
    #     samplesheet=True,
    # ),
    "samples": LatchParameter(
        display_name="Samples", samplesheet=True, batch_table_column=True
    ),
    "output_folder_name": LatchParameter(
        display_name="Name of Output Folder",
        description="Use only letters, underscores, or dashes (no spaces)",
    ),
    "pooling_param": LatchParameter(
        display_name="Pooling parameter",
        description="Pooling parameter setting for dada(...). Pseudo pooling is recommended for quicker modeling across all samples. Check the about page for more information on this parameter.",
    ),
    "omega_a": LatchParameter(
        display_name="OMEGA_A Sensitivity Parameter",
        description="Sets the level of “statistical evidence” (think p-value) required for inferences of a new ASV for dada(...)",
    ),
    # "just_forward": LatchParameter(display_name="Only Forward Reads", description="Only use forward reads for DADA2"),
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
    # "path_to_results_files": LatchParameter(
    #     display_name="Previous Results Folder",
    #     description="Used to find top Shannon Diversity score from past runs",
    # ),
    "minLen": LatchParameter(
        display_name="minLen",
        description="Remove reads with length less than this value",
    ),
    # "maxN": LatchParameter(
    #     display_name="maxN",
    #     description="After truncation, sequences with more than this number of Ns will be discarded",
    # ),
    # "minQ": LatchParameter(
    #     display_name="minQ",
    #     description="Reads that contain a quality score less than this will be discarded",
    # ),
    # "maxEE": LatchParameter(
    #     display_name="maxEE",
    #     description="After truncation, reads with higher expected errors than this will be discarded.",
    # ),
    # "truncQ": LatchParameter(
    #     display_name="truncQ",
    #     description="Truncate reads at the first instance of a quality score less than or equal to this",
    # ),
    "truncLen": LatchParameter(
        display_name="truncLen",
        description="Truncate reads to a this length",
    ),
    # "trimLeft": LatchParameter(
    #     display_name="trimLeft",
    #     description="The number of nucleotides to remove from the start of each read.",
    # ),
    # "trimRight": LatchParameter(
    #     display_name="trimRight",
    #     description="The number of nucleotides to remove from the end of each read.",
    # ),
    # "customFiltering": LatchParameter(
    #     display_name="Custom Filtering",
    #     description="Strip N characters from reads without conducting any other filtering/trimming. If this is True, filterAndTrim will not be conducted.",
    # ),
    "maxMismatch": LatchParameter(
        display_name="maxMismatch",
        description="The maxMismatch parameter in mergePairs.",
    ),
    "mergeReads": LatchParameter(
        display_name="Merge Reads",
        description="Use both forward and backward reads, and merge them. If false, only forward reads will be used.",
    ),
    "minOverlap": LatchParameter(
        display_name="minOverlap",
        description="The maxMismatch parameter in minOverlap.",
    ),
    "derep": LatchParameter(
        display_name="dereplication",
        description="Conduct dereplication",
    ),
    "detect_singletons": LatchParameter(
        display_name="DETECT_SINGLETONS",
        description="Detect singletons during sample composition inference",
    ),
    "removeChimeric": LatchParameter(
        display_name="removeChimeric",
        description="Remove chimeric sequences",
    ),
    "results_output_directory": LatchParameter(
        display_name="Output Directory",
        description="Select a directory to save results to",
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
        "Output Format",
        # Text(
        #     "Directory containing FASTQ files"
        # ),
        # Fork(
        #     "input_fork",
        #     "Input samples",
        #     # input_dir=ForkBranch("Input Folder", Params("input_dir")),
        #     samples=ForkBranch("Reads", Params("samples")),
        # ),
        Params(
            # "samples",
            "output_folder_name",
            "mergeReads",
            # "path_to_results_files",
            "results_output_directory",
        ),
    ),
    Section(
        "Reference Files",
        Text(
            " The taxonomy reference should be a training set of reference sequences with known taxonomy"
            " which will be used by the `assignTaxonomy` function. The species assignment reference will be used"
            " as input to the `addSpecies` function. See further explanation on the"
            " [DADA2 documentation](http://benjjneb.github.io/dada2/tutorial.html#assign-taxonomy)."
        ),
        Text(
            "We already provide [SILVA v138.1](https://doi.org/10.5281/zenodo.4587954),"
            " [RDP 11.5](https://zenodo.org/record/4310151), [UNITE 9.0](https://unite.ut.ee/repository.php)"
            "and [GREENGENES v13.8](https://zenodo.org/records/158955)."
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
            "Choose a species assignment database. If GREENGENES is selected for the taxonomy reference database, it is not necessary to select a species assignment. \
            GREENGENES already comes with a species assignment. If you would like to use a different species database other than GREENGENES for data that used GREENGENES \
            as a taxonomy database, this is still possible.",
            species_assignment=ForkBranch("Database", Params("species_assignment")),
            species_assignment_fasta=ForkBranch(
                "FASTA Reference", Params("species_assignment_fasta")
            ),
        ),
    ),
    Section(
        "Processing Parameters",
        Text("Parameters provided to the DADA2 workflow."),
        Params(
            # "customFiltering",
            "derep",
            "removeChimeric",
            # "just_forward",
        ),
    ),
    Section(
        "filterAndTrim Parameters",
        Text(
            "Parameters provided to the `filterAndTrim` function"
            " in the DADA2 workflow"
        ),
        Params(
            "minLen",
            # "maxN",
            # "minQ",
            # "maxEE",
            # "truncQ",
            # "trimLeft",
            # "trimRight",
            # "maxMismatch",
            # "minOverlap",
            "truncLen",
            # "maxLen",
        ),
    ),
    Section(
        "mergePairs Parameters",
        Text(
            "Parameters provided to the `mergePairs` function" " in the DADA2 workflow"
        ),
        Params(
            "maxMismatch",
            "minOverlap",
        ),
    ),
    Section(
        "dada() Parameters",
        Text(
            "Parameters provided to the `dada(...)` function" " in the DADA2 workflow"
        ),
        Params(
            "pooling_param",
            "omega_a",
            "detect_singletons",
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
