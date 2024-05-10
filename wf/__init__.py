import csv
import glob
import math
import os
import random
import re
import shutil
import subprocess
import warnings
from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import List, Optional
from urllib.parse import urlparse

import numpy as np
import pandas as pd
from dataclasses_json import dataclass_json
from ete3 import NCBITaxa
from latch import medium_task, workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import LatchDir, LatchFile

from wf.docs import wf_docs
from wf.types import (
    SPECIES_ASSIGN_DICT,
    TAXONOMY_REF_DICT,
    SpeciesAssignmentReference,
    TaxonomyReference,
)


@dataclass
class Sample:
    read1: LatchFile
    read2: LatchFile


class PoolingSetting(Enum):
    pseudo = "pseudo"
    TRUE = "TRUE"
    FALSE = "FALSE"


def convert_to_csv(file_path):
    # Read the text file and extract the data
    data = []
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if line:
                key, value = line.split(":")
                data.append([key.strip(), value.strip()])

    # Replace the text file with a CSV file
    csv_file_path = file_path.replace(".txt", ".csv")
    with open(csv_file_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Category", "Value"])  # Write the header
        writer.writerows(data)  # Write the data rows


def aws_cp(src: str, dst: str, show_progress: bool = False):
    cmd = ["aws", "s3", "cp", src, dst]
    if not show_progress:
        cmd.append("--no-progress")
    subprocess.run(cmd, check=True)


def interpret_range_1(row):
    if row["% hit"] < row["Range_low (%)"]:
        return "low"
    elif row["% hit"] > row["Range_high (%)"]:
        return "high"
    else:
        return "optimal"


# Function to find percentile based on sample name
def get_percentile_by_sample_name(df, sample_name):
    result = df[df["Sample"] == sample_name]["Percentile"]
    if not result.empty:
        return result.iloc[0]
    else:
        return "Sample not found"


def interpret_score_1(row):
    if row["Interpretation"] == "low":
        return row["% hit"] / row["Range_low (%)"] * 70
    elif row["Interpretation"] == "high":
        return row["Range_high (%)"] / row["% hit"] * 70
    else:
        return 70


def parse_file(file_path):
    scores = {}

    # Open the file and read the lines
    with open(file_path, "r") as file:
        for line in file:
            key, value = line.split(":")
            scores[key.strip()] = float(value.strip())

    return scores


def interpret_range_2(row):
    if row["% hit"] > 0.1:
        return "high"
    elif (row["% hit"] < 0.1) & (row["% hit"] > 0.0):
        return "low"
    else:
        return "absent"


def interpret_score_2(row):
    if (row["Interpretation"] == "low") | (row["Interpretation"] == "absent"):
        return 0
    else:
        return -2


# def get_tax_name(tax_name, ncbi):
#     try:
#         res = ncbi.get_name_translator([tax_name])[tax_name]
#         if(len(res) > 1):
#             return int(res[0])
#         else:
#             return int(res[0])
#     except:
#         return np.nan


def replace_with_random_nucleotide(sequence):
    def replace_match(match):
        return random.choice(["A", "C", "G", "T"])

    return re.sub(r"[^ACGT\n]", replace_match, sequence)


def clean_sequences(input_file, output_file):
    with gzip.open(input_file, "rt") as f_in, gzip.open(output_file, "wt") as f_out:
        for line in f_in:
            if line.startswith("@") or line.startswith("+") or line.strip() == "":
                # Write identifier, plus line, or empty line as is
                f_out.write(line)
            else:
                cleaned_line = replace_with_random_nucleotide(line)
                f_out.write(cleaned_line)


def get_tax_name(tax_name, tax_rank, ncbi):
    try:
        if (tax_name == "Prevotella copri") or (tax_name == "Segatella copri"):
            return 165179
        res = ncbi.get_name_translator([tax_name])[tax_name]
        if len(res) > 1:
            ids = ncbi.get_name_translator([tax_name])[tax_name]
            info = ncbi.get_rank(ids)
            for key, value in info.items():
                if value == ((tax_rank).lower()):
                    return key
        else:
            return int(res[0])
    except:
        return np.nan


def get_parent(tax_id, tax_rank, ncbi):
    try:
        lineage = ncbi.get_lineage(tax_id)
        annotated_lineage = ncbi.get_rank(lineage)

        if tax_rank == "Kingdom":
            return lineage[-2]
        elif tax_rank == "Phylum":
            return_val = np.nan
            for vals in annotated_lineage:
                if annotated_lineage[vals] == "superkingdom":
                    return_val = vals
            return return_val
        elif tax_rank == "Class":
            return_val = np.nan
            for vals in annotated_lineage:
                if annotated_lineage[vals] == "phylum":
                    return_val = vals
            return return_val
        elif tax_rank == "Order":
            return_val = np.nan
            for vals in annotated_lineage:
                if annotated_lineage[vals] == "class":
                    return_val = vals
            return return_val
        elif tax_rank == "Family":
            return_val = np.nan
            for vals in annotated_lineage:
                if annotated_lineage[vals] == "order":
                    return_val = vals
            return return_val
        elif tax_rank == "Genus":
            return_val = np.nan
            for vals in annotated_lineage:
                if annotated_lineage[vals] == "family":
                    return_val = vals
            return return_val
        elif tax_rank == "Species":
            return_val = np.nan
            for vals in annotated_lineage:
                if annotated_lineage[vals] == "genus":
                    return_val = vals
            return return_val
        else:
            return np.nan
    except:
        return np.nan


# Function to assign score based on percentile
def assign_shannon_score(percentile):
    if 0 <= percentile <= 50:
        return 10
    elif 50 < percentile <= 80:
        return 20
    elif 80 < percentile <= 95:
        return 25
    elif 95 < percentile <= 100:
        return 30
    else:
        print("ERROR, undefined percentile. ")
        print(percentile)
        return "Undefined percentile range"  # Handle unexpected case


def shannon_val(row, denominator, blanks):
    if blanks:
        p = row["num_hits"] / denominator
        return p * math.log(p)
    else:
        if len(row["tax_name"]) > 2:
            p = row["num_hits"] / denominator
            return p * math.log(p)
        else:
            return 0


def fetch_species_with_taxid(family_name, ncbi):
    family_id = ncbi.get_name_translator([family_name])[family_name][0]
    descendant_ids = ncbi.get_descendant_taxa(family_id, collapse_subspecies=True)
    species_ids = [
        taxid for taxid in descendant_ids if ncbi.get_rank([taxid])[taxid] == "species"
    ]
    species_names = ncbi.get_taxid_translator(species_ids)
    return [(taxid, species_names[taxid]) for taxid in species_ids]


def fill_missing_ranks(row, ncbi):
    ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    for i in range(len(ranks) - 1):
        if (pd.isna(row[ranks[i]]) or (row[ranks[i]] == "")) and not (
            pd.isna(row[ranks[i + 1]]) or (row[ranks[i + 1]] == "")
        ):
            if ranks[i + 1] == "Species":
                # Search for species using substring
                species_list = fetch_species_with_taxid(row["Family"], ncbi)
                substring = (
                    row["Species"].split()[1]
                    if " " in row["Species"]
                    else row["Species"]
                )
                matching_species = [
                    name for taxid, name in species_list if substring in name
                ]
                if matching_species:
                    row[ranks[i]] = matching_species[0].split()[0]  # Fill in the Genus
            else:
                # Use NCBI taxonomy to fill in the missing rank
                try:
                    name_to_id = ncbi.get_name_translator([row[ranks[i + 1]]])
                    lineage = ncbi.get_lineage(name_to_id[row[ranks[i + 1]]][0])
                    rank_dict = ncbi.get_rank(lineage)
                except:
                    break
                for taxid in lineage:
                    if rank_dict[taxid] == ranks[i].lower():
                        row[ranks[i]] = ncbi.get_taxid_translator([taxid])[taxid]
                        break
    return row


@medium_task
def run_dada2(
    samples: List[Sample],
    # input_dir: Optional[LatchDir],
    results_output_directory: LatchDir,
    taxonomy_reference: TaxonomyReference,
    species_assignment: Optional[SpeciesAssignmentReference],
    taxonomy_ref_fasta: Optional[LatchFile],
    species_assignment_fasta: Optional[LatchFile],
    # maxN: int,
    # minQ: int,
    # maxEE: int,
    # truncQ: int,
    # trimLeft: int,
    # trimRight: int,
    maxMismatch: int,
    minOverlap: int,
    derep: bool,
    removeChimeric: bool,
    omega_a: str,
    output_folder_name: str,
    truncLen: int = 150,
    mergeReads: bool = False,
    detect_singletons: bool = True,
    pooling_param: PoolingSetting = PoolingSetting.pseudo,
    # path_to_results_files: LatchDir = LatchDir("latch:///dada2_results/"),
    minLen: int = 20,
) -> LatchDir:
    """Task to run dada2"""

    # if (samples is None) and (input_dir is None):
    #     print("Error: A list of samples or an input directory must be included.")
    #     return

    if taxonomy_ref_fasta:
        taxonomy_loc = taxonomy_ref_fasta.local_path
    else:
        remote_path = TAXONOMY_REF_DICT[taxonomy_reference.value]
        local_path = Path.cwd() / Path(urlparse(remote_path).path).name
        taxonomy_loc = str(local_path.resolve())
        aws_cp(remote_path, taxonomy_loc)

    if species_assignment:
        remote_path = SPECIES_ASSIGN_DICT[species_assignment.value]
        local_path = Path.cwd() / Path(urlparse(remote_path).path).name
        species_assign_loc = str(local_path.resolve())
        aws_cp(remote_path, species_assign_loc)
    elif species_assignment_fasta:
        species_assign_loc = species_assignment_fasta.local_path
    else:
        species_assign_loc = None

    output_dir = f"/root/outputs/{output_folder_name}"
    output_dirpath = Path(output_dir).resolve()
    output_dirpath.mkdir(parents=True, exist_ok=True)
    print(str(output_dirpath))

    # Move local files:
    # if (samples is not None) and (input_dir is None):
    #     sample_input_dir = "/root/inputs"
    #     sample_input_dirpath = Path(sample_input_dir).resolve()
    #     sample_input_dirpath.mkdir(parents=True, exist_ok=True)

    #     for sample in samples:
    #         shutil.move(sample.read1.local_path, str(sample_input_dirpath))
    #         shutil.move(sample.read2.local_path, str(sample_input_dirpath))
    #     input_dir_final = sample_input_dirpath
    # else:
    #     input_dir_final = input_dir.local_path

    sample_input_dir = "/root/inputs"
    sample_input_dirpath = Path(sample_input_dir).resolve()
    sample_input_dirpath.mkdir(parents=True, exist_ok=True)

    for sample in samples:
        shutil.move(sample.read1.local_path, str(sample_input_dirpath))
        shutil.move(sample.read2.local_path, str(sample_input_dirpath))
    input_dir_final = sample_input_dirpath

    if derep:
        derep_setting = "derep"
    else:
        derep_setting = "no_derep"

    if removeChimeric:
        chimeric_setting = "remove"
    else:
        chimeric_setting = "dont_remove"

    if detect_singletons:
        singletons_setting = "detect_singletons"
    else:
        singletons_setting = "dont_detect"

    if mergeReads:
        paired_end_setting = "merge"
    else:
        paired_end_setting = "dont_merge"

    _run_cmd = [
        "Rscript",
        "/root/dada2.R",
        str(input_dir_final),
        str(output_dirpath),
        taxonomy_loc,
        str(truncLen),
        str(omega_a),
        str(singletons_setting),
        str(derep_setting),
        str(chimeric_setting),
        str(pooling_param.value),
        str(paired_end_setting),
        str(maxMismatch),
        str(minOverlap),
        str(minLen),
    ]

    print(_run_cmd)
    if species_assign_loc:
        _run_cmd.append(species_assign_loc)

    subprocess.run(_run_cmd)

    # Get highest Shannon Score found:

    if (species_assignment or species_assignment_fasta) or (
        taxonomy_reference.value == "GREENGENES"
    ):
        ncbi = NCBITaxa()
        ncbi.update_taxonomy_database()

        stat_sample = []

        paths = [x[0] for x in os.walk(output_dirpath)]
        for index, curr_path in enumerate(paths):
            print(index, curr_path)
            if index != 0:
                stat_sample.append(curr_path)

                asv_table = pd.read_csv(f"{curr_path}/asv_table.csv")
                species_table = pd.read_csv(f"{curr_path}/species_table.csv")

                if taxonomy_reference.value == "GREENGENES" and (
                    not species_assignment
                ):
                    species_table = species_table.rename(
                        columns={"Unnamed: 0": "sequence"}
                    )
                    asv_table = asv_table.T
                    asv_table = asv_table.reset_index().rename(
                        columns={"index": "sequence"}
                    )
                    asv_table.drop(index=asv_table.index[0], axis=0, inplace=True)
                    asv_table["num_hits"] = asv_table.drop("sequence", axis=1).sum(
                        axis=1
                    )
                    asv_table = asv_table[["sequence", "num_hits"]]
                    asv_table.reset_index(inplace=True, drop=True)
                    merged_table = species_table.merge(
                        asv_table, on="sequence", how="left"
                    )
                    sum_val = merged_table["num_hits"].sum()
                    merged_table["percent_hits"] = (
                        merged_table["num_hits"] / sum_val * 100
                    )

                    if merged_table["Kingdom"].str.contains("__").any():
                        merged_table["Kingdom"] = merged_table["Kingdom"].str[3:]
                        merged_table["Phylum"] = merged_table["Phylum"].str[3:]
                        merged_table["Class"] = merged_table["Class"].str[3:]
                        merged_table["Order"] = merged_table["Order"].str[3:]
                        merged_table["Family"] = merged_table["Family"].str[3:]
                        merged_table["Genus"] = merged_table["Genus"].str[3:]
                        merged_table["Species"] = merged_table["Species"].str[3:]

                    merged_table = merged_table.replace(
                        to_replace=r"\[.*?\]", value="", regex=True
                    )

                    merged_path = f"{curr_path}/sequence_summary_table.csv"
                    merged_table.to_csv(merged_path, index=False)

                else:
                    if taxonomy_reference.value == "GREENGENES":
                        species_table = species_table.rename(
                            columns={"Unnamed: 0": "sequence"}
                        )
                        species_table = species_table.drop(["Species"], axis=1)
                        species_table = species_table.reset_index().rename(
                            columns={"Species.1": "Species"}
                        )
                    else:
                        species_table = species_table.rename(
                            columns={"Unnamed: 0": "sequence"}
                        )

                    # Read asv table
                    asv_table = asv_table.T
                    asv_table = asv_table.reset_index().rename(
                        columns={"index": "sequence"}
                    )
                    asv_table.drop(index=asv_table.index[0], axis=0, inplace=True)
                    asv_table["num_hits"] = asv_table.drop("sequence", axis=1).sum(
                        axis=1
                    )
                    asv_table = asv_table[["sequence", "num_hits"]]
                    asv_table.reset_index(inplace=True, drop=True)

                    # Merge tables
                    merged_table = species_table.merge(
                        asv_table, on="sequence", how="left"
                    )
                    sum_val = merged_table["num_hits"].sum()
                    merged_table["percent_hits"] = (
                        merged_table["num_hits"] / sum_val * 100
                    )

                    if merged_table["Kingdom"].str.contains("__").any():
                        merged_table["Kingdom"] = merged_table["Kingdom"].str[3:]
                        merged_table["Phylum"] = merged_table["Phylum"].str[3:]
                        merged_table["Class"] = merged_table["Class"].str[3:]
                        merged_table["Order"] = merged_table["Order"].str[3:]
                        merged_table["Family"] = merged_table["Family"].str[3:]
                        merged_table["Genus"] = merged_table["Genus"].str[3:]

                    merged_table = merged_table.replace(
                        to_replace=r"\[.*?\]", value="", regex=True
                    )

                    merged_path = f"{curr_path}/sequence_summary_table.csv"
                    merged_table.to_csv(merged_path, index=False)

                merged_table = merged_table.apply(fill_missing_ranks, axis=1, ncbi=ncbi)

                # Summary Table Converged
                sum_table = merged_table.copy()
                sum_table = sum_table[
                    [
                        "Kingdom",
                        "Phylum",
                        "Class",
                        "Order",
                        "Family",
                        "Genus",
                        "Species",
                        "num_hits",
                    ]
                ]
                sum_table = sum_table.groupby(
                    [
                        "Kingdom",
                        "Phylum",
                        "Class",
                        "Order",
                        "Family",
                        "Genus",
                        "Species",
                    ],
                    as_index=False,
                    dropna=False,
                ).agg({"num_hits": "sum"})
                sum_val = sum_table["num_hits"].sum()
                sum_table["percent_hits"] = sum_table["num_hits"] / sum_val * 100
                sum_table = sum_table.sort_values("percent_hits", ascending=False)
                sum_table_path = f"{curr_path}/summary_table.csv"
                sum_table.to_csv(sum_table_path, index=False)

                tax_table = merged_table.copy()

                # Create taxonomy table
                columns = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]

                tax_table.loc[
                    (tax_table["Species"].notnull())
                    & (tax_table["Species"] != "")
                    & (tax_table["Genus"].notnull())
                    & (tax_table["Genus"] != ""),
                    "Species",
                ] = (
                    tax_table["Genus"] + " " + tax_table["Species"]
                )

                tax_table_v1 = (
                    tax_table.groupby(["Kingdom"])["num_hits"].sum().reset_index()
                )
                tax_table_v1 = tax_table_v1.rename(columns={"Kingdom": "tax_name"})
                tax_table_v1 = tax_table_v1.sort_values("num_hits", ascending=False)
                tax_table_v1["tax_rank"] = "Kingdom"

                for col in columns:
                    tax_table_temp = (
                        tax_table.groupby([col])["num_hits"].sum().reset_index()
                    )
                    tax_table_temp = tax_table_temp.sort_values(
                        "num_hits", ascending=False
                    )
                    tax_table_temp = tax_table_temp.rename(columns={col: "tax_name"})
                    tax_table_temp["tax_rank"] = col
                    tax_table_v1 = pd.concat(
                        [tax_table_v1, tax_table_temp], ignore_index=True, axis=0
                    )

                tax_table_v1 = tax_table_v1[["tax_name", "tax_rank", "num_hits"]]

                # Get taxonomy IDs
                tax_table_v1["taxon"] = tax_table_v1.apply(
                    lambda x: get_tax_name(x["tax_name"], x["tax_rank"], ncbi), axis=1
                )

                # Get parent IDs
                tax_table_v1["parent"] = tax_table_v1.apply(
                    lambda x: get_parent(x["taxon"], x["tax_rank"], ncbi), axis=1
                )

                sum_val = tax_table_v1.loc[
                    tax_table_v1["tax_rank"] == "Kingdom", "num_hits"
                ].sum()
                tax_table_v1["percent_hits"] = tax_table_v1["num_hits"] / sum_val * 100

                tax_table_v1["taxon"] = tax_table_v1["taxon"].astype("Int64")
                tax_table_v1["parent"] = tax_table_v1["parent"].astype("Int64")

                # tax_table_v1["percentile_rank"] = tax_table_v1["num_hits"].rank(
                #     pct=True
                # )

                tax_table_v2 = tax_table_v1.copy()
                tax_table_v2["Sample"] = curr_path

                taxonomy_path = f"{curr_path}/taxonomy_table.csv"
                tax_table_v2.to_csv(taxonomy_path, index=False)

    # return LatchDir(str(output_dirpath), f"latch:///dada2_results/{output_dir}")
    return LatchDir("/root/outputs/", results_output_directory.remote_path)


@workflow(wf_docs)
def dada2(
    samples: List[Sample],
    results_output_directory: LatchDir,
    tax_ref_fork: str,
    output_folder_name: str,
    taxonomy_reference: TaxonomyReference,
    species_assign_fork: str = "none",
    species_assignment: Optional[SpeciesAssignmentReference] = None,
    taxonomy_ref_fasta: Optional[LatchFile] = None,
    species_assignment_fasta: Optional[LatchFile] = None,
    minLen: int = 20,
    maxMismatch: int = 5,
    minOverlap: int = 10,
    derep: bool = False,
    removeChimeric: bool = False,
    pooling_param: PoolingSetting = PoolingSetting.pseudo,
    omega_a: str = "1e-30",
    truncLen: int = 150,
    mergeReads: bool = False,
    detect_singletons: bool = True,
    # input_dir: Optional[LatchDir] = None,
    # path_to_results_files: LatchDir = LatchDir("latch:///dada2_results/"),
) -> LatchDir:
    """A workflow for fast and accurate sample inference from amplicon data

    DADA2
    ------

    This is a workflow to run the DADA2[^1] (v1.26.0) software for amplicon data.
    It currently runs the following DADA2 functions under the hood:

    - `filterAndTrim` for filtering and trimming sequence files
    - `derepFastq` for dereplication
    - `learnErrors` to learn the error rates from the data
    - `dada` the main DADA2 function, to infer sample composition
    - `mergePairs` to merge paired-end reads
    - `makeSequenceTable` to build the per-sample sequence table
    - `removeBimeraDenovo` to remove sequencing chimeras
    - `assignTaxonomy` and `addSpecies` for taxonomic classification of ASVs

    # Workflow Output Files

    The workflow outputs two data files, one corresponding to the per-sequence table after chimera removal
    (`asv_table.csv`) and one for the taxonomic classification (`species_table.csv`). More files are also generated in addition to these if GREENGENES is used or if
    a species assignment database is selected.

    *For every run (default output of DADA2):*

    - species_table.csv: Taxonomy assigned to each ASV sequence and the number of times it was observed in the sample
    - asv_table.csv: Number of times each ASV was observed in the sample

    *If species assignment database is selected:*

    - summary_table.csv: Number of hits and percent hits for each taxa (overall summary without sequence information)
    - sequence_summary_table.csv: Taxonomy assigned to each ASV sequence and the number of times it was observed in the sample. Percent of total observations is included for each sequence.
    - taxonomy_table.csv: Number of observations for every taxonomy level (Kingdom to Species) with percentages. NCBI Taxonomy database was used to map taxa and parent IDs.

    # Note on Pooling
    To learn more about the settings for the pooling parameter, you can check out [this DADA2 resource](https://benjjneb.github.io/dada2/pseudo.html#pseudo-pooling).
    This is a very important parameter that increasing the sensitivity of DADA2, especially across pooled samples.

    See the [DADA2 tutorial](http://benjjneb.github.io/dada2/tutorial.html) for further explanation
    on the general workflow.

    [^1]: Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016).
    “DADA2: High-resolution sample inference from Illumina amplicon data.”
    Nature Methods, 13, 581-583. https://doi.org/10.1038/nmeth.3869.
    """

    return run_dada2(
        samples=samples,
        # input_dir=input_dir,
        taxonomy_reference=taxonomy_reference,
        species_assignment=species_assignment,
        taxonomy_ref_fasta=taxonomy_ref_fasta,
        species_assignment_fasta=species_assignment_fasta,
        minLen=minLen,
        # maxN=maxN,
        # minQ=minQ,
        # maxEE=maxEE,
        # truncQ=truncQ,
        # trimLeft=trimLeft,
        # trimRight=trimRight,
        maxMismatch=maxMismatch,
        detect_singletons=detect_singletons,
        mergeReads=mergeReads,
        minOverlap=minOverlap,
        derep=derep,
        removeChimeric=removeChimeric,
        pooling_param=pooling_param,
        omega_a=omega_a,
        truncLen=truncLen,
        output_folder_name=output_folder_name,
        # path_to_results_files=path_to_results_files,
        results_output_directory=results_output_directory,
    )


LaunchPlan(
    dada2,
    "Fungus metagenome",
    {
        "samples": [
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
        "taxonomy_reference": TaxonomyReference.UNITE,
    },
)

LaunchPlan(
    dada2,
    "MiSeq SOP",
    {
        "samples": [
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D0_S188_L001_R1_001.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D0_S188_L001_R2_001.fastq"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D1_S189_L001_R1_001.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D1_S189_L001_R2_001.fastq"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D2_S190_L001_R1_001.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D2_S190_L001_R2_001.fastq"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D3_S191_L001_R1_001.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D3_S191_L001_R2_001.fastq"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D5_S193_L001_R1_001.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D5_S193_L001_R2_001.fastq"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D6_S194_L001_R1_001.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D6_S194_L001_R2_001.fastq"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D7_S195_L001_R1_001.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D7_S195_L001_R2_001.fastq"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D8_S196_L001_R1_001.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D8_S196_L001_R2_001.fastq"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D9_S197_L001_R1_001.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D9_S197_L001_R2_001.fastq"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D141_S207_L001_R1_001.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D141_S207_L001_R2_001.fastq"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D142_S208_L001_R1_001.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D142_S208_L001_R2_001.fastq"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D143_S209_L001_R1_001.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D143_S209_L001_R2_001.fastq"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D144_S210_L001_R1_001.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D144_S210_L001_R2_001.fastq"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D145_S211_L001_R1_001.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D145_S211_L001_R2_001.fastq"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D146_S212_L001_R1_001.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D146_S212_L001_R2_001.fastq"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D147_S213_L001_R1_001.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D147_S213_L001_R2_001.fastq"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D148_S214_L001_R1_001.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D148_S214_L001_R2_001.fastq"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D149_S215_L001_R1_001.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D149_S215_L001_R2_001.fastq"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D150_S216_L001_R1_001.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/F3D150_S216_L001_R2_001.fastq"
                ),
            ),
            Sample(
                read1=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/Mock_S280_L001_R1_001.fastq"
                ),
                read2=LatchFile(
                    "s3://latch-public/test-data/4318/MiSeq_SOP/Mock_S280_L001_R2_001.fastq"
                ),
            ),
        ],
        "taxonomy_reference": TaxonomyReference.SILVA,
        "species_assignment": SpeciesAssignmentReference.SILVA,
    },
)
