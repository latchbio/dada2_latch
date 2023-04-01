import subprocess
from pathlib import Path

from latch import medium_task, workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import LatchDir, LatchFile, file_glob
from typing import List

from .docs import wf_docs
from .types import Sample


@medium_task
def run_dada2(sample_name: str) -> List[LatchFile]:
    """Task to run a software"""

    sample_name = sample_name
    results_path = "dada2_results"
    output_dir = Path(results_path).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    _run_cmd = ["Rscript", "/root/dada2.R", sample_name]

    subprocess.run(_run_cmd)

    return file_glob(f"{sample_name}*tsv", f"latch:///{results_path}")


@workflow(wf_docs)
def dada2(sample_name: str) -> List[LatchFile]:
    """Workflow to do X

    Header
    ------

    This is a workflow that does X.

    """
    return run_dada2(sample_name=sample_name)


LaunchPlan(
    dada2,
    "Test Data",
    {
        "sample_name": "SRR579292",
    },
)
