"""Testing strategy for HTSinfer."""

import os
import time
import subprocess as sp
import datetime as dt
from pathlib import Path
from typing import List
import pandas as pd

# Generate run_id and set params for run
JOB_ID = "remove_falbicollis"

RUN_ID = str(dt.datetime.now().strftime("%m%d_%H%M%S_") + f"{JOB_ID}")
RECORDS = 1000000
READ_MIN_MATCH = 0.1
READ_MIN_FREQ = 2
LIB_MIN_MATCH = 2
LIB_MIN_FREQ = 2


# test files
ZARP_DIR = (Path(__file__).resolve().parents[3] /
            "zarp")
RESULTS_SRA_DIR = (Path(__file__).resolve().parent /
                   "results_sra_downloads")
RESULTS_HTS_DIR = (Path(__file__).resolve().parent /
                   "results_htsinfer")
MINED_DATA = (Path(__file__).resolve().parent /
              "mined_test_data_all.tsv")


RUN_DIR = "/".join([str(RESULTS_HTS_DIR), RUN_ID])
os.makedirs(RUN_DIR, exist_ok=True)

result_data = {}
graph_data = {}


def compile_sra_command() -> List[str]:
    """Compile command for SRA download workflow."""
    cmd_sra = ["conda", "run", "-n", "zarp", "snakemake"]
    cmd_sra.extend(["--snakefile",
                    "".join([str(ZARP_DIR),
                            "/workflow/rules/sra_download.smk"])])
    cmd_sra.extend(["--profile",
                    "".join([str(ZARP_DIR),
                             "/profiles/local-conda"])])
    cmd_sra.extend(["--config", "".join(["samples=", str(MINED_DATA)]),
                    "".join(["outdir=", str(RESULTS_SRA_DIR)]),
                    "".join(["samples_out=", str(RESULTS_SRA_DIR),
                             "/mined_data.out.tsv"]),
                    "".join(["log_dir=", "".join([str(ZARP_DIR), "/logs"])]),
                    "".join(["cluster_log_dir=", str(ZARP_DIR),
                             "/logs/cluster_log"])])
    return cmd_sra


class TestDownloadSRASamples:
    """Download the SRA Samples."""

    def test_zarp_download(self):
        """Download SRA samples."""
        cmd_sra = compile_sra_command()
        sp.Popen(cmd_sra).communicate()
