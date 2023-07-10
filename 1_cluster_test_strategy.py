"""Testing strategy for HTSinfer."""

import os
import subprocess as sp
import pandas as pd
import datetime as dt
from pathlib import Path
from typing import List

# Generate run_id and set params for run
job_id = "rec_100mil"

run_id = str(dt.datetime.now().strftime("%m%d_%H%M%S_") + f"{job_id}")
records = 100000000
read_min_match = 0.1
read_min_freq = 2
lib_min_match = 2
lib_min_freq = 2


# test files
ZARP_DIR = (Path(__file__).resolve().parents[3] /
            "zarp")
RESULTS_SRA_DIR = (Path(__file__).resolve().parents[1] /
                   "htsinfer/tests/cluster_tests/results_sra_downloads")
RESULTS_HTS_DIR = (Path(__file__).resolve().parent /
                   "results_htsinfer")
MINED_DATA = (Path(__file__).resolve().parent /
              "mined_test_data.tsv")


RUN_DIR = "/".join([str(RESULTS_HTS_DIR), run_id])
os.makedirs(RUN_DIR, exist_ok=True)

result_data = {}
graph_data = {}


def compile_sra_command() -> List[str]:
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


def compile_hts_command() -> List[str]:
    cmd_hts = ["htsinfer"]
    cmd_hts.extend(["--output-directory", RUN_DIR])
    cmd_hts.extend(["--temporary-directory", "$TMPDIR"])
    cmd_hts.extend(["--cleanup-regime", "KEEP_ALL"])
    cmd_hts.extend(["--verbosity", "DEBUG"])
    cmd_hts.extend(["--read-layout-min-match-percentage",
                    str(read_min_match)])
    cmd_hts.extend(["--read-layout-min-frequency-ratio",
                    str(read_min_freq)])
    cmd_hts.extend(["--library-source-min-match-percentage",
                    str(lib_min_match)])
    cmd_hts.extend(["--library-source-min-frequency-ratio",
                    str(lib_min_freq)])
    cmd_hts.extend(["--records", str(records)])
    return cmd_hts


class TestDownloadSRASamples:
    """Download the SRA Samples."""
    def test_zarp_download(self):
        cmd_sra = compile_sra_command()
        sp.Popen(cmd_sra).communicate()


class TestHTSinfer:
    """Test HTSinfer on downloaded samples."""

    def test_htsinfer_se_pe(self):
        # Read in sample CSV file using pandas
        source = pd.read_csv(MINED_DATA, sep='\t')
        # Group data into chunks of 15 samples
        chunks = [source[i:i+5] for i in range(1, len(source) + 1, 5)]

        # Create run directory to store job scripts
        JOBS_DIR = "/".join([str(RUN_DIR), '_job_scripts'])
        os.makedirs(JOBS_DIR, exist_ok=True)
        RUN_RESULT_DIR = "/".join([str(RUN_DIR), '_results'])
        os.makedirs(RUN_RESULT_DIR, exist_ok=True)

        # Iterate over each chunk and create a Slurm job script
        for i, chunk in enumerate(chunks):
            # Create a job script filename
            job_script_filename = f'{JOBS_DIR}/{run_id}_job_{i}.sh'

            # Write Slurm job script
            with open(job_script_filename, 'w') as f:
                f.write('#!/bin/bash\n\n')
                f.write(f'#SBATCH --job-name=htsinfer_{job_id}_{i}\n')
                f.write('#SBATCH --cpus-per-task=8\n')
                f.write('#SBATCH --mem-per-cpu=4G\n')
                f.write('#SBATCH --time=1-00:00:00\n')
                f.write('#SBATCH --qos=1day\n')
                f.write('#SBATCH --output=/dev/null\n')
                f.write('#SBATCH --error=/dev/null\n\n')
                f.write('source /scicore/home/zavolan/${USER}/.bashrc\n')
                f.write('conda activate htsinfer\n')

                # Iterate over each row in the chunk
                for _, row in chunk.iterrows():
                    # Determine the command based on the layout
                    if row['layout'] == 'SE':
                        if 'fq1' in row and not pd.isnull(row['fq1']):
                            sample = row['sample']
                            sample_se = row['fq1']
                        else:
                            sample = row['sample']
                            sample_se = "/".join([
                                str(RESULTS_SRA_DIR),
                                sample, sample
                                ]) + '.fastq.gz'
                        cmd_hts = compile_hts_command()
                        cmd_hts.extend([sample_se])
                    elif row['layout'] == 'PE':
                        if 'fq1' in row and not pd.isnull(row['fq1']):
                            sample = row['sample']
                            sample_1 = row['fq1']
                            sample_2 = row['fq2']
                        else:
                            sample = row['sample']
                            sample_1 = "/".join([
                                str(RESULTS_SRA_DIR),
                                sample, sample
                                ]) + '_1.fastq.gz'
                            sample_2 = "/".join([
                                str(RESULTS_SRA_DIR),
                                sample, sample
                                ]) + '_2.fastq.gz'
                        cmd_hts = compile_hts_command()
                        cmd_hts.extend([sample_1, sample_2])
                    else:
                        continue

                    cmd_hts.extend([
                        f'> {RUN_RESULT_DIR}/{sample}_result.json'
                        ])
                    cmd_hts.extend([
                        f'2> {RUN_RESULT_DIR}/{sample}_error.txt'
                        ])
                    cmd_hts = " ".join(cmd_hts)

                    f.write(f'{cmd_hts}\n')

            # Submit the job script using subprocess
            sp.run(['sbatch', job_script_filename])
