# generate_slurm_script.py
import os
import subprocess as sp
from pathlib import Path
from typing import List
import pandas as pd

def generate_slurm_script(run_id, sample_table_path, blocks=1):

    # Create §htsinfer directory in the parent directory
    results_dir = os.path.abspath(os.path.join(Path.cwd(), "results_htsinfer"))
    os.makedirs(results_dir, exist_ok=True)

    # Create a directory within results_htsinfer with the given run_id
    run_dir = os.path.join(results_dir, run_id)
    os.makedirs(run_dir, exist_ok=True)

    # Create a slurm_scripts directory within the run_id directory
    slurm_scripts_dir = os.path.join(run_dir, "slurm_scripts")
    os.makedirs(slurm_scripts_dir, exist_ok=True)

    with open(sample_table_path, "r", encoding="utf-8") as csv_file:
        source = pd.read_csv(csv_file, sep=',')
        # Group data into chunks of 15 samples
        chunks = [source[i:i+blocks] for i in range(0, len(source), blocks)]

        for i, chunk in enumerate(chunks):
            # Create a job script filename
            job_script_filename = os.path.join(slurm_scripts_dir, f'{run_id}_job_{i}.sh')

            # Write Slurm job script
            with open(job_script_filename, 'w', encoding="utf-8") as file:
                file.write('#!/bin/bash\n\n')
                file.write(f'#SBATCH --job-name={run_id}_{i}\n')
                file.write('#SBATCH --cpus-per-task=4\n')
                file.write('#SBATCH --mem-per-cpu=4G\n')
                file.write('#SBATCH --time=1-00:00:00\n')
                file.write('#SBATCH --qos=1day\n')
                file.write('#SBATCH --output=/dev/null\n')
                file.write('#SBATCH --error=/dev/null\n\n')
                file.write('source /scicore/home/zavolan/${USER}/.bashrc\n')
                file.write('conda activate htsinfer\n')

                # Iterate over each row in the chunk
                for _, row in chunk.iterrows():
                    sample = row['sample']
                    if row['fq2'] is None:
                        sample_path = row['fq1']
                        cmd_hts = compile_hts_command(outdir=run_dir)
                        cmd_hts.extend([str(sample_path)])
                    else:
                        sample_path_1 = row['fq1']
                        sample_path_2 = row['fq2']
                        cmd_hts = compile_hts_command(outdir=run_dir)
                        cmd_hts.extend([str(sample_path_1), str(sample_path_2)])

                    cmd_hts.extend([
                        f'> {run_dir}/{sample}_result.json'
                        ])
                    cmd_hts.extend([
                        f'2> {run_dir}/{sample}_error.txt'
                        ])
                    cmd_hts = " ".join(cmd_hts)

                    file.write(f'{cmd_hts}\n')


def compile_hts_command(
        outdir,
        read_layout_min_match=0.1,
        read_layout_min_freq=2,
        lib_source_min_match=5,
        lib_source_min_freq=2,
        lib_type_max=1000,
        lib_type_cutoff=0.85,
        read_orient_min=20,
        read_orient_freq=0.75,
        verbosity="DEBUG",
        cleanup="KEEP_NONE",
        records=1000000,
        threads=8,
) -> List[str]:
    """Compile command for HTSinfer."""
    cmd_hts = ["htsinfer"]
    cmd_hts.extend(["--output-directory", outdir])
    cmd_hts.extend(["--temporary-directory", "$TMPDIR"])
    cmd_hts.extend(["--threads", str(threads)])
    cmd_hts.extend(["--cleanup-regime", cleanup])
    cmd_hts.extend(["--verbosity", verbosity])
    cmd_hts.extend(["--library-type-mates-cutoff", str(lib_type_cutoff)])
    cmd_hts.extend(["--read-layout-min-match-percentage", str(read_layout_min_match)])
    cmd_hts.extend(["--read-layout-min-frequency-ratio", str(read_layout_min_freq)])
    cmd_hts.extend(["--library-source-min-match-percentage", str(lib_source_min_match)])
    cmd_hts.extend(["--library-source-min-frequency-ratio", str(lib_source_min_freq)])
    cmd_hts.extend(["--records", str(records)])
    cmd_hts.extend(["--library-type-max-distance", str(lib_type_max)])
    cmd_hts.extend(["--library-type-mates-cutoff", str(lib_type_cutoff)])
    cmd_hts.extend(["--read-orientation-min-mapped-reads", str(read_orient_min)])
    cmd_hts.extend(["--read-orientation-min-fraction", str(read_orient_freq)])
    return cmd_hts





generate_slurm_script("test_id", "test_mined_data.csv", 20)

