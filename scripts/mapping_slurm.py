# mapping_slurm.py
import os
import time
import subprocess as sp
from pathlib import Path
from typing import List
import pandas as pd

def generate_mapping_slurm_script(run_id, ref_genomes_dir, sample_table_path, blocks=1):

    # Create §htsinfer directory in the parent directory
    results_dir = os.path.abspath(os.path.join(Path.cwd(), "results_mapping"))
    os.makedirs(results_dir, exist_ok=True)

    # Create a directory within results_htsinfer with the given run_id
    run_dir = os.path.join(results_dir, run_id)
    os.makedirs(run_dir, exist_ok=True)

    # Create a slurm_scripts directory within the run_id directory
    slurm_scripts_dir = os.path.join(run_dir, "slurm_scripts")
    os.makedirs(slurm_scripts_dir, exist_ok=True)

    with open(sample_table_path, "r", encoding="utf-8") as csv_file:
        source = pd.read_csv(csv_file, sep=',')
        # Group data into chunks of n samples
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
                file.write(f'#SBATCH --output={run_id}_{i}.o%j\n')
                file.write(f'#SBATCH --error={run_id}_{i}.e%j\n\n')
                file.write('source /scicore/home/zavolan/${USER}/.bashrc\n')
                file.write('conda activate htsinfer\n')

                # Iterate over each row in the chunk
                for _, row in chunk.iterrows():
                    sample = row['sample']
                    print(f"Current sample: {sample}")
                    sample_result_dir = os.path.join(run_dir, f"{sample}")
                    os.makedirs(sample_result_dir, exist_ok=True)

                    if pd.isna(row['fq2']):
                        sample_path = row['fq1']

                        org_id_ref_genome = os.path.join(ref_genomes_dir, str(row['org_id']))
                        filename = os.listdir(org_id_ref_genome)[0]
                        org_id_ref_file = os.path.join(org_id_ref_genome, filename)
                        org_id_ref_index = os.path.join(sample_result_dir, str(row['org_id']), "index")

                        org_1_id_ref_genome = os.path.join(ref_genomes_dir, str(row['1_org_1_id']))
                        filename = os.listdir(org_1_id_ref_genome)[0]
                        org_1_id_ref_file = os.path.join(org_1_id_ref_genome, filename)
                        org_1_id_ref_index = os.path.join(sample_result_dir, str(row['1_org_1_id']), "index")

                        org_2_id_ref_genome = os.path.join(ref_genomes_dir, str(row['1_org_2_id']))
                        filename = os.listdir(org_2_id_ref_genome)[0]
                        org_2_id_ref_file = os.path.join(org_2_id_ref_genome, filename)
                        org_2_id_ref_index = os.path.join(sample_result_dir, str(row['1_org_2_id']), "index")

                        org_3_id_ref_genome = os.path.join(ref_genomes_dir, str(row['1_org_3_id']))
                        filename = os.listdir(org_3_id_ref_genome)[0]
                        org_3_id_ref_file = os.path.join(org_3_id_ref_genome, filename)
                        org_3_id_ref_index = os.path.join(sample_result_dir, str(row['1_org_3_id']), "index")

                        cmd_index_org = compile_create_star_index(
                            genome_fasta=org_id_ref_file,
                            index_dir=org_id_ref_index
                        )
                        cmd_index_1 = compile_create_star_index(
                            genome_fasta=org_1_id_ref_file,
                            index_dir=org_1_id_ref_index
                        )
                        cmd_index_2 = compile_create_star_index(
                            genome_fasta=org_2_id_ref_file,
                            index_dir=org_2_id_ref_index
                        )
                        cmd_index_3 = compile_create_star_index(
                            genome_fasta=org_3_id_ref_file,
                            index_dir=org_3_id_ref_index
                        )

                        # STAR align commands
                        cmd_align_org = compile_star_align_cmd(
                            sample=sample_path,
                            outdir=sample_result_dir,
                            genomeDir=org_id_ref_index
                        )
                        cmd_align_org_1 = compile_star_align_cmd(
                            sample=sample_path,
                            outdir=sample_result_dir,
                            genomeDir=org_1_id_ref_index
                        )
                        cmd_align_org_2 = compile_star_align_cmd(
                            sample=sample_path,
                            outdir=sample_result_dir,
                            genomeDir=org_2_id_ref_index
                        )
                        cmd_align_org_3 = compile_star_align_cmd(
                            sample=sample_path,
                            outdir=sample_result_dir,
                            genomeDir=org_3_id_ref_index
                        )

                        file.write(" ".join(cmd_index_org) + "\n")
                        file.write(" ".join(cmd_index_1) + "\n")
                        file.write(" ".join(cmd_index_2) + "\n")
                        file.write(" ".join(cmd_index_3) + "\n")
                        file.write(" ".join(cmd_align_org) + "\n")
                        file.write(" ".join(cmd_align_org_1) + "\n")
                        file.write(" ".join(cmd_align_org_2) + "\n")
                        file.write(" ".join(cmd_align_org_3) + "\n")

                    else:
                        sample_path_1 = row['fq1']
                        sample_path_2 = row['fq2']

                        # compile indexing commands
                        org_id_ref_genome = os.path.join(ref_genomes_dir, str(row['org_id']))
                        filename = os.listdir(org_id_ref_genome)[0]
                        org_id_ref_file = os.path.join(org_id_ref_genome, filename)
                        org_id_ref_index = os.path.join(sample_result_dir, str(row['org_id']), "index")

                        org_1_id_ref_genome_1 = os.path.join(ref_genomes_dir, str(row['1_org_1_id']))
                        filename = os.listdir(org_1_id_ref_genome_1)[0]
                        org_1_id_ref_file_1 = os.path.join(org_1_id_ref_genome_1, filename)
                        org_1_id_ref_index_1 = os.path.join(sample_result_dir, str(row['1_org_1_id']), "index")

                        org_2_id_ref_genome_1 = os.path.join(ref_genomes_dir, str(row['1_org_2_id']))
                        filename = os.listdir(org_2_id_ref_genome_1)[0]
                        org_2_id_ref_file_1 = os.path.join(org_2_id_ref_genome_1, filename)
                        org_2_id_ref_index_1 = os.path.join(sample_result_dir, str(row['1_org_2_id']), "index")

                        org_3_id_ref_genome_1 = os.path.join(ref_genomes_dir, str(row['1_org_3_id']))
                        filename = os.listdir(org_3_id_ref_genome_1)[0]
                        org_3_id_ref_file_1 = os.path.join(org_3_id_ref_genome_1, filename)
                        org_3_id_ref_index_1 = os.path.join(sample_result_dir, str(row['1_org_3_id']), "index")

                        org_1_id_ref_genome_2 = os.path.join(ref_genomes_dir, str(int(row['2_org_1_id'])))
                        filename = os.listdir(org_1_id_ref_genome_1)[0]
                        org_1_id_ref_file_2 = os.path.join(org_1_id_ref_genome_2, filename)
                        org_1_id_ref_index_2 = os.path.join(sample_result_dir, str(int(row['2_org_1_id'])), "index")

                        org_2_id_ref_genome_2 = os.path.join(ref_genomes_dir, str(int(row['2_org_2_id'])))
                        filename = os.listdir(org_2_id_ref_genome_2)[0]
                        org_2_id_ref_file_2 = os.path.join(org_2_id_ref_genome_2, filename)
                        org_2_id_ref_index_2 = os.path.join(sample_result_dir, str(int(row['2_org_2_id'])), "index")

                        org_3_id_ref_genome_2 = os.path.join(ref_genomes_dir, str(int(row['2_org_3_id'])))
                        filename = os.listdir(org_3_id_ref_genome_2)[0]
                        org_3_id_ref_file_2 = os.path.join(org_3_id_ref_genome_2, filename)
                        org_3_id_ref_index_2 = os.path.join(sample_result_dir, str(int(row['2_org_3_id'])), "index")

                        cmd_index_org = compile_create_star_index(
                            genome_fasta=org_id_ref_file,
                            index_dir=org_id_ref_index
                        )
                        cmd_index_1_1 = compile_create_star_index(
                            genome_fasta=org_1_id_ref_file_1,
                            index_dir=org_1_id_ref_index_1
                        )
                        cmd_index_2_1= compile_create_star_index(
                            genome_fasta=org_2_id_ref_file_1,
                            index_dir=org_2_id_ref_index_1
                        )
                        cmd_index_3_1 = compile_create_star_index(
                            genome_fasta=org_3_id_ref_file_1,
                            index_dir=org_3_id_ref_index_1
                        )
                        cmd_index_1_2 = compile_create_star_index(
                            genome_fasta=org_1_id_ref_file_2,
                            index_dir=org_1_id_ref_index_2
                        )
                        cmd_index_2_2= compile_create_star_index(
                            genome_fasta=org_2_id_ref_file_2,
                            index_dir=org_2_id_ref_index_2
                        )
                        cmd_index_3_2 = compile_create_star_index(
                            genome_fasta=org_3_id_ref_file_2,
                            index_dir=org_3_id_ref_index_2
                        )
                        

                        # STAR align commands
                        cmd_align_org_1 = compile_star_align_cmd(
                            sample=sample_path_1,
                            outdir=sample_result_dir,
                            genomeDir=org_id_ref_index
                        )
                        cmd_align_org_1_1 = compile_star_align_cmd(
                            sample=sample_path_1,
                            outdir=sample_result_dir,
                            genomeDir=org_1_id_ref_index_1
                        )
                        cmd_align_org_2_1 = compile_star_align_cmd(
                            sample=sample_path_1,
                            outdir=sample_result_dir,
                            genomeDir=org_2_id_ref_index_1
                        )
                        cmd_align_org_3_1 = compile_star_align_cmd(
                            sample=sample_path_1,
                            outdir=sample_result_dir,
                            genomeDir=org_3_id_ref_index_1
                        )
                        cmd_align_org_2 = compile_star_align_cmd(
                            sample=sample_path_2,
                            outdir=sample_result_dir,
                            genomeDir=org_id_ref_index
                        )
                        cmd_align_org_1_2 = compile_star_align_cmd(
                            sample=sample_path_2,
                            outdir=sample_result_dir,
                            genomeDir=org_1_id_ref_index_2
                        )
                        cmd_align_org_2_2 = compile_star_align_cmd(
                            sample=sample_path_2,
                            outdir=sample_result_dir,
                            genomeDir=org_2_id_ref_index_2
                        )
                        cmd_align_org_3_2 = compile_star_align_cmd(
                            sample=sample_path_2,
                            outdir=sample_result_dir,
                            genomeDir=org_3_id_ref_index_2
                        )

                        file.write(" ".join(cmd_index_org) + "\n")
                        file.write(" ".join(cmd_index_1_1) + "\n")
                        file.write(" ".join(cmd_index_2_1) + "\n")
                        file.write(" ".join(cmd_index_3_1) + "\n")
                        file.write(" ".join(cmd_index_1_2) + "\n")
                        file.write(" ".join(cmd_index_2_2) + "\n")
                        file.write(" ".join(cmd_index_3_2) + "\n")

                        file.write(" ".join(cmd_align_org_1) + "\n")
                        file.write(" ".join(cmd_align_org_1_1) + "\n")
                        file.write(" ".join(cmd_align_org_2_1) + "\n")
                        file.write(" ".join(cmd_align_org_3_1) + "\n")
                        file.write(" ".join(cmd_align_org_2) + "\n")
                        file.write(" ".join(cmd_align_org_1_2) + "\n")
                        file.write(" ".join(cmd_align_org_2_2) + "\n")
                        file.write(" ".join(cmd_align_org_3_2) + "\n")
            
            sp.run(['sbatch', job_script_filename])
            time.sleep(2)

def compile_create_star_index(
        genome_fasta: Path,
        index_dir,
        chr_bin_bits: int = 18,
        index_string_size: int = 5,
        threads: int = 8,
) -> List[str]:
    """Prepare STAR index."""

    cmd = [
        "STAR",
        "--runMode", "genomeGenerate",
        "--genomeSAindexNbases", f"{str(index_string_size)}",
        "--genomeChrBinNbits", f"{str(chr_bin_bits)}",
        "--runThreadN", f"{str(threads)}",
        "--genomeDir", f"{str(index_dir)}",
        "--genomeFastaFiles", f"{str(genome_fasta)}",
    ]
    return cmd

def compile_star_align_cmd(
        sample,
        outdir,
        genomeDir,
        threads: int = 8,
) -> List[str]:
    """Compile command for STAR alignment."""
    cmd_hts = ["STAR"]
    cmd_hts.extend(["--alignIntronMax", "1"])
    cmd_hts.extend(["--alignEndsType", "Local"])
    cmd_hts.extend(["--runThreadN", str(threads)])
    cmd_hts.extend(["--genomeDir", str(genomeDir)])
    cmd_hts.extend(["--outFilterMultimapNmax", "50"])
    cmd_hts.extend(["--outSAMorder", "PairedKeepInputOrder"])
    cmd_hts.extend(["--outSAMunmapped", "Within", "KeepPairs"])
    cmd_hts.extend(["--readFilesCommand", "zcat"])
    cmd_hts.extend(["--readFilesIn", str(sample)])
    cmd_hts.extend(["--outFileNamePrefix", str(outdir)])
    return cmd_hts

generate_mapping_slurm_script("test1", "mapping_tests/ref_genomes", "mismatch_mappings.csv", blocks=1)
