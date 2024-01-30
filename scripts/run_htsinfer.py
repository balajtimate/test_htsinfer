# run_htsinfer.py
import csv
import subprocess

def run_htsinfer(sample_table_path, verbosity, min_match_percentage, lib_type):
    with open(sample_table_path, "r") as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            sample = row["Sample"]
            fq1 = row["fq1"]
            fq2 = row["fq2"] if "fq2" in row else None
            script_path = generate_slurm_script(sample, fq1, fq2, verbosity, min_match_percentage, lib_type)
            subprocess.run(["sbatch", script_path])

# Example usage:
run_htsinfer("sample_tables/sample_table.csv", "DEBUG", "0.5", "SINGLE")
