"""Extract metadata for samples"""

import pandas as pd
import csv
import subprocess
import os

samples_path = '/scicore/projects/openbis/userstore/biozentrum_zavolan'
result_path = '~/Work/test_htsinfer/zavolan_rnaseq_samples.tsv'

# List of column names to extract
columns_to_extract = [
    "EXTERNAL_SAMPLE_NAME",
    "KIT",
    "NCBI_ORGANISM_TAXONOMY",
    "SAMPLE_CODE",
    "SAMPLE_KIND",
    "SAMPLE_TYPE",
    "SEQUENCING_APPLICATION",
    "CYCLES_REQUESTED_BY_CUSTOMER",
    "END_TYPE",
]

# Initialize an empty dictionary to store extracted data
extracted_data = {col: [] for col in columns_to_extract}
fastq_cols_list = []
subdirectory_list = []

# Find TSV files containing "RNA_SEQ" from 2023/2022/2021
find_command = f'find {samples_path} \( -path "{samples_path}/2023*" ' \
               f'-o -path "{samples_path}/2022*" ' \
               f'-o -path "{samples_path}/2021*" \) ' \
               f'-type f -name "*.tsv" -exec grep -l "RNA_SEQ" {{}} \;'
output = subprocess.check_output(find_command, shell=True)
tsv_files = output.decode().splitlines()

for i, file in enumerate(tsv_files):
    # Extract subdirectory name
    subdirectory = os.path.dirname(file).replace(samples_path + "/", "")
    subdirectory_list.append(subdirectory)

    # Read the TSV file and transpose the DataFrame
    df = pd.read_csv(
        file, sep='\t', header=None,
        index_col=0, on_bad_lines='skip'
        ).transpose()

    # Extract the required columns
    for col in columns_to_extract:
        if col in df.columns:
            extracted_values = df[col].values
            extracted_data[col].extend(extracted_values)
        else:
            extracted_data[col].extend([None] * len(df))

    # Extract the fastq file paths
    fastq_cols = []
    if "FASTQ_FILES" in df.columns and "DATASET PROPERTIES" in df.columns:
        start_idx = df.columns.get_loc("FASTQ_FILES") + 1
        end_idx = df.columns.get_loc("DATASET PROPERTIES")
        if start_idx < end_idx:
            fastq_cols = df.iloc[:, start_idx:end_idx].columns.tolist()

    # Append FASTQ column headers to the list
    fastq_cols_list.append(fastq_cols)

# Create separate columns for FASTQ file headers
max_cols = max(len(cols) for cols in fastq_cols_list)
for j in range(1, max_cols + 1):
    extracted_data[f"FASTQ_FILE{j}"] = [
        cols[j-1] if j <= len(cols) else None for cols in fastq_cols_list
        ]

# Add subdirectory column
extracted_data["SUBDIRECTORY"] = subdirectory_list

# Create a DataFrame from the extracted data
data = pd.DataFrame(extracted_data)

# Write the DataFrame to a CSV file
data.to_csv(f'{result_path}', index=False, sep='\t', quoting=csv.QUOTE_ALL)
