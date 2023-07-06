"""Get results from HTSinfer run."""

import os
import json
import pandas as pd
import numpy as np
import subprocess as sp
from datetime import datetime
from pathlib import Path

run_id = str(input("Run ID: "))

RESULTS_HTS_DIR = (Path(__file__).resolve().parent /
                   "results_htsinfer")
MINED_DATA = (Path(__file__).resolve().parent /
              "mined_test_data.tsv")

# Read in the tsv file
data = pd.read_csv(MINED_DATA, sep='\t')

table_data = {}

# Construct results dataframe from json files
for index, row in data.iterrows():
    try:
        sample = row['sample']
        layout = row['layout']
        results_folder = run_id

        if layout == 'SE':
            # 1. Read in final results from _results
            with open(
                os.path.join(RESULTS_HTS_DIR, results_folder, "_results",
                             f"{sample}_result.json"), encoding="utf-8"
                        ) as json_file:
                result_model = json.load(json_file)
                table_data[index] = {
                    'pred_org': result_model['library_source']['file_1']['short_name'],
                    'pred_orient': result_model['read_orientation']['file_1'],
                    'pred_adapter': result_model['read_layout']['file_1']['adapt_3'],
                    'pred_length_min_1': result_model['library_stats']['file_1']['read_length']['min'],
                    'pred_length_max_1': result_model['library_stats']['file_1']['read_length']['max'],
                    'pred_length_min_2': np.nan,
                    'pred_length_max_2': np.nan
                }
            # 2. Read in library_source results
            with open(
                os.path.join(
                    RESULTS_HTS_DIR, results_folder, "_results_library_source",
                    f"library_source_{sample}.fastq.json"
                ), encoding="utf-8"
            ) as lib_file:
                lib_model = json.load(lib_file)
                table_data[index].update({
                    '1_tpm_1': lib_model['data'][0][0],
                    '1_org_1': lib_model['data'][0][1][0],
                    '1_tpm_2': lib_model['data'][1][0],
                    '1_org_2': lib_model['data'][1][1][0],
                    '2_tpm_1': np.nan,
                    '2_org_1': np.nan,
                    '2_tpm_2': np.nan,
                    '2_org_2': np.nan
                })
            # 3. Read in read_layout results
            with open(
                os.path.join(
                    RESULTS_HTS_DIR, results_folder, "_results_read_layout",
                    f"read_layout_{sample}.fastq.json"
                ), encoding="utf-8"
            ) as read_file:
                read_model = json.load(read_file)
                table_data[index].update({
                    '1_adapt_1': read_model['data'][0][0],
                    '1_percent_1': read_model['data'][0][1],
                    '1_adapt_2': read_model['data'][1][0],
                    '1_percent_2': read_model['data'][1][1],
                    '2_adapt_1': np.nan,
                    '2_percent_1': np.nan,
                    '2_adapt_2': np.nan,
                    '2_percent_2': np.nan
                })

        elif layout == 'PE':
            # 1. Read in final results from _results
            with open(
                os.path.join(
                    RESULTS_HTS_DIR, results_folder, "_results",
                    f"{sample}_result.json"
                ), encoding="utf-8"
            ) as json_file:
                result_model = json.load(json_file)
                table_data[index] = {
                    'pred_org': result_model['library_source']['file_1']['short_name'],
                    'pred_orient': result_model['read_orientation']['relationship'],
                    'pred_adapter': result_model['read_layout']['file_1']['adapt_3'],
                    'pred_length_min_1': result_model['library_stats']['file_1']['read_length']['min'],
                    'pred_length_max_1': result_model['library_stats']['file_1']['read_length']['max'],
                    'pred_length_min_2': result_model['library_stats']['file_2']['read_length']['min'],
                    'pred_length_max_2': result_model['library_stats']['file_2']['read_length']['max']
                }
            # 2. Read in library_source results
            with open(
                os.path.join(
                    RESULTS_HTS_DIR, results_folder, "_results_library_source",
                    f"library_source_{sample}_1.fastq.json"
                ), encoding="utf-8"
            ) as json_file_1, open(
                os.path.join(
                    RESULTS_HTS_DIR, results_folder, "_results_library_source",
                    f"library_source_{sample}_2.fastq.json"
                ), encoding="utf-8"
            ) as json_file_2:
                lib_model_1 = json.load(json_file_1)
                lib_model_2 = json.load(json_file_2)
                table_data[index].update({
                    '1_tpm_1': lib_model_1['data'][0][0],
                    '1_org_1': lib_model_1['data'][0][1][0],
                    '1_tpm_2': lib_model_1['data'][1][0],
                    '1_org_2': lib_model_1['data'][1][1][0],
                    '2_tpm_1': lib_model_2['data'][0][0],
                    '2_org_1': lib_model_2['data'][0][1][0],
                    '2_tpm_2': lib_model_2['data'][1][0],
                    '2_org_2': lib_model_2['data'][1][1][0]
                })
            # 3. Read in read_layout results
            with open(
                os.path.join(
                    RESULTS_HTS_DIR, results_folder, "_results_read_layout",
                    f"read_layout_{sample}_1.fastq.json"
                ), encoding="utf-8"
            ) as read_file_1, open(
                os.path.join(
                    RESULTS_HTS_DIR, results_folder, "_results_read_layout",
                    f"read_layout_{sample}_2.fastq.json"
                ), encoding="utf-8"
            ) as read_file_2:
                read_model_1 = json.load(read_file_1)
                read_model_2 = json.load(read_file_2)
                table_data[index].update({
                    '1_adapt_1': read_model_1['data'][0][0],
                    '1_percent_1': read_model_1['data'][0][1],
                    '1_adapt_2': read_model_1['data'][1][0],
                    '1_percent_2': read_model_1['data'][1][1],
                    '2_adapt_1': read_model_2['data'][0][0],
                    '2_percent_1': read_model_2['data'][0][1],
                    '2_adapt_2': read_model_2['data'][1][0],
                    '2_percent_2': read_model_2['data'][1][1]
                })

        # Additional code to read error file
        error_file = os.path.join(
            RESULTS_HTS_DIR, results_folder, "_results",
            f"{sample}_error.txt"
        )

        # 1. Processing reads
        proc_start = datetime.strptime(sp.check_output(
            f"grep 'Processing read file 1:' {error_file} | grep -oE '[0-9]{{4}}-[0-9]{{2}}-[0-9]{{2}} [0-9]{{2}}:[0-9]{{2}}:[0-9]{{2}}'",
            shell=True
        ).decode("utf-8").strip().splitlines()[0], '%Y-%m-%d %H:%M:%S')
        proc_end = datetime.strptime(sp.check_output(
            f"awk -v n=2 '/Processing read file 1:/ {{ for (i = 1; i <= n; i++) getline; print }}' {error_file} | grep -oE '[0-9]{{4}}-[0-9]{{2}}-[0-9]{{2}} [0-9]{{2}}:[0-9]{{2}}:[0-9]{{2}}'",
            shell=True
        ).decode("utf-8").strip().splitlines()[0], '%Y-%m-%d %H:%M:%S')
        table_data[index].update({
            'processing_time': (proc_end - proc_start).total_seconds(),
        })

        # 2. Extract read length
        extract_start = datetime.strptime(sp.check_output(
            f"awk -v n=2 '/Determining library statistics/ {{ for (i = 1; i <= n; i++) getline; print }}' {error_file} | grep -oE '[0-9]{{4}}-[0-9]{{2}}-[0-9]{{2}} [0-9]{{2}}:[0-9]{{2}}:[0-9]{{2}}'",
            shell=True
        ).decode("utf-8").strip().splitlines()[0], '%Y-%m-%d %H:%M:%S')
        extract_end = datetime.strptime(sp.check_output(
            f"awk -v n=3 '/Determining library statistics/ {{ for (i = 1; i <= n; i++) getline; print }}' {error_file} | grep -oE '[0-9]{{4}}-[0-9]{{2}}-[0-9]{{2}} [0-9]{{2}}:[0-9]{{2}}:[0-9]{{2}}'",
            shell=True
        ).decode("utf-8").strip().splitlines()[0], '%Y-%m-%d %H:%M:%S')
        table_data[index].update({
            'extract_time': (extract_end - extract_start).total_seconds(),
        })

        # 3. Kallisto quantification
        kallisto_start = datetime.strptime(sp.check_output(
            f"awk -v n=3 '/Creating kallisto index for:/ {{ for (i = 1; i <= n; i++) getline; print }}' {error_file} | grep -oE '[0-9]{{4}}-[0-9]{{2}}-[0-9]{{2}} [0-9]{{2}}:[0-9]{{2}}:[0-9]{{2}}'",
            shell=True
        ).decode("utf-8").strip().splitlines()[0], '%Y-%m-%d %H:%M:%S')
        kallisto_end = datetime.strptime(sp.check_output(
            f"awk -v n=4 '/Creating kallisto index for:/ {{ for (i = 1; i <= n; i++) getline; print }}' {error_file} | grep -oE '[0-9]{{4}}-[0-9]{{2}}-[0-9]{{2}} [0-9]{{2}}:[0-9]{{2}}:[0-9]{{2}}'",
            shell=True
        ).decode("utf-8").strip().splitlines()[0], '%Y-%m-%d %H:%M:%S')
        table_data[index].update({
            'kallisto_time': (kallisto_end - kallisto_start).total_seconds(),
        })

        # 4. Align reads with STAR
        align_start = datetime.strptime(sp.check_output(
            f"grep 'Aligning reads with STAR' {error_file} | grep -oE '[0-9]{{4}}-[0-9]{{2}}-[0-9]{{2}} [0-9]{{2}}:[0-9]{{2}}:[0-9]{{2}}'",
            shell=True
        ).decode("utf-8").strip().splitlines()[0], '%Y-%m-%d %H:%M:%S')
        align_end = datetime.strptime(sp.check_output(
            f"awk -v n=1 '/Aligning reads with STAR/ {{ for (i = 1; i <= n; i++) getline; print }}' {error_file} | grep -oE '[0-9]{{4}}-[0-9]{{2}}-[0-9]{{2}} [0-9]{{2}}:[0-9]{{2}}:[0-9]{{2}}'",
            shell=True
        ).decode("utf-8").strip().splitlines()[0], '%Y-%m-%d %H:%M:%S')
        table_data[index].update({
            'align_time': (align_end - align_start).total_seconds(),
        })

        # 5. Parse with Cutadapt
        cutadapt_start = datetime.strptime(sp.check_output(
            f"awk -v n=2 '/Determining read layout/ {{ for (i = 1; i <= n; i++) getline; print }}' {error_file} | grep -oE '[0-9]{{4}}-[0-9]{{2}}-[0-9]{{2}} [0-9]{{2}}:[0-9]{{2}}:[0-9]{{2}}'",
            shell=True
        ).decode("utf-8").strip().splitlines()[0], '%Y-%m-%d %H:%M:%S')
        cutadapt_end = datetime.strptime(sp.check_output(
            f"awk -v n=3 '/Determining read layout/ {{ for (i = 1; i <= n; i++) getline; print }}' {error_file} | grep -oE '[0-9]{{4}}-[0-9]{{2}}-[0-9]{{2}} [0-9]{{2}}:[0-9]{{2}}:[0-9]{{2}}'",
            shell=True
        ).decode("utf-8").strip().splitlines()[0], '%Y-%m-%d %H:%M:%S')
        table_data[index].update({
            'cutadapt_time': (cutadapt_end - cutadapt_start).total_seconds(),
        })

        # 6. Total time
        total_start = datetime.strptime(sp.check_output(
            f"grep 'Started HTSinfer' {error_file} | grep -oE '[0-9]{{4}}-[0-9]{{2}}-[0-9]{{2}} [0-9]{{2}}:[0-9]{{2}}:[0-9]{{2}}'",
            shell=True
        ).decode("utf-8").strip().splitlines()[0], '%Y-%m-%d %H:%M:%S')
        total_end = datetime.strptime(sp.check_output(
            f"grep 'INFO] Done' {error_file} | grep -oE '[0-9]{{4}}-[0-9]{{2}}-[0-9]{{2}} [0-9]{{2}}:[0-9]{{2}}:[0-9]{{2}}'",
            shell=True
        ).decode("utf-8").strip().splitlines()[0], '%Y-%m-%d %H:%M:%S')
        table_data[index].update({
            'total_time': (total_end - total_start).total_seconds(),
        })

    except (
        FileNotFoundError, json.decoder.JSONDecodeError, sp.CalledProcessError
        ):
        table_data[index] = {
            key: np.nan for key in [
                'pred_org', 'pred_orient', 'pred_adapter',
                'pred_length_min_1', 'pred_length_max_1',
                '1_tpm_1', '1_org_1', '1_tpm_2', '1_org_2',
                '2_tpm_1', '2_org_1', '2_tpm_2', '2_org_2',
                '1_adapt_1', '1_percent_1', '1_adapt_2', '1_percent_2',
                '2_adapt_1', '2_percent_1', '2_adapt_2', '2_percent_2',
                'processing_time', 'extract_time', 'kallisto_time',
                'align_time',  'cutadapt_time', 'total_time'
                ]
                }
    continue

print(pd.DataFrame.from_dict(table_data, orient='index'))

# Concatenate the original data with the graph data
final_result = pd.concat([data, pd.DataFrame.from_dict(table_data, orient='index')], axis=1, join="inner")

# Comparison of results
final_result['match_org'] = np.where(
                final_result['org'] == final_result['pred_org'],
                True,
                np.where(
                    pd.isna(final_result['pred_org']),
                    'Undecided',
                    False
                )
            )
final_result['match_orient'] = np.where(
    final_result['pred_orient'] ==
    final_result['orient'], True, False)
final_result['match_adapter'] = final_result.apply(
    lambda x: str(x.pred_adapter) in str(x.adapter), axis=1)
final_result['match_length'] = np.where(
    final_result['pred_length_max_1'] ==
    final_result['length_max_1'], True, False)

# Write result to csv file
pd.DataFrame.to_csv(final_result, f"{run_id}_result.csv",
                    index=False)
