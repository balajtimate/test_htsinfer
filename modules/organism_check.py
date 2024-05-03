import pandas as pd
import subprocess

# Read the TSV file into a pandas DataFrame
tsv_file_path = 'org_check_data.tsv'  # Replace with the actual path to your TSV file
df = pd.read_csv(tsv_file_path, sep='\t')

# Define the function to run the bash command and extract information
def run_bash_command(sample):
    command = f"esearch -db sra -query '{sample}' | esummary | xtract -pattern DocumentSummary -element Experiment@name,LIBRARY_STRATEGY,Organism@taxid,Organism@ScientificName"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    return result.stdout.strip()

# Create a new DataFrame to store the results
result_df = pd.DataFrame(columns=['Sample', 'Experiment_name', 'Library_strategy', 'Tax_id', 'ScientificName'])

# Iterate over each sample and run the bash command
for index, row in df.iterrows():
    sample = row['sample']
    result = run_bash_command(sample)
    print(result)
    
    # Parse the result and add it to the result DataFrame
    result_list = result.split()
    if len(result_list) == 4:
        result_df = result_df.append(
            {'Sample': sample, 'Experiment_name': result_list[0], 'Library_strategy': result_list[1], 'Tax_id': result_list[2], 'ScientificName': result_list[3]},
            ignore_index=True
        )

# Save the result DataFrame to a CSV file
result_csv_path = 'organism_check_mismatch.csv'  # Replace with the desired output CSV file path
result_df.to_csv(result_csv_path, index=False)

print(f"Results saved to {result_csv_path}")
