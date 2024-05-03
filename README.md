# HTSinfer testing

This the repository for the testing [HTSinfer](https://github.com/zavolanlab/htsinfer). 
It contains the data and the scripts used to generate the plots.

## Installation

### 1. Clone the repository

Go to the desired directory/folder on your file system, then clone/get the 
repository and move into the respective directory with:

```bash
git clone https://github.com/zavolanlab/test_htsinfer
cd test_htsinfer
```

### 2. Conda and Mamba installation

Workflow dependencies can be conveniently installed with the [Conda](https://docs.conda.io/projects/conda/en/stable/)
package manager. We recommend that you install [Miniconda](https://docs.anaconda.com/free/miniconda/miniconda-install/)
for your system (Linux). Be sure to select the Python 3 option. 

```bash
conda install -y mamba -n base -c conda-forge
```

### 3. Create environment

Install the remaining dependencies with:
```bash
mamba env create -f environment.yml
```

### 4. Activate environment

Activate the Conda environment with:

```bash
conda activate test_htsinfer
```

## Running the workflow

To start the workflow, in your activated `test_htsinfer` conda environment, run

```bash
nextflow main.nf -profile conda
```
> Currently, running the workflow is only supported with conda

For running on SLURM:
```bash
nextflow main.nf -profile slurm,conda
```