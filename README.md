# HTSinfer testing

This the repository for the testing [HTSinfer](https://github.com/zavolanlab/htsinfer). 
It contains the data and the scripts used to generate the plots.

### Structure

- 0_zavolan_rnaseq_samples_extract.py - Extract the metadata from the Zavolan lab sequenced RNA-Seq samples
- 1_cluster_test_strategy.py - Run HTSinfer on downloaded samples on the cluster
- 2_get_results.py - Obtain the results for predicted org/adapter/orientation etc. as well as performance metrics
- 3_create_stat_plots.ipynb - Notebook collection of the generated stats/plots.

The data used:

- mined_test_data.tsv - Collection of 770 Illumina RNA-Seq samples downloaded from SRA, from various organisms.
- zavolan_rnaseq_samples.tsv - Extracted from sequencing samples used in the Zavolan lab.
- zavolan_rnaseq_samples_filtered.tsv - A filtered table of 101 samples where the metadata was correct, used for testing.

A GoogleDoc with the progress and ideas for the tool [can be found here.](https://docs.google.com/document/d/12vBwjZ7N6aS9bBJq-hs5RKqWOneRg3D8Q4-41elksR8/edit?usp=sharing)

This repo is continously updated to keep track with the publication.