# core-stress-transcriptome
Project Title: Exploring core abiotic stress responses in *Zea mays* through meta-analysis of publicly available RNA-seq data.

Background: Abiotic stressors such as drought, salt, flooding, cold, heat, and low nitrogen are major constraints on agricultural productivity, including that of maize (*Zea mays*), one of the most important global crops. Gene expression changes under these stress conditions have been extensively studied, leading to a large amount of publicly available RNA-sequencing data for stress-treated maize. In this project, I seek to leverage these public data to explore conserved gene expression responses to the six stressors listed above: the core stress response.

### Repository Contents
All code, scripts, and Jupyter notebooks are contained in the `scripts` directory. Please note this is an active project and changes will be made frequently.  
Within `scripts` are the following directories:
- `exploratory_analyses`: Jupyter notebooks for exploratory principal component analysis and hierarchical clustering, just after RNA-seq data processing.
- `mapping_rate_exploration`: Python script to assemble a file with mapping rates from a Salmon output directory (`find_map_rates.py`) and Jupyter notebooks for finding which samples have low mapping rates.
- `tximport`: R script for running tximport to calculate transcripts per million (TPM) on the Michigan State University HPCC.
