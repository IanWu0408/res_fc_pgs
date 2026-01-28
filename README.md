# res_fc_pgs

This repository contains the full analysis pipeline for studying the association between
**polygenic scores (PGS)** and **resting-state functional connectivity (rsFC)**, with a focus on
resilience.

The codebase is organized to clearly separate:
1. **core algorithms**,  
2. **parameterized execution scripts**, and
3. **data preprocessing**
to facilitate transparency, reproducibility, and reuse.

---

## Repository Structure
### 01_core_analysis
Core analytical functions and computational algorithms used in the paper.  
This folder contains **model definitions and statistical procedures only**, without hard-coded
dataset-specific parameters.

> Think of this folder as the “methods section in code”.

---

### 02_run_analysis
High-level scripts that **execute the core algorithms** using predefined parameters.  
These scripts specify:
- phenotype definitions  
- covariates  
- cross-validation settings  
- model configurations  

and are intended to reproduce the main results reported in the paper.

---

### 03_genotype_processing
Code for preprocessing **genotype data** from the **UK Biobank (UKB)**, including:
- Genotype quality control (QC)
- Genome wide-association study (GWAS)
- preparation of inputs for polygenic score construction  

---

### 04_fc_processing
Code for preprocessing **resting-state functional connectivity (rsFC)** data from UKB and Human Connectome Project (HCP).

This folder includes scripts for:
- functional connectivity matrix construction  
- quality control and harmonization across datasets  

---

## Notes
- This repository does **not** contain raw UK Biobank or HCP data.
- Access to these datasets must comply with their respective data use agreements.
- File paths and environment-specific settings should be adapted by the user.

---

## Citation
If you use this code, please cite the corresponding paper:

> *[Citation to be added]*

