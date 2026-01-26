# UKB RAP Genotype QC Flow

Author: Tung-Yen Wu  
Platform: UK Biobank Research Analysis Platform (RAP)  
Tool: Swiss Army Knife (v5.0.0)

Genome-wide association analysis was performed using **BOLT-LMM**, a
linear mixed model (LMM) approach.

## Prepare modelSnps for BOLT-LMM (ukb21007_pruned.prune.in.txt)
```bash
plink --bfile ukb21007_QC_done --indep-pairwise 50 5 0.2 --out ukb21007_pruned
```

## Phenotype and covariate file (PhenoCov.txt)

The file **PhenoCov.txt** is a tab-delimited table containing the following
required columns:

- **FID / IID**  
  Family and individual identifiers matching the genotype files.

- **INTresilience**  
  Inverse-normal–transformed resilience phenotype.  
  Resilience was defined as the **sum score of six BRS items**, and the
  inverse normal transformation was performed **within the GWAS cohort**
  prior to association analysis.

- **Age**  
  Age at the time of completing the BRS questionnaire.

- **Age2**  
  Squared age term (Age²), included to model non-linear age effects.

- **Sex**  
  Sex recorded at **UK Biobank instance 0**.

- **Batch**  
  **Genotype batch**, included as a categorical covariate to account for
  batch effects in genotyping.

- **PC1–PC20**  
  The first 20 genetic principal components provided by UK Biobank,
  used to control for population stratification.

  
## Run BOLT-LMM
```bash
bolt \
  --bfile ukb21007_QC_done \
  --bgenFile ukb21007_QC_done.bgen \
  --sampleFile ukb21007_QC_done.sample \
  --phenoFile PhenoCov.txt \
  --phenoCol INTresilience \
  --covarFile PhenoCov.txt \
  --covarCol Sex \
  --covarCol Batch \
  --qCovarCol Age \
  --qCovarCol Age2 \
  --qCovarCol PC1 \
  --qCovarCol PC2 \
  --qCovarCol PC3 \
  --qCovarCol PC4 \
  --qCovarCol PC5 \
  --qCovarCol PC6 \
  --qCovarCol PC7 \
  --qCovarCol PC8 \
  --qCovarCol PC9 \
  --qCovarCol PC10 \
  --qCovarCol PC11 \
  --qCovarCol PC12 \
  --qCovarCol PC13 \
  --qCovarCol PC14 \
  --qCovarCol PC15 \
  --qCovarCol PC16 \
  --qCovarCol PC17 \
  --qCovarCol PC18 \
  --qCovarCol PC19 \
  --qCovarCol PC20 \
  --modelSnps ukb21007_pruned.prune.in.txt \
  --LDscoresFile LDSCORE.1000G_EUR.GRCh38.tab.gz \
  --LDscoresMatchBp \
  --geneticMapFile genetic_map_hg38_withX.txt \
  --lmmForceNonInf \
  --numThreads 16 \
  --statsFile resilience_array_white_hg38.stats.gz \
  --statsFileBgenSnps resilience_dosage.stats.gz
```
