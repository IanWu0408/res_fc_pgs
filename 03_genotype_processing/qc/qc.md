# UKB RAP Genotype QC Flow

Author: Tung-Yen Wu  
Platform: UK Biobank Research Analysis Platform (RAP)  
Tool: Swiss Army Knife (PLINK 2.0, v5.0.0)

This file records the genotype QC steps performed on UK Biobank data.
All commands were executed on UKB RAP. Paths and project-specific IDs
are omitted or generalized for portability.

## Merge chromosomes (example: chr1–chr3)

```bash
# Create merge list
for chr in {2..3}; do
  echo "ukb21007_c${chr}_b0_v1.cleaned" >> mergelist.txt
done

# Merge chromosomes
plink2 \
  --pfile ukb21007_c1_b0_v1.cleaned \
  --set-all-var-ids '@:#:$r:$a' \
  --new-id-max-allele-len 70 \
  --max-alleles 2 \
  --pmerge-list mergelist.txt \
  --make-pgen \
  --threads 4 \
  --memory 30000 \
  --out ukb21007_merged_chr_1_3
```

## Sample QC
```bash
plink2 \
  --pfile ukb21007_merged_chr_1_22 \
  --remove Sexcheck_Het_Missing_Ethnicity.txt \
  --threads 32 \
  --memory 220000 \
  --make-pgen \
  --out ukb21007_sampleQC
```
## SNP QC
```bash
# Identify variants with imputation INFO/R2 < 0.8
for chr in {1..22}; do
  bcftools query -f '%ID\t%INFO/R2\n' \
    ukb21007_c${chr}_b0_v1.sites.vcf.gz \
  | awk '$2 < 0.8 { print $1 }' \
  >> all_snps_r2_lt_08.txt
done

# Exclude low-quality variants
plink2 \
  --pfile ukb21007_sampleQC \
  --exclude all_snps_r2_lt_08.txt \
  --make-pgen \
  --out ukb21007_filtered_R2_08

# Standard SNP-level QC
plink2 \
  --pfile ukb21007_filtered_R2_08 \
  --maf 0.01 \
  --hwe 1e-6 \
  --geno 0.02 \
  --threads 32 \
  --memory 220000 \
  --make-pgen \
  --out ukb21007_QC_done
```
## Export to PLINK BED format
```bash
plink2 \
  --pfile ukb21007_QC_done \
  --make-bed \
  --threads 32 \
  --memory 220000 \
  --out ukb21007_QC_done
```

