
# Generate SNP Helper Map on UKB RAP (Posit Workbench)

Author: Tung-Yen Wu  
Platform: **UK Biobank Research Analysis Platform (RAP)**  
Environment: **Posit Workbench (RStudio)**  
Language: **R**

This script generates a SNP helper mapping file (`map.txt`) from UK Biobank
imputed genotype *sites* VCF files.  
The resulting map is used for downstream GWAS summary statistics
harmonization (e.g. aligning `CHR:BP` to rsID, REF/ALT, and allele frequency).

The workflow is designed to be executed **within the UKB RAP Posit Workbench**
environment.

---

## R script: Merge and extract helper files
```r
if (!requireNamespace("data.table", quietly = TRUE)) {
install.packages("data.table")
}
library(data.table)

# Base directory containing helper VCF files
base_dir <- file.path(
"/mnt/project", "Bulk", "Imputation",
"Imputation from genotype (TOPmed)", "helper_files"
)

# Utility function to extract INFO tags (robust to missing trailing semicolons)
extract_tag <- function(x, tag) {
m <- regexpr(paste0("(^|;)", tag, "=([^;]+)"), x, perl = TRUE)
out <- rep(NA_character_, length(x))
hit <- m > 0
if (any(hit)) {
  s <- regmatches(x, m)
  out[hit] <- sub(paste0(".*", tag, "=([^;]+).*"), "\\1", s)
}
out
}

process_one_chr <- function(chr) {
vcf_path <- file.path(
  base_dir,
  sprintf("ukb21007_c%d_b0_v1.sites.vcf.gz", chr)
)
if (!file.exists(vcf_path)) {
  message(sprintf("[warn] file not found, skip: %s", vcf_path))
  return(NULL)
}

# Read VCF (remove '##' meta lines)
cmd <- sprintf("zcat %s | awk '!/^##/'", shQuote(vcf_path))
DT <- fread(
  cmd = cmd,
  sep = "\t",
  header = TRUE,
  data.table = TRUE,
  showProgress = FALSE
)

# Rename chromosome column if needed
if ("#CHROM" %in% names(DT)) {
  setnames(DT, "#CHROM", "CHR")
}

# Remove 'chr' prefix
DT[, CHR := sub("^chr", "", CHR)]

# Extract allele frequency and filter
DT[, AF := as.numeric(extract_tag(INFO, "AF"))]
DT <- DT[!is.na(AF) & AF >= 0.001]

# Construct CHR:BP identifier
DT[, original_name := paste0(CHR, ":", POS)]

# Select output columns
out <- DT[, .(
  ID,
  original_name,
  CHR,
  BP = POS,
  REF,
  ALT,
  AF
)]

# Sort by chromosome and position
suppressWarnings(out[, CHR_num := as.integer(CHR)])
setorder(out, CHR_num, BP)
out[, CHR_num := NULL]

message(sprintf("[ok] chr%d: %s rows",
                chr, format(nrow(out), big.mark=",")))
out
}

# Process chromosomes 1–22
res_list <- lapply(1:22, process_one_chr)
res <- rbindlist(res_list, use.names = TRUE, fill = TRUE)

# Write output map file
out_path <- "/home/rstudio-server/map.txt"
fwrite(res, file = out_path, sep = "\t", quote = FALSE)

message(sprintf(
"[done] wrote %s rows to %s",
format(nrow(res), big.mark=","),
out_path
))