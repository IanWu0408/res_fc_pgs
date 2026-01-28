

#!/usr/bin/env Rscript
# ============================================================
# Prepare GWAS summary statistics for LDSC / LDpred2
#
# This script:
#   1) Reads GWAS results (e.g. BOLT-LMM output)
#   2) Harmonizes effect direction and allele frequency
#   3) Aligns SNP IDs and alleles using a helper map file
#   4) Computes Z-scores
#   5) Exports LDSC-ready summary statistics (.gz)
#
# IMPORTANT:
#   - map.txt is generated from helper files
#     (e.g. reference SNP annotation / liftover helpers)
#
# Author: Tung-Yen Wu
# ============================================================


# ============================================================
# 0) USER CONFIG — EDIT THIS SECTION ONLY
# ============================================================

CFG <- list(
  
  ## ---- Input GWAS ----
  gwas_path = "/PATH/TO/resilience_dosage.stats.gz",
  
  ## ---- Helper map file ----
  # map.txt is derived from helper files
  # (e.g. SNP ID mapping / genome build harmonization helpers)
  map_path  = "/PATH/TO/map.txt",
  
  ## ---- Output ----
  out_name  = "resilience_ldsc_ready.sumstats.gz",
  
  ## ---- Sample size ----
  # Used to compute effective N = N_total * INFO
  n_total   = 130394
)


# ============================================================
# 1) Packages
# ============================================================

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table")
  library(data.table)
})


# ============================================================
# 2) Load GWAS summary statistics
# ============================================================

dat <- fread(CFG$gwas_path)
setDT(dat)

cat("[info] GWAS rows:", format(nrow(dat), big.mark=","), "\n")


# ============================================================
# 3) Harmonize key fields
# ============================================================

# P-value (BOLT-LMM infinitesimal model)
dat[, P := P_BOLT_LMM_INF]

# Effective sample size
dat[, N := as.integer(CFG$n_total * INFO)]

# Flip effect direction (project-specific convention)
dat[, BETA := -BETA]

# Convert A1 frequency to ALT frequency
dat[, FREQ := 1 - A1FREQ]

# Alignment key: CHR:BP
dat[, CHR_BP := paste(CHR, BP, sep = ":")]


# ============================================================
# 4) Load helper map file (SNP ID & allele reference)
# ============================================================

map <- fread(CFG$map_path)

needed <- c("original_name", "ID", "REF", "ALT")
stopifnot(all(needed %in% names(map)))

cat("[info] Helper map rows:",
    format(nrow(map), big.mark=","), "\n")


# ============================================================
# 5) Align SNP IDs using CHR:BP
# ============================================================

idx <- match(dat$CHR_BP, map$original_name)

message(sprintf(
  "[align] GWAS rows: %s | matched by CHR:BP: %s",
  format(nrow(dat), big.mark=","), format(sum(!is.na(idx)), big.mark=",")
))

# Overwrite SNP ID for matched variants
dat[!is.na(idx), SNP := map$ID[idx]]


# ============================================================
# 6) Overwrite alleles using helper map
# ============================================================

dat[!is.na(idx), `:=`(
  ALLELE0 = map$REF[idx],
  ALLELE1 = map$ALT[idx]
)]

setnames(dat,
         old = c("ALLELE0", "ALLELE1"),
         new = c("REF", "ALT"))


# ============================================================
# 7) Compute Z-scores
# ============================================================

dat[, Z := fifelse(
  is.na(SE) | SE == 0,
  NA_real_,
  BETA / SE
)]


# ============================================================
# 8) Select LDSC-ready columns & export
# ============================================================

out <- dat[, .(
  SNP, CHR, BP,
  REF, ALT,
  N, BETA, SE, Z,
  FREQ, INFO, P
)]

out_path <- file.path(dirname(CFG$gwas_path), CFG$out_name)

fwrite(
  out,
  file = out_path,
  sep = "\t",
  quote = FALSE,
  compress = "gzip"
)

message(sprintf(
  "[done] wrote %s | rows: %s | cols: %s",
  out_path,
  format(nrow(out), big.mark=","),
  ncol(out)
))
