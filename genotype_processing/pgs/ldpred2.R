
#!/usr/bin/env Rscript
# ============================================================
# LDpred2 pipeline 
#
# This script:
#   1) Matches GWAS summary statistics to HapMap3+ SNPs
#   2) Builds genome-wide LD matrix from precomputed LDref
#   3) Runs LDpred2-inf and LDpred2-auto
#   4) Computes polygenic scores (PGS) in target genotypes
#
# Requirements:
#   - GWAS summary statistics (hg38 recommended)
#   - Target genotypes in PLINK BED format
#   - HapMap3+ LD reference (LDpred2 format)
#
# Author: Tung-Yen (Ian) Wu
# ============================================================


# ============================================================
# 0) USER SETTINGS — EDIT THIS SECTION ONLY
# ============================================================

CFG <- list(
  
  ## ---- Input data ----
  sumstats_path = "/PATH/TO/gwas.sumstats.gz",       # GWAS summary statistics
  bed_prefix    = "/PATH/TO/target_genotypes",       # without .bed/.bim/.fam
  map_path      = "/PATH/TO/map_hm3_plus.rds",       # HapMap3+ SNP map
  ldref_zip     = "/PATH/TO/ldref_hm3_plus.zip",     # LD reference (zip)
  
  ## ---- Working directories ----
  work_dir      = "/PATH/TO/workdir",
  ldref_dir     = "/PATH/TO/workdir/ldref_hm3_plus",
  bigsnp_prefix = "/PATH/TO/workdir/target_bigSNP",
  
  ## ---- Output ----
  out_inf  = "/PATH/TO/output/PGS_LDpred2_inf.sscore",
  out_auto = "/PATH/TO/output/PGS_LDpred2_auto.sscore",
  ##### Recommend output to "/home/rstudio-server/" and download manually.
  
  ## ---- Model parameters ----
  h2_estimate = 0.0732,    # SNP-heritability (from LDSC)
  burn_in     = 500,
  num_iter    = 1000
)

dir.create(dirname(CFG$out_inf), recursive = TRUE, showWarnings = FALSE)
dir.create(CFG$work_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 1) Packages
# ============================================================

suppressPackageStartupMessages({
  pkgs <- c("bigsnpr","bigreadr","bigstatsr",
            "bigsparser","data.table","dplyr","R.utils")
  for (p in pkgs)
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  
  library(bigsnpr)
  library(bigreadr)
  library(bigstatsr)
  library(bigsparser)
  library(data.table)
  library(dplyr)
  library(R.utils)
})


# ============================================================
# 2) Load GWAS summary statistics
# ============================================================

sumstats <- fread2(CFG$sumstats_path)
stopifnot(all(c("SNP","CHR","BP","REF","ALT","BETA","SE","P") %in% names(sumstats)))

sumstats <- sumstats %>%
  rename(
    rsid = SNP,
    chr  = CHR,
    pos  = BP,
    a0   = REF,
    a1   = ALT,
    beta = BETA,
    beta_se = SE,
    p = P
  )

cat("[info] GWAS variants:", nrow(sumstats), "\n")


# ============================================================
# 3) Match GWAS to HapMap3+ map
# ============================================================

map_ld  <- readRDS(CFG$map_path)
info_snp <- snp_match(sumstats, map_ld,
                      join_by_pos = FALSE,
                      strand_flip = TRUE)

cat("[info] SNPs after GWAS–HM3+ matching:", nrow(info_snp), "\n")
stopifnot("pos_hg38" %in% names(info_snp))


# ============================================================
# 4) Load target genotypes (PLINK BED → bigSNP)
# ============================================================

if (!file.exists(paste0(CFG$bigsnp_prefix, ".rds"))) {
  snp_readBed(paste0(CFG$bed_prefix, ".bed"),
              backingfile = CFG$bigsnp_prefix)
}

obj <- snp_attach(paste0(CFG$bigsnp_prefix, ".rds"))
G   <- obj$genotypes
CHR <- obj$map$chromosome
POS <- obj$map$physical.pos

cat("[info] Target samples:", nrow(G),
    " | SNPs:", ncol(G), "\n")


# ============================================================
# 5) Align genotypes to GWAS SNPs (hg38 positions)
# ============================================================

bim <- fread2(paste0(CFG$bed_prefix, ".bim"), header = FALSE)
setnames(bim, c("CHR","SNP","CM","BP","A1","A2"))

matched_idx <- match(
  paste(bim$CHR, bim$BP),
  paste(info_snp$chr, info_snp$pos_hg38)
)

ind.keep <- which(!is.na(matched_idx))
info_sub <- info_snp[matched_idx[ind.keep], ]

cat("[info] SNPs used for PGS:", length(ind.keep), "\n")


# ============================================================
# 6) Build genome-wide LD matrix (subset HM3+ LDref)
# ============================================================

if (!dir.exists(CFG$ldref_dir)) {
  unzip(CFG$ldref_zip, exdir = CFG$ldref_dir)
}

ld_files <- list.files(
  CFG$ldref_dir,
  pattern = "^LD_with_blocks_chr\\d+\\.rds$",
  full.names = TRUE
)
stopifnot(length(ld_files) == 22)

keep_ids <- sort(unique(info_sub$`_NUM_ID_`))
first <- TRUE

for (chr in 1:22) {
  
  idx_map_chr <- which(map_ld$chr == chr)
  keep_chr    <- intersect(idx_map_chr, keep_ids)
  if (!length(keep_chr)) next
  
  sel_local <- match(keep_chr, idx_map_chr)
  M <- readRDS(ld_files[chr])
  M_sub <- M[sel_local, sel_local, drop = FALSE]
  
  if (first) {
    corr <- as_SFBM(M_sub,
                    backingfile = file.path(CFG$work_dir, "corr_ldpred2"),
                    compact = TRUE)
    first <- FALSE
  } else {
    corr$add_columns(M_sub, nrow(corr))
  }
  
  cat(sprintf("[info] chr%02d: kept %d SNPs\n",
              chr, length(sel_local)))
}

stopifnot(ncol(corr) == length(keep_ids))


# ============================================================
# 7) LDpred2-inf
# ============================================================

df_beta <- data.frame(
  beta    = info_sub$beta,
  beta_se = info_sub$beta_se,
  n_eff   = info_sub$N
)

beta_inf <- snp_ldpred2_inf(
  corr = corr,
  df_beta = df_beta,
  h2 = CFG$h2_estimate
)


# ============================================================
# 8) Compute PGS (LDpred2-inf)
# ============================================================

G2 <- snp_fastImputeSimple(G, method = "mean2")
pos_in_keep <- match(info_sub$`_NUM_ID_`, keep_ids)

scores_inf <- big_prodVec(
  G2,
  beta_inf[pos_in_keep],
  ind.col = ind.keep
)

pgs_inf <- data.frame(
  FID = obj$fam$family.ID,
  IID = obj$fam$sample.ID,
  PGS_LDpred2_inf = scores_inf
)

fwrite(pgs_inf, CFG$out_inf, sep = "\t")
cat("[done] LDpred2-inf PGS written\n")


# ============================================================
# 9) LDpred2-auto
# ============================================================

df_beta_auto <- data.frame(
  beta    = rep(0, length(keep_ids)),
  beta_se = rep(1e9, length(keep_ids)),
  n_eff   = rep(NA_real_, length(keep_ids))
)

pos <- match(info_snp$`_NUM_ID_`, keep_ids)
ok  <- !is.na(pos)

df_beta_auto$beta[pos[ok]]    <- info_snp$beta[ok]
df_beta_auto$beta_se[pos[ok]] <- info_snp$beta_se[ok]
if ("N" %in% names(info_snp))
  df_beta_auto$n_eff[pos[ok]] <- info_snp$N[ok]

beta_auto <- snp_ldpred2_auto(
  corr = corr,
  df_beta = df_beta_auto,
  h2_init = CFG$h2_estimate,
  burn_in = CFG$burn_in,
  num_iter = CFG$num_iter
)

beta_auto_mat  <- do.call(cbind, lapply(beta_auto, `[[`, "beta_est"))
beta_auto_mean <- rowMeans(beta_auto_mat, na.rm = TRUE)
beta_auto_mean[is.na(beta_auto_mean)] <- 0


# ============================================================
# 10) Compute PGS (LDpred2-auto)
# ============================================================

w_auto <- beta_auto_mean[pos_in_keep]

scores_auto <- big_prodVec(
  G2,
  w_auto,
  ind.col = ind.keep
)

pgs_auto <- data.frame(
  FID = obj$fam$family.ID,
  IID = obj$fam$sample.ID,
  PGS_LDpred2_auto = scores_auto
)

fwrite(pgs_auto, CFG$out_auto, sep = "\t")
cat("[done] LDpred2-auto PGS written\n")