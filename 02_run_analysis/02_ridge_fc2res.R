
#!/usr/bin/env Rscript
# ------------------------------------------------------------
# ONLY configuration:
#   - paths / parameters
#   - ROI groups
#   - model definitions (rules)
# No pipeline logic here.
#
# Author: Tung-Yen Wu
# ------------------------------------------------------------

# ------------------------------------------------------------
# 1) Fix working directory (this script's folder)
# ------------------------------------------------------------
get_script_dir <- function() {
  # Case 1: Rscript
  if (!is.null(sys.frame(1)$ofile)) {
    return(dirname(normalizePath(sys.frame(1)$ofile)))
  }
  
  # Case 2: RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    if (rstudioapi::isAvailable()) {
      ctx <- rstudioapi::getSourceEditorContext()
      if (!is.null(ctx$path) && nzchar(ctx$path)) {
        return(dirname(normalizePath(ctx$path)))
      }
    }
  }
}

RUN_DIR <- get_script_dir()
setwd(RUN_DIR)
message("[INFO] Working dir: ", getwd())


# ------------------------------------------------------------
# 2) Load core pipeline + config
# ------------------------------------------------------------
source("../01_define_analysis/02_ridge.R")

# ------------------------------------------------------------
# 3) Run
# ------------------------------------------------------------

options(stringsAsFactors = FALSE)

CFG <- list(
  # base_dir is 02_run_analysis/ because we setwd(RUN_DIR) in runner
  base_dir = ".",
  
  # results go under 02_run_analysis/results/...
  out_parent_dir = file.path("results", "RIDGE_GLMNET"),
  run_tag = "ROI_v2",
  
  # inputs are under 02_run_analysis/data/...
  brain_cohort_file = file.path("data", "brain_cohort.txt"),
  cov_file          = file.path("data", "cov", "cov.csv"),
  fc_file_gz        = file.path("data", "RAW_FC", "UKB_FC_Glasser_Tian_S2_allSubjects_10_60_merged.txt.gz"),
  fc_file_txt       = file.path("data", "RAW_FC", "UKB_FC_Glasser_Tian_S2_allSubjects_10_60_merged.txt"),
  mapping_file      = file.path("data", "RAW_FC", "UKB_FC_Glasser_Tian_S2_feature_mapping.csv"),
  
  # Column names in your files
  brain_id_col  = "IID",
  cov_id_col    = "eid",
  outcome_var   = "INTresilience",
  
  # Brain cohort covariates to include (model.matrix will handle factors)
  brain_cov = c(paste0("PC", 1:20), "Age", "Age2", "Sex"),
  
  # Extra covariates from cov/cov.csv
  cov_cont = c(
    "p25000_i2","p25756_i2","p25757_i2","p25758_i2","p25759_i2",
    "p24432_i2","p24434_i2","p25739_i2","p25744_i2","p24438_i2"
  ),
  cov_cat  = c("p24410_i2"),
  
  # CV
  seed   = 20260117,
  nfolds = 10,
  
  # Ridge (glmnet)
  ridge_scale_X = TRUE,              # scale FC within outer fold using TRAIN stats only
  ridge_tune_lambda = TRUE,          # nested inner CV to choose lambda
  ridge_inner_folds = 5,
  ridge_lambda_list = 10^seq(-4, 2, length.out = 13),
  ridge_lambda_fixed = 1e-2,
  ridge_report_lambda = "lambda.min", # "lambda.min" or "lambda.1se"
  ridge_type_measure  = "mse",
  
  # Global FC complete-case policy
  # - "cc_models": union of FC features from models with cc_include=TRUE
  # - "allFC": all mapping features (required if you enable a type="allFC" model)
  global_fc_cc_scope = "cc_models",
  
  # FC big file reading
  fc_scan_chunk_n = 25000,
  fc_max_chunk_cells = 5e9,
  fc_scan_show_progress = TRUE,
  fc_na_strings = c("NA","NaN","nan","Inf","-Inf"),
  
  # Cache FC matrix for ROI-based models (pool features) in RAM once
  fc_cache_pool_matrix = TRUE,
  
  # Decompress .gz -> .txt if fread cannot read gz directly
  decompress_if_needed = TRUE,
  
  # Dev
  test_mode = FALSE,
  test_n    = 100,
  
  # Threads
  dt_threads = 16,
  
  # Saving
  save_inner_cv_table_per_fold = TRUE,
  save_ridge_fold_rds = TRUE
)

# ------------------------------------------------------------
# ROI GROUPS
# kind:
#   - "glasser": tokens -> Glasser_L_*_ROI + Glasser_R_*_ROI
#   - "tian":    labels -> Tian_<label>
#   - "union":   union of other ROI_GROUPS (by name)
# ------------------------------------------------------------
ROI_GROUPS <- list(
  # -------- Glasser cortical tokens --------
  INS = list(kind = "glasser", tokens = c("AAIC", "AVI", "MI", "PI", "PoI1", "PoI2")),
  ACC = list(kind = "glasser", tokens = c("s32","a24","p24","p32","d32","a32pr","a24pr","33pr","p32pr","p24pr","24dd","24dv")),
  
  DMN_vmPFC     = list(kind="glasser", tokens=c("10v","10r","25")),
  DMN_dmPFC     = list(kind="glasser", tokens=c("10d","9m","8BM","10pp","a10p","p10p")),
  DMN_PCC       = list(kind="glasser", tokens=c("31pv","31pd","d23ab","v23ab")),
  DMN_Precuneus = list(kind="glasser", tokens=c("7m","7Pm","7PC")),
  DMN_Angular   = list(kind="glasser", tokens=c("PGi","PGp","PGs")),
  DMN           = list(kind="union", members=c("DMN_vmPFC","DMN_dmPFC","DMN_PCC","DMN_Precuneus","DMN_Angular")),
  
  CEN_dlPFC = list(kind="glasser", tokens=c("46","9-46d","a9-46v","p9-46v","9a","9p","8Ad","8Av","8BL","8C","i6-8","s6-8","SFL")),
  CEN_IPS   = list(kind="glasser", tokens=c("IP0","IP1","IP2","AIP","MIP","VIP","LIPd","LIPv")),
  CEN_SPL   = list(kind="glasser", tokens=c("7AL","7Am","7PL")),
  CEN       = list(kind="union", members=c("CEN_dlPFC","CEN_IPS","CEN_SPL")),
  
  # PFC (re-using some DMN/CEN parts)
  PFC_OFC   = list(kind="glasser", tokens=c("47s","47m","11l","13l","OFC","pOFC","47l","a47r","p47r")),
  PFC_vmPFC = list(kind="union", members=c("DMN_vmPFC")),
  PFC_dmPFC = list(kind="union", members=c("DMN_dmPFC")),
  PFC_dlPFC = list(kind="union", members=c("CEN_dlPFC")),
  PFC_vlPFC = list(kind="glasser", tokens=c("44","45","IFJa","IFJp","IFSa","IFSp")),
  PFC_ACC   = list(kind="union", members=c("ACC")),
  PFC       = list(kind="union", members=c("PFC_OFC","PFC_vmPFC","PFC_dmPFC","PFC_dlPFC","PFC_vlPFC","PFC_ACC")),
  
  PCC       = list(kind="union", members=c("DMN_PCC")),
  Precuneus = list(kind="union", members=c("DMN_Precuneus")),
  
  # Convenience super-set used by "Triple-Network"
  TRIPLE = list(kind="union", members=c("INS","ACC","DMN","CEN")),
  
  # -------- Tian subcortical labels --------
  AMYG = list(kind="tian", labels=c("lAMY-rh","mAMY-rh","lAMY-lh","mAMY-lh")),
  NAc  = list(kind="tian", labels=c("NAc-shell-rh","NAc-core-rh","NAc-shell-lh","NAc-core-lh")),
  HIP  = list(kind="tian", labels=c("aHIP-lh","pHIP-lh","aHIP-rh","pHIP-rh")),
  THA  = list(kind="tian", labels=c("THA-DP-lh","THA-VP-lh","THA-VA-lh","THA-DA-lh",
                                    "THA-DP-rh","THA-VP-rh","THA-VA-rh","THA-DA-rh")),
  CAU  = list(kind="tian", labels=c("aCAU-rh","pCAU-rh","aCAU-lh","pCAU-lh"))
)

# ------------------------------------------------------------
# MODEL_SPECS (edit here to add/remove models)
#
# type:
#   - "cov_only" : baseline OLS covariates only
#   - "ridge_fc" : covariates + FC (glmnet ridge); FC defined by rules
#   - "union"    : union of other models' FC features (still ridge)
#   - "allFC"    : covariates + ALL FC features (memory heavy)
#
# rules:
#   - within: op="within", set="ROI_GROUP_NAME"
#   - cross : op="cross",  setA="ROI_GROUP_A", setB="ROI_GROUP_B"
#
# cc_include:
#   if TRUE, this model's FC features contribute to GLOBAL FC complete-case
# ------------------------------------------------------------
MODEL_SPECS <- list(
  cov_only = list(
    enabled = TRUE,
    type = "cov_only",
    cc_include = FALSE
  ),
  
  cov_triple = list(
    enabled = TRUE,
    type = "ridge_fc",
    cc_include = TRUE,
    rules = list(
      list(op="within", set="TRIPLE")
    )
  ),
  
  cov_appraisals = list(
    enabled = TRUE,
    type = "ridge_fc",
    cc_include = TRUE,
    rules = list(
      list(op="cross", setA="PFC",       setB="AMYG"),
      list(op="cross", setA="INS",       setB="AMYG"),
      list(op="cross", setA="ACC",       setB="AMYG"),
      list(op="cross", setA="PCC",       setB="AMYG"),
      list(op="cross", setA="Precuneus", setB="AMYG"),
      list(op="within", set="AMYG"),
      list(op="cross", setA="HIP", setB="AMYG"),
      list(op="cross", setA="THA", setB="AMYG")
    )
  ),
  
  cov_coping = list(
    enabled = TRUE,
    type = "ridge_fc",
    cc_include = TRUE,
    rules = list(
      list(op="cross", setA="PFC", setB="NAc"),
      list(op="cross", setA="PFC", setB="CAU"),
      list(op="cross", setA="INS", setB="NAc"),
      list(op="cross", setA="INS", setB="CAU"),
      list(op="cross", setA="ACC", setB="NAc"),
      list(op="cross", setA="ACC", setB="CAU"),
      list(op="cross", setA="NAc", setB="CAU"),
      list(op="within", set="NAc"),
      list(op="within", set="CAU"),
      list(op="cross", setA="HIP", setB="NAc"),
      list(op="cross", setA="HIP", setB="CAU"),
      list(op="cross", setA="THA", setB="NAc"),
      list(op="cross", setA="THA", setB="CAU")
    )
  ),
  
  cov_appraisals_coping = list(
    enabled = TRUE,
    type = "union",
    cc_include = TRUE,
    members = c("cov_appraisals", "cov_coping")
  ),
  
  cov_triple_appraisals_coping = list(
    enabled = TRUE,
    type = "union",
    cc_include = TRUE,
    members = c("cov_triple", "cov_appraisals", "cov_coping")
  ),
  
  # VERY MEMORY HEAVY. If enabled=TRUE, you must set CFG$global_fc_cc_scope="allFC"
  cov_allFC = list(
    enabled = FALSE,
    type = "allFC",
    cc_include = FALSE
  )
)